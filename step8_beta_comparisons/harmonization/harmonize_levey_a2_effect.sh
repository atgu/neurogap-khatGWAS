#!/bin/bash
set -euo pipefail

cat > /tmp/harmonize_levey_chunked.py <<'EOF'
import pandas as pd
from pathlib import Path

# -------------------------
# Input and output paths
# -------------------------

levey_dir = Path(
    "~/sud/levey-2022/liftover_only"
).expanduser()

khat_dir = Path(
    "~/sud/khat/a2_effect"
).expanduser()

out_dir = Path(
    "~/rsid_match/levey_lifted_harmonized/a2_effect"
).expanduser()

out_dir.mkdir(parents=True, exist_ok=True)


ku_file = khat_dir / "ku_gp08_meta_cleaned_with_N.txt"
kuf_file = khat_dir / "kuf_gp08_meta_cleaned_with_N.tsv"

levey_files = [
    "AFR_w_rsid_GRCh37.cleaned.tsv.lifted_hg38.tsv",
    "EUR_w_rsid_GRCh37.cleaned.tsv.lifted_hg38.tsv",
]

# -------------------------
# Validation
# -------------------------

required_files = [ku_file, kuf_file]

for filename in levey_files:
    required_files.append(levey_dir / filename)

for path in required_files:
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")

# -------------------------
# Allele helper functions
# -------------------------

def is_snv(allele):
    return (
        isinstance(allele, str)
        and len(allele) == 1
        and allele.upper() in {"A", "C", "G", "T"}
    )


def is_ambiguous(allele1, allele2):
    pair = {str(allele1).upper(), str(allele2).upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def complement_allele(allele):
    return {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }.get(str(allele).upper())


# -------------------------
# Prepare KU/KUF results
# -------------------------

def prep_khat(path, beta_col, p_col):
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={"MarkerName": str},
        low_memory=False,
    )

    required_columns = [
        "MarkerName",
        "CHR",
        "BP",
        "Allele1",
        "Allele2",
        "Effect",
        "P-value",
    ]

    missing_columns = [
        column for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing columns in {path}: "
            f"{', '.join(missing_columns)}"
        )

    df = df[required_columns].rename(
        columns={
            "MarkerName": "khat_marker",
            "Allele1": "khat_A1",
            "Allele2": "khat_A2",
            "Effect": beta_col,
            "P-value": p_col,
        }
    )

    # In the final METAL output:
    # khat_A1 = METAL Allele1 = effect allele
    # khat_A2 = METAL Allele2 = other allele

    df["CHR"] = (
        df["CHR"]
        .astype(str)
        .str.replace(r"\.0$", "", regex=True)
        .str.replace("^chr", "", case=False, regex=True)
    )

    df["BP"] = pd.to_numeric(
        df["BP"],
        errors="coerce",
    )

    df[beta_col] = pd.to_numeric(
        df[beta_col],
        errors="coerce",
    )

    df[p_col] = pd.to_numeric(
        df[p_col],
        errors="coerce",
    )

    df = df.dropna(
        subset=["CHR", "BP", beta_col, p_col]
    ).copy()

    df["BP"] = df["BP"].astype(int)

    df["khat_A1"] = (
        df["khat_A1"]
        .astype(str)
        .str.upper()
    )

    df["khat_A2"] = (
        df["khat_A2"]
        .astype(str)
        .str.upper()
    )

    # Retain only biallelic SNVs.
    df = df[
        df["khat_A1"].apply(is_snv)
        & df["khat_A2"].apply(is_snv)
    ].copy()

    # Remove strand-ambiguous A/T and C/G variants.
    df = df[
        ~df.apply(
            lambda row: is_ambiguous(
                row["khat_A1"],
                row["khat_A2"],
            ),
            axis=1,
        )
    ].copy()

    df = df.drop_duplicates(
        subset=[
            "CHR",
            "BP",
            "khat_A1",
            "khat_A2",
        ]
    )

    return df


# -------------------------
# Harmonize one Levey chunk
# -------------------------

def harmonize_chunk(
    khat,
    levey,
    trait,
    khat_beta_col,
    khat_p_col,
):
    merged = khat.merge(
        levey,
        on=["CHR", "BP"],
        how="inner",
    )

    if merged.empty:
        return merged

    merged["levey_A1_comp"] = (
        merged["levey_A1"]
        .apply(complement_allele)
    )

    merged["levey_A2_comp"] = (
        merged["levey_A2"]
        .apply(complement_allele)
    )

    exact_match = (
        (merged["khat_A1"] == merged["levey_A1"])
        & (merged["khat_A2"] == merged["levey_A2"])
    )

    swapped_match = (
        (merged["khat_A1"] == merged["levey_A2"])
        & (merged["khat_A2"] == merged["levey_A1"])
    )

    complement_match = (
        (
            merged["khat_A1"]
            == merged["levey_A1_comp"]
        )
        & (
            merged["khat_A2"]
            == merged["levey_A2_comp"]
        )
    )

    complement_swapped_match = (
        (
            merged["khat_A1"]
            == merged["levey_A2_comp"]
        )
        & (
            merged["khat_A2"]
            == merged["levey_A1_comp"]
        )
    )

    keep = (
        exact_match
        | swapped_match
        | complement_match
        | complement_swapped_match
    )

    merged = merged.loc[keep].copy()

    if merged.empty:
        return merged

    exact_match = exact_match.loc[merged.index]
    swapped_match = swapped_match.loc[merged.index]
    complement_match = complement_match.loc[merged.index]
    complement_swapped_match = (
        complement_swapped_match.loc[merged.index]
    )

    beta_raw_col = f"beta_{trait}"
    beta_harmonized_col = (
        f"beta_{trait}_to_khat_A1"
    )

    merged[beta_raw_col] = pd.to_numeric(
        merged["BETA"],
        errors="coerce",
    )

    merged[beta_harmonized_col] = (
        merged[beta_raw_col]
    )

    # Flip Levey beta when its alleles are reversed relative
    # to the Khat METAL effect allele.
    flip_mask = (
        swapped_match
        | complement_swapped_match
    )

    merged.loc[
        flip_mask,
        beta_harmonized_col,
    ] = (
        -1
        * merged.loc[
            flip_mask,
            beta_raw_col,
        ]
    )

    merged["allele_match_status"] = "exact"

    merged.loc[
        swapped_match,
        "allele_match_status",
    ] = "swapped"

    merged.loc[
        complement_match,
        "allele_match_status",
    ] = "complement"

    merged.loc[
        complement_swapped_match,
        "allele_match_status",
    ] = "complement_swapped"

    keep_cols = [
        "khat_marker",
        "CHR",
        "BP",
        "khat_A1",
        "khat_A2",
        khat_beta_col,
        khat_p_col,
        "levey_marker",
        "levey_A1",
        "levey_A2",
        beta_raw_col,
        beta_harmonized_col,
        "pval",
        "SE",
        "ancestry",
        "allele_match_status",
    ]

    for column in [
        "Weight",
        "Zscore",
        "SNP",
        "rsid",
    ]:
        if column in merged.columns:
            keep_cols.append(column)

    return merged[keep_cols]


# -------------------------
# Load KU and KUF results
# -------------------------

print(f"Reading KU results: {ku_file}", flush=True)
ku = prep_khat(
    ku_file,
    "beta_ku",
    "p_ku",
)

print(f"Reading KUF results: {kuf_file}", flush=True)
kuf = prep_khat(
    kuf_file,
    "beta_kuf",
    "p_kuf",
)

print(
    f"Prepared {len(ku):,} KU variants and "
    f"{len(kuf):,} KUF variants",
    flush=True,
)

# -------------------------
# Process Levey ancestries
# -------------------------

for filename in levey_files:
    trait = (
        filename
        .replace(".lifted_hg38.tsv", "")
        .replace(".cleaned.tsv", "")
        .replace("_w_rsid_GRCh37", "")
    )

    input_path = levey_dir / filename

    print(f"Processing {trait}", flush=True)
    print(f"  Input: {input_path}", flush=True)

    ku_out = (
        out_dir
        / f"ku_{trait}_harmonized.txt"
    )

    kuf_out = (
        out_dir
        / f"kuf_{trait}_harmonized.txt"
    )

    for output_file in [ku_out, kuf_out]:
        if output_file.exists():
            output_file.unlink()

    first_ku = True
    first_kuf = True

    total_ku = 0
    total_kuf = 0

    reader = pd.read_csv(
        input_path,
        sep="\t",
        dtype={
            "rsid": str,
            "SNP": str,
        },
        chunksize=500000,
        low_memory=False,
    )

    for chunk_i, df in enumerate(reader, start=1):
        required_levey_columns = [
            "CHR_hg38",
            "BP_hg38",
            "rsid",
            "a1",
            "a2",
            "BETA",
            "pval",
            "SE",
            "ancestry",
        ]

        missing_columns = [
            column for column in required_levey_columns
            if column not in df.columns
        ]

        if missing_columns:
            raise ValueError(
                f"Missing columns in {input_path}: "
                f"{', '.join(missing_columns)}"
            )

        df["CHR"] = (
            df["CHR_hg38"]
            .astype(str)
            .str.replace(r"\.0$", "", regex=True)
            .str.replace("^chr", "", case=False, regex=True)
        )

        df["BP"] = pd.to_numeric(
            df["BP_hg38"],
            errors="coerce",
        )

        df = df.dropna(
            subset=["CHR", "BP"]
        ).copy()

        df["BP"] = df["BP"].astype(int)

        df["levey_marker"] = (
            df["rsid"]
            .astype(str)
        )

        df["levey_A1"] = (
            df["a1"]
            .astype(str)
            .str.upper()
        )

        df["levey_A2"] = (
            df["a2"]
            .astype(str)
            .str.upper()
        )

        # Retain only biallelic SNVs.
        df = df[
            df["levey_A1"].apply(is_snv)
            & df["levey_A2"].apply(is_snv)
        ].copy()

        # Remove strand-ambiguous A/T and C/G variants.
        df = df[
            ~df.apply(
                lambda row: is_ambiguous(
                    row["levey_A1"],
                    row["levey_A2"],
                ),
                axis=1,
            )
        ].copy()

        df = df.drop_duplicates(
            subset=[
                "CHR",
                "BP",
                "levey_A1",
                "levey_A2",
            ]
        )

        ku_match = harmonize_chunk(
            ku,
            df,
            trait,
            "beta_ku",
            "p_ku",
        )

        kuf_match = harmonize_chunk(
            kuf,
            df,
            trait,
            "beta_kuf",
            "p_kuf",
        )

        if not ku_match.empty:
            ku_match.to_csv(
                ku_out,
                sep="\t",
                index=False,
                mode="w" if first_ku else "a",
                header=first_ku,
            )

            first_ku = False
            total_ku += len(ku_match)

        if not kuf_match.empty:
            kuf_match.to_csv(
                kuf_out,
                sep="\t",
                index=False,
                mode="w" if first_kuf else "a",
                header=first_kuf,
            )

            first_kuf = False
            total_kuf += len(kuf_match)

        print(
            f"  Chunk {chunk_i}: "
            f"KU {len(ku_match):,}, "
            f"KUF {len(kuf_match):,}",
            flush=True,
        )

    print(
        f"Finished {trait}: "
        f"KU {total_ku:,}, "
        f"KUF {total_kuf:,}",
        flush=True,
    )

    print(f"  Wrote: {ku_out}", flush=True)
    print(f"  Wrote: {kuf_out}", flush=True)

print("Done.", flush=True)
EOF

python3 /tmp/harmonize_levey_chunked.py
