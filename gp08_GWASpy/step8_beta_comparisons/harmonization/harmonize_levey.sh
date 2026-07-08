#!/bin/bash
set -euo pipefail

cat > /tmp/harmonize_levey_chunked.py << 'EOF'
import pandas as pd
from pathlib import Path

levey_dir = Path("~/sud/levey-2022/liftover_only").expanduser()
khat_dir = Path("~/sud/khat").expanduser()
out_dir = Path("~/rsid_match/levey_lifted_harmonized").expanduser()
out_dir.mkdir(parents=True, exist_ok=True)

ku_file = khat_dir / "ku_gp08_meta_cleaned_with_N.txt"
kuf_file = khat_dir / "kuf_se_meta.txt"

levey_files = [
    "AFR_w_rsid_GRCh37.cleaned.tsv.lifted_hg38.tsv",
    "EUR_w_rsid_GRCh37.cleaned.tsv.lifted_hg38.tsv"
]

def is_snv(a):
    return isinstance(a, str) and len(a) == 1 and a.upper() in {"A", "C", "G", "T"}

def is_ambiguous(a1, a2):
    pair = {str(a1).upper(), str(a2).upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}

def complement_allele(a):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(str(a).upper(), None)

def prep_khat(path, beta_col, p_col):
    df = pd.read_csv(path, sep="\t", dtype={"MarkerName": str}, low_memory=False)

    df = df[
        ["MarkerName", "CHR", "BP", "Allele1", "Allele2", "Effect", "P-value"]
    ].rename(columns={
        "MarkerName": "khat_marker",
        "Allele1": "khat_A1",
        "Allele2": "khat_A2",
        "Effect": beta_col,
        "P-value": p_col
    })

    df["CHR"] = df["CHR"].astype(str)
    df["BP"] = df["BP"].astype(int)
    df["khat_A1"] = df["khat_A1"].astype(str).str.upper()
    df["khat_A2"] = df["khat_A2"].astype(str).str.upper()

    df = df[
        df["khat_A1"].apply(is_snv) &
        df["khat_A2"].apply(is_snv)
    ].copy()

    df = df[
        ~df.apply(lambda r: is_ambiguous(r["khat_A1"], r["khat_A2"]), axis=1)
    ].copy()

    df = df.drop_duplicates(subset=["CHR", "BP", "khat_A1", "khat_A2"])

    return df

def harmonize_chunk(khat, levey, trait, khat_beta_col, khat_p_col):
    merged = khat.merge(levey, on=["CHR", "BP"], how="inner")

    if len(merged) == 0:
        return merged

    merged["levey_A1_comp"] = merged["levey_A1"].apply(complement_allele)
    merged["levey_A2_comp"] = merged["levey_A2"].apply(complement_allele)

    exact_match = (
        (merged["khat_A1"] == merged["levey_A1"]) &
        (merged["khat_A2"] == merged["levey_A2"])
    )

    swapped_match = (
        (merged["khat_A1"] == merged["levey_A2"]) &
        (merged["khat_A2"] == merged["levey_A1"])
    )

    complement_match = (
        (merged["khat_A1"] == merged["levey_A1_comp"]) &
        (merged["khat_A2"] == merged["levey_A2_comp"])
    )

    complement_swapped_match = (
        (merged["khat_A1"] == merged["levey_A2_comp"]) &
        (merged["khat_A2"] == merged["levey_A1_comp"])
    )

    keep = exact_match | swapped_match | complement_match | complement_swapped_match

    merged = merged.loc[keep].copy()

    if len(merged) == 0:
        return merged

    exact_match = exact_match.loc[merged.index]
    swapped_match = swapped_match.loc[merged.index]
    complement_match = complement_match.loc[merged.index]
    complement_swapped_match = complement_swapped_match.loc[merged.index]

    beta_raw_col = f"beta_{trait}"
    beta_harmonized_col = f"beta_{trait}_to_khat_A1"

    merged[beta_raw_col] = merged["BETA"]
    merged[beta_harmonized_col] = merged["BETA"]

    flip_mask = swapped_match | complement_swapped_match
    merged.loc[flip_mask, beta_harmonized_col] = -1 * merged.loc[flip_mask, "BETA"]

    merged["allele_match_status"] = "exact"
    merged.loc[swapped_match, "allele_match_status"] = "swapped"
    merged.loc[complement_match, "allele_match_status"] = "complement"
    merged.loc[complement_swapped_match, "allele_match_status"] = "complement_swapped"

    keep_cols = [
        "khat_marker", "CHR", "BP", "khat_A1", "khat_A2",
        khat_beta_col, khat_p_col,
        "levey_marker", "levey_A1", "levey_A2",
        beta_raw_col, beta_harmonized_col,
        "pval", "SE", "ancestry", "allele_match_status"
    ]

    for col in ["Weight", "Zscore", "SNP", "rsid"]:
        if col in merged.columns:
            keep_cols.append(col)

    return merged[keep_cols]

ku = prep_khat(ku_file, "beta_ku", "p_ku")
kuf = prep_khat(kuf_file, "beta_kuf", "p_kuf")

for filename in levey_files:
    trait = filename.replace(".lifted_hg38.tsv", "").replace(".cleaned.tsv", "").replace("_w_rsid_GRCh37", "")

    print(f"Processing {trait}", flush=True)

    ku_out = out_dir / f"ku_{trait}_harmonized.txt"
    kuf_out = out_dir / f"kuf_{trait}_harmonized.txt"

    for f in [ku_out, kuf_out]:
        if f.exists():
            f.unlink()

    first_ku = True
    first_kuf = True
    total_ku = 0
    total_kuf = 0

    reader = pd.read_csv(
        levey_dir / filename,
        sep="\t",
        dtype={"rsid": str, "SNP": str},
        chunksize=500000,
        low_memory=False
    )

    for chunk_i, df in enumerate(reader):
        df["CHR"] = df["CHR_hg38"].astype(str)
        df["BP"] = df["BP_hg38"].astype(int)
        df["levey_marker"] = df["rsid"].astype(str)
        df["levey_A1"] = df["a1"].astype(str).str.upper()
        df["levey_A2"] = df["a2"].astype(str).str.upper()

        df = df[
            df["levey_A1"].apply(is_snv) &
            df["levey_A2"].apply(is_snv)
        ].copy()

        df = df[
            ~df.apply(lambda r: is_ambiguous(r["levey_A1"], r["levey_A2"]), axis=1)
        ].copy()

        df = df.drop_duplicates(subset=["CHR", "BP", "levey_A1", "levey_A2"])

        ku_match = harmonize_chunk(ku, df, trait, "beta_ku", "p_ku")
        kuf_match = harmonize_chunk(kuf, df, trait, "beta_kuf", "p_kuf")

        if len(ku_match) > 0:
            ku_match.to_csv(ku_out, sep="\t", index=False, mode="w" if first_ku else "a", header=first_ku)
            first_ku = False
            total_ku += len(ku_match)

        if len(kuf_match) > 0:
            kuf_match.to_csv(kuf_out, sep="\t", index=False, mode="w" if first_kuf else "a", header=first_kuf)
            first_kuf = False
            total_kuf += len(kuf_match)

        print(f"  chunk {chunk_i}: KU {len(ku_match):,}, KUF {len(kuf_match):,}", flush=True)

    print(f"Finished {trait}: KU {total_ku:,}, KUF {total_kuf:,}", flush=True)
    print(f"  Wrote: {ku_out}", flush=True)
    print(f"  Wrote: {kuf_out}", flush=True)

print("Done.", flush=True)
EOF

python3 /tmp/harmonize_levey_chunked.py
