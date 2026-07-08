#!/bin/bash
set -euo pipefail

cat > /tmp/harmonize_hatoum_european_lifted.py << 'EOF'
import pandas as pd
from pathlib import Path

hatoum_file = Path("~/sud/hatoum-2023/liftover_only/Hatoum2023AddictionEuropean.txt.lifted_hg38.tsv").expanduser()
khat_dir = Path("~/sud/khat").expanduser()
out_dir = Path("~/rsid_match/hatoum_lifted_harmonized").expanduser()

out_dir.mkdir(parents=True, exist_ok=True)

ku_file = khat_dir / "ku_gp08_meta_cleaned_with_N.txt"
kuf_file = khat_dir / "kuf_se_meta.txt"

trait = "Hatoum2023AddictionEuropean"

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
    ].rename(
        columns={
            "MarkerName": "khat_marker",
            "Allele1": "khat_A1",
            "Allele2": "khat_A2",
            "Effect": beta_col,
            "P-value": p_col
        }
    )

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

def prep_hatoum(path):
    df = pd.read_csv(path, sep="\t", dtype={"SNP": str}, low_memory=False)

    df["CHR"] = df["CHR_hg38"].astype(str)
    df["BP"] = df["BP_hg38"].astype(int)

    df["hatoum_marker"] = df["SNP"].astype(str)
    df["hatoum_A1"] = df["A1"].astype(str).str.upper()
    df["hatoum_A2"] = df["A2"].astype(str).str.upper()

    df = df[
        df["hatoum_A1"].apply(is_snv) &
        df["hatoum_A2"].apply(is_snv)
    ].copy()

    df = df[
        ~df.apply(lambda r: is_ambiguous(r["hatoum_A1"], r["hatoum_A2"]), axis=1)
    ].copy()

    df = df.drop_duplicates(subset=["CHR", "BP", "hatoum_A1", "hatoum_A2"])

    return df

def harmonize_to_khat(khat, hatoum, khat_beta_col, khat_p_col):
    merged = khat.merge(hatoum, on=["CHR", "BP"], how="inner")

    print("  Raw lifted CHR:BP overlap before allele harmonization: {:,}".format(len(merged)))

    merged["hatoum_A1_comp"] = merged["hatoum_A1"].apply(complement_allele)
    merged["hatoum_A2_comp"] = merged["hatoum_A2"].apply(complement_allele)

    exact_match = (
        (merged["khat_A1"] == merged["hatoum_A1"]) &
        (merged["khat_A2"] == merged["hatoum_A2"])
    )

    swapped_match = (
        (merged["khat_A1"] == merged["hatoum_A2"]) &
        (merged["khat_A2"] == merged["hatoum_A1"])
    )

    complement_match = (
        (merged["khat_A1"] == merged["hatoum_A1_comp"]) &
        (merged["khat_A2"] == merged["hatoum_A2_comp"])
    )

    complement_swapped_match = (
        (merged["khat_A1"] == merged["hatoum_A2_comp"]) &
        (merged["khat_A2"] == merged["hatoum_A1_comp"])
    )

    keep = exact_match | swapped_match | complement_match | complement_swapped_match

    print("  Allele harmonization before dropping unmatched:")
    print("    exact: {:,}".format(exact_match.sum()))
    print("    swapped: {:,}".format(swapped_match.sum()))
    print("    complement: {:,}".format(complement_match.sum()))
    print("    complement_swapped: {:,}".format(complement_swapped_match.sum()))
    print("    unmatched: {:,}".format((~keep).sum()))

    merged = merged.loc[keep].copy()

    exact_match = exact_match.loc[merged.index]
    swapped_match = swapped_match.loc[merged.index]
    complement_match = complement_match.loc[merged.index]
    complement_swapped_match = complement_swapped_match.loc[merged.index]

    beta_col = "Beta"
    beta_harmonized_col = "beta_Hatoum2023AddictionEuropean_to_khat_A1"

    merged["beta_Hatoum2023AddictionEuropean"] = merged[beta_col]
    merged[beta_harmonized_col] = merged[beta_col]

    flip_mask = swapped_match | complement_swapped_match

    merged.loc[flip_mask, beta_harmonized_col] = (
        -1 * merged.loc[flip_mask, beta_col]
    )

    merged["allele_match_status"] = "exact"
    merged.loc[swapped_match, "allele_match_status"] = "swapped"
    merged.loc[complement_match, "allele_match_status"] = "complement"
    merged.loc[complement_swapped_match, "allele_match_status"] = "complement_swapped"

    keep_cols = [
        "khat_marker",
        "CHR",
        "BP",
        "khat_A1",
        "khat_A2",
        khat_beta_col,
        khat_p_col,
        "hatoum_marker",
        "hatoum_A1",
        "hatoum_A2",
        "beta_Hatoum2023AddictionEuropean",
        beta_harmonized_col,
        "P",
        "allele_match_status"
    ]

    return merged[keep_cols]

ku = prep_khat(ku_file, "beta_ku", "p_ku")
kuf = prep_khat(kuf_file, "beta_kuf", "p_kuf")

hatoum = prep_hatoum(hatoum_file)

print("Hatoum European variants after filtering: {:,}".format(len(hatoum)))

print("Processing KU")
ku_match = harmonize_to_khat(
    ku,
    hatoum,
    "beta_ku",
    "p_ku"
)

print("Processing KUF")
kuf_match = harmonize_to_khat(
    kuf,
    hatoum,
    "beta_kuf",
    "p_kuf"
)

ku_out = out_dir / "ku_Hatoum2023AddictionEuropean_harmonized.txt"
kuf_out = out_dir / "kuf_Hatoum2023AddictionEuropean_harmonized.txt"

ku_match.to_csv(ku_out, sep="\t", index=False)
kuf_match.to_csv(kuf_out, sep="\t", index=False)

print("KU harmonized matches: {:,}".format(len(ku_match)))
print("KU flipped: {:,}".format(
    (
        (ku_match["allele_match_status"] == "swapped") |
        (ku_match["allele_match_status"] == "complement_swapped")
    ).sum()
))

print("KUF harmonized matches: {:,}".format(len(kuf_match)))
print("KUF flipped: {:,}".format(
    (
        (kuf_match["allele_match_status"] == "swapped") |
        (kuf_match["allele_match_status"] == "complement_swapped")
    ).sum()
))

print("Wrote:", ku_out)
print("Wrote:", kuf_out)
print("Done.")
EOF

python3 /tmp/harmonize_hatoum_european_lifted.py
