#!/bin/bash
set -euo pipefail

cat > /tmp/harmonize_walters_beta.py << 'EOF'
import pandas as pd
from pathlib import Path

walters_file = Path("~/sud/walters-2018/liftover_only/pgc_alcdep.trans_fe_unrel_geno.aug2018_release.txt.lifted_hg38.tsv").expanduser()
khat_dir = Path("~/sud/khat").expanduser()
out_dir = Path("~/rsid_match/walters_lifted_harmonized").expanduser()
out_dir.mkdir(parents=True, exist_ok=True)

ku_file = khat_dir / "ku_gp08_meta_cleaned_with_N.txt"
kuf_file = khat_dir / "kuf_se_meta.txt"

trait = "pgc_alcdep_trans_fe_unrel_geno"

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

def prep_walters(path):
    df = pd.read_csv(path, sep="\t", dtype={"SNP": str}, low_memory=False)

    df["CHR"] = df["CHR_hg38"].astype(str)
    df["BP"] = df["BP_hg38"].astype(int)

    df["walters_marker"] = df["SNP"].astype(str)
    df["walters_A1"] = df["A1"].astype(str).str.upper()
    df["walters_A2"] = df["A2"].astype(str).str.upper()

    df = df[
        df["walters_A1"].apply(is_snv) &
        df["walters_A2"].apply(is_snv)
    ].copy()

    df = df[
        ~df.apply(lambda r: is_ambiguous(r["walters_A1"], r["walters_A2"]), axis=1)
    ].copy()

    df = df.drop_duplicates(subset=["CHR", "BP", "walters_A1", "walters_A2"])
    return df

def harmonize_to_khat(khat, walters, khat_beta_col, khat_p_col):
    merged = khat.merge(walters, on=["CHR", "BP"], how="inner")

    print("  Raw lifted CHR:BP overlap before allele harmonization: {:,}".format(len(merged)))

    merged["walters_A1_comp"] = merged["walters_A1"].apply(complement_allele)
    merged["walters_A2_comp"] = merged["walters_A2"].apply(complement_allele)

    exact_match = (
        (merged["khat_A1"] == merged["walters_A1"]) &
        (merged["khat_A2"] == merged["walters_A2"])
    )

    swapped_match = (
        (merged["khat_A1"] == merged["walters_A2"]) &
        (merged["khat_A2"] == merged["walters_A1"])
    )

    complement_match = (
        (merged["khat_A1"] == merged["walters_A1_comp"]) &
        (merged["khat_A2"] == merged["walters_A2_comp"])
    )

    complement_swapped_match = (
        (merged["khat_A1"] == merged["walters_A2_comp"]) &
        (merged["khat_A2"] == merged["walters_A1_comp"])
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
        "walters_marker", "walters_A1", "walters_A2",
        beta_raw_col, beta_harmonized_col,
        "SE", "P", "Neff", "INFO", "allele_match_status"
    ]

    return merged[keep_cols]

ku = prep_khat(ku_file, "beta_ku", "p_ku")
kuf = prep_khat(kuf_file, "beta_kuf", "p_kuf")
walters = prep_walters(walters_file)

print("Walters BETA variants after filtering: {:,}".format(len(walters)))

ku_match = harmonize_to_khat(ku, walters, "beta_ku", "p_ku")
kuf_match = harmonize_to_khat(kuf, walters, "beta_kuf", "p_kuf")

ku_out = out_dir / f"ku_{trait}_harmonized.txt"
kuf_out = out_dir / f"kuf_{trait}_harmonized.txt"

ku_match.to_csv(ku_out, sep="\t", index=False)
kuf_match.to_csv(kuf_out, sep="\t", index=False)

print("KU harmonized matches: {:,}".format(len(ku_match)))
print("KUF harmonized matches: {:,}".format(len(kuf_match)))
print("Wrote:", ku_out)
print("Wrote:", kuf_out)
print("Done.")
EOF

python3 /tmp/harmonize_walters_beta.py
