#!/bin/bash
set -euo pipefail

# Requirements
#
# Python >= 3.8
#
# Install required Python package:
#
#   pip install pandas
#
# or:
#
#   conda install pandas
#
# Standard library modules used:
#   pathlib
#
# If you have an environment for LDSC already, you can just use that instead 

mkdir -p ~/rsid_match/walters

cat > /tmp/match_walters_rsids.py << 'EOF'
import pandas as pd
from pathlib import Path

walters_dir = Path("~/sud/walters-2018").expanduser()
khat_dir = Path("~/sud/khat").expanduser()
out_dir = Path("~/rsid_match/walters").expanduser()

out_dir.mkdir(parents=True, exist_ok=True)

ku_file = khat_dir / "ku_gp08_meta_cleaned_with_N.txt"
kuf_file = khat_dir / "kuf_se_meta.txt"


def is_snv(a):
    return isinstance(a, str) and len(a) == 1 and a.upper() in {"A", "C", "G", "T"}


def is_ambiguous(a1, a2):
    pair = {str(a1).upper(), str(a2).upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def complement_allele(a):
    comp = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    return comp.get(str(a).upper(), None)


def prep_khat(path, beta_col, p_col):
    df = pd.read_csv(path, sep="\t", dtype={"MarkerName": str})
    df = df.drop_duplicates(subset="MarkerName")

    df = df[["MarkerName", "Allele1", "Allele2", "Effect", "P-value"]].rename(
        columns={
            "MarkerName": "RSID",
            "Allele1": "khat_A1",
            "Allele2": "khat_A2",
            "Effect": beta_col,
            "P-value": p_col
        }
    )

    df["khat_A1"] = df["khat_A1"].str.upper()
    df["khat_A2"] = df["khat_A2"].str.upper()

    df = df[
        df["khat_A1"].apply(is_snv) &
        df["khat_A2"].apply(is_snv)
    ]

    df = df[
        ~df.apply(
            lambda r: is_ambiguous(r["khat_A1"], r["khat_A2"]),
            axis=1
        )
    ]

    return df


def harmonize_to_khat(khat, walters, trait, khat_beta_col, khat_p_col):

    merged = khat.merge(walters, on="RSID", how="inner")

    print("  Raw RSID overlap before allele harmonization: {:,}".format(len(merged)))

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

    keep = (
        exact_match |
        swapped_match |
        complement_match |
        complement_swapped_match
    )

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

    beta_col = "beta_{}".format(trait)
    p_col = "p_{}".format(trait)
    beta_harmonized_col = "beta_{}_to_khat_A1".format(trait)

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
        "RSID",
        "khat_A1",
        "khat_A2",
        khat_beta_col,
        khat_p_col,
        "walters_A1",
        "walters_A2",
        beta_col,
        beta_harmonized_col,
        p_col,
        "allele_match_status"
    ]

    return merged[keep_cols]


ku = prep_khat(ku_file, "beta_ku", "p_ku")
kuf = prep_khat(kuf_file, "beta_kuf", "p_kuf")

walters_files = [
    (
        walters_dir / "pgc_alcdep.trans_fe_unrel_geno.aug2018_release.txt",
        "WaltersTransFE",
        "STAT1_FE"
    )
]

for f, trait, stat_col in walters_files:

    print("Processing {}".format(trait))

    walters = pd.read_csv(
        f,
        sep=r"\s+",
        dtype={"SNP": str}
    )

    walters.columns = walters.columns.str.strip()

    if stat_col not in walters.columns:
        raise ValueError(
            "Could not find {} in {}. Columns found: {}".format(
                stat_col,
                f.name,
                list(walters.columns)
            )
        )

    walters = walters[
        ["SNP", "A1", "A2", stat_col, "P"]
    ].rename(
        columns={
            "SNP": "RSID",
            "A1": "walters_A1",
            "A2": "walters_A2",
            stat_col: "beta_{}".format(trait),
            "P": "p_{}".format(trait)
        }
    )

    walters["walters_A1"] = walters["walters_A1"].str.upper()
    walters["walters_A2"] = walters["walters_A2"].str.upper()

    walters = walters.drop_duplicates(subset="RSID")

    walters = walters[
        walters["RSID"].notna() &
        (walters["RSID"] != ".") &
        (walters["RSID"] != "NA")
    ]

    walters = walters[
        walters["walters_A1"].apply(is_snv) &
        walters["walters_A2"].apply(is_snv)
    ]

    walters = walters[
        ~walters.apply(
            lambda r: is_ambiguous(r["walters_A1"], r["walters_A2"]),
            axis=1
        )
    ]

    print("  Walters rows after filtering: {:,}".format(len(walters)))

    ku_match = harmonize_to_khat(
        ku,
        walters,
        trait,
        "beta_ku",
        "p_ku"
    )

    kuf_match = harmonize_to_khat(
        kuf,
        walters,
        trait,
        "beta_kuf",
        "p_kuf"
    )

    ku_match.to_csv(
        out_dir / "ku_{}.txt".format(trait),
        sep="\t",
        index=False
    )

    kuf_match.to_csv(
        out_dir / "kuf_{}.txt".format(trait),
        sep="\t",
        index=False
    )

    print("  KU harmonized matches: {:,}".format(len(ku_match)))
    print("  KU flipped: {:,}".format(
        (
            (ku_match["allele_match_status"] == "swapped") |
            (ku_match["allele_match_status"] == "complement_swapped")
        ).sum()
    ))

    print("  KUF harmonized matches: {:,}".format(len(kuf_match)))
    print("  KUF flipped: {:,}".format(
        (
            (kuf_match["allele_match_status"] == "swapped") |
            (kuf_match["allele_match_status"] == "complement_swapped")
        ).sum()
    ))

print("Done.")
EOF

python3 /tmp/match_walters_rsids.py
