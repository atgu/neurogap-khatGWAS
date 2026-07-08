#!/usr/bin/env bash
set -euo pipefail

# If the var exists at a respective site, we will add that site's respective overall sample size to that var's N
BIM_DIR="$HOME/neurogap_passed_all_qc"

# -------------------------
# 1) ku_se_meta.txt
# Already has: N_case N_ctrl N
# Standardize name/format only
# -------------------------
cp ku_se_meta.txt ku_gp08_meta_cleaned_with_N.txt


# -------------------------
# 2) kuf_se_meta.txt
# Add N based on whether MarkerName exists in each site's .bim
# Site Ns:
# KEMRI=2587
# AAU=10127
# Moi=3724
# Uganda=8558
# -------------------------

awk '
BEGIN { OFS="\t" }

FNR==NR {
    if (FILENAME ~ /KEMRI/) site["KEMRI"][$2]=1
    else if (FILENAME ~ /AAU/) site["AAU"][$2]=1
    else if (FILENAME ~ /Moi/) site["Moi"][$2]=1
    else if (FILENAME ~ /Uganda/) site["Uganda"][$2]=1
    next
}

FNR==1 {
    print $0, "N"
    next
}

{
    marker=$1
    n=0

    if (marker in site["KEMRI"])  n += 2587
    if (marker in site["AAU"])    n += 10127
    if (marker in site["Moi"])    n += 3724
    if (marker in site["Uganda"]) n += 8558

    print $0, n
}
' \
"$BIM_DIR/KEMRI_passed_all_qc.bim" \
"$BIM_DIR/AAU_passed_all_qc.bim" \
"$BIM_DIR/Moi_passed_all_qc.bim" \
"$BIM_DIR/Uganda_passed_all_qc.bim" \
kuf_se_meta.txt > kuf_gp08_meta_cleaned_with_N.tsv
