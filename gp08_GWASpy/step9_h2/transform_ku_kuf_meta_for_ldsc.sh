#!/usr/bin/env bash
set -euo pipefail

BIM_DIR="$HOME/neurogap_passed_all_qc"

# -------------------------
# 1) ku_se_meta.txt
# Add/recalculate N as effective sample size:
# N_eff = 4 / ((1/N_case) + (1/N_ctrl)) --> N_eff is N
# -------------------------

awk '
BEGIN { OFS="\t" }

NR==1 {
    print "MarkerName","CHR","BP","Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P-value","Direction","N_case","N_ctrl","N"
    next
}

{
    neff = 4 / ((1/$14) + (1/$15))

    for (i=1; i<=15; i++) {
        printf "%s%s", $i, OFS
    }

    printf "%g\n", neff
}
' ku_se_meta.txt > ku_gp08_meta_cleaned_with_N.txt


# -------------------------
# 2) kuf_se_meta.txt
# Add N based on whether MarkerName exists in each site's .bim
# -------------------------

awk '
BEGIN { OFS="\t" }

FNR==NR {
    site["KEMRI", $2] = 1
    next
}

ARGIND==2 {
    site["AAU", $2] = 1
    next
}

ARGIND==3 {
    site["Moi", $2] = 1
    next
}

ARGIND==4 {
    site["Uganda", $2] = 1
    next
}

ARGIND==5 && FNR==1 {
    print $0, "N"
    next
}

ARGIND==5 {
    marker=$1
    n=0

    if (("KEMRI", marker) in site)  n += 2587
    if (("AAU", marker) in site)    n += 10127
    if (("Moi", marker) in site)    n += 3724
    if (("Uganda", marker) in site) n += 8558

    print $0, n
}
' \
"$BIM_DIR/KEMRI_passed_all_qc.bim" \
"$BIM_DIR/AAU_passed_all_qc.bim" \
"$BIM_DIR/Moi_passed_all_qc.bim" \
"$BIM_DIR/Uganda_passed_all_qc.bim" \
kuf_se_meta.txt > kuf_gp08_meta_cleaned_with_N.tsv

echo "Done."
echo "Wrote:"
echo "  ku_gp08_meta_cleaned_with_N.txt"
echo "  kuf_gp08_meta_cleaned_with_N.tsv"
