#!/usr/bin/env bash
set -euo pipefail

# Python versions get a bit weird with ldsc. Starting off, use this one: https://github.com/cbiit/ldsc
# this makes a conda environment called ldsc39 that we can initially use, we later transition to a different ldsc environment
# conda activate ldsc39

# NOTE: Because of negative h2 estimates, we decided to not go forward with calculating rg using LDSC or popcorn

# ----------------------------
# Paths
# ----------------------------

SUMSTATS_DIR="${HOME}/sud/khat"
REF_DIR="${HOME}/hgdp1kg_ref"

KU_SUMSTATS="${SUMSTATS_DIR}/ku_gp08_meta_cleaned_with_N.txt"
KUF_SUMSTATS="${SUMSTATS_DIR}/kuf_gp08_meta_cleaned_with_N.tsv"

OUT_DIR="${SUMSTATS_DIR}/ldsc_h2_gp08"
mkdir -p "${OUT_DIR}"

# Change this if your LDSC scripts are elsewhere
LDSC_DIR="${HOME}/ldsc"

MUNGE="${LDSC_DIR}/munge_sumstats.py"
LDSC="${LDSC_DIR}/ldsc.py"

LD_PREFIX="${REF_DIR}/HGDP_1KG_AFR_hapmap"

# Use the same panel as regression weights unless you have separate weights
W_PREFIX="${LD_PREFIX}"

# ----------------------------
# Check inputs
# ----------------------------

echo "Checking files..."

for f in \
    "${KU_SUMSTATS}" \
    "${KUF_SUMSTATS}" \
    "${LD_PREFIX}.l2.ldscore.gz" \
    "${LD_PREFIX}.l2.M" \
    "${LD_PREFIX}.l2.M_5_50"
do
    if [[ ! -f "$f" ]]; then
        echo "Missing required file: $f"
        exit 1
    fi
done

# ----------------------------
# Munge summary statistics
# ----------------------------

echo "Munging KU summary stats..."

python "${MUNGE}" \
    --sumstats "${KU_SUMSTATS}" \
    --snp MarkerName \
    --a1 Allele1 \
    --a2 Allele2 \
    --p "P-value" \
    --signed-sumstats Effect,0 \
    --N-col N \
    --out "${OUT_DIR}/ku_gp08"

echo "Munging KUF summary stats..."

python "${MUNGE}" \
    --sumstats "${KUF_SUMSTATS}" \
    --snp MarkerName \
    --a1 Allele1 \
    --a2 Allele2 \
    --p "P-value" \
    --signed-sumstats Effect,0 \
    --N-col N \
    --out "${OUT_DIR}/kuf_gp08"

# ----------------------------
# Estimate SNP heritability
# ----------------------------

echo "Running LDSC h2 for KU..."

python "${LDSC}" \
    --h2 "${OUT_DIR}/ku_gp08.sumstats.gz" \
    --ref-ld "${LD_PREFIX}" \
    --w-ld "${W_PREFIX}" \
    --out "${OUT_DIR}/ku_gp08_h2"

echo "Running LDSC h2 for KUF..."

python "${LDSC}" \
    --h2 "${OUT_DIR}/kuf_gp08.sumstats.gz" \
    --ref-ld "${LD_PREFIX}" \
    --w-ld "${W_PREFIX}" \
    --out "${OUT_DIR}/kuf_gp08_h2"

echo "Done."
echo "Results:"
echo "${OUT_DIR}/ku_gp08_h2.log"
echo "${OUT_DIR}/kuf_gp08_h2.log"
