#!/bin/bash

# Exit on error
set -e

# Variables
PHENO="assist_khat_amt"
SITE="AAU" # AAU, KEMRI, Moi, UCT, Uganda
PLINK_PREFIX="${SITE}_passed_all_qc"
NULL_MODEL="${SITE}_kufNas0_polmmNull.rds"
OUTPUT_DIR="output/Nas0_MAF1e-2"
WORKDIR="$HOME/polmm"

mkdir -p "${WORKDIR}/${OUTPUT_DIR}"
cd "${WORKDIR}"

echo "Using PLINK files and null model in the working directory..."
echo "Starting POLMM association test for ${SITE}..."

# Run the Docker container using the fixed image
sudo docker run --rm \
    -v "${WORKDIR}:/data" \
    -w /data \
    nsvalencia/polmm:v6 \
    Rscript run_polmm_assoc.R ${NULL_MODEL} ${PLINK_PREFIX} "${OUTPUT_DIR}/assoc_${SITE}_Nas0_MAF1e-2"

echo "POLMM association test completed for ${SITE}. Output: ${OUTPUT_DIR}/assoc_${SITE}_Nas0_MAF1e-2"
