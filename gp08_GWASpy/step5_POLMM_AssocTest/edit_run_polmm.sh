#!/bin/bash
set -euo pipefail

# --- CONFIG ---
PHENO="assist_khat_amt"
SITE="KEMRI"   # AAU, KEMRI, Moi, UCT, Uganda

# Host paths (on your machine)
WORKDIR="$HOME/polmm"
HOST_GP_DIR="$HOME/gp08_data"

# Paths visible inside the container
PLINK_PREFIX="/gp08_data/adjusted/NeuroGAP_assoc_test_file"   # no ~ inside container
NULL_MODEL="w_site_${SITE}_pheno_polmm_null_dense.rds"

# Output (relative to WORKDIR)
OUTPUT_DIR="output/assoc_${SITE}_adjusted_filters_test3"

# Optional: allow passing a single chromosome from the command line (e.g., ./run_polmm.sh 2)
CHR_ARG="${1:-}"

# --- PREP ---
mkdir -p "${WORKDIR}/${OUTPUT_DIR}"
cd "${WORKDIR}"

# Sanity checks on host
[[ -f "${WORKDIR}/${NULL_MODEL}" ]] || { echo "Missing null model: ${WORKDIR}/${NULL_MODEL}"; exit 1; }
[[ -f "${HOST_GP_DIR}/adjusted/NeuroGAP_assoc_test_file.bim" ]] || { echo "Missing PLINK .bim in ${HOST_GP_DIR}/adjusted"; exit 1; }
[[ -f "${HOST_GP_DIR}/adjusted/NeuroGAP_assoc_test_file.bed" ]] || { echo "Missing PLINK .bed in ${HOST_GP_DIR}/adjusted"; exit 1; }
[[ -f "${HOST_GP_DIR}/adjusted/NeuroGAP_assoc_test_file.fam" ]] || { echo "Missing PLINK .fam in ${HOST_GP_DIR}/adjusted"; exit 1; }

echo "Starting POLMM association test for ${SITE}..."
echo "Null model: ${NULL_MODEL}"
echo "PLINK prefix (container): ${PLINK_PREFIX}"
echo "Output dir: ${OUTPUT_DIR}"
if [[ -n "${CHR_ARG}" ]]; then
  echo "Restricting to chromosome: ${CHR_ARG}"
fi

# --- RUN (non-interactive: no -t/-i flags) ---
if [[ -n "${CHR_ARG}" ]]; then
  # wrapper supports optional 4th arg = chromosome
  sudo docker run --rm \
    -v "${WORKDIR}:/data" \
    -v "${HOST_GP_DIR}:/gp08_data" \
    -w /data \
    nsvalencia/polmm:v6 \
    Rscript run_polmm_assoc_with_loco.R "${NULL_MODEL}" "${PLINK_PREFIX}" "${OUTPUT_DIR}" "${CHR_ARG}"
else
  sudo docker run --rm \
    -v "${WORKDIR}:/data" \
    -v "${HOST_GP_DIR}:/gp08_data" \
    -w /data \
    nsvalencia/polmm:v6 \
    Rscript run_polmm_assoc_with_loco.R "${NULL_MODEL}" "${PLINK_PREFIX}" "${OUTPUT_DIR}"
fi

echo "Done. Output in ${OUTPUT_DIR}"
