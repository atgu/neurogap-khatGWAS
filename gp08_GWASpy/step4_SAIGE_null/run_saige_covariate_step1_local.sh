# Made by Nico Matthew S. Valencia

#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="$(pwd)"
OUT_DIR="$DATA_DIR/saige_step1" 

# Docker image & resources
SAIGE_IMAGE="wzhou88/saige:1.3.0"
THREADS=16

# Phenotype settings
PHENO="assist_khat"          # Only run binary phenotype
TRAIT_TYPE="binary"
INV_NORMALIZE="FALSE"        # No inverse-normalization for binary

# Phenotype files for each site
declare -A SITE_TO_PHENO=(
  ["AAU"]="w_site_AAU_pheno_noNA.txt"
  ["Uganda"]="w_site_Uganda_pheno_noNA.txt"
  ["UCT"]="w_site_UCT_pheno.txt"
  ["KEMRI"]="w_site_KEMRI_pheno.txt"
  ["Moi"]="w_site_Moi_pheno.txt"
)

mkdir -p "$OUT_DIR"

for SITE in "${!SITE_TO_PHENO[@]}"; do
  PHENO_FILE="${DATA_DIR}/${SITE_TO_PHENO[$SITE]}"
  PLINK_PREFIX="${DATA_DIR}/gp08_${SITE}_ldpruned"
  OUTPUT_PREFIX="${OUT_DIR}/${SITE}_${PHENO}_saige_step1"

  if [[ ! -f "${PLINK_PREFIX}.bed" ]]; then
    echo "WARNING: Skipping $SITE — missing PLINK files with prefix: $PLINK_PREFIX"
    continue
  fi
  if [[ ! -f "$PHENO_FILE" ]]; then
    echo "WARNING: Skipping $SITE — missing phenotype file: $PHENO_FILE"
    continue
  fi

  # Include study_site only for AAU and Uganda
  if [[ "$SITE" == "AAU" || "$SITE" == "Uganda" ]]; then
    COVARS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age,sex,study_site"
  else
    COVARS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age,sex"
  fi


  docker run --rm \
    -v "$DATA_DIR":"$DATA_DIR" \
    "$SAIGE_IMAGE" \
    step1_fitNULLGLMM.R \
      --plinkFile="$PLINK_PREFIX" \
      --phenoFile="$PHENO_FILE" \
      --phenoCol="$PHENO" \
      --covarColList="$COVARS" \
      --sampleIDColinphenoFile="IID" \
      --traitType="$TRAIT_TYPE" \
      --outputPrefix="$OUTPUT_PREFIX" \
      --nThreads="$THREADS" \
      --LOCO=TRUE \
      --invNormalize="$INV_NORMALIZE"

  echo "Done: $SITE"
  echo
done

echo "All requested SAIGE step1 runs completed."
