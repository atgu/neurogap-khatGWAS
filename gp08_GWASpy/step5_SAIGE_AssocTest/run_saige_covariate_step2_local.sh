#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$HOME/gp08_data/gwaspy_filtered/pre_ldprune/post_ldprune"
STEP1_DIR="$BASE_DIR/saige_step1"
STEP2_DIR="$STEP1_DIR/saige_step2"
mkdir -p "$STEP2_DIR"

SAIGE_IMAGE="wzhou88/saige:1.3.0"

phenotypes=("assist_khat") 
# sites=("AAU" "KEMRI" "Moi" "Uganda" "UCT")
sites=("Moi")

# Set this to TRUE only if step1 was run with LOCO
LOCO_FLAG=TRUE

for phenotype in "${phenotypes[@]}"; do
  for site in "${sites[@]}"; do
    echo "Running SAIGE step2 for SITE=$site, PHENO=$phenotype"

    rda_file="${STEP1_DIR}/${site}_${phenotype}_saige_step1.rda"
    variance_file="${STEP1_DIR}/${site}_${phenotype}_saige_step1.varianceRatio.txt"

    if [[ ! -f "$rda_file" || ! -f "$variance_file" ]]; then
      echo "WARNING: Missing step1 outputs for $site/$phenotype. Skipping."
      echo "  Expected: $rda_file"
      echo "            $variance_file"
      continue
    fi

    # Use $HOME so the path expands; no quotes around ~
    PLINK_PREFIX="$HOME/gp08_data/gwaspy_filtered/${site}_passed_all_qc"
    if [[ ! -f "${PLINK_PREFIX}.bed" || ! -f "${PLINK_PREFIX}.bim" || ! -f "${PLINK_PREFIX}.fam" ]]; then
      echo "WARNING: Missing PLINK trio for $site at prefix: $PLINK_PREFIX. Skipping."
      continue
    fi

    for chr in {1..22}; do
      output_file="${STEP2_DIR}/test/${site}_chr${chr}_${phenotype}_saige_step2.txt"

      docker run --rm \
        -v "$BASE_DIR":"$BASE_DIR" \
        -v "$HOME/gp08_data":"$HOME/gp08_data" \
        "$SAIGE_IMAGE" \
        step2_SPAtests.R \
          --bedFile="${PLINK_PREFIX}.bed" \
          --bimFile="${PLINK_PREFIX}.bim" \
          --famFile="${PLINK_PREFIX}.fam" \
          --AlleleOrder=alt-first \
          --GMMATmodelFile="$rda_file" \
          --varianceRatioFile="$variance_file" \
          --SAIGEOutputFile="$output_file" \
          --LOCO=$LOCO_FLAG \
          --chrom "$chr" \
          --is_output_moreDetails=TRUE \
          --is_Firth=TRUE \
          --pCutoffforFirth=0.05
         #--minMAF=0.01

      echo "Wrote: $output_file"
    done

    echo "Done SITE=$site PHENO=$phenotype"
    echo
  done
done

echo "All SAIGE step2 jobs completed."
