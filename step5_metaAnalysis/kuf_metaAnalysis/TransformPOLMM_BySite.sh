#!/bin/bash
set -euo pipefail

# Transform raw POLMM association files and annotate alleles from matching .bim files.
# No PLINK command, .bed file, or .fam file is used.
#
# POLMM orientation:
#   switch.allele = TRUE  -> beta/Stat are for A1, so negate them to return to A2.
#   switch.allele = FALSE -> beta/Stat are already for A2, so leave them unchanged.
#
# The output beta and Stat are therefore oriented to the annotated A2 allele.

declare -A BIM_FILES=(
  [AAU]="AAU_passed_all_qc.bim"
  [KEMRI]="KEMRI_passed_all_qc.bim"
  [Moi]="Moi_passed_all_qc.bim"
  [Uganda]="Uganda_passed_all_qc.bim"
)

declare -A COHORT_NS=(
  [AAU]=10127
  [KEMRI]=2587
  [Moi]=3724
  [Uganda]=8558
)

for site in AAU KEMRI Moi Uganda; do
  input_file="assoc_${site}_studysite_covar_all"
  bim_file="${BIM_FILES[$site]}"
  output_file="A2_oriented_annotated_${input_file}.tsv"

  if [[ ! -f "${input_file}" ]]; then
    echo "[${site}] Skipping missing file: ${input_file}" >&2
    continue
  fi

  if [[ ! -f "${bim_file}" ]]; then
    echo "[${site}] Skipping missing allele file: ${bim_file}" >&2
    continue
  fi

  echo "[${site}] Transforming ${input_file} with ${bim_file}; N_SITE=${COHORT_NS[$site]}"

  awk -v OFS='\t' \
      -v BIM="${bim_file}" \
      -v N_SITE="${COHORT_NS[$site]}" '
    BEGIN {
      # PLINK .bim columns: CHR, SNPID, CM, BP, A1, A2.
      while ((getline < BIM) > 0) {
        snp = tolower($2)
        a1[snp] = $5
        a2[snp] = $6
      }
      close(BIM)

      print "SNPID", "A1", "A2", "chr", "MAF", "missing.rate", "Stat", "VarW", "VarP", \
            "beta", "pval.norm", "pval.spa", "original.switch.allele", "SE", "N_SITE"
    }

    NR == 1 { next }

    {
      gsub(/"/, "", $0)

      snp = tolower($1)
      if (!(snp in a1)) {
        unmatched++
        next
      }

      switch_value = toupper($11)
      stat_a2 = $5 + 0
      beta_a2 = $8 + 0

      if (switch_value == "TRUE" || switch_value == "T" || switch_value == "1") {
        stat_a2 = -stat_a2
        beta_a2 = -beta_a2
        flipped++
      } else if (switch_value == "FALSE" || switch_value == "F" || switch_value == "0") {
        unchanged++
      } else {
        unexpected_switch++
        next
      }

      if (($7 + 0) <= 0) {
        invalid_varp++
        next
      }

      # POLMM beta = Stat/VarP, so its model-based SE is 1/sqrt(VarP).
      se = 1 / sqrt($7 + 0)

      print snp, a1[snp], a2[snp], $2, $3, $4, stat_a2, $6, $7, beta_a2, \
            $9, $10, switch_value, se, N_SITE
      written++
    }

    END {
      printf "[%s] Written: %d; TRUE rows flipped: %d; FALSE rows unchanged: %d\n", \
             FILENAME, written + 0, flipped + 0, unchanged + 0 > "/dev/stderr"

      if (unexpected_switch > 0)
        printf "[%s] Skipped unexpected switch values: %d\n", \
               FILENAME, unexpected_switch > "/dev/stderr"

      if (invalid_varp > 0)
        printf "[%s] Skipped non-positive VarP values: %d\n", \
               FILENAME, invalid_varp > "/dev/stderr"

      if (unmatched > 0)
        printf "[%s] Skipped SNPIDs absent from .bim: %d\n", \
               FILENAME, unmatched > "/dev/stderr"
    }
  ' "${input_file}" > "${output_file}"

  echo "[${site}] Done -> ${output_file}"
done

echo "All available raw POLMM files processed."
