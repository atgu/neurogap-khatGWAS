#!/bin/bash
set -euo pipefail

for input in \
  ku_gp08_meta_cleaned_with_N.txt \
  kuf_gp08_meta_cleaned_with_N_and_freq.tsv
do
  output="${input%.*}_SNP_A1_A2_Z.tsv"

  awk '
    BEGIN { OFS = "\t" }

    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "MarkerName") snp_col = i
        else if ($i == "Allele1") a1_col = i
        else if ($i == "Allele2") a2_col = i
        else if ($i == "Effect") effect_col = i
        else if ($i == "StdErr") se_col = i
      }

      if (!snp_col || !a1_col || !a2_col || !effect_col || !se_col) {
        print "ERROR: Required columns were not found in " FILENAME > "/dev/stderr"
        exit 1
      }

      print "SNP", "A1", "A2", "Z"
      next
    }

    {
      if ($effect_col == "" || $effect_col == "NA" ||
          $se_col == "" || $se_col == "NA" ||
          ($se_col + 0) == 0) {
        z = "NA"
      } else {
        z = ($effect_col + 0) / ($se_col + 0)
      }

      print $snp_col, toupper($a1_col), toupper($a2_col), z
    }
  ' "$input" > "$output"

  echo "Created: $output"
done
