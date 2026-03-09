#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Inputs / outputs
# -----------------------------
# BIM file to use for allele lookup (SNP -> A1/A2)
BIM_FILE="all_sites_all_phenos.bim"

# Assoc outputs directory (set to where your assoc_* files live)
ASSOC_DIR="."

# -----------------------------
# Cohort sample sizes (NO UCT):
# Reported from Null Model outputs
# -----------------------------
declare -A cohort_sizes
cohort_sizes=(
  [KEMRI]=2587
  [AAU]=10142
  [Moi]=3724
  [Uganda]=9416
)

# -----------------------------
# Build allele lookup from BIM
# Output format: snpid_lower \t A1 \t A2
# -----------------------------
if [[ ! -s "${BIM_FILE}" ]]; then
  echo "ERROR: BIM file not found or empty: ${BIM_FILE}" >&2
  exit 1
fi

ALLELE_LOOKUP="snp_alleles_lookup.sorted.txt"
awk 'BEGIN{OFS="\t"} {print tolower($2), $5, $6}' "${BIM_FILE}" | sort -k1,1 > "${ALLELE_LOOKUP}"

# -----------------------------
# Process each POLMM assoc output
# -----------------------------
shopt -s nullglob
for file in "${ASSOC_DIR}"/assoc_*_Nas0_MAF1e-2*; do
  [[ -f "$file" ]] || continue

  base="$(basename "$file")"

  # Expect: assoc_<COHORT>_Nas0_MAF1e-2...
  cohort="${base#assoc_}"
  cohort="${cohort%%_Nas0_MAF1e-2*}"

  N_SITE="${cohort_sizes[$cohort]:-}"
  if [[ -z "${N_SITE}" ]]; then
    echo "Skipping '$base' (parsed cohort='$cohort' not in cohort_sizes: ${!cohort_sizes[*]})" >&2
    continue
  fi

  temp_file="temp_${base}"
  output_file="metal_ready_${base}"

  # Intermediate TSV: add A1/A2, carry switch.allele, append N_SITE
  {
    echo -e "SNPID\tA1\tA2\tFREQ\tmissing.rate\tStat\tVarW\tVarP\tbeta\tpval.norm\tpval.spa\tswitch.allele\tN_SITE"

    awk -v N_SITE="${N_SITE}" 'BEGIN {
      OFS="\t"
      while ((getline < "snp_alleles_lookup.sorted.txt") > 0) {
        a1[$1] = $2
        a2[$1] = $3
      }
    }
    NR > 1 {
      gsub(/\"/, "", $0)
      snp = tolower($1)

      # Keep only rows with non-empty VarP (as in your prior script: $7)
      if ($7 != "") {
        aa1 = (snp in a1 ? a1[snp] : "NA")
        aa2 = (snp in a2 ? a2[snp] : "NA")
        # $3 FREQ, $4 missing.rate, $5 Stat, $6 VarW, $7 VarP, $8 beta, $9 pval.norm, $10 pval.spa, $14 switch.allele
        print snp, aa1, aa2, $3, $4, $5, $6, $7, $8, $9, $10, $14, N_SITE
      }
    }' "$file"
  } > "$temp_file"

  # Compute SE and create METAL-ready EA/NEA based on switch.allele
  python3 - <<EOF
import numpy as np
import pandas as pd
from scipy.stats import norm

input_file = "$temp_file"
output_file = "$output_file"

df = pd.read_csv(input_file, sep="\\t")

# --- SE from pval.spa (your prior approach) ---
p = df["pval.spa"].astype(float)
p = p.clip(lower=np.nextafter(0, 1), upper=1 - np.finfo(float).eps)
z = norm.ppf(1 - p / 2)
df["SE"] = (df["beta"].astype(float).abs() / np.abs(z))

# --- POLMM.plink effect allele logic ---
# switch.allele can be TRUE/FALSE, 1/0, Yes/No, etc.
sw = df["switch.allele"].astype(str).str.strip().str.lower().isin(["true", "t", "1", "yes", "y"])

# If switch=FALSE: beta corresponds to A2 -> EA=A2, NEA=A1
df["EA"]  = df["A2"]
df["NEA"] = df["A1"]

# If switch=TRUE: beta corresponds to A1 -> EA=A1, NEA=A2
df.loc[sw, "EA"]  = df.loc[sw, "A1"]
df.loc[sw, "NEA"] = df.loc[sw, "A2"]

# Optional normalization
df["EA"] = df["EA"].astype(str).str.upper()
df["NEA"] = df["NEA"].astype(str).str.upper()

# Final column order for METAL (keep pval.spa as reference/QC)
df = df[[
    "SNPID",
    "EA",
    "NEA",
    "FREQ",
    "missing.rate",
    "beta",
    "SE",
    "pval.spa",
    "N_SITE"
]]

df.to_csv(output_file, sep="\\t", index=False)
EOF

  rm -f "$temp_file"
  echo "Wrote: $output_file (cohort=$cohort, N_SITE=$N_SITE)"
done
