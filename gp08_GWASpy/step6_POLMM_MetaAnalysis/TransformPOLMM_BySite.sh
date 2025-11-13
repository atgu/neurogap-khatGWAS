#!/bin/bash
set -euo pipefail

PLINK2="/home/tboltz/apps/plink2"

# Map site -> plink bfile prefix (exact filenames in your folder)
declare -A BFILES=(
  [AAU]="AAU_passed_all_qc"
  [KEMRI]="KEMRI_passed_all_qc"
  [Moi]="Moi_passed_all_qc"
  [Uganda]="Uganda_passed_all_qc"
  [UCT]="UCT_passed_all_qc"
)

# hard-code cohort Ns (from assoc logs)
declare -A COHORT_NS=(
  [KEMRI]=2587
  [AAU]=10127
  [Moi]=3724
  [Uganda]=8558
  [UCT]=6963
)


# Process each site
for site in AAU KEMRI Moi Uganda UCT; do
  bfile="${BFILES[$site]}"
  assoc_file="assoc_${site}_Nas0_MAF1e-2"

  if [[ ! -f "${bfile}.bed" || ! -f "${bfile}.bim" || ! -f "${bfile}.fam" ]]; then
    echo "[${site}] Skipping: missing PLINK files for ${bfile} (.bed/.bim/.fam)."
    continue
  fi

  if [[ ! -f "${assoc_file}" ]]; then
    echo "[${site}] Skipping: missing ${assoc_file}."
    continue
  fi

  echo "[${site}] Step 0: per-site allele freq for N-per-SNP"
  ${PLINK2} --bfile "${bfile}" --freq --out "snp_n_sizes_${site}"

  echo "[${site}] Step 1: SNP -> (CHR BP A1 A2) from site .bim"
  awk 'BEGIN{OFS="\t"} {print tolower($2), $1, $4, $5, $6}' "${bfile}.bim" \
    | sort -k1,1 > "snp_info_lookup_${site}.sorted.txt"

  echo "[${site}] Step 2: SNP -> N from site .afreq"
  # Using your original $7 column; keep consistent with your plink2 version/output
  awk 'NR>1 {print tolower($2), $7}' "snp_n_sizes_${site}.afreq" \
    | sort -k1,1 > "snp_n_lookup_${site}.sorted.txt"

  echo "[${site}] Step 3: resolve N_SITE"
  N_SITE="${COHORT_NS[$site]}"
  echo "[${site}] N_SITE=${N_SITE}"

  echo "[${site}] Step 4: join lookups with ${assoc_file}"
  tmp="temp_${assoc_file}_${site}.tsv"
  {
    echo -e "SNPID\tA1\tA2\tFREQ\tmissing.rate\tStat\tVarW\tVarP\tbeta\tpval.norm\tpval.spa\tswitch.allele\tN_SNP\tN_SITE"
    awk -v N_SITE="${N_SITE}" -v INFO="snp_info_lookup_${site}.sorted.txt" -v NLK="snp_n_lookup_${site}.sorted.txt" 'BEGIN{
        OFS="\t"
        while((getline<INFO)>0){a1[$1]=$4; a2[$1]=$5}
        while((getline<NLK)>0){n[$1]=$2}
      }
      NR>1{
        gsub(/\"/,"",$0)
        snp=tolower($1)
        # require N_SNP present and (keep your original field guards)
        if (snp in n && $7 != ""){
          print snp, a1[snp], a2[snp], $3, $4, $5, $6, $7, $8, $9, $10, $14, n[snp], N_SITE
        }
      }' "${assoc_file}"
  } > "${tmp}"

  echo "[${site}] Step 5: compute SE from beta and pval.spa"
  python3 - <<EOF
import pandas as pd
from scipy.stats import norm

# Read in the intermediate file
df = pd.read_csv("${tmp}", sep="\t")

# Compute standard error from beta and pval.spa
df["SE"] = abs(df["beta"] / norm.ppf(1 - df["pval.spa"] / 2))

# Reorder and export columns
outcols = [
    "SNPID","A1","A2","FREQ","missing.rate","Stat","VarW","VarP",
    "beta","pval.norm","pval.spa","switch.allele","SE","N_SNP","N_SITE"
]
df[outcols].to_csv(f"cleaned_${site}.tsv", sep="\t", index=False)
EOF

  rm -f "${tmp}"
  echo "[${site}] Done -> cleaned_${site}.tsv"
done

echo "All sites processed."
