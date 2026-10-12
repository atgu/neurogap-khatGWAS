#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
shopt -s nullglob

SAUNDERS_DIR="${HOME}/sud/saunders-2022"
KHAT_DIR="${HOME}/sud/khat/a2_effect"
OUT_DIR="${HOME}/rsid_match/saunders_lifted_harmonized/a2_effect"

KU_FILE="${KHAT_DIR}/ku_gp08_meta_cleaned_with_N.txt"
KUF_FILE="${KHAT_DIR}/kuf_gp08_meta_cleaned_with_N.tsv"

mkdir -p "${OUT_DIR}"

for file in "${KU_FILE}" "${KUF_FILE}"; do
  [[ -f "${file}" ]] || {
    echo "ERROR: required Khat file not found: ${file}" >&2
    exit 1
  }
done

# These patterns select only Saunders files with directly comparable BETA,
# EFFECT_ALLELE, OTHER_ALLELE, SE, and P columns. MULTI meta-regression,
# MDS, and PRS files are deliberately excluded.
SAUNDERS_FILES=(
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_AFR.txt
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_AFR.txt.gz
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_AMR.txt
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_AMR.txt.gz
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_EAS.txt
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_EAS.txt.gz
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_EUR.txt
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_EUR.txt.gz
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_EUR_LOO_UKB_23andMe.txt
  "${SAUNDERS_DIR}"/GSCAN_*_2022_GWAS_SUMMARY_STATS_EUR_LOO_UKB_23andMe.txt.gz
)

if (( ${#SAUNDERS_FILES[@]} == 0 )); then
  echo "ERROR: no Saunders ancestry-stratified BETA files found in ${SAUNDERS_DIR}" >&2
  exit 1
fi

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/saunders_harmonize.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
echo "Temporary workspace: ${WORK_DIR}"

prepare_khat() {
  local input=$1 output=$2

  # Output: key marker CHR BP effect_allele other_allele beta pvalue
  awk 'BEGIN{FS=OFS="\t"}
    function clean(x){gsub(/\r|"/,"",x); return x}
    function normchr(x){x=clean(x); sub(/^chr/,"",x); sub(/\.0$/,"",x); return x}
    function snv(x){return x~/^[ACGT]$/}
    function amb(a,b){return a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC"}
    NR==1{
      for(i=1;i<=NF;i++){
        h=clean($i)
        if(h=="MarkerName")m=i; else if(h=="CHR")c=i; else if(h=="BP")pos=i
        else if(h=="Allele1")a1=i; else if(h=="Allele2")a2=i
        else if(h=="Effect")b=i; else if(h=="P-value")p=i
      }
      if(!m||!c||!pos||!a1||!a2||!b||!p){
        print "ERROR: Khat header missing required columns in " FILENAME > "/dev/stderr"; exit 2
      }
      next
    }
    {
      marker=clean($(m)); chr=normchr($(c)); bp=int($(pos))
      effect=toupper(clean($(a1))); other=toupper(clean($(a2)))
      beta=clean($(b)); pv=clean($(p))
      if(marker!=""&&chr!=""&&bp>0&&beta!=""&&pv!=""&&snv(effect)&&snv(other)&&!amb(effect,other))
        print chr ":" bp,marker,chr,bp,effect,other,beta,pv
    }
  ' "${input}" | sort -t $'\t' -k1,1 -k5,6 -u > "${output}"
}

stream_file() {
  local input=$1
  if [[ "${input}" == *.gz ]]; then
    gzip -cd -- "${input}"
  else
    command cat -- "${input}"
  fi
}

prepare_saunders() {
  local input=$1 output=$2

  # Output: key rsid effect_allele other_allele effect_AF beta SE P N
  stream_file "${input}" |
  awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"}
    function clean(x){gsub(/\r|"/,"",x); return x}
    function normchr(x){x=clean(x); sub(/^chr/,"",x); sub(/\.0$/,"",x); return x}
    function snv(x){return x~/^[ACGT]$/}
    function amb(a,b){return a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC"}
    NR==1{
      for(i=1;i<=NF;i++){
        h=toupper(clean($i))
        if(h=="CHR")c=i; else if(h=="POS")pos=i; else if(h=="RSID")rsid=i
        else if(h=="EFFECT_ALLELE")ea=i; else if(h=="OTHER_ALLELE")oa=i
        else if(h=="AF_1000G")af=i; else if(h=="BETA")beta=i
        else if(h=="SE")se=i; else if(h=="P")p=i; else if(h=="N")n=i
      }
      if(!c||!pos||!rsid||!ea||!oa||!af||!beta||!se||!p||!n){
        print "ERROR: Saunders header missing required columns" > "/dev/stderr"; exit 2
      }
      next
    }
    {
      chr=normchr($(c)); bp=int($(pos)); marker=clean($(rsid))
      effect=toupper(clean($(ea))); other=toupper(clean($(oa)))
      b=clean($(beta)); pv=clean($(p))
      if(chr!=""&&bp>0&&marker!=""&&marker!="."&&marker!="NA"&&b!=""&&pv!=""&&snv(effect)&&snv(other)&&!amb(effect,other))
        print chr ":" bp,marker,effect,other,clean($(af)),b,clean($(se)),pv,clean($(n))
    }
  ' | sort -t $'\t' -k1,1 -k3,4 -u > "${output}"
}

harmonize() {
  local khat_sorted=$1 saunders_sorted=$2 phenotype=$3 trait=$4 output
  output="${OUT_DIR}/${phenotype}_${trait}_harmonized.txt"

  printf 'khat_marker\tCHR\tBP\tkhat_A1\tkhat_A2\tbeta_%s\tp_%s\tsaunders_marker\tsaunders_effect_allele\tsaunders_other_allele\tsaunders_effect_allele_AF_1000G\tbeta_%s_raw\tbeta_%s_to_khat_A1\tsaunders_SE\tsaunders_P\tsaunders_N\tallele_match_status\n' \
    "${phenotype}" "${phenotype}" "${trait}" "${trait}" > "${output}"

  # Joined fields: 1 key; 2-8 Khat; 9-16 Saunders.
  join -t $'\t' -1 1 -2 1 "${khat_sorted}" "${saunders_sorted}" |
  awk 'BEGIN{FS=OFS="\t"}
    function comp(x){return x=="A"?"T":x=="T"?"A":x=="C"?"G":x=="G"?"C":""}
    {
      k1=$5; k2=$6; effect=$10; other=$11; status=""; flip=0
      # Saunders BETA is based on EFFECT_ALLELE. Align it to Khat Allele1.
      if(k1==effect&&k2==other) status="effect_match"
      else if(k1==other&&k2==effect){status="effect_swapped";flip=1}
      else if(k1==comp(effect)&&k2==comp(other)) status="effect_complement"
      else if(k1==comp(other)&&k2==comp(effect)){status="effect_complement_swapped";flip=1}
      else next

      raw=$13+0; aligned=(flip?-raw:raw)
      print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,raw,aligned,$14,$15,$16,status
    }
  ' >> "${output}"

  awk 'BEGIN{FS="\t"} NR>1{n++; if($17~/swapped$/)f++}
       END{printf "  Wrote %s: %d matches; %d effects flipped\n",FILENAME,n+0,f+0}' "${output}"
}

echo "Preparing KU"
prepare_khat "${KU_FILE}" "${WORK_DIR}/ku.sorted.tsv"
echo "Preparing KUF"
prepare_khat "${KUF_FILE}" "${WORK_DIR}/kuf.sorted.tsv"

# Avoid processing both compressed and decompressed copies of the same file.
declare -A SEEN
for input in "${SAUNDERS_FILES[@]}"; do
  basename_no_gz="$(basename "${input}")"
  basename_no_gz="${basename_no_gz%.gz}"
  if [[ -n "${SEEN[${basename_no_gz}]:-}" ]]; then
    echo "Skipping duplicate representation: ${input}"
    continue
  fi
  SEEN[${basename_no_gz}]=1

  trait="${basename_no_gz%.txt}"
  prepared="${WORK_DIR}/${trait}.sorted.tsv"

  echo "Processing ${trait}"
  prepare_saunders "${input}" "${prepared}"
  harmonize "${WORK_DIR}/ku.sorted.tsv" "${prepared}" ku "${trait}"
  harmonize "${WORK_DIR}/kuf.sorted.tsv" "${prepared}" kuf "${trait}"
  rm -f -- "${prepared}"
done

echo "Done. Outputs are in ${OUT_DIR}"
