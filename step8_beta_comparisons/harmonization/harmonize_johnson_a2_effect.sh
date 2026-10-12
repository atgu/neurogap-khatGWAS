#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

JOHNSON_DIR="${HOME}/sud/johnson-2020/liftover_only"
KHAT_DIR="${HOME}/sud/khat/a2_effect"
OUT_DIR="${HOME}/rsid_match/johnson_lifted_harmonized/a2_effect"

KU_FILE="${KHAT_DIR}/ku_gp08_meta_cleaned_with_N.txt"
KUF_FILE="${KHAT_DIR}/kuf_gp08_meta_cleaned_with_N.tsv"

JOHNSON_FILES=(
  "CUD_AFR_full_public_11.14.2020.lifted_hg38.tsv"
  "CUD_EUR_casecontrol_public_11.14.2020.lifted_hg38.tsv"
  "CUD_EUR_full_public_11.14.2020.lifted_hg38.tsv"
)

mkdir -p "${OUT_DIR}"
for file in "${KU_FILE}" "${KUF_FILE}"; do
  [[ -f "${file}" ]] || { echo "ERROR: required file not found: ${file}" >&2; exit 1; }
done
for filename in "${JOHNSON_FILES[@]}"; do
  [[ -f "${JOHNSON_DIR}/${filename}" ]] || {
    echo "ERROR: required Johnson file not found: ${JOHNSON_DIR}/${filename}" >&2
    exit 1
  }
done

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/johnson_harmonize.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
echo "Temporary workspace: ${WORK_DIR}"

prepare_khat() {
  local input=$1 output=$2
  # key marker CHR BP effect_allele other_allele beta pvalue
  awk 'BEGIN{FS=OFS="\t"}
    function normchr(x){gsub(/\r|"/,"",x); sub(/^chr/,"",x); sub(/\.0$/,"",x); return x}
    function snv(x){return x~/^[ACGT]$/}
    function amb(a,b){return a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC"}
    NR==1{
      for(i=1;i<=NF;i++){
        h=$i; gsub(/\r|"/,"",h)
        if(h=="MarkerName")m=i; else if(h=="CHR")c=i; else if(h=="BP")pos=i
        else if(h=="Allele1")a1=i; else if(h=="Allele2")a2=i
        else if(h=="Effect")b=i; else if(h=="P-value")p=i
      }
      if(!m||!c||!pos||!a1||!a2||!b||!p){print "ERROR: bad Khat header in " FILENAME > "/dev/stderr"; exit 2}
      next
    }
    {
      marker=$(m); gsub(/\r|"/,"",marker); chr=normchr($(c)); bp=int($(pos))
      x=toupper($(a1)); y=toupper($(a2)); gsub(/\r|"/,"",x); gsub(/\r|"/,"",y)
      if(marker!=""&&chr!=""&&bp>0&&snv(x)&&snv(y)&&!amb(x,y)&&$(b)!=""&&$(p)!="")
        print chr ":" bp,marker,chr,bp,x,y,$(b),$(p)
    }' "${input}" | sort -t $'\t' -k1,1 -k5,6 -u > "${output}"
}

prepare_johnson() {
  local input=$1 output=$2
  # key marker effect_A1 other_A2 raw_effect P SE N N_CAS N_CON effect_type
  awk 'BEGIN{FS=OFS="\t"}
    function normchr(x){gsub(/\r|"/,"",x); sub(/^chr/,"",x); sub(/\.0$/,"",x); return x}
    function snv(x){return x~/^[ACGT]$/}
    function amb(a,b){return a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC"}
    NR==1{
      for(i=1;i<=NF;i++){
        h=$i; gsub(/\r|"/,"",h); hu=toupper(h)
        if(hu=="CHR_HG38")c=i; else if(hu=="BP_HG38")pos=i; else if(hu=="SNP")m=i
        else if(hu=="A1")a1=i; else if(hu=="A2")a2=i
        else if(hu=="BETA")beta=i; else if(hu=="Z")z=i; else if(hu=="P")p=i
        else if(hu=="SE")se=i; else if(hu=="N")n=i; else if(hu=="N_CAS")ncas=i; else if(hu=="N_CON")ncon=i
      }
      if(!c||!pos||!m||!a1||!a2||(!beta&&!z)||!p){print "ERROR: bad Johnson header in " FILENAME > "/dev/stderr"; exit 2}
      effect=(beta?beta:z); etype=(beta?"beta":"z"); next
    }
    {
      chr=normchr($(c)); bp=int($(pos)); marker=$(m); gsub(/\r|"/,"",marker)
      ea=toupper($(a1)); oa=toupper($(a2)); gsub(/\r|"/,"",ea); gsub(/\r|"/,"",oa)
      if(chr!=""&&bp>0&&marker!=""&&snv(ea)&&snv(oa)&&!amb(ea,oa)&&$(effect)!="")
        print chr ":" bp,marker,ea,oa,$(effect),$(p),(se?$(se):""),(n?$(n):""),(ncas?$(ncas):""),(ncon?$(ncon):""),etype
    }' "${input}" | sort -t $'\t' -k1,1 -k3,4 -u > "${output}"
}

harmonize() {
  local khat=$1 johnson=$2 phenotype=$3 trait=$4 output
  output="${OUT_DIR}/${phenotype}_${trait}_harmonized.txt"

  printf 'khat_marker\tCHR\tBP\tkhat_A1\tkhat_A2\tbeta_%s\tp_%s\tjohnson_marker\tjohnson_effect_allele_A1\tjohnson_other_allele_A2\tjohnson_effect_type\teffect_%s\teffect_%s_to_khat_A1\tP\tSE\tN\tN_CAS\tN_CON\tallele_match_status\n' \
    "${phenotype}" "${phenotype}" "${trait}" "${trait}" > "${output}"

  # Joined fields: key; Khat fields 2-8; Johnson fields 9-18.
  join -t $'\t' -1 1 -2 1 "${khat}" "${johnson}" |
  awk 'BEGIN{FS=OFS="\t"}
    function comp(x){return x=="A"?"T":x=="T"?"A":x=="C"?"G":x=="G"?"C":""}
    {
      k1=$5; k2=$6; effect=$10; other=$11; status=""; flip=0
      # Johnson A1 is its effect allele. Align that effect to Khat Allele1.
      if(k1==effect&&k2==other) status="effect_match"
      else if(k1==other&&k2==effect){status="effect_swapped";flip=1}
      else if(k1==comp(effect)&&k2==comp(other)) status="effect_complement"
      else if(k1==comp(other)&&k2==comp(effect)){status="effect_complement_swapped";flip=1}
      else next
      raw=$12+0; aligned=(flip?-raw:raw)
      print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$18,raw,aligned,$13,$14,$15,$16,$17,status
    }' >> "${output}"

  awk 'BEGIN{FS="\t"} NR>1{n++; if($19~/swapped$/)f++} END{printf "  Wrote %s: %d matches; %d flipped\n",FILENAME,n+0,f+0}' "${output}"
}

echo "Preparing KU"
prepare_khat "${KU_FILE}" "${WORK_DIR}/ku.tsv"
echo "Preparing KUF"
prepare_khat "${KUF_FILE}" "${WORK_DIR}/kuf.tsv"

for filename in "${JOHNSON_FILES[@]}"; do
  trait="${filename%.lifted_hg38.tsv}"
  prepared="${WORK_DIR}/${trait}.tsv"
  echo "Processing ${trait}"
  prepare_johnson "${JOHNSON_DIR}/${filename}" "${prepared}"
  harmonize "${WORK_DIR}/ku.tsv" "${prepared}" ku "${trait}"
  harmonize "${WORK_DIR}/kuf.tsv" "${prepared}" kuf "${trait}"
  rm -f -- "${prepared}"
done

echo "Done. Outputs are in ${OUT_DIR}"
