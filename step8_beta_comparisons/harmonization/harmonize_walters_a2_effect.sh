#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
shopt -s nullglob

WALTERS_DIR="${HOME}/sud/walters-2018/liftover_only"
KHAT_DIR="${HOME}/sud/khat/a2_effect"
OUT_DIR="${HOME}/rsid_match/walters_lifted_harmonized/a2_effect"

KU_FILE="${KHAT_DIR}/ku_gp08_meta_cleaned_with_N.txt"
KUF_FILE="${KHAT_DIR}/kuf_gp08_meta_cleaned_with_N.tsv"

mkdir -p "${OUT_DIR}"
for file in "${KU_FILE}" "${KUF_FILE}"; do
  [[ -f "${file}" ]] || { echo "ERROR: required Khat file not found: ${file}" >&2; exit 1; }
done

WALTERS_FILES=("${WALTERS_DIR}"/pgc_alcdep.*.aug2018_release.txt.lifted_hg38.tsv)
(( ${#WALTERS_FILES[@]} > 0 )) || { echo "ERROR: no lifted Walters files found in ${WALTERS_DIR}" >&2; exit 1; }

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/walters_harmonize.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
echo "Temporary workspace: ${WORK_DIR}"

prepare_khat() {
  local input=$1 output=$2
  awk 'BEGIN{FS=OFS="\t"}
    function clean(x){gsub(/\r|"/,"",x);return x}
    function normchr(x){x=clean(x);sub(/^chr/,"",x);sub(/\.0$/,"",x);return x}
    function snv(x){return x~/^[ACGT]$/}
    function amb(a,b){return a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC"}
    NR==1{
      for(i=1;i<=NF;i++){
        h=clean($i)
        if(h=="MarkerName")m=i;else if(h=="CHR")c=i;else if(h=="BP")pos=i
        else if(h=="Allele1")a1=i;else if(h=="Allele2")a2=i
        else if(h=="Effect")b=i;else if(h=="P-value")p=i
      }
      if(!m||!c||!pos||!a1||!a2||!b||!p){print "ERROR: bad Khat header in " FILENAME > "/dev/stderr";exit 2}
      next
    }
    {
      marker=clean($(m));chr=normchr($(c));bp=int($(pos));effect=toupper(clean($(a1)));other=toupper(clean($(a2)))
      if(marker!=""&&chr!=""&&bp>0&&$(b)!=""&&$(p)!=""&&snv(effect)&&snv(other)&&!amb(effect,other))
        print chr ":" bp,marker,chr,bp,effect,other,$(b),$(p)
    }' "${input}" | sort -t $'\t' -k1,1 -k5,6 -u > "${output}"
}

prepare_walters() {
  local input=$1 output=$2
  # Output: key marker effect_A1 other_A2 raw P INFO SE Neff Weight type
  awk 'BEGIN{FS=OFS="\t"}
    function clean(x){gsub(/\r|"/,"",x);return x}
    function normchr(x){x=clean(x);sub(/^chr/,"",x);sub(/\.0$/,"",x);return x}
    function snv(x){return x~/^[ACGT]$/}
    function amb(a,b){return a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC"}
    NR==1{
      for(i=1;i<=NF;i++){
        h=toupper(clean($i))
        if(h=="CHR_HG38")c=i;else if(h=="BP_HG38")pos=i;else if(h=="SNP")m=i
        else if(h=="A1")a1=i;else if(h=="A2")a2=i;else if(h=="BETA")beta=i
        else if(h=="Z")z=i;else if(h=="OR")orcol=i;else if(h=="P")p=i
        else if(h=="INFO")info=i;else if(h=="SE")se=i;else if(h=="NEFF")neff=i;else if(h=="WEIGHT")weight=i
      }
      if(!c||!pos||!m||!a1||!a2||!p){print "ERROR: bad Walters header in " FILENAME > "/dev/stderr";exit 2}
      if(beta){effectcol=beta;type="BETA"}else if(z){effectcol=z;type="Z"}else if(orcol){effectcol=orcol;type="OR"}
      else exit 3
      next
    }
    {
      chr=normchr($(c));bp=int($(pos));marker=clean($(m));effect=toupper(clean($(a1)));other=toupper(clean($(a2)));raw=clean($(effectcol))
      if(chr!=""&&bp>0&&marker!=""&&raw!=""&&$(p)!=""&&snv(effect)&&snv(other)&&!amb(effect,other))
        print chr ":" bp,marker,effect,other,raw,$(p),(info?$(info):""),(se?$(se):""),(neff?$(neff):""),(weight?$(weight):""),type
    }' "${input}" | sort -t $'\t' -k1,1 -k3,4 -u > "${output}"
}

harmonize() {
  local khat=$1 walters=$2 phenotype=$3 trait=$4 type=$5 output
  output="${OUT_DIR}/${phenotype}_${trait}_harmonized.txt"

  printf 'khat_marker\tCHR\tBP\tkhat_A1\tkhat_A2\tbeta_%s\tp_%s\twalters_marker\twalters_effect_allele_A1\twalters_other_allele_A2\twalters_effect_type\teffect_%s_raw\teffect_%s_to_khat_A1\tlogOR_%s_to_khat_A1\twalters_P\twalters_INFO\twalters_SE\twalters_Neff\twalters_Weight\tallele_match_status\n' \
    "${phenotype}" "${phenotype}" "${trait}" "${trait}" "${trait}" > "${output}"

  # Joined fields: 1 key; 2-8 Khat; 9-18 Walters.
  join -t $'\t' -1 1 -2 1 "${khat}" "${walters}" |
  awk -v TYPE="${type}" 'BEGIN{FS=OFS="\t"}
    function comp(x){return x=="A"?"T":x=="T"?"A":x=="C"?"G":x=="G"?"C":""}
    {
      k1=$5;k2=$6;effect=$10;other=$11;status="";flip=0
      if(k1==effect&&k2==other)status="effect_match"
      else if(k1==other&&k2==effect){status="effect_swapped";flip=1}
      else if(k1==comp(effect)&&k2==comp(other))status="effect_complement"
      else if(k1==comp(other)&&k2==comp(effect)){status="effect_complement_swapped";flip=1}
      else next

      raw=$12+0;logor=""
      if(TYPE=="OR"){
        if(raw<=0)next
        aligned=(flip?1/raw:raw)
        logor=(flip?-log(raw):log(raw))
      }else aligned=(flip?-raw:raw)

      print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,TYPE,raw,aligned,logor,$13,$14,$15,$16,$17,status
    }' >> "${output}"

  awk 'BEGIN{FS="\t"}NR>1{n++;if($20~/swapped$/)f++}END{printf "  Wrote %s: %d matches; %d reoriented\n",FILENAME,n+0,f+0}' "${output}"
}

echo "Preparing KU"
prepare_khat "${KU_FILE}" "${WORK_DIR}/ku.tsv"
echo "Preparing KUF"
prepare_khat "${KUF_FILE}" "${WORK_DIR}/kuf.tsv"

for input in "${WALTERS_FILES[@]}"; do
  filename="$(basename "${input}")"
  trait="${filename%.txt.lifted_hg38.tsv}"

  header="$(head -n 1 "${input}" | tr '\t' ' ')"
  if [[ " ${header^^} " == *" BETA "* ]]; then type="BETA"
  elif [[ " ${header^^} " == *" Z "* ]]; then type="Z"
  elif [[ " ${header^^} " == *" OR "* ]]; then type="OR"
  else
    echo "Skipping ${trait}: no signed BETA/Z or directional OR field"
    continue
  fi

  prepared="${WORK_DIR}/${trait}.tsv"
  echo "Processing ${trait} (${type})"
  prepare_walters "${input}" "${prepared}"
  harmonize "${WORK_DIR}/ku.tsv" "${prepared}" ku "${trait}" "${type}"
  harmonize "${WORK_DIR}/kuf.tsv" "${prepared}" kuf "${trait}" "${type}"
  rm -f -- "${prepared}"
done

echo "Done. Outputs are in ${OUT_DIR}"
