#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

CANLU_DIR="${HOME}/sud/bright-2026"
KHAT_DIR="${HOME}/sud/khat/a2_effect"
OUT_DIR="${HOME}/rsid_match/canlu_harmonized/a2_effect"

KU_FILE="${KHAT_DIR}/ku_gp08_meta_cleaned_with_N.txt"
KUF_FILE="${KHAT_DIR}/kuf_gp08_meta_cleaned_with_N.tsv"
EUR_FILE="${CANLU_DIR}/canlu_eur_meta.txt"
CROSS_FILE="${CANLU_DIR}/canlu_cross_ancestry"

mkdir -p "${OUT_DIR}"

for file in "${KU_FILE}" "${KUF_FILE}" "${EUR_FILE}" "${CROSS_FILE}"; do
  if [[ ! -f "${file}" ]]; then
    echo "ERROR: required file not found: ${file}" >&2
    exit 1
  fi
done

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/canlu_harmonize.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT

echo "Temporary workspace: ${WORK_DIR}"
echo "Preparing CanLU files"

# Output columns:
# key rsid A1 A2 A1_FREQ beta SE Z P N
awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"}
  function norm_chr(x) {sub(/^chr/, "", x); sub(/\.0$/, "", x); return x}
  function snv(x) {return x ~ /^[ACGT]$/}
  function ambiguous(a,b) {return (a b=="AT" || a b=="TA" || a b=="CG" || a b=="GC")}
  NR==1 {
    for(i=1;i<=NF;i++) {
      if($i=="RSID") rsid=i; else if($i=="CHR") chr=i; else if($i=="BP") bp=i
      else if($i=="A1") a1=i; else if($i=="A2") a2=i; else if($i=="A1_FREQ") freq=i
      else if($i=="BETA") beta=i; else if($i=="SE") se=i; else if($i=="Z") z=i
      else if($i=="P") p=i; else if($i=="N") n=i
    }
    if(!rsid||!chr||!bp||!a1||!a2||!freq||!beta||!se||!z||!p||!n) {
      print "ERROR: CanLU EUR header is missing required columns" > "/dev/stderr"; exit 2
    }
    next
  }
  {
    c=norm_chr($(chr)); pos=int($(bp)); x=toupper($(a1)); y=toupper($(a2))
    if(c!="" && pos>0 && snv(x) && snv(y) && !ambiguous(x,y))
      print c ":" pos, $(rsid), x, y, $(freq), $(beta), $(se), $(z), $(p), $(n)
  }
' "${EUR_FILE}" | sort -t $'\t' -k1,1 -k3,4 -u > "${WORK_DIR}/eur.sorted.tsv"

# Output columns:
# key rsid marker A1 A2 Z P weight direction
awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"}
  function norm_chr(x) {sub(/^chr/, "", x); sub(/\.0$/, "", x); return x}
  function snv(x) {return x ~ /^[ACGT]$/}
  function ambiguous(a,b) {return (a b=="AT" || a b=="TA" || a b=="CG" || a b=="GC")}
  NR==1 {
    for(i=1;i<=NF;i++) {
      if($i=="rsID") rsid=i; else if($i=="MarkerName") marker=i
      else if($i=="CHR") chr=i; else if($i=="BP") bp=i
      else if($i=="Allele1") a1=i; else if($i=="Allele2") a2=i
      else if($i=="Zscore") z=i; else if($i=="P-value") p=i
      else if($i=="Weight") weight=i; else if($i=="Direction") direction=i
    }
    if(!rsid||!marker||!chr||!bp||!a1||!a2||!z||!p||!weight||!direction) {
      print "ERROR: CanLU cross-ancestry header is missing required columns" > "/dev/stderr"; exit 2
    }
    next
  }
  {
    c=norm_chr($(chr)); pos=int($(bp)); x=toupper($(a1)); y=toupper($(a2))
    if(c!="" && pos>0 && snv(x) && snv(y) && !ambiguous(x,y))
      print c ":" pos, $(rsid), $(marker), x, y, $(z), $(p), $(weight), $(direction)
  }
' "${CROSS_FILE}" | sort -t $'\t' -k1,1 -k4,5 -u > "${WORK_DIR}/cross.sorted.tsv"

prepare_khat() {
  local input=$1
  local output=$2

  # Output columns:
  # key khat_marker CHR BP khat_A1 khat_A2 khat_beta khat_p
  awk 'BEGIN{FS="\t"; OFS="\t"}
    function norm_chr(x) {sub(/^chr/, "", x); sub(/\.0$/, "", x); return x}
    function snv(x) {return x ~ /^[ACGT]$/}
    function ambiguous(a,b) {return (a b=="AT" || a b=="TA" || a b=="CG" || a b=="GC")}
    NR==1 {
      for(i=1;i<=NF;i++) {
        h=$i; gsub(/\r|"/, "", h)
        if(h=="MarkerName") marker=i; else if(h=="CHR") chr=i; else if(h=="BP") bp=i
        else if(h=="Allele1") a1=i; else if(h=="Allele2") a2=i
        else if(h=="Effect") beta=i; else if(h=="P-value") p=i
      }
      if(!marker||!chr||!bp||!a1||!a2||!beta||!p) {
        print "ERROR: Khat header is missing required columns in " FILENAME > "/dev/stderr"; exit 2
      }
      next
    }
    {
      m=$(marker); gsub(/"/, "", m); c=norm_chr($(chr)); pos=int($(bp))
      x=toupper($(a1)); y=toupper($(a2)); b=$(beta); pv=$(p)
      if(m!="" && c!="" && pos>0 && b!="" && pv!="" && snv(x) && snv(y) && !ambiguous(x,y))
        print c ":" pos, m, c, pos, x, y, b, pv
    }
  ' "${input}" | sort -t $'\t' -k1,1 -k5,6 -u > "${output}"
}

harmonize_eur() {
  local khat_sorted=$1
  local label=$2
  local output="${OUT_DIR}/${label}_canlu_eur_harmonized.txt"

  printf 'khat_marker\tCHR\tBP\tkhat_A1\tkhat_A2\tbeta_%s\tp_%s\trsid\tcanlu_A1\tcanlu_A2\tcanlu_A1_FREQ\tbeta_canlu_raw\tbeta_canlu_to_khat_A1\tcanlu_SE\tcanlu_Z\tcanlu_P\tcanlu_N\tallele_match_status\n' "${label}" "${label}" > "${output}"

  # Joined fields:
  # 1 key; 2-8 Khat; 9-17 CanLU EUR
  join -t $'\t' -1 1 -2 1 "${khat_sorted}" "${WORK_DIR}/eur.sorted.tsv" |
  awk 'BEGIN{FS=OFS="\t"}
    function comp(x) {return x=="A"?"T":x=="T"?"A":x=="C"?"G":x=="G"?"C":""}
    {
      k1=$5; k2=$6; c1=$10; c2=$11; status=""; flip=0
      if(k1==c1 && k2==c2) status="exact"
      else if(k1==c2 && k2==c1) {status="swapped"; flip=1}
      else if(k1==comp(c1) && k2==comp(c2)) status="complement"
      else if(k1==comp(c2) && k2==comp(c1)) {status="complement_swapped"; flip=1}
      else next
      raw=$13+0; aligned=(flip ? -raw : raw)
      print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,raw,aligned,$14,$15,$16,$17,status
    }
  ' >> "${output}"
  echo "Wrote ${output}"
}

harmonize_cross() {
  local khat_sorted=$1
  local label=$2
  local output="${OUT_DIR}/${label}_canlu_cross_ancestry_harmonized.txt"

  printf 'khat_marker\tCHR\tBP\tkhat_A1\tkhat_A2\tbeta_%s\tp_%s\trsid\tcanlu_marker\tcanlu_A1\tcanlu_A2\tz_canlu_raw\tz_canlu_to_khat_A1\tcanlu_P\tcanlu_weight\tcanlu_direction\tallele_match_status\n' "${label}" "${label}" > "${output}"

  # Joined fields:
  # 1 key; 2-8 Khat; 9-16 CanLU cross-ancestry
  join -t $'\t' -1 1 -2 1 "${khat_sorted}" "${WORK_DIR}/cross.sorted.tsv" |
  awk 'BEGIN{FS=OFS="\t"}
    function comp(x) {return x=="A"?"T":x=="T"?"A":x=="C"?"G":x=="G"?"C":""}
    {
      k1=$5; k2=$6; c1=$11; c2=$12; status=""; flip=0
      if(k1==c1 && k2==c2) status="exact"
      else if(k1==c2 && k2==c1) {status="swapped"; flip=1}
      else if(k1==comp(c1) && k2==comp(c2)) status="complement"
      else if(k1==comp(c2) && k2==comp(c1)) {status="complement_swapped"; flip=1}
      else next
      raw=$13+0; aligned=(flip ? -raw : raw)
      print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,raw,aligned,$14,$15,$16,status
    }
  ' >> "${output}"
  echo "Wrote ${output}"
}

process_khat() {
  local input=$1
  local label=$2
  local khat_sorted="${WORK_DIR}/${label}.sorted.tsv"

  echo "Preparing ${label^^}"
  prepare_khat "${input}" "${khat_sorted}"
  echo "Harmonizing ${label^^} with CanLU EUR"
  harmonize_eur "${khat_sorted}" "${label}"
  echo "Harmonizing ${label^^} with CanLU cross-ancestry"
  harmonize_cross "${khat_sorted}" "${label}"
  rm -f -- "${khat_sorted}"
}

process_khat "${KU_FILE}" ku
process_khat "${KUF_FILE}" kuf

echo "Done. Outputs are in ${OUT_DIR}"
