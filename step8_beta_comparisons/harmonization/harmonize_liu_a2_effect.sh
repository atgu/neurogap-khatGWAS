#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

LIU_DIR="${HOME}/sud/liu-2019/liftover_only"
KHAT_DIR="${HOME}/sud/khat/a2_effect"
OUT_DIR="${HOME}/rsid_match/liu_lifted_harmonized/a2_effect"

KU_FILE="${KHAT_DIR}/ku_gp08_meta_cleaned_with_N.txt"
KUF_FILE="${KHAT_DIR}/kuf_gp08_meta_cleaned_with_N.tsv"

LIU_FILES=(
  "AgeOfInitiation.WithoutUKB.txt.lifted_hg38.tsv"
  "AgeofInitiation.txt.lifted_hg38.tsv"
  "CigarettesPerDay.WithoutUKB.txt.lifted_hg38.tsv"
  "CigarettesPerDay.txt.lifted_hg38.tsv"
  "DrinksPerWeek.WithoutUKB.txt.lifted_hg38.tsv"
  "DrinksPerWeek.txt.lifted_hg38.tsv"
  "SmokingCessation.WithoutUKB.txt.lifted_hg38.tsv"
  "SmokingCessation.txt.lifted_hg38.tsv"
  "SmokingInitiation.WithoutUKB.txt.lifted_hg38.tsv"
  "SmokingInitiation.txt.lifted_hg38.tsv"
)

mkdir -p "${OUT_DIR}"

for file in "${KU_FILE}" "${KUF_FILE}"; do
  if [[ ! -f "${file}" ]]; then
    echo "ERROR: required file not found: ${file}" >&2
    exit 1
  fi
done

for filename in "${LIU_FILES[@]}"; do
  if [[ ! -f "${LIU_DIR}/${filename}" ]]; then
    echo "ERROR: required Liu file not found: ${LIU_DIR}/${filename}" >&2
    exit 1
  fi
done

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/liu_harmonize.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT

echo "Temporary workspace: ${WORK_DIR}"

prepare_khat() {
  local input=$1
  local output=$2

  # Output fields:
  # key marker CHR BP A1 A2 beta pvalue
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
        print "ERROR: Khat header missing required columns in " FILENAME > "/dev/stderr"; exit 2
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

prepare_liu() {
  local input=$1
  local output=$2

  # Liu documentation states that Beta is based on ALT.
  # Fixed output fields:
  # key marker effect_allele(ALT) other_allele(REF) Beta PValue SE N Effective_N AF STAT studies ANNO ANNOFULL
  awk 'BEGIN{FS="\t"; OFS="\t"}
    function norm_chr(x) {sub(/^chr/, "", x); sub(/\.0$/, "", x); return x}
    function snv(x) {return x ~ /^[ACGT]$/}
    function ambiguous(a,b) {return (a b=="AT" || a b=="TA" || a b=="CG" || a b=="GC")}
    NR==1 {
      for(i=1;i<=NF;i++) {
        h=$i; gsub(/\r|"/, "", h)
        hu=toupper(h)
        if(hu=="CHR_HG38") chr=i; else if(hu=="BP_HG38") bp=i; else if(hu=="RSID") marker=i
        else if(hu=="REF") ref=i; else if(hu=="ALT") alt=i; else if(hu=="BETA") beta=i
        else if(hu=="PVALUE") p=i; else if(hu=="SE") se=i; else if(hu=="N") n=i
        else if(hu=="EFFECTIVE_N") neff=i; else if(hu=="AF") af=i; else if(hu=="STAT") stat=i
        else if(hu=="NUMBER_OF_STUDIES") studies=i; else if(hu=="ANNO") anno=i; else if(hu=="ANNOFULL") annofull=i
      }
      if(!chr||!bp||!marker||!ref||!alt||!beta||!p||!se||!n||!neff) {
        print "ERROR: Liu header missing required columns in " FILENAME > "/dev/stderr"; exit 2
      }
      next
    }
    {
      c=norm_chr($(chr)); pos=int($(bp)); m=$(marker)
      effect=toupper($(alt)); other=toupper($(ref))
      gsub(/"/, "", m)
      if(c!="" && pos>0 && m!="" && snv(effect) && snv(other) && !ambiguous(effect,other)) {
        afv=(af ? $(af) : ""); stv=(stat ? $(stat) : ""); nsv=(studies ? $(studies) : "")
        av=(anno ? $(anno) : ""); afull=(annofull ? $(annofull) : "")
        print c ":" pos, m, effect, other, $(beta), $(p), $(se), $(n), $(neff), afv, stv, nsv, av, afull
      }
    }
  ' "${input}" | sort -t $'\t' -k1,1 -k3,4 -u > "${output}"
}

harmonize() {
  local khat_sorted=$1
  local liu_sorted=$2
  local phenotype=$3
  local trait=$4
  local output="${OUT_DIR}/${phenotype}_${trait}_harmonized.txt"

  printf 'khat_marker\tCHR\tBP\tkhat_A1\tkhat_A2\tbeta_%s\tp_%s\tliu_marker\tliu_effect_allele_ALT\tliu_other_allele_REF\tbeta_%s\tbeta_%s_to_khat_A1\tPVALUE\tSE\tN\tEFFECTIVE_N\tallele_match_status\tAF\tSTAT\tNumber_of_Studies\tANNO\tANNOFULL\n' \
    "${phenotype}" "${phenotype}" "${trait}" "${trait}" > "${output}"

  # Joined fields:
  # 1 key
  # 2-8  Khat: marker CHR BP A1 A2 beta pvalue
  # 9-21 Liu: marker A1 A2 BETA PVALUE SE N EFFECTIVE_N AF STAT studies ANNO ANNOFULL
  join -t $'\t' -1 1 -2 1 "${khat_sorted}" "${liu_sorted}" |
  awk 'BEGIN{FS=OFS="\t"}
    function comp(x) {return x=="A"?"T":x=="T"?"A":x=="C"?"G":x=="G"?"C":""}
    {
      k1=$5; k2=$6; effect=$10; other=$11; status=""; flip=0
      # Liu Beta is for ALT. Keep it when ALT is Khat Allele1;
      # negate it when ALT is Khat Allele2.
      if(k1==effect && k2==other) status="effect_match"
      else if(k1==other && k2==effect) {status="effect_swapped"; flip=1}
      else if(k1==comp(effect) && k2==comp(other)) status="effect_complement"
      else if(k1==comp(other) && k2==comp(effect)) {status="effect_complement_swapped"; flip=1}
      else next

      raw=$12+0; aligned=(flip ? -raw : raw)
      print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,raw,aligned,$13,$14,$15,$16,status,$17,$18,$19,$20,$21
    }
  ' >> "${output}"

  echo "Wrote ${output}"
}

echo "Preparing and sorting KU"
prepare_khat "${KU_FILE}" "${WORK_DIR}/ku.sorted.tsv"

echo "Preparing and sorting KUF"
prepare_khat "${KUF_FILE}" "${WORK_DIR}/kuf.sorted.tsv"

for filename in "${LIU_FILES[@]}"; do
  trait="${filename%.txt.lifted_hg38.tsv}"
  trait="${trait/.WithoutUKB/_WithoutUKB}"
  liu_sorted="${WORK_DIR}/liu.${trait}.sorted.tsv"

  echo "Preparing Liu trait: ${trait}"
  prepare_liu "${LIU_DIR}/${filename}" "${liu_sorted}"

  echo "Harmonizing KU with ${trait}"
  harmonize "${WORK_DIR}/ku.sorted.tsv" "${liu_sorted}" ku "${trait}"

  echo "Harmonizing KUF with ${trait}"
  harmonize "${WORK_DIR}/kuf.sorted.tsv" "${liu_sorted}" kuf "${trait}"

  rm -f -- "${liu_sorted}"
done

echo "Done. Outputs are in ${OUT_DIR}"
