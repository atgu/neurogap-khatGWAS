#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

WALTERS_DIR="${HOME}/sud/walters-2018"
OUT_DIR="${HOME}/sud/walters-2018/liftover_only"
CHAIN="${HOME}/liftover_files/hg19ToHg38.over.chain.gz"
LIFTOVER_BIN="${HOME}/liftOver"

FILES=(
  "pgc_alcdep.discovery.aug2018_release.txt"
  "pgc_alcdep.afr_discovery.aug2018_release.txt"
  "pgc_alcdep.eur_discovery.aug2018_release.txt"
  "pgc_alcdep.afr_unrelated.aug2018_release.txt"
  "pgc_alcdep.eur_unrelated.aug2018_release.txt"
  "pgc_alcdep.afr_unrel_genotyped.aug2018_release.txt"
  "pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt"
  "pgc_alcdep.trans_fe_unrel_geno.aug2018_release.txt"
  "pgc_alcdep.trans_re2_unrel_geno.aug2018_release.txt"
  "pgc_alcdep.trans_mantra_unrel_geno.aug2018_release.txt"
)

mkdir -p "${OUT_DIR}"
[[ -x "${LIFTOVER_BIN}" ]] || { echo "ERROR: liftOver not executable: ${LIFTOVER_BIN}" >&2; exit 1; }
[[ -f "${CHAIN}" ]] || { echo "ERROR: chain file not found: ${CHAIN}" >&2; exit 1; }

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/walters_liftover.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
echo "Temporary workspace: ${WORK_DIR}"

stream_file() {
  local input=$1
  if [[ "${input}" == *.gz ]]; then gzip -cd -- "${input}"; else command cat -- "${input}"; fi
}

processed=0
for filename in "${FILES[@]}"; do
  if [[ -f "${WALTERS_DIR}/${filename}" ]]; then
    input="${WALTERS_DIR}/${filename}"
  elif [[ -f "${WALTERS_DIR}/${filename}.gz" ]]; then
    input="${WALTERS_DIR}/${filename}.gz"
  else
    echo "WARNING: skipping missing Walters file: ${filename}[.gz]" >&2
    continue
  fi

  echo "Processing ${filename}"
  prefix="${OUT_DIR}/${filename}"
  header="${WORK_DIR}/${filename}.header"
  pre="${WORK_DIR}/${filename}.pre"
  pre_sorted="${WORK_DIR}/${filename}.pre.sorted"
  lookup="${WORK_DIR}/${filename}.lookup"
  lookup_sorted="${WORK_DIR}/${filename}.lookup.sorted"
  bed_hg19="${prefix}.hg19.bed"
  bed_hg38="${prefix}.hg38.bed"
  unmapped="${prefix}.unmapped.bed"
  output="${prefix}.lifted_hg38.tsv"

  stream_file "${input}" |
  awk -v HEADER="${header}" -v PRE="${pre}" -v BED="${bed_hg19}" '
    BEGIN{FS="[[:space:]]+"; OFS="\t"}
    function snv(x){return x~/^[ACGT]$/}
    function amb(a,b){return a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC"}
    NR==1{
      for(i=1;i<=NF;i++){
        h=toupper($i)
        if(h=="CHR")chr=i; else if(h=="SNP")snp=i; else if(h=="BP")bp=i
        else if(h=="A1")a1=i; else if(h=="A2")a2=i
        printf "%s%s",$i,(i<NF?OFS:ORS) > HEADER
      }
      close(HEADER)
      if(!chr||!snp||!bp||!a1||!a2){print "ERROR: required Walters columns missing" > "/dev/stderr"; exit 2}
      next
    }
    {
      marker=$(snp); c=$(chr); sub(/^chr/,"",c); pos=int($(bp))
      effect=toupper($(a1)); other=toupper($(a2))
      if(marker==""||marker=="."||marker=="NA"||c==""||pos<=0||!snv(effect)||!snv(other)||amb(effect,other)) next
      $(chr)=c; $(bp)=pos; $(a1)=effect; $(a2)=other
      id=sprintf("var%012d",++kept)
      printf "%s",id > PRE; for(i=1;i<=NF;i++) printf "%s%s%s",OFS,$i,(i<NF?"":ORS) > PRE
      print "chr" c,pos-1,pos,id > BED
    }
  '

  "${LIFTOVER_BIN}" "${bed_hg19}" "${CHAIN}" "${bed_hg38}" "${unmapped}"

  sort -t $'\t' -k1,1 "${pre}" > "${pre_sorted}"
  awk 'BEGIN{FS=OFS="\t"}{c=$1;sub(/^chr/,"",c);print $4,c,$3}' "${bed_hg38}" > "${lookup}"
  sort -t $'\t' -k1,1 "${lookup}" > "${lookup_sorted}"

  awk 'BEGIN{FS=OFS="\t"}{print $0,"row_id","CHR_hg38","BP_hg38"}' "${header}" > "${output}"
  join -t $'\t' -1 1 -2 1 "${pre_sorted}" "${lookup_sorted}" |
  awk 'BEGIN{FS=OFS="\t"}{for(i=2;i<=NF-2;i++)printf "%s%s",$i,OFS; print $1,$(NF-1),$NF}' >> "${output}"

  input_n="$(wc -l < "${pre}")"
  lifted_n="$(wc -l < "${bed_hg38}")"
  echo "  Filtered nonambiguous SNVs: ${input_n}"
  echo "  Lifted records:             ${lifted_n}"
  echo "  Wrote: ${output}"
  processed=$((processed+1))
done

(( processed > 0 )) || { echo "ERROR: no Walters files were processed" >&2; exit 1; }
echo "Done. Lifted files are in ${OUT_DIR}"
