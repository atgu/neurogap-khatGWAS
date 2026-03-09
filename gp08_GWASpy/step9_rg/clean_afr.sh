#!/usr/bin/env bash
set -euo pipefail

SUD=/home/nvalenci/sud

# Explicit KU/KUF ldsc.tsv (edit if these move)
KU_LDSC="/home/nvalenci/sud/khat/ku_gp08_meta_cleaned.txt.ldsc.tsv"
KUF_LDSC="/home/nvalenci/sud/khat/kuf.ldsc.tsv"

OUTDIR="${SUD}/munged_inputs_afr"
mkdir -p "$OUTDIR"

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

for f in "$KU_LDSC" "$KUF_LDSC"; do
  [[ -f "$f" ]] || die "Missing required file: $f"
done

log "Collecting AFR *.ldsc.tsv files under $SUD"
mapfile -t AFR_LDSC < <(find "$SUD" -type f -name "*.ldsc.tsv" | grep -i "afr" | sort)

ALL_LDSC=("${AFR_LDSC[@]}" "$KU_LDSC" "$KUF_LDSC")
log "Total files to clean: ${#ALL_LDSC[@]}"
log "Output directory: $OUTDIR"

standardize_one () {
  local in_tsv="$1"
  local base out_tsv header ncol

  base="$(basename "$in_tsv" .ldsc.tsv)"
  out_tsv="${OUTDIR}/${base}.munge.tsv"
  header="$(head -n 1 "$in_tsv")"

  if echo "$header" | grep -qw "N_eff"; then
    ncol="N_eff"
  elif echo "$header" | grep -qw "N_eff_overall"; then
    ncol="N_eff_overall"
  elif echo "$header" | grep -qw "N"; then
    ncol="N"
  else
    die "No usable N column in: $in_tsv"
  fi

  log "Cleaning: $(basename "$in_tsv") -> $(basename "$out_tsv") (N <- $ncol)"

  awk -v OFS="\t" -v NCOL="$ncol" '
    NR==1{
      for(i=1;i<=NF;i++){h[$i]=i}
      print "SNP","CHR","BP","A1","A2","Z","P","N"
      next
    }
    { print $h["SNP"],$h["CHR"],$h["BP"],$h["A1"],$h["A2"],$h["Z"],$h["P"],$h[NCOL] }
  ' "$in_tsv" > "$out_tsv"
}

for f in "${ALL_LDSC[@]}"; do
  standardize_one "$f"
done

log "Done cleaning."
log "Check outputs in: $OUTDIR"
