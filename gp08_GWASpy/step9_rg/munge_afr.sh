#!/usr/bin/env bash
set -euo pipefail

SUD=/home/nvalenci/sud
LDSC_DIR=/home/nvalenci/ldsc

MUNGE="${LDSC_DIR}/munge_sumstats.py"
HM3="${LDSC_DIR}/w_hm3.snplist"

INDIR="${SUD}/munged_inputs_afr"
OUTDIR="${SUD}/munged_afr"
mkdir -p "$OUTDIR"

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

[[ -f "$MUNGE" ]] || die "Missing: $MUNGE"
[[ -f "$HM3" ]] || die "Missing: $HM3"
[[ -d "$INDIR" ]] || die "Missing directory: $INDIR"

mapfile -t INPUTS < <(ls -1 "$INDIR"/*.munge.tsv 2>/dev/null | sort)
[[ ${#INPUTS[@]} -gt 0 ]] || die "No inputs found in $INDIR (run clean_afr_plus_ku_kuf.sh first)"

log "Munging ${#INPUTS[@]} files from: $INDIR"
log "Output directory: $OUTDIR"

for f in "${INPUTS[@]}"; do
  base="$(basename "$f" .munge.tsv)"
  outpref="${OUTDIR}/${base}"

  log "Munging: $base"
  python "$MUNGE" \
    --sumstats "$f" \
    --out "$outpref" \
    --merge-alleles "$HM3" \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --N-col N \
    --signed-sumstats Z,0
done

log "Done munging."
log "Munged outputs in: $OUTDIR"
