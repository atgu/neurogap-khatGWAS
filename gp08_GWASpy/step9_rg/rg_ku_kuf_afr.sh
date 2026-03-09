#!/usr/bin/env bash
set -euo pipefail

LDSC=~/ldsc/ldsc.py
REF=/home/nvalenci/sud/ref_panel/HGDP_1KG_AFR_hapmap

AFR_DIR=/home/nvalenci/sud/munged_afr
OUTDIR=/home/nvalenci/sud/ldsc_rg_afr
mkdir -p "$OUTDIR"

# KU + KUF are in munged_afr (per your ls)
KU="${AFR_DIR}/ku_gp08_meta_cleaned.txt.sumstats.gz"
KUF="${AFR_DIR}/kuf.sumstats.gz"

timestamp(){ date "+%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(timestamp)] $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }

# --- checks ---
[[ -f "$LDSC" ]] || die "Missing LDSC script: $LDSC"
[[ -f "${REF}.l2.ldscore.gz" ]] || die "Missing ref panel: ${REF}.l2.ldscore.gz"
[[ -f "${REF}.l2.M" ]] || die "Missing ref panel: ${REF}.l2.M"

[[ -f "$KU" ]] || die "Missing KU sumstats: $KU"
[[ -f "$KUF" ]] || die "Missing KUF sumstats: $KUF"
[[ -d "$AFR_DIR" ]] || die "Missing AFR munged dir: $AFR_DIR"

# Collect AFR sumstats
mapfile -t AFR_SUMSTATS < <(ls -1 "$AFR_DIR"/*.sumstats.gz 2>/dev/null | sort)
[[ ${#AFR_SUMSTATS[@]} -gt 0 ]] || die "No AFR sumstats found in $AFR_DIR"

# Remove KU and KUF from the “other AFR traits” list
AFR_FILTERED=()
for f in "${AFR_SUMSTATS[@]}"; do
  bn="$(basename "$f")"
  if [[ "$f" == "$KU" || "$f" == "$KUF" ]]; then
    continue
  fi
  AFR_FILTERED+=("$f")
done

SUMMARY="${OUTDIR}/rg_summary.tsv"
echo -e "trait1\ttrait2\trg\trg_se\tp\tintercept\tintercept_se" > "$SUMMARY"

run_rg () {
  local a="$1"
  local b="$2"
  local label="$3"
  local outpref="${OUTDIR}/${label}"

  log "rg: $(basename "$a") vs $(basename "$b") -> ${outpref}"

  python "$LDSC" \
    --rg "${a},${b}" \
    --ref-ld "$REF" \
    --w-ld "$REF" \
    --out "$outpref" >/dev/null

  local logfile="${outpref}.log"
  if [[ ! -f "$logfile" ]]; then
    log "WARNING: missing log file: $logfile"
    return 0
  fi

  # Parse values from LDSC log
  local rg rgse p int intse
  rg="$(awk -F': ' '/Genetic Correlation:/ {print $2}' "$logfile" | head -n1 | awk '{print $1}')"
  rgse="$(awk -F'[()]' '/Genetic Correlation:/ {print $2}' "$logfile" | head -n1 | awk '{print $1}')"
  p="$(awk -F': ' '/P:/ {print $2}' "$logfile" | head -n1 | awk '{print $1}')"
  int="$(awk -F': ' '/Intercept:/ {print $2}' "$logfile" | head -n1 | awk '{print $1}')"
  intse="$(awk -F'[()]' '/Intercept:/ {print $2}' "$logfile" | head -n1 | awk '{print $1}')"

  echo -e "$(basename "$a" .sumstats.gz)\t$(basename "$b" .sumstats.gz)\t${rg:-NA}\t${rgse:-NA}\t${p:-NA}\t${int:-NA}\t${intse:-NA}" >> "$SUMMARY"
}

log "AFR sumstats found: ${#AFR_SUMSTATS[@]} total"
log "AFR traits excluding KU/KUF: ${#AFR_FILTERED[@]}"
log "Results dir: $OUTDIR"

# 1) KU vs KUF
run_rg "$KU" "$KUF" "KU_vs_KUF"

# 2) KU vs each other AFR trait
for f in "${AFR_FILTERED[@]}"; do
  run_rg "$KU" "$f" "KU_vs_$(basename "$f" .sumstats.gz)"
done

# 3) KUF vs each other AFR trait
for f in "${AFR_FILTERED[@]}"; do
  run_rg "$KUF" "$f" "KUF_vs_$(basename "$f" .sumstats.gz)"
done

log "Done."
log "Summary (top 30 lines):"
column -t "$SUMMARY" | head -n 30
