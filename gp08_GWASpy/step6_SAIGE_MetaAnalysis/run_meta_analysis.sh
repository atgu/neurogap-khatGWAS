#!/usr/bin/env bash
set -euo pipefail

IMAGE="hailgenetics/genetics:0.2.37"

# Input dir: 
HOST_STEP2_DIR="$HOME/gp08_data/gwaspy_filtered/pre_ldprune/post_ldprune/saige_step1/saige_step2"

# Output dir:
HOST_OUTDIR="$HOME/gp08_data/gwaspy_filtered/pre_ldprune/post_ldprune/saige_step1/saige_step2/meta_analysis"


# Container paths
CTR_STEP2_DIR="/in"
CTR_OUTDIR="/out"
CTR_TMPDIR="$CTR_OUTDIR/tmp"

# --- Prep ---
docker pull "$IMAGE" >/dev/null
mkdir -p "$HOST_OUTDIR"

# --- Run inside container ---
docker run --rm -it \
  -u "$(id -u)":"$(id -g)" \
  -v "$HOST_STEP2_DIR":"$CTR_STEP2_DIR" \
  -v "$HOST_OUTDIR":"$CTR_OUTDIR" \
  "$IMAGE" bash -lc '
    set -euo pipefail

    CTR_STEP2_DIR="'"$CTR_STEP2_DIR"'"
    CTR_OUTDIR="'"$CTR_OUTDIR"'"
    CTR_TMPDIR="'"$CTR_TMPDIR"'"

    echo "Container step2 dir: $CTR_STEP2_DIR"
    echo "Container out dir  : $CTR_OUTDIR"
    echo "Container tmp dir  : $CTR_TMPDIR"
    mkdir -p "$CTR_OUTDIR" "$CTR_TMPDIR"

    echo "Sanity check: first few entries under $CTR_STEP2_DIR"
    ls -lah "$CTR_STEP2_DIR" | head -n 30 || true
    echo

    # Ensure PLINK 1.9 (plink2 lacks --meta-analysis)
    if ! command -v plink >/dev/null 2>&1; then
      echo "ERROR: plink (1.9) not found in container." >&2
      exit 1
    fi
    if plink --version 2>/dev/null | head -1 | grep -qi "plink 2"; then
      echo "ERROR: Detected plink2; need PLINK 1.9 for --meta-analysis." >&2
      exit 1
    fi

    # SAIGE (binary) yields log-odds betas → use logscale
    declare -A TRAIT_ARG=(
      ["assist_khat"]="logscale"
    )

    SITES=("AAU" "KEMRI" "Moi" "UCT" "Uganda")

    # Print missing files nicely if any
    require_complete_site () {
      local site="$1"
      local pheno="$2"
      local missing=()
      for chr in {1..22}; do
        local f="$CTR_STEP2_DIR/${site}_chr${chr}_${pheno}_saige_step2.txt"
        [[ -s "$f" ]] || missing+=("$f")
      done
      if (( ${#missing[@]} )); then
        echo "ERROR: Missing files for site=${site}, phenotype=${pheno}:"
        printf "  %s\n" "${missing[@]}"
        return 1
      fi
    }

    # Concatenate chr1..22 into one per-site file (keep header once)
    concat_site_files () {
      local site="$1"
      local pheno="$2"
      local outfile="$CTR_TMPDIR/${site}_allchr_${pheno}_saige_step2.txt"
      : > "$outfile"

      local first=1
      for chr in {1..22}; do
        local f="$CTR_STEP2_DIR/${site}_chr${chr}_${pheno}_saige_step2.txt"
        if (( first )); then
          cat "$f" >> "$outfile"            # keep header on first file
          first=0
        else
          tail -n +2 "$f" >> "$outfile"     # drop header for subsequent files
        fi
      done
      echo "$outfile"
    }

    # ---- Run meta-analysis for the phenotype(s) you care about ----
    for PHENO in "assist_khat"; do
      echo
      echo "================  META-ANALYSIS for phenotype: ${PHENO}  ==============="

      # Discover sites with ≥1 matching file, using robust nullglob matching
      AVAILABLE=()
      shopt -s nullglob
      for site in "${SITES[@]}"; do
        matches=( "$CTR_STEP2_DIR/${site}_chr"*"_${PHENO}_saige_step2.txt" )
        if (( ${#matches[@]} )); then
          AVAILABLE+=("$site")
        fi
      done
      shopt -u nullglob

      if (( ${#AVAILABLE[@]} < 2 )); then
        echo "ERROR: Need ≥2 sites; found: ${AVAILABLE[*]:-none} for ${PHENO}" >&2
        exit 1
      fi
      echo "Candidate sites (found ≥1 chr file): ${AVAILABLE[*]}"

      # Ensure all chromosomes exist per included site
      for site in "${AVAILABLE[@]}"; do
        require_complete_site "$site" "$PHENO" || {
          echo "Stopping: missing chromosome(s) for site=${site}, phenotype=${PHENO}" >&2
          exit 2
        }
      done
      echo "All candidate sites have complete chr1..22 for ${PHENO}."

      # Build per-site meta input files
      META_INPUTS=()
      for site in "${AVAILABLE[@]}"; do
        out=$(concat_site_files "$site" "$PHENO")
        echo "Concatenated -> $out"
        META_INPUTS+=("$out")
      done

      TRAIT="${TRAIT_ARG[$PHENO]}"
      BASENAME="saige_all_sites_${PHENO}"

      echo "Running PLINK meta-analysis (${TRAIT})..."
      plink \
        --meta-analysis "${META_INPUTS[@]}" + "$TRAIT" \
        --meta-analysis-a1-field Allele1 \
        --meta-analysis-a2-field Allele2 \
        --meta-analysis-snp-field MarkerID \
        --meta-analysis-chr-field CHR \
        --meta-analysis-bp-field POS \
        --out "$CTR_OUTDIR/$BASENAME"

      echo "Done: $CTR_OUTDIR/$BASENAME.*"
    done
  '

echo
echo "Meta-analysis outputs are in: $HOST_OUTDIR"

