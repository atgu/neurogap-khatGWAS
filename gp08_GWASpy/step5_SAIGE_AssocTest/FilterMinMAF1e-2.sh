#!/usr/bin/env bash
set -euo pipefail

IN_DIR="${1:-$PWD}"
OUT_DIR="${2:-$IN_DIR/minMAF1e-2}"
THRESH="${3:-0.01}"   # want to match SAIGE outputs to POLMM (minMAF 1e-2)

mkdir -p "$OUT_DIR"

# find files with given name pattern
mapfile -d '' FILES < <(find "$IN_DIR" -maxdepth 1 -type f -name '*_chr*_assist_khat_saige_step2.txt' -print0)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No SAIGE step2 .txt files found in: $IN_DIR"
  exit 1
fi

for f in "${FILES[@]}"; do
  bn="$(basename "$f")"
  echo "Processing $bn → $OUT_DIR/$bn"

  awk -v FS='[[:space:]]+' -v OFS='\t' -v thresh="$THRESH" '
    NR==1{
      for(i=1;i<=NF;i++) h[$i]=i
      # association/stat columns to blank when AF_Allele2 < thresh
      stats["BETA"]; stats["SE"]; stats["Tstat"]; stats["var"]; stats["p.value"]; stats["p.value.NA"]; stats["Is.SPA"];
      # add these if you want them blanked too:
      # stats["AF_case"]; stats["AF_ctrl"]; stats["N_case"]; stats["N_ctrl"];
      # stats["N_case_hom"]; stats["N_case_het"]; stats["N_ctrl_hom"]; stats["N_ctrl_het"];
      print; next
    }
    {
      af = $(h["AF_Allele2"])+0
      if (af < thresh) {
        for (s in stats) if (s in h) $(h[s])="NA"
      }
      print
    }' "$f" > "$OUT_DIR/$bn"
done

echo "Done. Output in: $OUT_DIR"
