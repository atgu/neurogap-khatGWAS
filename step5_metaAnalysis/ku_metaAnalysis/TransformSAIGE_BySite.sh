#!/bin/bash
set -euo pipefail

# Add total analyzed sample size N = N_case + N_ctrl to KU SAIGE outputs.
# Source files are not overwritten. Each output is prefixed with metal_ready_.

# input file names
files=(
  "gp08_KEMRI_assist_khat_saige_assoc_allchr.txt"
  "gp08_AAU_assist_khat_saige_assoc_allchr.txt"
  "gp08_Moi_assist_khat_saige_assoc_allchr.txt"
  "gp08_UCT_assist_khat_saige_assoc_allchr.txt"
  "gp08_Uganda_assist_khat_saige_assoc_allchr.txt"
)

for input_file in "${files[@]}"; do
  if [[ ! -f "${input_file}" ]]; then
    echo "Skipping missing file: ${input_file}" >&2
    continue
  fi

  output_file="metal_ready_${input_file}"

  INPUT_FILE="${input_file}" OUTPUT_FILE="${output_file}" python3 - <<'PY'
import csv
import math
import os

input_file = os.environ["INPUT_FILE"]
output_file = os.environ["OUTPUT_FILE"]

required = {
    "MarkerID",
    "Allele1",
    "Allele2",
    "AF_Allele2",
    "BETA",
    "SE",
    "p.value",
    "N_case",
    "N_ctrl",
}

with open(input_file, "r", newline="") as source:
    reader = csv.DictReader(source, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"{input_file}: missing header")

    missing = sorted(required.difference(reader.fieldnames))
    if missing:
        raise SystemExit(
            f"{input_file}: missing required column(s): {', '.join(missing)}"
        )

    fieldnames = list(reader.fieldnames)
    if "N" not in fieldnames:
        fieldnames.append("N")

    row_count = 0
    n_values = set()

    with open(output_file, "w", newline="") as destination:
        writer = csv.DictWriter(
            destination,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row_number, row in enumerate(reader, start=2):
            try:
                n_case = float(row["N_case"])
                n_ctrl = float(row["N_ctrl"])
            except (TypeError, ValueError):
                raise SystemExit(
                    f"{input_file}:{row_number}: invalid N_case or N_ctrl"
                )

            if not math.isfinite(n_case) or not math.isfinite(n_ctrl):
                raise SystemExit(
                    f"{input_file}:{row_number}: non-finite N_case or N_ctrl"
                )

            total_n = n_case + n_ctrl
            row["N"] = str(int(total_n)) if total_n.is_integer() else str(total_n)
            n_values.add(row["N"])
            writer.writerow(row)
            row_count += 1

print(
    f"Created {output_file}: {row_count:,} variants; "
    f"observed N value(s): {', '.join(sorted(n_values))}"
)
PY
done


