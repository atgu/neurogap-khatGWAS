#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import shutil

# Site BIM files and their sample sizes.
site_files = {
    "KEMRI_passed_all_qc.bim": 2587,
    "AAU_passed_all_qc.bim": 10142,
    "Moi_passed_all_qc.bim": 3724,
    "Uganda_passed_all_qc.bim": 9416,
}

# Input and output files.
input_file = Path("kuf_se_meta.txt")
output_file = Path("kuf_gp08_meta_cleaned_with_N.tsv")
backup_file = Path("kuf_gp08_meta_cleaned_with_N.tsv.backup")
temporary_file = Path("kuf_gp08_meta_cleaned_with_N.tsv.tmp")

# Check that all required files exist.
for bim_file in site_files:
    if not Path(bim_file).is_file():
        raise FileNotFoundError(f"Missing BIM file: {bim_file}")

if not input_file.is_file():
    raise FileNotFoundError(f"Missing meta-analysis input: {input_file}")

if output_file.exists():
    raise FileExistsError(
        f"Output already exists: {output_file}. "
        "Rename or remove it before running to avoid overwriting it."
    )

# Variant -> summed sample size across sites containing that variant.
variant_n = defaultdict(int)

# Calculate N for each variant based on which site BIM files contain it.
for bim_file, sample_size in site_files.items():
    variants = set()

    with open(bim_file) as infile:
        for line in infile:
            fields = line.split()

            if len(fields) >= 2:
                variants.add(fields[1])

    for variant in variants:
        variant_n[variant] += sample_size

    print(
        f"{bim_file}: {len(variants):,} unique variants, "
        f"N = {sample_size}"
    )

# Write annotated output without modifying the original input.
matched = 0
unmatched = 0

try:
    with input_file.open() as infile, temporary_file.open("w") as outfile:
        header = infile.readline().split()

        if not header:
            raise ValueError(f"{input_file}: missing header")

        if "MarkerName" not in header:
            raise ValueError(
                f"{input_file}: expected MarkerName column"
            )

        marker_index = header.index("MarkerName")

        # Add N if it does not already exist.
        if "N" not in header:
            header.append("N")

        n_index = header.index("N")

        outfile.write("\t".join(header) + "\n")

        for line_number, line in enumerate(infile, start=2):
            fields = line.split()

            if not fields:
                continue

            if len(fields) > len(header):
                raise ValueError(
                    f"{input_file}, line {line_number}: "
                    f"more fields than header ({len(fields)})"
                )

            fields.extend([""] * (len(header) - len(fields)))

            marker = fields[marker_index]
            new_n = variant_n.get(marker, 0)

            # Replace N with the total sample size across contributing sites.
            fields[n_index] = str(new_n)

            if new_n > 0:
                matched += 1
            else:
                unmatched += 1

            outfile.write("\t".join(fields) + "\n")

    # Back up an existing output before replacing it, if applicable.
    if output_file.exists():
        shutil.copy2(output_file, backup_file)

    temporary_file.replace(output_file)

except Exception:
    if temporary_file.exists():
        temporary_file.unlink()
    raise

print(f"\nInput:  {input_file}")
print(f"Output: {output_file}")
print(f"Matched variants:   {matched:,}")
print(f"Unmatched variants: {unmatched:,}")
