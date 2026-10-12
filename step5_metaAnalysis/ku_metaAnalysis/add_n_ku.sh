#!/usr/bin/env python3

# calculates effective sample size but labels it N
from pathlib import Path
from collections import defaultdict
import shutil

site_files = [
    "gp08_AAU_assist_khat_saige_assoc_allchr.txt",
    "gp08_KEMRI_assist_khat_saige_assoc_allchr.txt",
    "gp08_Moi_assist_khat_saige_assoc_allchr.txt",
    "gp08_UCT_assist_khat_saige_assoc_allchr.txt",
    "gp08_Uganda_assist_khat_saige_assoc_allchr.txt",
]

# Input and output files
input_file = Path("ku_se_meta.txt")
output_file = Path("ku_gp08_meta_cleaned_with_N.txt")
backup_file = Path("ku_gp08_meta_cleaned_with_N.txt.backup")
temporary_file = Path("ku_gp08_meta_cleaned_with_N.txt.tmp")

# Check that all required files exist.
for filename in site_files:
    if not Path(filename).is_file():
        raise FileNotFoundError(f"Missing SAIGE input file: {filename}")

if not input_file.is_file():
    raise FileNotFoundError(f"Missing meta-analysis input: {input_file}")

if output_file.exists():
    raise FileExistsError(
        f"Output already exists: {output_file}. "
        "Rename or remove it before running to avoid overwriting it."
    )

# marker -> [summed N_case, summed N_ctrl, summed effective N]
variant_counts = defaultdict(lambda: [0, 0, 0.0])

for site_file in site_files:
    site_variants = {}

    with open(site_file) as infile:
        header = infile.readline().split()

        if not header:
            raise ValueError(f"{site_file}: missing header")

        if "MarkerID" not in header or "N_case" not in header:
            raise ValueError(
                f"{site_file}: expected MarkerID and N_case columns"
            )

        marker_index = header.index("MarkerID")
        n_case_index = header.index("N_case")

        # N_ctrl is assumed to be immediately after N_case in data rows.
        n_ctrl_index = n_case_index + 1

        for line_number, line in enumerate(infile, start=2):
            fields = line.split()

            if not fields:
                continue

            required_index = max(
                marker_index, n_case_index, n_ctrl_index
            )

            if len(fields) <= required_index:
                raise ValueError(
                    f"{site_file}, line {line_number}: "
                    f"not enough columns ({len(fields)})"
                )

            marker = fields[marker_index]

            try:
                n_case_float = float(fields[n_case_index])
                n_ctrl_float = float(fields[n_ctrl_index])
            except ValueError:
                raise ValueError(
                    f"{site_file}, line {line_number}: "
                    "invalid N_case or N_ctrl"
                )

            if (
                not n_case_float.is_integer()
                or not n_ctrl_float.is_integer()
                or n_case_float < 0
                or n_ctrl_float < 0
            ):
                raise ValueError(
                    f"{site_file}, line {line_number}: "
                    "sample counts must be non-negative integers"
                )

            n_case = int(n_case_float)
            n_ctrl = int(n_ctrl_float)

            # Store once per site so duplicate rows cannot double-count.
            if marker in site_variants:
                if site_variants[marker] != (n_case, n_ctrl):
                    raise ValueError(
                        f"{site_file}: duplicate {marker} "
                        "has conflicting counts"
                    )
            else:
                site_variants[marker] = (n_case, n_ctrl)

    # Sum counts and effective sample size across contributing sites.
    for marker, (n_case, n_ctrl) in site_variants.items():
        effective_n = (
            4 * n_case * n_ctrl / (n_case + n_ctrl)
            if n_case + n_ctrl > 0
            else 0.0
        )

        variant_counts[marker][0] += n_case
        variant_counts[marker][1] += n_ctrl
        variant_counts[marker][2] += effective_n

    print(f"{site_file}: {len(site_variants):,} unique variants")

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

        # Add sample-size columns if they do not already exist.
        for column in ["N_case", "N_ctrl", "N"]:
            if column not in header:
                header.append(column)

        n_case_index = header.index("N_case")
        n_ctrl_index = header.index("N_ctrl")
        n_index = header.index("N")

        outfile.write("\t".join(header) + "\n")

        for line_number, line in enumerate(infile, start=2):
            fields = line.split()

            if not fields:
                continue

            original_length = len(fields)
            expected_original_columns = len(header) - sum(
                1 for column in ["N_case", "N_ctrl", "N"]
                if column not in fields[:0]
            )

            # Ensure all original fields are present.
            if original_length > len(header):
                raise ValueError(
                    f"{input_file}, line {line_number}: "
                    f"more fields than header ({original_length})"
                )

            fields.extend([""] * (len(header) - len(fields)))

            marker = fields[marker_index]

            if marker in variant_counts:
                n_case, n_ctrl, effective_n = variant_counts[marker]
                matched += 1
                fields[n_case_index] = str(n_case)
                fields[n_ctrl_index] = str(n_ctrl)
                fields[n_index] = (
                    f"{effective_n:.10f}".rstrip("0").rstrip(".")
                )
            else:
                unmatched += 1
                fields[n_case_index] = ""
                fields[n_ctrl_index] = ""
                fields[n_index] = ""

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
