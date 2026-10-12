#!/usr/bin/env bash
set -euo pipefail

# Change this to the name/path of your GWAS results file.
INPUT_FILE="ku_gp08_meta_cleaned_with_N.txt"

# Output file.
OUTPUT_FILE="selected_variants_with_eaf_se_ku.txt"

if [[ ! -s "$INPUT_FILE" ]]; then
    echo "ERROR: Input file is missing or empty:"
    echo "  $INPUT_FILE"
    exit 1
fi

marker_file=$(mktemp)
trap 'rm -f "$marker_file"' EXIT

cat > "$marker_file" <<'MARKERS'
rs114077623
rs138362315
rs140323147
rs151311894
rs189085944
rs189679719
rs192966269
rs193228507
rs200034328
rs28478222
rs34334893
rs535972653
rs544591434
rs564111892
rs568695780
rs60863078
rs74088774
rs79163444
rs937215889
MARKERS

# Allow either a plain-text or gzip-compressed GWAS file.
if [[ "$INPUT_FILE" == *.gz ]]; then
    input_command=(gzip -cd "$INPUT_FILE")
else
    input_command=(cat "$INPUT_FILE")
fi

awk '
BEGIN {
    FS = "[[:space:]]+"
    OFS = "\t"
}

# Read the requested marker list.
NR == FNR {
    wanted[$1] = 1
    marker_order[++number_markers] = $1
    next
}

# Locate columns from the GWAS header.
FNR == 1 {
    for (i = 1; i <= NF; i++) {
        column[$i] = i
    }

    required[1] = "MarkerName"
    required[2] = "CHR"
    required[3] = "BP"
    required[4] = "Allele1"
    required[5] = "Allele2"
    required[6] = "Freq1"
    required[7] = "Effect"
    required[8] = "StdErr"
    required[9] = "P-value"
    required[10] = "Direction"
    required[11] = "N_case"
    required[12] = "N_ctrl"
    required[13] = "N"

    for (i = 1; i <= 13; i++) {
        if (!(required[i] in column)) {
            print "ERROR: Missing required column: " required[i] \
                > "/dev/stderr"
            exit 1
        }
    }

    marker_col = column["MarkerName"]

    print \
        "Chromosome", \
        "BP", \
        "MarkerName", \
        "Effect_Allele", \
        "Non_effect_Allele", \
        "EAF", \
        "Effect", \
        "SE", \
        "P_value", \
        "Direction", \
        "N_case", \
        "N_ctrl", \
        "N"

    next
}

# Save requested variants.
($marker_col in wanted) {
    marker = $marker_col

    result[marker] = \
        $(column["CHR"]) OFS \
        $(column["BP"]) OFS \
        $(column["MarkerName"]) OFS \
        $(column["Allele1"]) OFS \
        $(column["Allele2"]) OFS \
        $(column["Freq1"]) OFS \
        $(column["Effect"]) OFS \
        $(column["StdErr"]) OFS \
        $(column["P-value"]) OFS \
        $(column["Direction"]) OFS \
        $(column["N_case"]) OFS \
        $(column["N_ctrl"]) OFS \
        $(column["N"])

    found[marker] = 1
}

END {
    for (i = 1; i <= number_markers; i++) {
        marker = marker_order[i]

        if (marker in found) {
            print result[marker]
        } else {
            print "WARNING: Marker not found: " marker > "/dev/stderr"
        }
    }
}
' "$marker_file" <("${input_command[@]}") > "$OUTPUT_FILE"

echo "Finished."
echo "Output: $OUTPUT_FILE"
echo
echo "Variants written:"
awk 'END {print NR - 1}' "$OUTPUT_FILE"
echo
column -t -s $'\t' "$OUTPUT_FILE"
