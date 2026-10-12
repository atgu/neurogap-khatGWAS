#!/usr/bin/env bash
set -euo pipefail

INPUT_FILE="kuf_gp08_meta_cleaned_with_N_and_freq.tsv"
OUTPUT_FILE="kuf_selected_variants_with_eaf_se.tsv"

if [[ ! -s "$INPUT_FILE" ]]; then
    echo "ERROR: Missing or empty input file: $INPUT_FILE" >&2
    exit 1
fi

marker_file=$(mktemp)
trap 'rm -f "$marker_file"' EXIT

cat > "$marker_file" <<'MARKERS'
rs111918073
rs11233220
rs138148405
rs370450693
rs7225138
rs944937859
MARKERS

awk '
BEGIN {
    FS = "[[:space:]]+"
    OFS = "\t"
}

# Read requested markers and preserve their order.
NR == FNR {
    wanted[$1] = 1
    marker_order[++n] = $1
    next
}

# Identify columns from the GWAS header.
FNR == 1 {
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }

    required[1]  = "MarkerName"
    required[2]  = "CHR"
    required[3]  = "BP"
    required[4]  = "Allele1"
    required[5]  = "Allele2"
    required[6]  = "Freq1"
    required[7]  = "Effect"
    required[8]  = "StdErr"
    required[9]  = "P-value"
    required[10] = "Direction"
    required[11] = "N"
    required[12] = "Frequency_Site_Count"
    required[13] = "Frequency_OBS_CT"

    for (i = 1; i <= 13; i++) {
        if (!(required[i] in col)) {
            print "ERROR: Missing required column: " required[i] \
                > "/dev/stderr"
            exit 1
        }
    }

    marker_col = col["MarkerName"]

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
        "N", \
        "Frequency_Site_Count", \
        "Frequency_OBS_CT"

    next
}

# Store matching variants.
($marker_col in wanted) {
    marker = $marker_col

    result[marker] = \
        $(col["CHR"]) OFS \
        $(col["BP"]) OFS \
        $(col["MarkerName"]) OFS \
        $(col["Allele1"]) OFS \
        $(col["Allele2"]) OFS \
        $(col["Freq1"]) OFS \
        $(col["Effect"]) OFS \
        $(col["StdErr"]) OFS \
        $(col["P-value"]) OFS \
        $(col["Direction"]) OFS \
        $(col["N"]) OFS \
        $(col["Frequency_Site_Count"]) OFS \
        $(col["Frequency_OBS_CT"])

    found[marker] = 1
}

END {
    for (i = 1; i <= n; i++) {
        marker = marker_order[i]

        if (marker in found) {
            print result[marker]
        } else {
            print "WARNING: Marker not found: " marker > "/dev/stderr"
        }
    }
}
' "$marker_file" "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Finished."
echo "Output: $OUTPUT_FILE"
echo
column -t -s $'\t' "$OUTPUT_FILE"
