#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

GENO_DIR="$HOME/neurogap_passed_all_qc/minmaf01"
ASSOC_DIR="$HOME/generic-metal"
KUF_FILE="$HOME/sud/khat/a2_effect/kuf_gp08_meta_cleaned_with_N.tsv"
OUTPUT_FILE="$HOME/sud/khat/a2_effect/kuf_gp08_meta_cleaned_with_N_and_freq.tsv"

sites=(AAU KEMRI Moi Uganda)

WORK_DIR="$(mktemp -d "${GENO_DIR}/kuf_eaf_awk.XXXXXX")"

cleanup() {
    rm -rf -- "$WORK_DIR"
}

trap cleanup EXIT

SORT_OPTIONS=(
    --temporary-directory="$WORK_DIR"
    --buffer-size=25%
)

# -------------------------
# Validate inputs
# -------------------------

if [[ ! -f "$KUF_FILE" ]]; then
    echo "ERROR: Missing KUF file: $KUF_FILE" >&2
    exit 1
fi

for site in "${sites[@]}"; do
    association_file="${ASSOC_DIR}/A2_oriented_annotated_assoc_${site}_studysite_covar_all.tsv"
    frequency_file="${GENO_DIR}/${site}_passed_all_qc_maf01.afreq"

    if [[ ! -f "$association_file" ]]; then
        echo "ERROR: Missing association file: $association_file" >&2
        exit 1
    fi

    if [[ ! -f "$frequency_file" ]]; then
        echo "ERROR: Missing frequency file: $frequency_file" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$OUTPUT_FILE")"

# -------------------------
# Extract final METAL effect alleles
# -------------------------

# Output:
# MarkerName  METAL_Allele1
#
# METAL Allele1 is the allele whose frequency we want.

META_TARGETS="${WORK_DIR}/meta_targets.sorted.tsv"

echo "Preparing final METAL effect-allele lookup"

awk -F'\t' '
    BEGIN {
        OFS = "\t"
    }

    NR == 1 {
        for (i = 1; i <= NF; i++) {
            header = $i
            gsub(/\r|"/, "", header)

            if (header == "MarkerName") marker_col = i
            if (header == "Allele1") allele_col = i
        }

        if (!marker_col || !allele_col) {
            print "ERROR: Missing MarkerName or Allele1 in " FILENAME > "/dev/stderr"
            exit 2
        }

        next
    }

    {
        marker = tolower($marker_col)
        effect_allele = toupper($allele_col)

        gsub(/"/, "", marker)
        gsub(/"/, "", effect_allele)
        gsub(/\r/, "", effect_allele)

        if (marker != "" && effect_allele != "") {
            print marker, effect_allele
        }
    }
' "$KUF_FILE" |
    sort "${SORT_OPTIONS[@]}" \
        -t $'\t' \
        -k1,1 \
        -k2,2 \
        -u \
        > "$META_TARGETS"

echo "Final METAL variants: $(wc -l < "$META_TARGETS")"

# -------------------------
# Calculate site frequencies
# -------------------------

ALL_CONTRIBUTIONS="${WORK_DIR}/all_site_contributions.tsv"
: > "$ALL_CONTRIBUTIONS"

for site in "${sites[@]}"; do
    echo
    echo "Processing $site"

    ASSOC_FILE="${ASSOC_DIR}/A2_oriented_annotated_assoc_${site}_studysite_covar_all.tsv"
    FREQ_FILE="${GENO_DIR}/${site}_passed_all_qc_maf01.afreq"

    ASSOC_IDS="${WORK_DIR}/${site}_assoc_ids.sorted.tsv"
    SITE_TARGETS="${WORK_DIR}/${site}_meta_targets.tsv"
    FREQ_SORTED="${WORK_DIR}/${site}_freq.sorted.tsv"
    MATCHED_RAW="${WORK_DIR}/${site}_matched.raw.tsv"
    MATCHED_SORTED="${WORK_DIR}/${site}_matched.sorted.tsv"
    SITE_CONTRIBUTIONS="${WORK_DIR}/${site}_contributions.tsv"

    # Association files are used only to establish that a
    # variant contributed to this site's KUF GWAS.
    echo "  Extracting contributing association variants"

    awk -F'\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                header = $i
                gsub(/\r|"/, "", header)

                if (header == "SNPID") snp_col = i
            }

            if (!snp_col) {
                print "ERROR: Missing SNPID in " FILENAME > "/dev/stderr"
                exit 2
            }

            next
        }

        {
            snp = tolower($snp_col)
            gsub(/"|\r/, "", snp)

            if (snp != "") print snp
        }
    ' "$ASSOC_FILE" |
        sort "${SORT_OPTIONS[@]}" \
            -u \
            > "$ASSOC_IDS"

    # Retain final METAL variants that were present in this site.
    #
    # Output:
    # MarkerName  final_METAL_Allele1
    join \
        -t $'\t' \
        -1 1 \
        -2 1 \
        "$ASSOC_IDS" \
        "$META_TARGETS" \
        > "$SITE_TARGETS"

    contributing_count=$(wc -l < "$SITE_TARGETS")

    echo "    Final meta-analysis variants contributed by $site: $contributing_count"

    # Prepare PLINK frequencies.
    #
    # Output:
    # ID REF ALT ALT_FREQS OBS_CT
    echo "  Sorting PLINK frequency records"

    awk '
        BEGIN {
            OFS = "\t"
        }

        NR == 1 {
            for (i = 1; i <= NF; i++) {
                header = $i
                gsub(/\r/, "", header)

                if (header == "ID") id_col = i
                if (header == "REF") ref_col = i
                if (header == "ALT") alt_col = i
                if (header == "ALT_FREQS") freq_col = i
                if (header == "OBS_CT") obs_col = i
            }

            if (!id_col || !ref_col || !alt_col || !freq_col || !obs_col) {
                print "ERROR: Missing required frequency columns in " FILENAME > "/dev/stderr"
                exit 2
            }

            next
        }

        {
            id = tolower($id_col)
            ref = toupper($ref_col)
            alt = toupper($alt_col)
            alt_freq = $freq_col
            obs = $obs_col

            gsub(/\r/, "", obs)

            if (id != "" && ref != "" && alt != "" && alt_freq != "" && obs != "") {
                print id, ref, alt, alt_freq, obs
            }
        }
    ' "$FREQ_FILE" |
        sort "${SORT_OPTIONS[@]}" \
            -t $'\t' \
            -k1,1 \
            -k2,2 \
            -k3,3 \
            -u \
            > "$FREQ_SORTED"

    # Join final METAL effect allele to PLINK frequencies.
    #
    # SITE_TARGETS:
    # ID METAL_Allele1
    #
    # FREQ_SORTED:
    # ID REF ALT ALT_FREQS OBS_CT
    #
    # join output:
    # ID METAL_Allele1 REF ALT ALT_FREQS OBS_CT
    echo "  Calculating frequency of final METAL Allele1"

    join \
        -t $'\t' \
        -1 1 \
        -2 1 \
        "$SITE_TARGETS" \
        "$FREQ_SORTED" |
        awk -F'\t' -v site="$site" '
            BEGIN {
                OFS = "\t"
            }

            {
                id = $1
                effect_allele = toupper($2)
                ref = toupper($3)
                alt = toupper($4)
                alt_freq = $5 + 0
                obs = $6 + 0

                if (obs <= 0) next

                if (effect_allele == alt) {
                    effect_frequency = alt_freq
                } else if (effect_allele == ref) {
                    effect_frequency = 1 - alt_freq
                } else {
                    next
                }

                if (effect_frequency < 0 || effect_frequency > 1) next

                key = id "|" effect_allele

                print key, effect_frequency, obs, site
            }
        ' > "$MATCHED_RAW"

    # Resolve duplicate PLINK records by retaining the one with
    # the greatest number of observed alleles.
    sort "${SORT_OPTIONS[@]}" \
        -t $'\t' \
        -k1,1 \
        -k3,3nr \
        "$MATCHED_RAW" \
        > "$MATCHED_SORTED"

    awk -F'\t' '
        BEGIN {
            OFS = "\t"
            previous_key = ""
        }

        $1 != previous_key {
            print
            previous_key = $1
        }
    ' "$MATCHED_SORTED" > "$SITE_CONTRIBUTIONS"

    matched_count=$(wc -l < "$SITE_CONTRIBUTIONS")
    unmatched_count=$((contributing_count - matched_count))

    echo "    Variants with final-effect-allele frequency: $matched_count"
    echo "    Variants without a usable PLINK frequency: $unmatched_count"

    cat "$SITE_CONTRIBUTIONS" >> "$ALL_CONTRIBUTIONS"

    rm -f \
        "$ASSOC_IDS" \
        "$SITE_TARGETS" \
        "$FREQ_SORTED" \
        "$MATCHED_RAW" \
        "$MATCHED_SORTED" \
        "$SITE_CONTRIBUTIONS"
done

# -------------------------
# Pool frequencies
# -------------------------

echo
echo "Pooling final effect-allele frequencies across sites"

SORTED_CONTRIBUTIONS="${WORK_DIR}/all_site_contributions.sorted.tsv"
POOLED_FREQ="${WORK_DIR}/kuf_pooled_frequencies.tsv"

sort "${SORT_OPTIONS[@]}" \
    -t $'\t' \
    -k1,1 \
    -k4,4 \
    "$ALL_CONTRIBUTIONS" \
    > "$SORTED_CONTRIBUTIONS"

# Input:
# key effect_frequency OBS_CT site
#
# Output:
# key Freq1 Frequency_Site_Count Frequency_OBS_CT

awk -F'\t' '
    BEGIN {
        OFS = "\t"
        previous_key = ""
    }

    function write_result() {
        if (previous_key != "" && total_obs > 0) {
            print previous_key, numerator / total_obs, site_count, total_obs
        }
    }

    {
        key = $1
        frequency = $2 + 0
        obs = $3 + 0
        site = $4

        if (key != previous_key) {
            write_result()

            previous_key = key
            numerator = frequency * obs
            total_obs = obs
            site_count = 1
            previous_site = site
        } else {
            numerator += frequency * obs
            total_obs += obs

            if (site != previous_site) {
                site_count++
                previous_site = site
            }
        }
    }

    END {
        write_result()
    }
' "$SORTED_CONTRIBUTIONS" > "$POOLED_FREQ"

echo "Variants with pooled frequency: $(wc -l < "$POOLED_FREQ")"

# -------------------------
# Annotate final KUF file
# -------------------------

echo
echo "Annotating final KUF file"

META_SORTED="${WORK_DIR}/kuf_meta.sorted.tsv"
JOINED_FILE="${WORK_DIR}/kuf_meta.joined.tsv"

# Add a temporary MarkerName|Allele1 key to every final METAL row.
awk -F'\t' '
    BEGIN {
        OFS = "\t"
    }

    NR == 1 {
        for (i = 1; i <= NF; i++) {
            header = $i
            gsub(/\r|"/, "", header)

            if (header == "MarkerName") marker_col = i
            if (header == "Allele1") allele_col = i
        }

        if (!marker_col || !allele_col) {
            print "ERROR: Missing MarkerName or Allele1 in " FILENAME > "/dev/stderr"
            exit 2
        }

        next
    }

    {
        marker = tolower($marker_col)
        effect_allele = toupper($allele_col)

        gsub(/"/, "", marker)
        gsub(/"/, "", effect_allele)
        gsub(/\r/, "", $NF)

        key = marker "|" effect_allele

        printf "%s", key

        for (i = 1; i <= NF; i++) {
            printf "%s%s", OFS, $i
        }

        printf "\n"
    }
' "$KUF_FILE" |
    sort "${SORT_OPTIONS[@]}" \
        -t $'\t' \
        -k1,1 \
        > "$META_SORTED"

join \
    -t $'\t' \
    -a 1 \
    -e "NA" \
    -o auto \
    "$META_SORTED" \
    "$POOLED_FREQ" \
    > "$JOINED_FILE"

# Write the header.
awk -F'\t' '
    BEGIN {
        OFS = "\t"
    }

    NR == 1 {
        gsub(/\r/, "", $NF)

        print \
            $0,
            "Freq1",
            "Frequency_Site_Count",
            "Frequency_OBS_CT"
    }
' "$KUF_FILE" > "$OUTPUT_FILE"

# Remove the temporary composite key.
cut -f2- "$JOINED_FILE" >> "$OUTPUT_FILE"

# -------------------------
# Summary
# -------------------------

total_variants=$(( $(wc -l < "$OUTPUT_FILE") - 1 ))

annotated_variants=$(
    awk -F'\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "Freq1") freq_col = i
            }

            next
        }

        freq_col && $freq_col != "NA" {
            count++
        }

        END {
            print count + 0
        }
    ' "$OUTPUT_FILE"
)

missing_variants=$((total_variants - annotated_variants))

echo
echo "Done."
echo "Final KUF variants: $total_variants"
echo "Variants annotated with Freq1: $annotated_variants"
echo "Variants without Freq1: $missing_variants"
echo "Wrote: $OUTPUT_FILE"
