#!/usr/bin/env bash
set -euo pipefail

PLINK="$HOME/apps/plink1.9/plink"

BFILE_DIR="$HOME/neurogap_passed_all_qc"
WORK_DIR="$BFILE_DIR/no_uct"

INITIAL="$WORK_DIR/initial_merge_no_uct"
FINAL="$WORK_DIR/NeuroGAP_passed_all_qc_no_UCT_no_multiallele"

SITES=(AAU KEMRI Moi Uganda)

mkdir -p "$WORK_DIR"

[[ -x "$PLINK" ]] || {
    echo "ERROR: PLINK 1.9 not found: $PLINK" >&2
    exit 1
}

# Confirm the original files exist.
for site in "${SITES[@]}"; do
    for ext in bed bim fam; do
        file="$BFILE_DIR/${site}_passed_all_qc.${ext}"

        [[ -s "$file" ]] || {
            echo "ERROR: missing original input: $file" >&2
            exit 1
        }
    done
done

# Remove only outputs from previous failed attempts.
rm -f \
    "$WORK_DIR"/*_no_merge_conflicts.{bed,bim,fam,log} \
    "$WORK_DIR"/initial_merge_no_uct* \
    "$WORK_DIR"/NeuroGAP_passed_all_qc_no_UCT_no_multiallele* \
    "$WORK_DIR"/merge_no_uct_cleaned.txt \
    "$WORK_DIR"/initial_merge_list.txt \
    "$WORK_DIR"/multiallelic_or_mismatching_variants.txt

# Create the initial merge list.
cat > "$WORK_DIR/initial_merge_list.txt" <<EOF
$BFILE_DIR/KEMRI_passed_all_qc
$BFILE_DIR/Moi_passed_all_qc
$BFILE_DIR/Uganda_passed_all_qc
EOF

echo "Detecting incompatible variants..."

set +e

"$PLINK" \
    --bfile "$BFILE_DIR/AAU_passed_all_qc" \
    --merge-list "$WORK_DIR/initial_merge_list.txt" \
    --make-bed \
    --out "$INITIAL"

initial_status=$?

set -e

MISSNP="${INITIAL}-merge.missnp"

if [[ "$initial_status" -eq 0 ]]; then
    echo "The original files merged without conflicts."

    mv "${INITIAL}.bed" "${FINAL}.bed"
    mv "${INITIAL}.bim" "${FINAL}.bim"
    mv "${INITIAL}.fam" "${FINAL}.fam"

else
    [[ -s "$MISSNP" ]] || {
        echo "ERROR: merge failed without creating: $MISSNP" >&2
        exit 1
    }

    sort -u "$MISSNP" \
        > "$WORK_DIR/multiallelic_or_mismatching_variants.txt"

    echo "Variants requiring exclusion:"
    wc -l "$WORK_DIR/multiallelic_or_mismatching_variants.txt"

    # Remove failed initial-merge genotype files.
    rm -f \
        "${INITIAL}.bed" \
        "${INITIAL}.bim" \
        "${INITIAL}.fam" \
        "${INITIAL}-merge.bed" \
        "${INITIAL}-merge.bim" \
        "${INITIAL}-merge.fam"

    # Create each cleaned cohort, verify it, and then remove its
    # original local BED to recover space. Original BIM/FAM stay.
    for site in "${SITES[@]}"; do
        CLEANED="$WORK_DIR/${site}_no_merge_conflicts"
        ORIGINAL="$BFILE_DIR/${site}_passed_all_qc"

        echo
        echo "Cleaning $site..."

        "$PLINK" \
            --bfile "$ORIGINAL" \
            --exclude "$WORK_DIR/multiallelic_or_mismatching_variants.txt" \
            --make-bed \
            --out "$CLEANED"

        for ext in bed bim fam; do
            [[ -s "${CLEANED}.${ext}" ]] || {
                echo "ERROR: missing cleaned file: ${CLEANED}.${ext}" >&2
                exit 1
            }
        done

        echo "Cleaned $site successfully."

        rm -f "${ORIGINAL}.bed"

        echo "Deleted original local BED: ${ORIGINAL}.bed"
        df -h "$BFILE_DIR" | tail -n 1
    done

    cat > "$WORK_DIR/merge_no_uct_cleaned.txt" <<EOF
$WORK_DIR/KEMRI_no_merge_conflicts
$WORK_DIR/Moi_no_merge_conflicts
$WORK_DIR/Uganda_no_merge_conflicts
EOF

    echo
    echo "Starting final merge..."

    "$PLINK" \
        --bfile "$WORK_DIR/AAU_no_merge_conflicts" \
        --merge-list "$WORK_DIR/merge_no_uct_cleaned.txt" \
        --make-bed \
        --out "$FINAL"
fi

# Validate final files.
for ext in bed bim fam; do
    [[ -s "${FINAL}.${ext}" ]] || {
        echo "ERROR: missing final file: ${FINAL}.${ext}" >&2
        exit 1
    }
done

variants=$(wc -l < "${FINAL}.bim")
samples=$(wc -l < "${FINAL}.fam")

echo
echo "Final variants: $variants"
echo "Final samples:  $samples"

[[ "$variants" -eq 31779760 ]] || {
    echo "ERROR: expected 31,779,760 variants." >&2
    exit 1
}

[[ "$samples" -eq 25873 ]] || {
    echo "ERROR: expected 25,873 samples." >&2
    exit 1
}

echo "Final merged dataset validated."

# Remove cleaned cohort copies only after the final files are valid.
for site in "${SITES[@]}"; do
    rm -f \
        "$WORK_DIR/${site}_no_merge_conflicts.bed" \
        "$WORK_DIR/${site}_no_merge_conflicts.bim" \
        "$WORK_DIR/${site}_no_merge_conflicts.fam" \
        "$WORK_DIR/${site}_no_merge_conflicts.log"
done

rm -f \
    "$WORK_DIR/initial_merge_list.txt" \
    "$WORK_DIR/merge_no_uct_cleaned.txt" \
    "$WORK_DIR/multiallelic_or_mismatching_variants.txt"

echo
echo "Complete:"
ls -lh "${FINAL}".{bed,bim,fam}

df -h "$BFILE_DIR"
