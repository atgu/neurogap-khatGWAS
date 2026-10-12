#!/usr/bin/env bash
set -euo pipefail

PLINK2="/home/tboltz/apps/plink2"

KU_BFILE="$HOME/neurogap_passed_all_qc/NeuroGAP_passed_all_qc_no_multi_allele"
KUF_BFILE="$HOME/neurogap_passed_all_qc/NeuroGAP_passed_all_qc_no_UCT_no_multiallele"

GWASDIR="$HOME/sud/khat/a2_effect"
CLUMPDIR="$HOME/clump"

BASE_OUTDIR="$HOME/sud/khat/a2_effect"
LD_CALCULATED_DIR="$BASE_OUTDIR/ld_calculated"
LD_ANNOTATED_DIR="$BASE_OUTDIR/ld_annotated"

KU_GWAS="$GWASDIR/ku_gp08_meta_cleaned_with_N.txt"
KUF_GWAS="$GWASDIR/kuf_gp08_meta_cleaned_with_N_and_freq.tsv"

KU_CLUMPS="$CLUMPDIR/ku_gwaspy_clump.clumps"
KUF_CLUMPS="$CLUMPDIR/kuf_gwaspy_clump.clumps"

WINDOW=500000
WINDOW_KB=500

mkdir -p \
    "$LD_CALCULATED_DIR" \
    "$LD_ANNOTATED_DIR"

check_file() {
    local file="$1"

    if [[ ! -s "$file" ]]; then
        echo "ERROR: missing or empty file: $file" >&2
        exit 1
    fi
}

check_bfile() {
    local prefix="$1"

    for extension in bed bim fam; do
        check_file "${prefix}.${extension}"
    done
}

check_file "$PLINK2"
check_bfile "$KU_BFILE"
check_bfile "$KUF_BFILE"

check_file "$KU_GWAS"
check_file "$KUF_GWAS"
check_file "$KU_CLUMPS"
check_file "$KUF_CLUMPS"

calculate_and_annotate_regions() {
    local gwas_file="$1"
    local clump_file="$2"
    local phenotype_prefix="$3"
    local bfile="$4"

    echo
    echo "============================================================"
    echo "Processing: $phenotype_prefix"
    echo "GWAS:       $gwas_file"
    echo "Clumps:     $clump_file"
    echo "LD ref:     $bfile"
    echo "============================================================"

    # Expected PLINK2 .clumps columns:
    # CHROM POS ID P ...
    tail -n +2 "$clump_file" |
    while read -r chr pos lead_snp p_value rest; do

        [[ -n "${chr:-}" ]] || continue
        [[ -n "${pos:-}" ]] || continue
        [[ -n "${lead_snp:-}" ]] || continue

        start=$((pos - WINDOW))
        end=$((pos + WINDOW))

        if [[ "$start" -lt 1 ]]; then
            start=1
        fi

        region_name="${phenotype_prefix}_${lead_snp}_chr${chr}_${start}_${end}_500kb"

        raw_region_file="$LD_CALCULATED_DIR/${region_name}.region.tsv"
        ld_prefix="$LD_CALCULATED_DIR/${region_name}"
        annotated_file="$LD_ANNOTATED_DIR/${region_name}.ld_annotated.tsv"

        echo
        echo "Lead SNP: $lead_snp"
        echo "Region:   chr${chr}:${start}-${end}"

        # ----------------------------------------------------
        # 1. Extract regional GWAS variants
        #
        # Expected GWAS columns:
        #   column 1 = MarkerName
        #   column 2 = CHR
        #   column 3 = BP
        # ----------------------------------------------------

        awk \
            -v target_chr="$chr" \
            -v region_start="$start" \
            -v region_end="$end" \
            '
            BEGIN {
                FS = OFS = "\t"
            }

            NR == 1 {
                print
                next
            }

            $2 == target_chr &&
            $3 >= region_start &&
            $3 <= region_end {
                print
            }
            ' \
            "$gwas_file" > "$raw_region_file"

        n_variants=$((
            $(wc -l < "$raw_region_file") - 1
        ))

        echo "Regional GWAS variants: $n_variants"

        if [[ "$n_variants" -lt 1 ]]; then
            echo "WARNING: no GWAS variants found; skipping $lead_snp" >&2
            continue
        fi

        # ----------------------------------------------------
        # 2. Confirm that the lead SNP is in the LD reference
        # ----------------------------------------------------

        if ! awk -v snp="$lead_snp" '$2 == snp { found=1; exit }
                                     END { exit !found }' \
            "${bfile}.bim"
        then
            echo "WARNING: lead SNP is absent from the NeuroGAP reference: $lead_snp" >&2

            awk -v lead="$lead_snp" '
                BEGIN {
                    FS = OFS = "\t"
                }

                NR == 1 {
                    print $0, "LEAD_SNP", "LD_R2"
                    next
                }

                {
                    print $0, lead, "NA"
                }
                ' \
                "$raw_region_file" > "$annotated_file"

            continue
        fi

        # ----------------------------------------------------
        # 3. Calculate unphased r² relative to the lead SNP
        # ----------------------------------------------------

        "$PLINK2" \
            --bfile "$bfile" \
            --chr "$chr" \
            --from-bp "$start" \
            --to-bp "$end" \
            --r2-unphased \
            --ld-snp "$lead_snp" \
            --ld-window-kb "$WINDOW_KB" \
            --ld-window 999999 \
            --ld-window-r2 0 \
            --out "$ld_prefix"

        # PLINK2 normally writes this as .vcor.
        # This also allows for version-specific .vcor extensions.
        ld_file=""

        for candidate in \
            "${ld_prefix}.vcor" \
            "${ld_prefix}.vcor1" \
            "${ld_prefix}.vcor2"
        do
            if [[ -s "$candidate" ]]; then
                ld_file="$candidate"
                break
            fi
        done

        if [[ -z "$ld_file" ]]; then
            ld_file=$(find "$LD_CALCULATED_DIR" \
                -maxdepth 1 \
                -type f \
                -name "$(basename "$ld_prefix").vcor*" \
                ! -name "*.zst" \
                -print \
                -quit)
        fi

        if [[ -z "${ld_file:-}" || ! -s "$ld_file" ]]; then
            echo "WARNING: no readable LD file created for $lead_snp" >&2

            awk -v lead="$lead_snp" '
                BEGIN {
                    FS = OFS = "\t"
                }

                NR == 1 {
                    print $0, "LEAD_SNP", "LD_R2"
                    next
                }

                {
                    r2 = ($1 == lead) ? 1 : "NA"
                    print $0, lead, r2
                }
                ' \
                "$raw_region_file" > "$annotated_file"

            continue
        fi

        echo "LD results: $ld_file"

        # ----------------------------------------------------
        # 4. Annotate regional GWAS data with lead-SNP r²
        #
        # This reads the PLINK2 header dynamically, allowing
        # UNPHASED_R2 or R2 as the r² column name.
        # ----------------------------------------------------

        awk \
            -v lead="$lead_snp" \
            '
            BEGIN {
                FS = OFS = "\t"
            }

            # Read the LD file first.
            NR == FNR {
                if (FNR == 1) {
                    for (i = 1; i <= NF; i++) {
                        header = $i
                        sub(/^#/, "", header)

                        if (header == "ID_A") {
                            id_a_col = i
                        }

                        if (header == "ID_B") {
                            id_b_col = i
                        }

                        if (header == "R2" ||
                            header == "UNPHASED_R2" ||
                            header ~ /R2$/) {
                            r2_col = i
                        }
                    }

                    next
                }

                if (!id_a_col || !id_b_col || !r2_col) {
                    next
                }

                id_a = $id_a_col
                id_b = $id_b_col
                r2 = $r2_col

                if (id_a == lead) {
                    ld_r2[id_b] = r2
                } else if (id_b == lead) {
                    ld_r2[id_a] = r2
                }

                next
            }

            # Read the regional GWAS file second.
            NR != FNR {
                if (FNR == 1) {
                    print $0, "LEAD_SNP", "LD_R2"
                    next
                }

                snp = $1

                if (snp == lead) {
                    r2 = 1
                } else if (snp in ld_r2) {
                    r2 = ld_r2[snp]
                } else {
                    r2 = "NA"
                }

                print $0, lead, r2
            }
            ' \
            "$ld_file" \
            "$raw_region_file" \
            > "$annotated_file"

        annotated_n=$((
            $(wc -l < "$annotated_file") - 1
        ))

        ld_matched_n=$(
            awk '
                BEGIN {
                    FS = "\t"
                    count = 0
                }

                NR > 1 && $NF != "NA" {
                    count++
                }

                END {
                    print count
                }
                ' \
                "$annotated_file"
        )

        echo "Annotated variants:       $annotated_n"
        echo "Variants with LD r²:      $ld_matched_n"
        echo "Annotated output:         $annotated_file"
    done
}

# KU
calculate_and_annotate_regions \
    "$KU_GWAS" \
    "$KU_CLUMPS" \
    "ku" \
    "$KU_BFILE"

# KUF
calculate_and_annotate_regions \
    "$KUF_GWAS" \
    "$KUF_CLUMPS" \
    "kuf" \
    "$KUF_BFILE"

echo
echo "============================================================"
echo "Finished"
echo "============================================================"
echo "Calculated LD files:"
echo "  $LD_CALCULATED_DIR"
echo
echo "LD-annotated GWAS files:"
echo "  $LD_ANNOTATED_DIR"
