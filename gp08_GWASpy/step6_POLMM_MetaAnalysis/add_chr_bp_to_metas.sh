#!/bin/bash

# Create lookup: RSID CHR BP
# using old NeuroGA file for this information since it all matches the same.
	# alternative is to make a script that matches the variant ID to the respective per-site plink file's rsID column and then important that respective information

awk 'BEGIN{OFS="\t"} {print $2, $1, $4}' NeuroGAP_GP08_no_dups.bim | sort -k1,1 > bim_lookup.txt

# Function to fix headers and content
process_meta_file() {
    local meta_file=$1
    local output_file=$2

    echo "Processing $meta_file -> $output_file"

    # Extract header
    head -n1 "$meta_file" > header.txt
    tail -n +2 "$meta_file" | sort -k1,1 > meta_sorted.txt

    # Join and fix column order: RSID CHR BP ...
    join -t $'\t' -1 1 -2 1 bim_lookup.txt meta_sorted.txt | \
      awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' > merged_body.txt

    # New header
    echo -e "MarkerName\tCHR\tBP\t$(cut -f2- header.txt)" > "$output_file"
    cat merged_body.txt >> "$output_file"
}

# Apply to each meta result
# process_meta_file se_minMAF1e-2_meta1.tbl KUF_MinMAF1e-2Meta_SE_cleaned.tbl
	# this above is the one for w/ UCT. One not commented out is the one we'll use 
process_meta_file se_MinMAF1e-2Meta_NoUCT1.tbl KUF_MinMAF1e-2Meta_NoUCT_SE_cleaned.tbl

# Clean up
rm -f header.txt meta_sorted.txt merged_body.txt bim_lookup.txt

echo " CHR and BP columns correctly aligned in meta result files."

