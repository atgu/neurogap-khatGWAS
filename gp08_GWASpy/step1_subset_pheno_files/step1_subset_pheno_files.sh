#!/bin/bash
set -euo pipefail

LERATO="lerato_pheno_study_site.txt"

# Build IID -> study_site map from lerato:
# IID = col 1
# study_site = second-to-last column (PDO is last)
awk 'BEGIN{FS="[ \t]+"; OFS="\t"}
NR==1{next}
{
  iid=$1
  site=$(NF-1)
  if(iid!=""){
    if(!(iid in seen)){
      print iid, site
      seen[iid]=1
    }
  }
}' "$LERATO" > lerato_iid_to_site.tsv


# Process each no_char file
for f in no_char_*_khat_pheno.txt; do
  [[ -f "$f" ]] || continue
  out="${f%.txt}_subset_lerato_with_studysite.tsv"

  awk -v OFS="\t" '
    BEGIN{FS="[ \t]+"}

    # Load IID -> study_site map
    NR==FNR { site[$1]=$2; next }

    # Header
    FNR==1 {
      amt_col=-1
      for(i=1;i<=NF;i++){
        if($i=="assist_khat_amt") amt_col=i
      }
      # Print header + study_site
      for(i=1;i<=NF;i++){
        printf "%s", $i
        if(i<NF) printf OFS
      }
      print OFS "study_site"
      next
    }

    # Data rows (subset to IIDs present in lerato)
    {
      iid=$1
      if(!(iid in site)) next

      # Convert NA -> 0 in assist_khat_amt
      if(amt_col>0 && $amt_col=="NA") $amt_col="0"

      # Print row as TSV
      for(i=1;i<=NF;i++){
        printf "%s", $i
        if(i<NF) printf OFS
      }
      print OFS site[iid]
    }
  ' lerato_iid_to_site.tsv "$f" > "$out"

  rows=$(awk 'END{print NR-1}' "$out")
  echo "Wrote: $out (data rows kept: $rows)"
done

echo "Done."
