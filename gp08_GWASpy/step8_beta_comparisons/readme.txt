Follow these steps:
1). Run liftover on respective substance use-adjacent files using the liftover scripts from liftovers
2). Run harmonization on those lifted over files to match information w/ KU/KUF GWAS
3). Run beta comparisons using beta_comparisons.Rmd locally

These scripts do a lot of things, and if you would like more info I believe I've commented most information into them. As a general summary:
1). Liftovers translate grch37 files to grch38
2). Harmonization performs the following steps to make sure alleles match before matching rsids:
  a) if alleles match exactly, keep the beta as is.
  b) if alleles are swapped, flip the beta sign.
  c) if alleles match after taking complements, keep or flip depending on whether the effect allele still aligns.
  d) if they still don’t match after checking complements, then drop the SNP.
  e) drop palindromic allele pairings (e.g., A1/A2 --> A/T --> drop)
3). Beta comparisons takes forver because there's so many files and variants to plot. First half plots beta comparisons between KU/KUF to all other substance use-adjacent GWASs
    Given certain pvalue thresholds, and also calculates their correlations. Second half plots beta comparisons between KU/KUF to all other substance use-adjacent GWASs, but now also
    highlights vars red if they are nominally significant (p<1e-5) in KU/KUF, as well as only plotting vars which are positive in KU/KUF while also calculating correlations
