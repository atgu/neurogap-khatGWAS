Follow these steps:
1). Run liftover on respective substance use-adjacent files using the liftover scripts from liftovers
2). Run harmonization on those lifted over files to match information w/ KU/KUF GWAS
3). Run beta comparisons using beta_comparisons.Rmd locally

As a general summary:
1). Liftovers translate grch37 files to grch38
2). Harmonization performs the following steps to make sure alleles match before matching rsids:
  a) if alleles match exactly, keep the beta as is.
  b) if alleles are swapped, flip the beta sign.
  c) if alleles match after taking complements, keep or flip depending on whether the effect allele still aligns.
  d) if they still don’t match after checking complements, then drop the SNP.
  e) drop palindromic allele pairings (e.g., A1/A2 --> A/T --> drop)
3). Beta comparisons takes a while because there's so many files and variants to plot. Functionally, we
	a) If effect is negative for KU or KUF, we multiply effect ffrom KU or KUF and their respective comparison file by -1 
		a1) Orients effects to only positive KU/KUF effect sizes
	b) Calculate Pearson correlation and pval
	c) Calculte Deming Regression and pval (full line)
	d) Create a line of best fit (dashed line)
	e) For Saunders Smoking Initiation, we plot only vars with GWAS pval < 0.05 as well
4) Create a table of Pearson Correlations, Deming Regressions, Variants in Figure, variant with GWAS pval < 0.05
