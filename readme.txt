Made these files using a covariates file that Lerato provided as well as the original phenotype files that Mary made.
Have been having difficulties matching between original NeuroGAP-release epi-level data's IIDs and other files, so I just 
carried over the study_site covariate from Lerato's file. If a study_site value had an NA for AAU or Uganda, I removed it.

Lerato also removed PC outliers in this phenotype file, so I removed those as well by copying the exact samples she used. 
For AAU and Uganda, there were some NA values at the study_site value. I just removed those as indicated below.

w_site_AAU_pheno.txt: removed 15 rows with study_site = NA
w_site_Uganda_pheno.txt: removed 861 rows with study_site = NA

Here are the samples before copying and after copying.

nvalenci@wm179-a4b neurogap-khatGWAS % gsutil cat gs://neurogap-bge-imputed-regional/nico/khat_gwas/phenotype_files/special_char_change/no_hyphen/w_site/lerato_pheno_study_site.txt | wc -l
   33211
nvalenci@khat:~/pheno/debug$ wc -l no_char_*
  10475 no_char_AAU_khat_pheno.txt
   2709 no_char_KEMRI_khat_pheno.txt
   4059 no_char_Moi_khat_pheno.txt
   7804 no_char_UCT_khat_pheno.txt
   9816 no_char_Uganda_khat_pheno.txt
  34863 total
nvalenci@khat:~/pheno/debug$ wc -l w_site_*
  10128 w_site_AAU_pheno_noNA.txt
   2588 w_site_KEMRI_pheno.txt
   3725 w_site_Moi_pheno.txt
   7337 w_site_UCT_pheno.txt
   8560 w_site_Uganda_pheno_noNA.txt
  32338 total

First file that is having rows measured is the phenotype file Lerato provided w/ all site info. 
Second file is the phenotype files before study_site is added.
Last file is after study_site is added. There's an extra 3 missing samples between this file and Lerato's study_site covariate file.
	This is just missing samples filtered out from other sites

SAIGE and POLMM drop samples w/ any missing data, so those NAs above would have been dropped irrespective.

I think that there may be some mismatched samples still? Will need to talk to everyone about it.
