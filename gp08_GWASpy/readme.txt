Following are notes on all the files that are here:

prestep1_GP_missingness_filter: Pre-filtering of data (GP0.8, sample relatedness filters). The most recent version is held by Lerato.
  Ran through Hailbatch clusters
  
step2_filter_snps_ldprune_by_site: Using PLINK to filter files by variant & sample-level missingness filters as well as LD Pruning.
  Ran on VM.

step3_fit_null_khat_binary: Creation of null model for Khat Use (KU) GWAS using SAIGE.
  Ran through Hailbatch.

step3_fit_null_khat_ordinal: Creation fo null model for Khat Use Frequency (KUF) GWAS using POLMM.
  Ran through Hailbatch.

step4_binary_association_test: Binary association testing for KU GWAS using SAIGE.
  Ran through Hailbatch.

step4_ordinal_association_test: Ordinal association testing for KUF GWAS using POLMM.
  Ran on spark cluster.

step5_binary_meta_analysis: Standard Error-based meta-analysis for KU GWAS using METAL.
  Ran on VM.

step5_ordinal_meta_analysis: Standard Error-based meta-analysis for KUF GWAS using METAL. 
  Ran on VM.

step6_manhattan_plots.Rmd: Creation of Manhattan Plots, QQ-Plots, and Lambda GC calculations for KU and KUF GWAS.
  Ran locally through Rstudio.

step7_clumping.txt: Clumping script used.
  Ran on VM.
  # Note, may need to rerun using *pass_all_qc.{bed/bed/bim} files
  # currently using NeuroGAP_passed_all_qc_no_multi_allele.{bed/bim/fam}

step8_beta_comparisons: Liftovers, harmonizations, and beta-comparisons script.
  Ran on VM - Liftovers & Harmonization
  Ran locally through Rstudio - beta_comparisons.Rmd

step9_h2: Using LDSC to calculate heritability estimate.
  Ran on VM

step10_twas_fusion: Using FUSION to run TWAS of KU and KUF variants.
  Ran on sparkcluster.

step11_twas_fusion_manhattan: Creation of Manhattan Plots for FUSION TWAS outputs as well as viewing of distribution of colocalization values.
  Ran locally through Rstudio

step12_twas_smr: TWAS using Summary Based Mendelian Randomization (SMR) done on tissues which had significant hits from KU and KUF FUSION TWAS outcomes.
  Ran through Hailbatch

step13_twas_smr_manhattan_plots: Creation of Manhattan plots for SMR TWAS outputs.
  Ran through Rstudio

To do:
  - Clarify clumping script
  - Locuszoom plots
  - GCTA GREML h2 estimate
  - PRS accuracy changes w/ KU and KUF
  - Add sparkcluster startup commands and approximate run times




