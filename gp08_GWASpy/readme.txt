Following are notes on all the files that are here:

prestep1_GP_missingness_filter: Pre-filtering of data (GP0.8, sample relatedness filters). The most recent version is held by Lerato.
  Ran through Hailbatch clusters.
  
step2_filter_snps_ldprune_by_site: Using PLINK to filter files by variant & sample-level missingness filters as well as LD Pruning.
  Ran on VM. Takes ~1-2 hrs in total? Maybe less.

step3_fit_null_khat_binary: Creation of null model for Khat Use (KU) GWAS using SAIGE.
  Ran through Hailbatch. Very fast, less than 5 minutes probably

step3_fit_null_khat_ordinal: Creation fo null model for Khat Use Frequency (KUF) GWAS using POLMM.
  Ran through Hailbatch. Takes ~3.5 hrs.

step4_binary_association_test: Binary association testing for KU GWAS using SAIGE.
  Ran through Hailbatch. Takes ~30 minutes.

step4_ordinal_association_test: Ordinal association testing for KUF GWAS using POLMM. Prefiltering of data (MAF 0.01 == POLMM minMAF 1e-2)
  Ran on spark cluster. Should take maybe ~2.5 - 3 days if ran all through one spark cluster. Can try to split it across multiple sparkclusters
hailctl dataproc start nv-fusion \
  --master-machine-type n1-highmem-16 \
  --worker-machine-type n1-highmem-16 \
  --master-boot-disk-size 450 \
  --num-workers 0 \
  --max-idle=100h \
  --requester-pays-allow-all \
  --big-executors \
  --no-off-heap-memory \
  --region us-central1 \
  --network broad-allow
If too small, bump up to  
  --master-machine-type n1-highmem-32 \
  --worker-machine-type n1-highmem-32 \

step5_binary_meta_analysis: Standard Error-based meta-analysis for KU GWAS using METAL.
  Ran on VM. Takes ~20-30 minutes.

step5_ordinal_meta_analysis: Standard Error-based meta-analysis for KUF GWAS using METAL. 
  Ran on VM. Takes ~20-30 minutes

step6_manhattan_plots.Rmd: Creation of Manhattan Plots, QQ-Plots, and Lambda GC calculations for KU and KUF GWAS.
  Ran locally through Rstudio. Takes ~30 minutes.

step7_clumping.txt: Clumping script used.
  Ran on VM. Takes ~10 minutes.
  # Note, may need to rerun using *pass_all_qc.{bed/bed/bim} files
  # currently using NeuroGAP_passed_all_qc_no_multi_allele.{bed/bim/fam}

step8_beta_comparisons: Liftovers, harmonizations, and beta-comparisons script.
  Ran on VM - Liftovers & Harmonization. Takes ~6 hrs .
  Ran locally through Rstudio - beta_comparisons.Rmd Takes a day.

step9_h2: Using LDSC to calculate heritability estimate.
  Ran on VM. Takes ~10 minutes

step10_twas_fusion: Using FUSION to run TWAS of KU and KUF variants.
  Ran on sparkcluster. Takes ~3-4 days (?).
  --master-machine-type n1-highmem-32 \
  --worker-machine-type n1-highmem-32 \
  --master-boot-disk-size 400 \
  --num-workers 0 \
  --max-idle=100h \
  --requester-pays-allow-all \
  --big-executors \
  --no-off-heap-memory \
  --region us-central1 \
I lost access to the setup notes for the sparkcluster. If n1-highmem-32 is too big, you can try to set MAX_JOBS to 5 (probably unstable).

step11_twas_fusion_manhattan: Creation of Manhattan Plots for FUSION TWAS outputs as well as viewing of distribution of colocalization values.
  Ran locally through Rstudio. Takes ~1-1.5 hrs.

step12_twas_smr: TWAS using Summary Based Mendelian Randomization (SMR) done on tissues which had significant hits from KU and KUF FUSION TWAS outcomes.
  Ran through Hailbatch. Takes ~1 hr.

step13_twas_smr_manhattan_plots: Creation of Manhattan plots for SMR TWAS outputs.
  Ran through Rstudio. Takes ~30 minutes.

To do:
  - Locuszoom plots
  - GCTA GREML h2 estimate
  - PRS accuracy changes w/ KU and KUF




