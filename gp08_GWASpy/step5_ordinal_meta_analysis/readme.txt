Ran on VM, so I got the Linux version of METAL from the UMich website for it: https://csg.sph.umich.edu/abecasis/Metal/download/.

METAL requires some edits before it can perform association testing using POLMM outputs. 
Refer tot he following steps

1) transform the association testing results from POLMM using transform_polmm.sh
	a) In this, you need the sample sizes that POLMM outputs at the end of its run.
           If you run the run_AllSites_polmm.sh version, it should output a log file for every output 
           from the VM. It says something like "X samples for X vars" and then next to it a 
           smaller sample size. Whatever that smaller sized one is, or if it is the same number
           is the one to use for this. This variation in number of samples comes from maybe some 
	   samples missing covariate or phenotype info when making the null model or performing associaiton testing. 
		Says something like "X samples out of X are being tested"
		a.1) I tried meta-analysis testing using number of replicaitons for each SNP and the data came out really inflated. 
		a.2) Sample size from individuals ONLY USED in association testing is necessary
	b) Add any new sites/change the sites included in that sample listing if necessary
2) Change the pval/se_submit_metal.txt scripts and run w/ ./metal pval/se_submit_metal.txt
	a) Change input file names is all
        b) pval meta analysis method was tested before but it gave bad data
           I added an old version of the pvalue meta-analysis script I made for other files. You can
           run it if you'd like, but the outputs are really inflated. SE is better.
3) Transform name of headers from METAL output for topr analysis using add_chr_bp_to_metas.sh



################################################################################################################

These are association testings done on files that have gone through the following adjusted processes:

1) gwaspy --> refer to pipeline on google slides/ask me or Lerato if you'd like to see
2) ldprune to Lerato's paramter (>10%)
3) minMAF 1e-2 --> appplied in POLMM, and manually applied to SAIGE outputs
4) Uganda and AAU null model fitting have study_site added as a covariate

