Ran on VM, so I got the Linux version of METAL from the UMich website for it: https://csg.sph.umich.edu/abecasis/Metal/download/.

METAL requires some edits before it can perform association testing using POLMM outputs. 
Refer tot he following steps

1) Transform the association testing results from POLMM using TransformPOLMM_BySite.sh
	a) In this, you need the sample sizes that POLMM outputs in the null model --> these are the samples that go into association testing
2) Change the pval/se_submit_metal.txt scripts and run w/ ./metal se_submit_metal.
3) Transform name of headers from METAL output for topr analysis using add_chr_bp_to_metas.sh
4) Add sample size using add_n_kuf.sh
	a) Requires files inputted into meta analysis
	b) NOTE: I have not run this script, the one used to create the working meta analysis scripts was run using the files which double reported sample size
		b1). Please let me know if this script works when you can

After this transformation, there are a few more steps required to get effect allele frequencies for running SMR:
1) Create neurogap file without uct for calculating effect allele freuqnecies for KUF using create_ngap_passed_all_qc_no_uct.sh 
2) Calculate effect allele frequency using calculate_kuf_eaf.sh
3) Annotate on those effects using annotate_kuf_eaf.sh
