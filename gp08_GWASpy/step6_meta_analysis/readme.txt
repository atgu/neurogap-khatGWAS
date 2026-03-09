Ran on VM, so I got the Linux version of METAL from the UMich website for it: https://csg.sph.umich.edu/abecasis/Metal/download/.

METAL requires some edits before it can perform association testing using POLMM outputs. 
Refer to the following steps

The steps can be broken down into the following for SAIGE and POLMM outputs:



SAIGE:
1) Change the pval/se_submit_metal.txt scripts and run w/ ./metal pval/se_submit_metal.txt
	a) Change input file names 
2) Transform name of headers from METAL output for topr analysis using add_chr_bp_to_metas.sh



POLMM:
1) Transform the association testing results from POLMM using transform_polmm.sh
        a) In this, you need the sample sizes that POLMM outputs from null model creation.
2) Change the pval/se_submit_metal.txt scripts and run w/ ./metal pval/se_submit_metal.txt
        a) Change input file names       
3) Transform name of headers from METAL output for topr analysis using add_chr_bp_to_metas.sh


You can run add_chr_bp_to_metas.sh after running meta-analysis.


