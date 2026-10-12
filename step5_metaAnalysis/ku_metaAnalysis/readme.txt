Ran on VM, so I got the Linux version of METAL from the UMich website for it: https://csg.sph.umich.edu/abecasis/Metal/download/.

METAL requires some edits before it can perform association testing using POLMM outputs. 
Refer to the following steps

1) Transform the association testing results from SAIGE using TransformSAIGE_BySite.sh
2) Change the pval/se_submit_metal.txt scripts and run w/ ./metal se_submit_metal.
3) Transform name of headers from METAL output for topr analysis using add_chr_bp_to_metas.sh
4) Annotate on effective sample size using add_n_ku.sh
	a) Requires files inputted into meta analysis
	b) NOTE: I have not run this script, the one used to create the working meta analysis scripts was run using the files which double reported sample size
                b1). Please let me know if this script works when you can
