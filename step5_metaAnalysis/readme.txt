NOTE: Cleaning these script still, will label more later

Ran on VM, so I got the Linux version of METAL from the UMich website for it: https://csg.sph.umich.edu/abecasis/Metal/download/.

METAL requires some edits before it can perform association testing using POLMM outputs. 
Refer to the following steps

The steps can be broken down into the following for SAIGE and POLMM outputs:

POLMM:
1) Transform the association testing results from POLMM/SAIGE using TransformSAIGE/POLMM_BySite.sh
        a) In this, you need the sample sizes that POLMM outputs from null model creation for POLMM
		a1) POLMM transformation script also orients all effects to A2 to match SAIGE
	b) SAIGE transformation script creates an N column that is sum of cases and controls for each association test
2) Change the pval/se_submit_metal.txt scripts and run w/ ./metal pval/se_submit_metal.txt
        a) Performs meta analysis using METAL's standard error approach (SE)
	b) Effect indicated by EFFECT A2 A1 for both scripts
		b1) first allele column indicated is your effect allele
		b2)  meta-analysis output is A1 (effect reported from provided A2), and A2 (non-effect allele reported from provided A1)       
3) Transform name of headers from METAL output for topr analysis using add_chr_bp_to_metas.sh
	a) This can be used on both files at once



