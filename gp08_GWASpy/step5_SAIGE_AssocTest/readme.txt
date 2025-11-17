This folder has the SAIGE step 2 scripts I used for association testing of the binary phenotype. 
I also have this other script, FilterMinMAF1e-2.sh, that automatically flags variants from the association test
output that have MAF < 0.01 (is equal to POLMM's step 2 association test minMAF 1e-2). This was
made for the purpose of comparison to the POLMM association test. I also commented it on the 
run_saige_covariate_step2_local.sh file the parameter if you'd like to run it.

Moving forward we agreed to not use this minMAF 0.01 filter because it removes too many variants.

