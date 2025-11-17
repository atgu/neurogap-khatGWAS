There are a few scripts here:

1) The site (*run_polmm.sh) scripts
2) run_AllSites_polmm.sh
3) run_polmm.sh
4) run_polmm_assoc.R

The *run_polmm.sh files is basically the run_polmm.sh file but the site has been specified to match the 
file name. These files were then ran subsequently using run_AllSites_polmm.sh

run_polmm.sh can be used to specify specific sites to run.

run_polmm_assoc.R is the R script that is used when you run. In it, you can specify the 
maxMissing filter (currently maxMissing 0.15, which matches with SAIGE) and the minMAF (currently
set to minMAF 1e-2, which is really strict and does not match with SAIGE).

I've been having a hard time running this through Hail even before the more recent issues
with reading in bigger plink files. Thats why we have it being run in a VM instead.

More CS-type of thoughts: run_AllSites_polmm.sh runs each *run_polmm.sh site file one by one.
POLMM's CPU usage for the current VM specs has it bounce between 20-30% CPU usage all the 
way up to 70-80%. This varies by site, with the smaller sites using less. It usually averages
around 40-60% usage though. This random large spiking in CPU usage causes the script to sometimes
crash due to using more than the VM's computational load. 

This run_AllSites_polmm.sh takes 5 days to run. It now takes very long because of the extra samples
we added by transforming assist_khat_amt NA values to 0. This increases the number of samples
almost by 8 times for each site and overall.

I have some thoughts on how this can be made more efficient, though.

1) You can breakup the runs into more chunks (refer to run_polmm_assoc.R for the information)
  --> with this, then you could probably run more sites at once
2) Upgrade the VM
3) Figuring out a way of instead flagging NAs, you instead drop them from the output files
  --> Unlike SAIGE which removes variants with its association testing filtering parameters,
	POLMM flags them as NA. 
  --> This might be an option in POLMM? But by default it just chooses to flag them
