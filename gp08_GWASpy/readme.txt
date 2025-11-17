These are the most recent version of the Khat GWAS using filters that match w/ Lerato's.
This version presumes the plink files are all split up by site.

However, I used MAC <= 20 in this. This is quite small and those variants that would be removed
through MAC <= 20 are then removed by MAF 0.05. MAF 0.05 removes MAC <= 20 and more

If you would like to access the original files that have gone through up to GWASpy, you can
find them here:
gs://neurogap-bge-imputed-regional/lerato/wave2/plink_files/*_passed_all_qc.bed
where * can just be replaced with the name of any of the sites.



What is below is notes on the scripts I made:

### Step 1 ###
Regarding step1 (making the phenotype files), I was having some issues with it in the beginning. I think this was an ID issue but I just moved 
forward with my part of the project without this. The only change I made was just adding on the study_site column for all sites.

### Step 2 ###
This version of step2 just transforms NA values of assist_khat_amt (khat use frequency in past 3 months) to 0 if assist_khat = 0 for that sample.
One sample was dropped in this @ UCT due to remaining NA after this. This script also creates 2 stacked barplots for the distribution of khat use and khat
use frequency across sites

### Step 3 ###
This filtering step just does MAF 0.05, MAC 20, and LD 0.10 using the data Lerato made going through the GWASpy pipeline. You'll need to get
that specific code she used to make the data that goes through this pipeline if you'd like to recreate this exactly. Otherwise, you can find the
data here: gs://neurogap-bge-imputed-regional/lerato/wave2/plink_files/*_passed_all_qc.bed
I talk about the filtering steps in these slides: https://docs.google.com/presentation/d/1ijAZ6LsRDVAbitzvqrd2rX2o3m_1rXD4s9Su2tjDRh8/edit?slide=id.g368ca55b1a0_0_14#slide=id.g368ca55b1a0_0_14

Regarding using MAF and MAC, usually one or the other is used. However, I just continued using what I thought was the current version that uses both.
A MAC of 20 is very small. MAF 0.05 removes MAC 20 vars and more so there's no need to worry there.

### Step 4 POLMM & SAIGE null ###
The biggest note here is that study_site was added as a covariate for AAU & Uganda w/ AAU, Uganda, KEMRI, Moi, and UCT all sharing PCs 1-10, age, & sex
I ran the POLMM null model (khat use frequency, ordinal data) creation steps locally on my laptop. I remember it working for hailbatch, but I haven't tried it recently since we've been
experiencing issues recently. This takes a while though (maybe a few hours for all of them). In the current script there is an option to run ALL sites at once.
I don't recommend doing that, just do it one at a time so you can see the outputs as they come.

SAIGE null model (khat use, binary data) was run through a VM. It was pretty fast from what I remember

### Step 5 Association Testing ###
There are filters for SAIGE & POLMM association testing. 
POLMM had variant level missingness filter of 0.15, and minMAF 1e-2 (removes so many vars)
SAIGE has the same variant level missingness filter, but instead uses minMAC 0.5 to ensure at least one minor allele is present in association testing.

Because of this discrepency, you could probably try filtering the association testing data MAC 1 then try doing step 4 then 5 for POLMM only if you
want to see if it results in more sites passing. I'm doubtful it will work as filtering out even more variants with a minMAF filter of 1e-3 or the
default 1e-4 resulted in so many sites still failing. In this case of testing, you could set minMAF to 1e-100000 or something extremely small to
let as many vars pass since that MAC filter done through PLINK makes the POLMM assoc test match SAIGE.

Another difference between POLMM & SAIGE is how they store association test outputs. POLMM will keep all variants but flag them as NA, and SAIGE just
drops them. I think this is what results in POLMM taking so long as compared to SAIGE alongside the more intensive calculation for ranked data.

Another note, the association testing step for POLMM is EXTREMELY slow. It took me maybe 5 days to run it. This could probably be
cut in half if you just run muliple association tests, but that entails making the vm a lot larger (maybe 2.5 times as large as it currently is).
POLMM has these large spikes where in the current VM's cpu specs it jumps between 20% - 80% and it will do this throughout association testing.
Current version of the script just runs things one after another. 


### Step 6 Meta-analysis ###
I did meta-analysis two ways:
1. SAIGE output (khat use, binary): PLINK --> this version was adapted from the original hailbatch script you made

2. POLMM output (khat use frequency, ordinal): METAL --> I talked to Zan about this since they were using it before and I saw it online for meta-analysis.
I used it because it was pretty easy to adapt. I made a readme.txt file in that folder with more information on how to this type of meta-analysis.
This part of the pipeline isn't that pretty and could definitely be made easier if you want, but following the readme.txt file results in working outputs.
	In this, I did meta-analysis through METAL using the standard error approach. There is a pvalue version, but everything comes out really inflated
		(e.g., the lambda gc of an SE meta-analysis was ~1.0 and the pval approach was ~3.8).
		Alicia and Toni just said to use that SE approach since it follows the same as the PLINK approach

### Step 7 Manhattan Plot ###
There are a few files I loaded in this script:
- Khat Use
- Khat Use with minMAF 1e-2 manually applied
- KUF w/out UCT
- KUF w/ UCT

Of note, you only want Khat Use and KUF w/out UCT Manhattan Plot outputs and the respective lambda GC values with their Q-Q plots. 
Those other two were tests.

### Step 8 Clumping ###

Pretty straightforward script, has more notes there if you'd like to know what each of the functions does. You can change the referenced
file to be anything, but you need the combined association testing file to capture LD. Tried combining the by-site files and it was a bit weird.
Try using this file: gs://neurogap-bge-imputed-regional/lerato/wave2/plink_files/all_sites_all_phenos.*


I'll try to cleanup things later on, but everything that I used to get to where the Khat GWAS currently is should be here.
