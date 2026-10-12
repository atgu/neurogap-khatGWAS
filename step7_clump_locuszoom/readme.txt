The order of running these scripts goes

1) Clumping using clumping.txt
2) Using the lead variants (of which I just passed directly as a vector of rsIDs), get all info for these lead vars for KU or KUF
3) Extract LD (r^2 calculated using plink) for vars +/-500kb of indicated lead variants
4) Plot locuszoom plots using ku_kuf_locuszoomR.Rmd (using GRCh38 files for recombination rates and gene names)
