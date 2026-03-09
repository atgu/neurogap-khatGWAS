This script is to perform rg calculation using LDSC

Ref panel for running rg_ku_kuf_afr.sh can be found here:
gs://neurogap-bge-imputed-regional/lerato/wave2/reference_panels/HGDP_1KG_AFR.bed

Please run in this order:

1) clean_afr.sh
	a) standardizes N, N_eff, Neff, etc. columns to be the same.  
	b) preprocessing of files was also done beforehand.
2) munge_afr.sh  
3) rg_ku_kuf_afr.sh
