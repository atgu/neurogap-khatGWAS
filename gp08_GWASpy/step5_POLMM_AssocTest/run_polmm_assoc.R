args <- commandArgs(trailingOnly = TRUE)
null_model_file <- args[1]
plink_prefix <- args[2]
output_prefix <- args[3]

library(POLMM)

objNull <- readRDS(null_model_file)
chrVec <- names(objNull$LOCOList)

# when working with more phenotypes or larger plink files, it's best to split the number of chunks to a smaller number. 
	# the equation generally goes (some large number)/(memory.chunk) = actual number of chunks job is split into
# this will take an eternity to run (~5 days for me)
# could look into parallelizing this by maybe making the VM larger and just running 2 sites at a time?
	# cpu usage bounced between ~30% - 40% (sometimes 60%)
		# going above that caused the job to crash and stop though, but that was before I split 
		# the memory chunks. Maybe it won't if you just run it as is like this 	
POLMM.plink(
  objNull       = objNull,
  PlinkFile     = plink_prefix,
  output.file   = output_prefix,
  chrVec.plink  = chrVec,
  memory.chunk  = 0.5,
  SPAcutoff     = 2,
  minMAF        = 1e-2,
  maxMissing    = 0.15,
  impute.method = "fixed",
  G.model       = "Add"
)
