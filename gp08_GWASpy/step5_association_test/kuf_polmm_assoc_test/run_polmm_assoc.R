args <- commandArgs(trailingOnly = TRUE)
null_model_file <- args[1]
plink_prefix <- args[2]
output_prefix <- args[3]

library(POLMM)

objNull <- readRDS(null_model_file)
chrVec <- names(objNull$LOCOList)

# when working with more phenotypes or larger plink files, with the current VM resources, it's required to split the job into many smaller chunks 
	# the equation generally goes (some large number)/(memory.chunk) = actual number of chunks job is split into
# this will take an eternity to run (~5 days for me)
# could look into parallelizing this by increasing the number of CPUs by however many jobs association tests you are going to run in parallel larger by (8 * num jobs).
# memory.chunk is default set to 4. With the current VM resources, I tried memory.chunk = 2 and 1 and both failed due to going over CPU resources. 
	# With the current number of chunks the job is split into, cpu usage bounced between ~30% - 40% (sometimes spiking ~60 - 80%)
		# the memory chunks this small. Maybe it won't if you just run it as is like this 	

# Sites with more samples took more resources (AAU being the most). Maybe you could run KEMRI and Moi at the same time if you increase num of CPUs by 50% to also accomodate 
# for spikes (?)

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
