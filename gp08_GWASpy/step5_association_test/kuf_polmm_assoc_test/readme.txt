There are a few scripts here:

1) The site (*run_polmm.sh) scripts
2) run_AllSites_polmm.sh
4) run_polmm_assoc.R

run_AllSites_polmm.sh runs all the *run_polmm.sh sites one after another. This takes a while (~1 day per site).

By making the chunks larger in run_polmm_assoc.R you can make this run faster but increase computational load.

Currently it is set to "memory.chunk  = 0.5".
