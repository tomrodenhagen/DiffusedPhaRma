


scripts = c("scripts/run_measure_diff_cycles.R",
	    "scripts/run_base.R",
	    "scripts/run_different_measurements.R",
	    "scripts/run_fixed_noise.R",
	    "scripts/join_and_compare.R")
for(script in scripts)
{ print(paste("Run", script))
  source(script)
}



