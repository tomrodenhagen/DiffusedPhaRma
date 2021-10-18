scripts = c("scripts/run_measure_diff_cycles.R FULL_VERY_EXT CONT",
	    "scripts/run_base.R FULL_VERY_EXT CONT",
	    "scripts/run_base_no_retry.R FULL_NO_RETRY CONT",
	    "scripts/run_different_measurements.R FULL_VERY_EXT CONT",
	    "scripts/run_fixed_noise.R FULL_VERY_EXT CONT",
	    "scripts/join_and_compare.R")
for(script in scripts)
{ print(paste("Run", script))
  system(paste("Rscript", script))
}
