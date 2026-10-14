# Execution-only entrypoint. Preserve the existing server's scientific config.R.
bundle_root <- Sys.getenv("ERA5_BUNDLE_ROOT",unset="/data/era5")
sys.source(file.path(bundle_root,"config.R"),envir=environment())
ERA5_CONFIG$run$scheduling <- "memory_lpt" # Alternatives: "lpt", "fifo" (no memory admission).
ERA5_CONFIG$run$workers <- 64L             # Pool ceiling; CPU/RAM/QC can reduce it.
ERA5_CONFIG$run$autotune <- FALSE
ERA5_CONFIG$run$memory <- list(
  reserve_fraction = 0.15,                 # At least 4 GiB, plus a separate master allowance.
  safety_factor = 1.30,                    # Applied to conservative worker peak estimates.
  backfill = TRUE
)
ERA5_CONFIG$run$resume <- "auto"           # Resume only this version's validated checkpoints.
# Never mix the old 14-tile benchmark with this changed execution fingerprint.
ERA5_CONFIG$remote$output_dir <- file.path(bundle_root,"results_memory_v2")
