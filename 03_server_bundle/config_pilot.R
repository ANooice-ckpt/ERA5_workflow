# Same scientific settings as config.R; isolated output prevents mixing with full runs.
bundle_root <- Sys.getenv("ERA5_BUNDLE_ROOT",unset="/data/era5")
sys.source(file.path(bundle_root,"config.R"),envir=environment())
ERA5_CONFIG$run$max_cells <- 300L
ERA5_CONFIG$remote$output_dir <- file.path(bundle_root,"results_pilot")
