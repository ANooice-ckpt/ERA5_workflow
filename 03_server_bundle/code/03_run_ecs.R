#!/usr/bin/env Rscript
# Production ECS runner: exact original brms M0/M1 -> fixed 2000 draw
# warm-boundary + streaming yearly CDD outputs. Uses only packed weather.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: Rscript --vanilla 03_run_ecs.R config.R")
cmd <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",grep("^--file=",cmd,value=TRUE)[1L])
base <- dirname(normalizePath(script,winslash="/",mustWork=TRUE))
source(file.path(base,"00_core_loader.R"),local=FALSE)
cfg <- stream_read_config(normalizePath(args[[1L]],mustWork=TRUE))

# Set environment before loading any BLAS/OpenMP-enabled package or PSOCK worker.
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
  VECLIB_MAXIMUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
stream_load_warm02(file.path(base,"warm02_reference.R"))

main <- function() {
  needed <- c("ncdf4","digest","Rcpp","brms","posterior","rstan","dlnm","FNN","splines")
  missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop("Install ECS R packages first: ",paste(missing,collapse=", "))
  if (!identical(as.integer(cfg$pack$years),2010:2019)) stop("Years must be 2010-2019.")
  r <- cfg$run
  if (!identical(as.integer(r$n_draws),2000L)) stop("This contract requires 2000 posterior draws.")
  scenarios <- fast_stream_validate_scenarios_preflight(cfg$cdd$scenarios,2000L)
  pack_dir <- cfg$remote$pack_dir
  index_path <- file.path(pack_dir,"pack_index.rds")
  if (!file.exists(index_path)) stop("Missing locally prepared pack_index.rds: ",index_path)
  index <- readRDS(index_path)
  if (!identical(as.integer(index$years),2010:2019) ||
      index$n_cells != cfg$pack$expected_selected_cells ||
      length(index$year_files)!=10L ||
      !identical(as.integer(cfg$pack$tile_side),as.integer(index$tile_side)))
    stop("Packed domain/years differ from configuration.")
  if (length(index$tiles$tile_id) != index$n_tiles ||
      any(diff(index$tiles$start_row) != head(index$tiles$n_cells,-1L)))
    stop("Packed tile continuity failed.")
  if (any(vapply(index$year_files,function(x) is.null(x$md5)||!nzchar(x$md5),logical(1))))
    stop("Packed input years have missing MD5 receipts.")
  if (!dir.exists(dirname(cfg$remote$output_dir)))
    dir.create(dirname(cfg$remote$output_dir),recursive=TRUE,showWarnings=FALSE)
  for (v in c("m0_file","m1_file","data_file","pred_file"))
    if (!file.exists(cfg$remote[[v]])) stop("Missing fitted-model support file: ",cfg$remote[[v]])

  # The reference's numerical tests run with unmodified original tile code.
  # They do not run a global projection and do not need original ERA5 files.
  message("Running original warm02 synthetic tile/checkpoint selftest...")
  fast_tile_selftest()
  fast_daily_cdd_selftest()
  source(file.path(base,"02_stream_overrides.R"),local=.GlobalEnv)
  source(file.path(base,"09_parallel_runtime.R"),local=.GlobalEnv)
  fast_stream_set_context(index)
  # Revalidate fully expanded parameters with exactly 2000 draws.
  stream_scenarios <- fast_stream_validate_scenarios(scenarios,2000L)

  # Explicit configuration: follow the original M0/M1 posterior-selection contract.
  config <- list(
    m1_file=cfg$remote$m1_file,m0_file=cfg$remote$m0_file,
    data_file=cfg$remote$data_file,pred_file=cfg$remote$pred_file,
    output_dir=cfg$remote$output_dir,
    fast=list(backend=r$backend,workers=as.integer(r$workers),
      worker_candidates=as.integer(r$worker_candidates),autotune=isTRUE(r$autotune),
      resume=r$resume,runtime_selftest=FALSE,target_total_hours=30,
      scheduling=r$scheduling,progress_seconds=r$progress_seconds),
    daily_cache=list(enabled=FALSE,temperature_mode="off"),
    n_draws=Inf,seed=20261002L,epred_qc=TRUE,strict_priors=TRUE,
    joint_support=isTRUE(r$joint_support),joint_k=20L,joint_quantile=0.95,
    joint_deduplicate=TRUE,scenario_weight_scheme=r$weight_scheme,
    figure_dpi=300L,
    global=list(input_dir=pack_dir,file_template="packed_%04d.nc",
      years=2010:2019,n_draws=2000L,spatial_stride=1L,
      tile_side=as.integer(cfg$pack$tile_side),day_chunk=as.integer(r$day_chunk),
      lat_range=cfg$pack$latitude_range,joint_support=isTRUE(r$joint_support),
      min_date_coverage=r$min_date_coverage,min_supported_days=as.integer(r$min_supported_days),
      min_supported_fraction=r$min_supported_fraction,
      support_fraction_sensitivity=r$support_fraction_sensitivity,
      weight_scheme=r$weight_scheme,save_annual_grid=isTRUE(r$save_annual_warm_grid),
      max_cells=r$max_cells,delta_threshold=0.20,preprocessing_acknowledged=TRUE,
      checksum_inputs=TRUE,
      pack_identity=index$identity,pack_index_md5=unname(tools::md5sum(index_path)),
      cdd_scenarios=stream_scenarios,
      engine_md5=unname(tools::md5sum(file.path(base,"02_stream_overrides.R"))),
      parallel_runtime_md5=unname(tools::md5sum(file.path(base,"09_parallel_runtime.R"))),
      original_reference_md5=unname(tools::md5sum(file.path(base,"warm02_reference.R"))),
      spatial_mask_rule="Natural Earth 10m any overlap AND 2015 GHSL population density >1")
  )
  # The scientific fingerprint includes cdd_scenarios, packed input hashes,
  # original posterior mapping and all overridden numerical function bodies.
  dest <- fast_prepare_run(config$output_dir,fast_defaults(config))
  numerical_done <- FALSE
  on.exit({
    if (!numerical_done)
      try(fast_stream_write_run_status(dest$path,"FAILED"),silent=TRUE)
    fast_release_run_lock(dest$path,dest$locktoken)
  },add=TRUE)
  fast_stream_write_run_status(dest$path,"RUNNING")
  config$output_root <- dest$root
  config$output_dir <- dest$path
  message("Run directory: ",config$output_dir)
  message("Packed land-population grid: ",index$n_cells," cells / ",index$n_tiles," spatial tiles")
  message("Fixed paired posterior draws: 2000; declared CDD scenarios: ",
    paste(vapply(stream_scenarios,`[[`,"","name"),collapse=", "))
  utils::write.csv(data.frame(scenario=vapply(stream_scenarios,`[[`,"","name"),
    draw_count=rep(2000L,length(stream_scenarios))),
    file.path(config$output_dir,"declared_cdd_scenarios.csv"),row.names=FALSE)
  saveRDS(list(config=config,packed_input_identity=index$identity,
    scenarios=stream_scenarios,started=Sys.time()),
    file.path(config$output_dir,"stream_run_plan.rds"))
  message("Loading already fitted M1/M0 models and running original posterior QC...")
  analysis <- load_analysis(config)
  if (nrow(analysis$draws_m1) < 2000L || nrow(analysis$draws_m0) < 2000L)
    stop("At least 2000 posterior samples are required in EACH fitted model.")
  saveRDS(list(paired_m1=analysis$draw_ids_m1,paired_m0=analysis$draw_ids_m0,
    pairing_explanation=analysis$draw_metadata),
    file.path(config$output_dir,"original_model_draw_ids.rds"))
  message("Running exact original warm prediction, annual CDD accumulation, and global summaries...")
  result <- run_global_era5(config,analysis)
  finish <- list(status="COMPLETED",finished=Sys.time(),
    fingerprint=result$metadata$fast_execution$fingerprint,
    cdd_file=file.path(result$output_dir,"cdd_global_posterior_draws.rds"),
    n_cells=index$n_cells,n_draws=2000L)
  fast_atomic_save_rds(finish,file.path(result$output_dir,"numerical_completion.rds"),replace=TRUE)
  numerical_done <- TRUE
  fast_stream_write_run_status(dest$path,"COMPLETED",finish$fingerprint)
  if (isTRUE(r$make_plots)) {
    message("Plotting warm-boundary maps from committed partitions...")
    fast_plot_preflight(config$output_dir)
    plot_global_maps(result)
    plot_global_annual_summary(result)
  }
  message("Completed ECS numeric projection. Retain all global/checkpoints/tile_*.rds files for later regional CDD posterior aggregation.")
  invisible(result)
}

# Parameter validation before loading any model/pack.
fast_stream_validate_scenarios_preflight <- function(scenarios,n_draws) {
  if (!is.list(scenarios) || !length(scenarios))
    stop("Fill cfg$cdd$scenarios with declared theta0/slope/kappa; cannot invent a scientific cooling scenario.")
  nm <- vapply(scenarios,function(x) if (is.list(x) && is.character(x$name) && length(x$name)==1L)
    x$name else "",character(1))
  if (any(!nzchar(nm)) || anyDuplicated(nm)) stop("Unique named scenarios required.")
  for (s in scenarios) for (key in c("theta0","slope","kappa")) {
    x <- s[[key]]
    if (!is.numeric(x) || !length(x) %in% c(1L,n_draws) || any(!is.finite(x)) ||
        (key=="slope" && any(x<=0)) || (key=="kappa" && any(x<0)))
      stop("Invalid scenario parameter ",s$name,"$",key)
  }
  scenarios
}

tryCatch(main(),error=function(e) {
  message("ECS ERROR: ",conditionMessage(e))
  stop(e)
})
