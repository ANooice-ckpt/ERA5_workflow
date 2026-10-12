#!/usr/bin/env Rscript
# Model-free validation after packing and/or ECS production.
# Usage: Rscript --vanilla 04_verify_outputs.R config.R [pack|run|both]
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) || length(args)>2L) stop("Usage: Rscript --vanilla 04_verify_outputs.R config.R [pack|run|both]")
mode <- if (length(args)==2L) args[[2L]] else "both"
if (!mode %in% c("pack","run","both")) stop("Unknown mode")
cmd <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",grep("^--file=",cmd,value=TRUE)[1L])
base <- dirname(normalizePath(script,mustWork=TRUE))
source(file.path(base,"00_core_loader.R"),local=FALSE)
cfg <- stream_read_config(normalizePath(args[[1L]],mustWork=TRUE))
packed <- if (dir.exists(cfg$remote$pack_dir)) cfg$remote$pack_dir else cfg$local$pack_dir
index_file <- file.path(packed,"pack_index.rds")
if (!file.exists(index_file)) stop("No pack index: ",index_file)
index <- readRDS(index_file)
stopifnot(identical(index$schema,"ERA5-indexed-weather-v1"),
  index$n_cells==cfg$pack$expected_selected_cells,
  index$n_cells==nrow(index$selected),index$n_tiles==nrow(index$tiles),
  !anyDuplicated(index$selected$cell_id),
  identical(index$selected$pack_row,seq_len(index$n_cells)))
if (mode %in% c("pack","both")) {
  if (!requireNamespace("ncdf4",quietly=TRUE)) stop("Install ncdf4 for packed file validation")
  for (i in seq_along(index$years)) {
    y<-index$years[i]; file<-file.path(packed,sprintf("packed_%04d.nc",y))
    if (!file.exists(file)) stop("Missing packed year: ",y)
    if (!identical(unname(tools::md5sum(file)),index$year_files[[i]]$md5))
      stop("Packed MD5 changed: ",file)
    nc<-ncdf4::nc_open(file)
    tryCatch({
      stopifnot(nc$dim$cell$len==index$n_cells,
        nc$dim$time$len==length(index$dates[[i]]),
        identical(as.integer(ncdf4::ncvar_get(nc,"cell_id")),as.integer(index$selected$cell_id)),
        all(c("B21_raw","TV_lag1_raw","TV_lag2_raw","TV_lag3_raw","Ta_C") %in% names(nc$var)))
    },finally=ncdf4::nc_close(nc))
    message("Verified packed year ",y," / ",index$year_files[[i]]$bytes/1e9," GB")
  }
  message("PACK VALIDATION PASS: ",index$n_cells," fixed grid cells; ",index$n_tiles," tiles")
}
if (mode %in% c("run","both")) {
  root<-cfg$remote$output_dir
  pointer<-file.path(root,"latest_run.txt")
  if (!file.exists(pointer)) stop("No completed run pointer: ",pointer)
  rundir<-trimws(readLines(pointer,warn=FALSE,encoding="UTF-8"))[1L]
  if (!dir.exists(rundir)) stop("Run directory does not exist: ",rundir)
  g<-file.path(rundir,"global")
  required<-c("fast_manifest.rds","warm_boundary_area_day_posterior_draws.rds",
    "cdd_global_posterior_draws.rds","numerical_completion.rds","global_run_metadata.rds")
  missing<-required[!file.exists(file.path(g,required))]
  if (length(missing)) stop("Missing final outputs: ",paste(missing,collapse=", "))
  manifest<-readRDS(file.path(g,"fast_manifest.rds"))
  complete<-readRDS(file.path(g,"numerical_completion.rds"))
  warm<-readRDS(file.path(g,"warm_boundary_area_day_posterior_draws.rds"))
  cool<-readRDS(file.path(g,"cdd_global_posterior_draws.rds"))
  if (!identical(complete$status,"COMPLETED") ||
      !identical(complete$fingerprint,manifest$fingerprint) ||
      !identical(cool$fingerprint,manifest$fingerprint))
    stop("Run fingerprint or completion status invalid.")
  if (!identical(warm$paired_draw_ids,cool$paired_draw_ids) ||
      nrow(cool$paired_draw_ids)!=2000L ||
      !identical(as.integer(dim(warm$estimates)[3L]),2000L) ||
      !identical(as.integer(dim(cool$area_km2_degree_days)[4L]),2000L))
    stop("Posterior draws/pairing lost or rearranged.")
  stream_load_warm02(file.path(base,"warm02_reference.R"))
  source(file.path(base,"02_stream_overrides.R"),local=.GlobalEnv)
  validation <- fast_stream_validate_global_cdd(cool)
  delta_err <- validation$max_delta_error
  checkpoints<-list.files(file.path(g,"checkpoints"),"^tile_[0-9]+[.]rds$",full.names=TRUE)
  expected_tiles <- if (isTRUE(cool$pilot_not_global)) manifest$total_tiles else index$n_tiles
  expected_cells <- if (isTRUE(cool$pilot_not_global)) manifest$total_land_points else index$n_cells
  if (length(checkpoints)!=expected_tiles)
    stop("Incomplete CDD tile checkpoints: ",length(checkpoints)," / ",expected_tiles)
  cat("RUN VALIDATION PASS (global level and checkpoint count): ",expected_cells,
      " cells; ",length(checkpoints)," checkpoints; 2000 fixed posterior pairs; max delta identity error ",
      format(delta_err,digits=12),"\n",sep="")
  message("Supported period/scenario pairs: ",validation$supported_period_scenarios,
    "; no eligible support (NA outputs): ",validation$empty_period_scenarios)
  message("For a full per-tile numerical integrity scan, use the original warm02 resume path plus the stream checkpoint validator.")
}
