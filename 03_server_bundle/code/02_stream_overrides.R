# Replace ONLY I/O/mask/cache interfaces and tile accumulator; inherit fitted-model
# prediction, posterior pairing, support checks and warm finalizer from warm02.
# This file must be sourced AFTER 00_core_loader.R has loaded warm02 definitions.

fast_stream_set_context <- function(index) {
  if (!identical(index$schema, "ERA5-indexed-weather-v1")) stop("Unexpected pack schema.")
  if (nrow(index$selected) != index$n_cells || anyDuplicated(index$selected$cell_id) ||
      anyDuplicated(index$selected$pack_row) ||
      !identical(as.integer(index$selected$pack_row), seq_len(index$n_cells)))
    stop("Invalid packed row/cell IDs.")
  lon <- index$lon_axis; lat <- index$lat_axis
  if (length(lon) != 1440L || !length(lat) || anyNA(lon) || anyNA(lat))
    stop("Invalid full source coordinate axes in pack index.")
  s <- index$selected
  lat_full_indices <- index$lat_indices[s$selected_lat_index]
  expected <- s$selected_lon_index + (lat_full_indices - 1L)*length(lon)
  if (anyNA(expected) || any(expected < 1L) || any(expected > length(lon)*length(lat)))
    stop("Selected grid indices invalid against original full axes.")
  lookup <- integer(length(lon)*length(lat))
  lookup[expected] <- s$pack_row
  ctx <- new.env(parent=emptyenv())
  ctx$index <- index
  ctx$lookup <- lookup
  assign(".fast_stream_ctx",ctx,envir=.GlobalEnv)
  invisible(TRUE)
}

fast_stream_context <- function() {
  if (!exists(".fast_stream_ctx", envir=.GlobalEnv, inherits=FALSE))
    stop("Pack context not initialized in this worker.")
  get(".fast_stream_ctx", envir=.GlobalEnv, inherits=FALSE)
}

# Retain global axes for original warm02 task construction and QC. Packed files
# use cell x time, and their original full coordinates are in pack_index.rds.
inspect_history_files <- function(g) {
  ctx <- fast_stream_context()$index
  if (!identical(as.integer(g$years),as.integer(ctx$years))) stop("Run years and packed years differ.")
  files <- vector("list",length(g$years))
  audit <- vector("list",length(g$years))
  for (i in seq_along(g$years)) {
    y <- g$years[i]
    path <- file.path(g$input_dir,sprintf("packed_%04d.nc",y))
    if (!file.exists(path)) stop("Missing packed input: ",path)
    saved <- ctx$year_files[[i]]
    hash <- unname(tools::md5sum(path))
    if (!identical(hash,saved$md5)) stop("Packed annual file MD5 mismatch: ",path)
    nc <- ncdf4::nc_open(path)
    tryCatch({
      if (nc$dim$cell$len != ctx$n_cells ||
          nc$dim$time$len != length(ctx$dates[[i]])) stop("Packed dimensions changed: ",path)
      att <- ncdf4::ncatt_get(nc,0,"pack_identity")
      if (!isTRUE(att$hasatt) || !identical(as.character(att$value),ctx$identity))
        stop("Packed identity differs from index: ",path)
      if (!all(c(era5_required_variables,"Ta_C") %in% names(nc$var)))
        stop("Packed weather variables missing: ",path)
      cell_ids <- ncdf4::ncvar_get(nc,"cell_id")
      if (!identical(as.integer(cell_ids),as.integer(ctx$selected$cell_id)))
        stop("Packed cell mapping changed: ",path)
    },finally=ncdf4::nc_close(nc))
    dates <- as.Date(ctx$dates[[i]])
    files[[i]] <- list(path=path,year=y,lon=ctx$lon_axis,lat=ctx$lat_axis,
      dates=dates,offsets=stats::setNames(rep(0,length(era5_required_variables)),era5_required_variables),
      provenance=list(source="locally packed from validated ERA5",time_zone="UTC",
        history_definition="Original raw thermal-history variables unchanged"))
    audit[[i]] <- data.frame(file=normalizePath(path,winslash="/",mustWork=TRUE),year=y,
      n_dates=length(dates),start=as.character(min(dates)),end=as.character(max(dates)),
      nlon=length(ctx$lon_axis),nlat=length(ctx$lat_axis),
      input_units="stored numeric after source unit conversions; see local prepack manifest",
      source="packed",md5=hash,bytes=unname(file.info(path)$size),
      pack_identity=ctx$identity,stringsAsFactors=FALSE)
  }
  list(files=files,audit=do.call(rbind,audit))
}

# Exact original geographic rectangle indices -> sparse packed cell indices.
# No interpolation, and nonselected cells remain NA in the temporary rectangle.
read_rectangular_history <- function(nc,variable,lon_indices,lat_indices,start_day,nday,offset=0) {
  if (abs(offset)>0) stop("Packed weather was already converted to its declared units.")
  ctx <- fast_stream_context()
  if (any(lon_indices < 1L | lon_indices > length(ctx$index$lon_axis)) ||
      any(lat_indices < 1L | lat_indices > length(ctx$index$lat_axis)))
    stop("Input grid indices outside packed original coordinate domain.")
  ids <- as.vector(outer(as.integer(lon_indices),
                         (as.integer(lat_indices)-1L)*length(ctx$index$lon_axis),"+"))
  rows <- ctx$lookup[ids]
  out <- matrix(NA_real_,nrow=length(ids),ncol=nday)
  present <- which(rows>0L)
  if (!length(present)) return(out)
  from <- min(rows[present]); to <- max(rows[present])
  if (to-from+1L > length(present) + 2L*ctx$index$tiles$n_cells[1L] + 512L)
    stop("Packed cell span unexpectedly fragmented; packing tile_side differs from ECS task tile_side.")
  if (!variable %in% c(era5_required_variables,"Ta_C")) stop("Unknown packed variable: ",variable)
  value <- ncdf4::ncvar_get(nc,variable,start=c(from,start_day),count=c(to-from+1L,nday),
    collapse_degen=FALSE)
  value <- matrix(as.numeric(value),nrow=to-from+1L,ncol=nday)
  out[present,] <- value[rows[present]-from+1L,,drop=FALSE]
  out
}

# Preserve the user's chosen land rule exactly. area_weight follows warm02's
# former binary-fraction * cos(latitude) weight, NOT an observed land fraction.
load_land_fraction <- function(g,lon,lat,selected_lon,selected_lat) {
  ctx <- fast_stream_context()$index
  if (!same_coordinate(ctx$lon_axis,lon) || !same_coordinate(ctx$lat_axis,lat) ||
      !identical(as.integer(selected_lat),as.integer(ctx$lat_indices)) ||
      !identical(as.integer(selected_lon),seq_along(lon)))
    stop("Mask/native-grid coordinates differ from validated local prepack.")
  points <- expand.grid(lon=lon[selected_lon],lat=lat[selected_lat],KEEP.OUT.ATTRS=FALSE)
  points$cell_id <- seq_len(nrow(points))
  points$lon_plot <- ((points$lon+180) %% 360)-180
  points$land_fraction <- 0
  points$land_fraction[ctx$selected$cell_id] <- 1
  points$cell_area_km2 <- NA_real_
  points$cell_area_km2[ctx$selected$cell_id] <- ctx$selected$cell_area_km2
  points$land <- points$land_fraction>0
  points$area_weight <- cos(points$lat*pi/180)*points$land_fraction
  list(points=points,method=paste("Natural Earth 10m cell-overlap, GHSL 2015 density >1,",
     "full-cell binary weight (not fractional land area), pack identity",ctx$identity))
}

# No daily-cache capacity scan: inference accumulates yearly per-draw values.
fast_daily_defaults <- function(config) {
  d <- list(enabled=FALSE,temperature_mode="off",max_abs_error=5e-5,
    available_space_gb=1,reserve_space_gb=0.5,verify_payload_on_resume=FALSE,
    temperature_checksum_inputs=FALSE)
  if (!is.null(config$daily_cache)) d <- utils::modifyList(d,config$daily_cache)
  if (isTRUE(d$enabled)) stop("The stream run must not enable TB-scale daily cache.")
  d
}
fast_daily_inspect_temperature <- function(config,inspected_history) {
  list(status="packed_Ta_C",mode="stream",files=list(),audit=data.frame(),
    reason="Daily temperature aligned and validated in packed annual NetCDF")
}
fast_daily_support_preflight <- function(d,tasks,files,support_analysis,g,input_audit,ndraw,out_dir) {
  list(supported_rows_by_tile=NULL,table=NULL,
       note="Full serial support scan removed; flags and counts computed once in each tile.")
}
fast_daily_snapshot <- function(config,analysis,bundle,tasks,fingerprint,out_dir) invisible(NULL)
fast_daily_finalize_inventory <- function(bundle,tasks,fingerprint,out_dir)
  list(status="NO_DAILY_CACHE",reason="Yearly grid posterior draws retained in tile checkpoints")

# Numerical completion is authoritative even if a process stopped between
# publishing the numerical receipt and the legacy run-status file. The inherited
# selector checks this function both before and after acquiring the run lock.
fast_existing_run_orig <- fast_existing_run
fast_existing_run <- function(path,root) {
  prior <- fast_existing_run_orig(path,root)
  global <- file.path(prior$path,"global")
  receipt <- file.path(global,"numerical_completion.rds")
  if (!file.exists(receipt)) return(prior)
  required <- c("fast_manifest.rds","warm_boundary_area_day_posterior_draws.rds",
                "cdd_global_posterior_draws.rds","global_run_metadata.rds")
  if (!all(file.exists(file.path(global,required))))
    stop("Numerical completion receipt exists but final outputs are missing: ",global)
  finish <- readRDS(receipt)
  manifest <- readRDS(file.path(global,"fast_manifest.rds"))
  fp <- manifest$fingerprint
  if (!is.list(finish) || !identical(finish$status,"COMPLETED") ||
      !is.character(fp) || length(fp)!=1L || is.na(fp) || !nzchar(fp) ||
      !identical(finish$fingerprint,fp))
    stop("Invalid stream numerical completion receipt/fingerprint: ",receipt)
  prior$status <- "COMPLETED"
  prior
}

fast_stream_write_run_status <- function(path,status,fingerprint=NULL) {
  if (!status %in% c("RUNNING","FAILED","COMPLETED")) stop("Invalid stream run status.")
  fast_atomic_save_rds(list(schema="ERA5-stream-run-status-v1",status=status,
    updated=Sys.time(),fingerprint=fingerprint),
    file.path(path,"standalone_run_status.rds"),replace=TRUE)
}

fast_stream_expand_parameter <- function(x,n_draws,name,positive=FALSE,nonnegative=FALSE) {
  if (!is.numeric(x) || !length(x) || !length(x) %in% c(1L,n_draws) ||
      any(!is.finite(x)) || (positive && any(x<=0)) || (nonnegative && any(x<0)))
    stop("Invalid CDD parameter ",name,": supply scalar or n_draws numeric values.")
  rep(as.numeric(x),length.out=n_draws)
}

fast_stream_validate_scenarios <- function(scenarios,n_draws) {
  if (!is.list(scenarios) || !length(scenarios)) stop("Declare at least one CDD scenario before ECS production.")
  ns <- vapply(scenarios,function(x) if(is.list(x) && is.character(x$name) && length(x$name)==1L)
    x$name else "",character(1))
  if (any(!nzchar(ns)) || anyDuplicated(ns)) stop("CDD scenario names must be unique and nonempty.")
  lapply(seq_along(scenarios),function(i) {
    s<-scenarios[[i]]
    list(name=s$name,
      theta0=fast_stream_expand_parameter(s$theta0,n_draws,paste0(s$name,".theta0")),
      slope=fast_stream_expand_parameter(s$slope,n_draws,paste0(s$name,".slope"),positive=TRUE),
      kappa=fast_stream_expand_parameter(s$kappa,n_draws,paste0(s$name,".kappa"),nonnegative=TRUE))
  })
}

# Worker initialization must reconstruct the local-only packed-row lookup.
fast_worker_initialize_orig <- fast_worker_initialize
fast_worker_initialize <- function(paths,cache_dir,backend,bundle,part_dir,
                                   checkpoint_dir,fingerprint,qc_history,qc_expected) {
  index <- readRDS(file.path(bundle$g$input_dir,"pack_index.rds"))
  fast_stream_set_context(index)
  fast_worker_initialize_orig(paths,cache_dir,backend,bundle,part_dir,
                              checkpoint_dir,fingerprint,qc_history,qc_expected)
}

# Strictly require the new yearly CDD draw arrays in committed checkpoints.
fast_read_checkpoint_orig <- fast_read_checkpoint
fast_read_checkpoint <- function(path,fingerprint,task,bundle) {
  cp <- fast_read_checkpoint_orig(path,fingerprint,task,bundle)
  ncells<-nrow(task$tile_points); ndraw<-as.integer(bundle$ndraw)
  ny<-length(bundle$files); ns<-length(bundle$g$cdd_scenarios)
  dim_expected <- list(cdd_year=c(ncells,ndraw,ny,2L,ns),
    cdd_all_weather=c(ncells,ndraw,ny,ns),
    cdd_valid_weather_days=c(ncells,ny),cdd_valid_supported_days=c(ncells,ny),
    cdd_history_finite_days=c(ncells,ny),cdd_primary_days=c(ncells,ny))
  for (nm in names(dim_expected)) {
    val<-cp[[nm]]
    if (!is.numeric(val) || !identical(as.integer(dim(val)),as.integer(dim_expected[[nm]])) ||
        any(!is.finite(val)) || any(val<0))
      stop("Invalid stream CDD checkpoint field: ",nm," in ",path)
  }
  if (!identical(cp$cdd_scenario_names,vapply(bundle$g$cdd_scenarios,`[[`,"", "name")) ||
      !identical(as.integer(cp$cdd_cell_ids),as.integer(task$tile_points$cell_id)) ||
      length(cp$cdd_cell_area_km2)!=ncells || any(!is.finite(cp$cdd_cell_area_km2) | cp$cdd_cell_area_km2<=0) ||
      !identical(as.numeric(cp$cdd_cell_area_km2),as.numeric(task$tile_points$cell_area_km2)))
    stop("CDD scenario/grid identity mismatch in ",path)
  planned <- matrix(rep(bundle$planned[seq_len(ny)],each=ncells),ncells,ny)
  count_fields <- c("cdd_valid_weather_days","cdd_valid_supported_days",
                    "cdd_history_finite_days","cdd_primary_days")
  for (nm in count_fields) {
    val <- cp[[nm]]
    if (any(val!=floor(val) | val>planned))
      stop("Invalid CDD integer/date count: ",nm," in ",path)
  }
  if (any(cp$cdd_valid_supported_days>cp$cdd_primary_days) ||
      any(cp$cdd_valid_supported_days>cp$cdd_valid_weather_days) ||
      any(cp$cdd_primary_days>cp$cdd_history_finite_days) ||
      any(cp$cdd_valid_supported_days < cp$cdd_primary_days+cp$cdd_valid_weather_days-planned))
    stop("Invalid CDD support/date counts in ",path)
  # Independent warm summaries must describe the same support in each saved period.
  for (period in unique(cp$summary_table$period)) {
    p <- match(period,bundle$periods)
    years <- if (p<=ny) p else seq_len(ny)
    tab <- cp$summary_table[cp$summary_table$period==period,,drop=FALSE]
    if (any(tab$n_valid_history_days!=rowSums(cp$cdd_history_finite_days[,years,drop=FALSE])) ||
        any(tab$n_common_support_days!=rowSums(cp$cdd_primary_days[,years,drop=FALSE])))
      stop("CDD and warm support counts differ in ",path)
  }
  # Validate per-year slices to avoid an extra full-size posterior-array copy.
  for (yy in seq_len(ny)) for (si in seq_len(ns)) {
    m0 <- matrix(cp$cdd_year[,,yy,1L,si],ncells,ndraw)
    m1 <- matrix(cp$cdd_year[,,yy,2L,si],ncells,ndraw)
    all_t <- matrix(cp$cdd_all_weather[,,yy,si],ncells,ndraw)
    empty_supported <- cp$cdd_valid_supported_days[,yy]==0
    empty_weather <- cp$cdd_valid_weather_days[,yy]==0
    same_support <- cp$cdd_valid_supported_days[,yy]==cp$cdd_valid_weather_days[,yy]
    tolerance <- 1e-10*pmax(all_t,1)
    if (any(m0[empty_supported,,drop=FALSE]!=0) || any(m1[empty_supported,,drop=FALSE]!=0) ||
        any(all_t[empty_weather,,drop=FALSE]!=0) || any(m0-all_t>tolerance) ||
        any(abs(m0[same_support,,drop=FALSE]-all_t[same_support,,drop=FALSE])>
            tolerance[same_support,,drop=FALSE]))
      stop("CDD accumulators contradict valid-day counts or baseline subset in ",path)
  }
  cp
}

# Elapsed (wall) seconds, never CPU seconds: stages are mutually exclusive.
# Timing sidecars are operational metadata and are not scientific checkpoints.
fast_stream_timer <- function() {
  e <- new.env(parent=emptyenv())
  e$started_at <- Sys.time()
  e$start <- e$last <- unname(proc.time()[["elapsed"]])
  e$stage <- "initialize_tile"
  stages <- c("initialize_tile","recover_checkpoint","validate_existing_checkpoint",
    "open_annual_history_files","read_history_chunk","evaluate_support","support_counts",
    "read_temperature_chunk","accumulate_all_weather_cdd","predict_supported_ab",
    "transform_warm_posterior","accumulate_supported_cdd","accumulate_warm_sums",
    "summarize_warm_tile","write_partition","checkpoint_write_validate","close_annual_files")
  e$seconds <- stats::setNames(numeric(length(stages)),stages)
  e
}

fast_stream_timer_switch <- function(timer,stage) {
  now <- unname(proc.time()[["elapsed"]])
  if (!timer$stage %in% names(timer$seconds)) timer$seconds[timer$stage] <- 0
  timer$seconds[timer$stage] <- timer$seconds[timer$stage] + now-timer$last
  timer$last <- now
  timer$stage <- stage
  invisible(NULL)
}

fast_stream_timer_snapshot <- function(timer) {
  fast_stream_timer_switch(timer,timer$stage)
  list(schema="ERA5-stream-timing-v1",started_at=timer$started_at,
    finished_at=Sys.time(),pid=Sys.getpid(),total_seconds=timer$last-timer$start,
    stage_seconds=timer$seconds)
}

# Adapted original tile accumulator; original warm-curve definitions unchanged.
fast_compute_tile <- function(task, bundle, part_dir, checkpoint_dir, fingerprint) {
  started <- Sys.time(); timer <- fast_stream_timer()
  context <- list(stage = "initialize_tile", tile = task$tile_id)
  set_context <- function(stage, ...) {
    fast_stream_timer_switch(timer,stage)
    context <<- c(list(stage = stage, tile = task$tile_id), list(...))
  }
  ncs <- list(); daily_state <- NULL
  on.exit(if (!is.null(daily_state)) fast_daily_cleanup(daily_state), add = TRUE)
  on.exit(invisible(lapply(ncs, function(nc) try(ncdf4::nc_close(nc), silent = TRUE))), add = TRUE)
  tryCatch({
    dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    tile_id <- as.integer(task$tile_id)
    checkpoint_path <- file.path(checkpoint_dir, sprintf("tile_%05d.rds", tile_id))
    set_context("recover_checkpoint")
    fast_daily_recover_warm(task, bundle, fingerprint, checkpoint_path, part_dir)
    if (file.exists(checkpoint_path)) {
      set_context("validate_existing_checkpoint")
      cp <- fast_read_checkpoint(checkpoint_path, fingerprint, task, bundle)
      filename <- fast_restore_partition(cp, part_dir)
      return(list(tile_id = tile_id, ncells = cp$ncells, elapsed_seconds = cp$elapsed_seconds,
        checkpoint = checkpoint_path, partition = filename, resumed = TRUE,
        timing = fast_stream_timer_snapshot(timer)))
    }
    set_context("initialize_tile")
    g <- bundle$g; fixed_set <- bundle$fixed_set; ndraw <- as.integer(bundle$ndraw)
    periods <- bundle$periods; planned <- bundle$planned; np <- length(periods)
    nsupport <- length(g$support_fraction_sensitivity); tile_points <- task$tile_points
    ncells <- nrow(tile_points); keep <- as.integer(task$keep); set_bounds <- bundle$set_bounds
    if (length(tile_id) != 1L || is.na(tile_id) || tile_id < 1L || ncells < 1L ||
        length(keep) != ncells || anyNA(keep) || any(keep < 1L) || anyDuplicated(keep) ||
        any(keep > length(task$lon_indices) * length(task$lat_indices))) stop("Invalid spatial tile.")
    if (np != length(bundle$files) + 1L || length(planned) != np ||
        any(!is.finite(planned) | planned <= 0) || ndraw < 1L ||
        length(fixed_set) != 1L || !identical(names(fixed_set), "P50") || !is.finite(fixed_set) ||
        length(set_bounds) != 2L || any(!is.finite(set_bounds)) || set_bounds[1L] >= set_bounds[2L] ||
        nsupport < 1L || any(!is.finite(g$support_fraction_sensitivity)) ||
        any(g$support_fraction_sensitivity < 0 | g$support_fraction_sensitivity > 1))
      stop("Invalid periods, draws, P50, SET bounds, or support thresholds.")
    if (any(!is.finite(tile_points$area_weight) | tile_points$area_weight <= 0)) stop("Invalid area weights.")
    metrics <- c("delta_M", "delta_N", "delta_U", "delta_TSV_P50", "U_M1", "U_M0",
                 "U_pair_in_SET_range", "absolute_delta_U")
    quality_names <- c("finite", "marginal", "joint", "primary", "B21_outside", "TV_lag1_outside",
                       "TV_lag2_outside", "TV_lag3_outside")
    aggregate_sums <- array(0, c(np, 8L, ndraw))
    aggregate_denominator <- numeric(np); aggregate_cell_counts <- integer(np)
    support_sums <- array(0, c(np, 8L, nsupport, ndraw))
    support_denominator <- matrix(0, np, nsupport); support_cell_counts <- matrix(0L, np, nsupport)
    sums <- array(0, c(ncells, ndraw, np, 8L))
    quality <- array(0L, c(ncells, np, 8L), dimnames = list(NULL, periods, quality_names))
    # STREAM ADDITION: yearly joint draw accumulators, never daily posterior cache.
    ny <- length(bundle$files)
    scenarios <- fast_stream_validate_scenarios(g$cdd_scenarios, ndraw)
    ns <- length(scenarios)
    cdd_year <- array(0, c(ncells, ndraw, ny, 2L, ns)) # CDD_M0 and CDD_M1
    cdd_all_weather <- array(0, c(ncells, ndraw, ny, ns)) # baseline on ALL finite t2m dates
    cdd_valid_weather_days <- matrix(0L, ncells, ny)
    cdd_valid_supported_days <- matrix(0L, ncells, ny)
    if (isTRUE(bundle$daily_cache$enabled))
      daily_state <- fast_daily_open(task, bundle, fingerprint, bundle$daily_cache$dir)
    set_context("open_annual_history_files")
    for (i in seq_along(bundle$files)) ncs[[i]] <- ncdf4::nc_open(bundle$files[[i]]$path)
    for (yy in seq_along(ncs)) {
      fi <- bundle$files[[yy]]; nt <- length(fi$dates)
      for (ds in seq.int(1L, nt, by = g$day_chunk)) {
        nday <- min(g$day_chunk, nt - ds + 1L)
        info <- list(year = fi$year, start_day = ds, n_days = nday, n_cells = ncells, n_draws = ndraw)
        set_context("read_history_chunk", chunk = info)
        mats <- lapply(era5_required_variables, function(v)
          read_rectangular_history(ncs[[yy]], v, task$lon_indices, task$lat_indices,
            ds, nday, fi$offsets[v])[keep, , drop = FALSE])
        history <- as.data.frame(setNames(lapply(mats, as.vector), era5_required_variables))
        cells <- rep(seq_len(ncells), times = nday)
        set_context("evaluate_support", chunk = info)
        support <- evaluate_support(bundle$support_analysis, history, joint = isTRUE(g$joint_support))
        set_context("support_counts", chunk = info)
        needed <- c("finite", "marginal_inside", "joint_inside", quality_names[5:8])
        if (!is.data.frame(support) || nrow(support) != nrow(history) || !all(needed %in% names(support)))
          stop("Support evaluator returned invalid rows or columns.")
        finite <- as.logical(support$finite); finite[is.na(finite)] <- FALSE
        marginal <- finite & !is.na(support$marginal_inside) & support$marginal_inside
        joint <- finite & !is.na(support$joint_inside) & support$joint_inside
        primary <- marginal & (if (isTRUE(g$joint_support)) joint else TRUE)
        flags <- cbind(finite, marginal, joint, primary, as.matrix(support[, quality_names[5:8], drop = FALSE]))
        flags[is.na(flags)] <- FALSE; flags[!finite, 5:8] <- FALSE
        for (q in seq_along(quality_names)) {
          counts <- tabulate(cells[as.logical(flags[, q])], nbins = ncells)
          quality[, yy, q] <- quality[, yy, q] + counts
          quality[, np, q] <- quality[, np, q] + counts
        }
        chunk_dates <- as.Date(fi$dates[seq.int(ds, length.out = nday)])
        temperature <- NULL
        # STREAM ADDITION: same-day Ta_C is already exactly matched in the
        # local packed annual input. Always accumulate all-weather baseline,
        # including dates for which the model history is unsupported.
        set_context("read_temperature_chunk", chunk = info)
        ta_mat <- read_rectangular_history(ncs[[yy]], "Ta_C", task$lon_indices,
          task$lat_indices, ds, nday, 0)[keep, , drop = FALSE]
        ta_all <- as.vector(ta_mat)
        if (length(ta_all) != length(primary)) stop("Packed temperature/history misaligned.")
        finite_t <- is.finite(ta_all)
        cdd_valid_weather_days[, yy] <- cdd_valid_weather_days[, yy] +
          tabulate(cells[finite_t], nbins = ncells)
        set_context("accumulate_all_weather_cdd", chunk = info)
        if (any(finite_t)) for (si in seq_len(ns)) {
          base_all <- pmax(outer(ta_all[finite_t], scenarios[[si]]$theta0, "-"), 0)
          cdd_all_weather[, , yy, si] <- matrix(cdd_all_weather[, , yy, si], ncells, ndraw) +
            sum_by_cell(base_all, rep(TRUE, sum(finite_t)), cells[finite_t], ncells)
        }
        if (!is.null(daily_state) && identical(bundle$temperature$status, "embedded")) {
          set_context("read_same_day_temperature", chunk = info)
          temperature <- fast_daily_temperature_chunk(bundle$temperature, task, chunk_dates)
        }
        # Record ALL dates/masks, including entirely unsupported chunks.
        if (!any(primary)) {
          if (!is.null(daily_state)) fast_daily_append(daily_state, chunk_dates,
            flags[, 1:4, drop = FALSE], NULL, temperature = temperature)
          next
        }
        set_context("predict_supported_ab", chunk = info, supported_rows = sum(primary))
        ab <- fast_predict_ab(bundle$cache, history[primary, , drop = FALSE])
        expected <- c(sum(primary), ndraw)
        if (!identical(as.integer(dim(ab$A1)), as.integer(expected)) ||
            !identical(as.integer(dim(ab$B1)), as.integer(expected)) ||
            length(ab$A0) != ndraw || length(ab$B0) != ndraw ||
            any(!is.finite(ab$A1)) || any(!is.finite(ab$B1)) ||
            any(!is.finite(ab$A0)) || any(!is.finite(ab$B0)) ||
            any(ab$A1 <= 0) || any(ab$A0 <= 0)) stop("Invalid positive-slope paired prediction output.")
        set_context("transform_warm_posterior", chunk = info)
        values <- fast_warm_metrics(ab, fixed_set, set_bounds)
        # STREAM ADDITION: nonlinear CDD transform BEFORE temporal aggregation.
        # Never average delta_U before positive-part truncation.
        set_context("accumulate_supported_cdd", chunk = info)
        supported_cells <- cells[primary]
        supported_t <- ta_all[primary]
        valid_t <- is.finite(supported_t)
        cdd_valid_supported_days[, yy] <- cdd_valid_supported_days[, yy] +
          tabulate(supported_cells[valid_t], nbins = ncells)
        if (any(valid_t)) for (si in seq_len(ns)) {
          scenario <- scenarios[[si]]
          raw <- outer(supported_t[valid_t], scenario$theta0, "-")
          v0 <- pmax(raw, 0)
          shift <- sweep(values$delta_U[valid_t, , drop = FALSE], 2L,
                         scenario$kappa / scenario$slope, "*")
          v1 <- pmax(raw - shift, 0)
          indices <- supported_cells[valid_t]
          cdd_year[, , yy, 1L, si] <- matrix(cdd_year[, , yy, 1L, si], ncells, ndraw) +
            sum_by_cell(v0, rep(TRUE, length(indices)), indices, ncells)
          cdd_year[, , yy, 2L, si] <- matrix(cdd_year[, , yy, 2L, si], ncells, ndraw) +
            sum_by_cell(v1, rep(TRUE, length(indices)), indices, ncells)
        }
        if (!is.null(daily_state)) {
          set_context("write_daily_energy_cache", chunk = info)
          fast_daily_append(daily_state, chunk_dates, flags[, 1:4, drop = FALSE],
            values$delta_U, temperature = temperature)
        }
        set_context("accumulate_warm_sums", chunk = info)
        selected <- rep(TRUE, sum(primary)); selected_cells <- cells[primary]
        for (metric in seq_along(metrics)) {
          if (any(!is.finite(values[[metric]]))) stop("Nonfinite transformed metric; coverage must not vary by draw.")
          increment <- sum_by_cell(values[[metric]], selected, selected_cells, ncells)
          sums[, , yy, metric] <- matrix(sums[, , yy, metric], ncells, ndraw) + increment
          sums[, , np, metric] <- matrix(sums[, , np, metric], ncells, ndraw) + increment
        }
      }
    }
    if (!is.null(bundle$capacity_supported_rows_by_tile)) {
      expected_supported <- bundle$capacity_supported_rows_by_tile[tile_id]
      actual_supported <- sum(quality[, np, "primary"])
      if (length(expected_supported) != 1L || !is.finite(expected_supported) ||
          actual_supported != expected_supported)
        stop("Exact capacity scan and production support counts disagree in tile ", tile_id,
          ": expected ", expected_supported, "; observed ", actual_supported,
          ". Halt before committing this tile; inspect input stability.")
    }
    set_context("summarize_warm_tile")
    rows <- list(); area <- tile_points$area_weight
    for (p in seq_along(periods)) {
      finite_n <- quality[, p, "finite"]; denominator <- quality[, p, "primary"]
      date_coverage <- finite_n / planned[p]; fraction <- denominator / planned[p]
      eligible <- date_coverage >= g$min_date_coverage & denominator >= g$min_supported_days &
                  fraction >= g$min_supported_fraction
      aggregate_denominator[p] <- sum(area[eligible] * denominator[eligible])
      aggregate_cell_counts[p] <- sum(eligible)
      save_grid <- p == np || isTRUE(g$save_annual_grid)
      qtables <- if (save_grid) vector("list", length(metrics)) else NULL
      for (metric in seq_along(metrics)) {
        mat <- matrix(sums[, , p, metric], ncells, ndraw)
        if (save_grid) qtables[[metric]] <- posterior_rows(mat, denominator)
        aggregate_sums[p, metric, ] <- colSums(mat[eligible, , drop = FALSE] * area[eligible])
        for (cut in seq_len(nsupport)) {
          chosen <- date_coverage >= g$min_date_coverage & denominator >= g$min_supported_days &
                    fraction >= g$support_fraction_sensitivity[cut]
          support_sums[p, metric, cut, ] <- colSums(mat[chosen, , drop = FALSE] * area[chosen])
          if (metric == 1L) {
            support_denominator[p, cut] <- sum(area[chosen] * denominator[chosen])
            support_cell_counts[p, cut] <- sum(chosen)
          }
        }
      }
      if (!save_grid) next
      row <- data.frame(cell_id = tile_points$cell_id, lon = tile_points$lon_plot, lat = tile_points$lat,
        land_fraction = tile_points$land_fraction, area_weight = area, period = periods[p],
        fixed_set_name = "P50", fixed_set = unname(fixed_set),
        fixed_SET_applies_only_to = "delta_TSV_P50; warm boundaries do not depend on P50",
        n_expected_days = planned[p], n_valid_history_days = finite_n,
        n_marginal_support_days = quality[, p, "marginal"], n_joint_support_days = quality[, p, "joint"],
        n_common_support_days = denominator, n_estimation_days = denominator,
        date_coverage = date_coverage, common_support_fraction = fraction,
        common_support_fraction_of_valid = ifelse(finite_n > 0, denominator / finite_n, NA_real_),
        eligible_for_area_day_summary = eligible, eligible_for_supported_day_map = eligible,
        required_common_support_fraction = g$min_supported_fraction,
        SET_range_lower = set_bounds[1L], SET_range_upper = set_bounds[2L],
        pilot_not_global = !is.null(g$max_cells))
      for (q in quality_names[5:8]) row[[paste0(q, "_fraction_of_valid")]] <-
        ifelse(finite_n > 0, quality[, p, q] / finite_n, NA_real_)
      for (metric in seq_along(metrics)) for (field in names(qtables[[metric]]))
        row[[paste0(metrics[metric], "_", field)]] <- qtables[[metric]][[field]]
      for (metric in metrics[1:4]) row[[paste0(metric, "_CrI_excludes_zero")]] <-
        row[[paste0(metric, "_lower")]] > 0 | row[[paste0(metric, "_upper")]] < 0
      rows[[length(rows) + 1L]] <- row
    }
    table <- do.call(rbind, rows); rownames(table) <- NULL
    joint_denominator <- as.numeric(quality[, np, "primary"])
    joint_draw_means <- array(NA_real_, c(ncells, ndraw, 4L), dimnames = list(NULL, NULL, metrics[1:4]))
    valid <- joint_denominator > 0
    for (metric in 1:4) if (any(valid)) {
      mat <- matrix(sums[, , np, metric], ncells, ndraw)
      joint_draw_means[valid, , metric] <- mat[valid, , drop = FALSE] / joint_denominator[valid]
    }
    # STREAM ADDITION: preserve annual eligibility/count denominators even
    # though warm02's default grid summary saves the whole decade only.
    cdd_history_finite_days <- matrix(quality[, seq_len(ny), "finite"], ncells, ny)
    cdd_primary_days <- matrix(quality[, seq_len(ny), "primary"], ncells, ny)
    if (any(cdd_valid_supported_days > cdd_primary_days))
      stop("CDD date coverage exceeds warm-model support.")
    set_context("write_partition")
    filename <- file.path(part_dir, sprintf("tile_%05d.rds", tile_id)); csv <- sub("[.]rds$", ".csv", filename)
    cp <- list(schema = "ERA5-warm02-tile-v1", fingerprint = fingerprint, tile_id = tile_id,
      task_signature = fast_tile_signature(task), ncells = ncells, ndraw = ndraw,
      periods = periods, planned = planned, metric_names = metrics,
      aggregate_sums = aggregate_sums, aggregate_denominator = aggregate_denominator,
      aggregate_cell_counts = aggregate_cell_counts, support_sums = support_sums,
      support_denominator = support_denominator, support_cell_counts = support_cell_counts,
      summary_table = table, joint_draw_means = joint_draw_means,
      joint_draw_denominator = joint_denominator, joint_draw_cell_ids = as.integer(tile_points$cell_id),
      joint_draw_note = "Full-period supported-date mean per cell and fixed paired draw; delta_U = delta_M + delta_N.",
      # Added fields do not change original warm02 checkpoint semantics.
      cdd_year = cdd_year, cdd_all_weather = cdd_all_weather,
      cdd_valid_weather_days = cdd_valid_weather_days,
      cdd_valid_supported_days = cdd_valid_supported_days,
      cdd_history_finite_days = cdd_history_finite_days,
      cdd_primary_days = cdd_primary_days,
      cdd_cell_area_km2 = as.numeric(tile_points$cell_area_km2),
      cdd_cell_ids = as.integer(tile_points$cell_id),
      cdd_scenario_names = vapply(scenarios, `[[`, "", "name"),
      cdd_units = "degree-C days per grid-cell per year, original daily Celsius t2m; not energy",
      history_row_counts = c(all = ncells * planned[np], finite = sum(quality[, np, "finite"]),
        predicted = sum(quality[, np, "primary"]),
        finite_skipped_before_prediction = sum(quality[, np, "finite"] - quality[, np, "primary"])),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), created = Sys.time())
    if (!is.null(daily_state)) {
      set_context("commit_daily_energy_cache")
      daily_state$recovery_checkpoint <- cp
      cp$daily_cache <- fast_daily_finish(daily_state)
      cp$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    }
    fast_atomic_save_rds(table, filename, replace = TRUE)
    fast_atomic_publish(csv, function(temporary) utils::write.csv(table, temporary, row.names = FALSE), replace = TRUE)
    set_context("checkpoint_write_validate")
    fast_atomic_publish(checkpoint_path, function(temporary) {
      saveRDS(cp, temporary, compress = "gzip")
      invisible(fast_read_checkpoint(temporary, fingerprint, task, bundle))
    })
    set_context("close_annual_files")
    invisible(lapply(ncs,ncdf4::nc_close)); ncs <- list()
    timing <- fast_stream_timer_snapshot(timer)
    timing$tile_id <- tile_id; timing$fingerprint <- fingerprint
    timing$ncells <- ncells; timing$history_row_counts <- cp$history_row_counts
    timing$checkpoint_bytes <- unname(file.info(checkpoint_path)$size)
    context <- list(stage="write_timing_report",tile=tile_id)
    timing_dir <- file.path(dirname(checkpoint_dir),"tile_timings")
    dir.create(timing_dir,recursive=TRUE,showWarnings=FALSE)
    fast_atomic_save_rds(timing,file.path(timing_dir,sprintf("tile_%05d.rds",tile_id)),replace=TRUE)
    list(tile_id = tile_id, ncells = ncells, elapsed_seconds = timing$total_seconds,
      checkpoint = checkpoint_path, partition = filename, resumed = FALSE, timing = timing)
  }, error = function(e) {
    report <- list(schema = "ERA5-warm02-tile-error-v1", fingerprint = fingerprint,
      tile_id = task$tile_id, time = Sys.time(), context = context, error = conditionMessage(e),
      timing = fast_stream_timer_snapshot(timer),
      call = paste(deparse(conditionCall(e), width.cutoff = 120L), collapse = " "))
    path <- file.path(checkpoint_dir, sprintf("tile_%05d_error.rds", as.integer(task$tile_id)))
    try(fast_atomic_save_rds(report, path, replace = TRUE), silent = TRUE)
    stop("Tile ", task$tile_id, " failed at ", context$stage, ": ", conditionMessage(e), call. = FALSE)
  })
}


# Shared model-free validation for the finalizer and 04_verify_outputs.R.
# An all-empty result is valid missing information, not a failed delta identity.
fast_stream_validate_global_cdd <- function(x) {
  fail <- function(why) stop("Invalid global CDD output: ",why,call.=FALSE)
  if (!identical(x$schema,"ERA5-stream-cdd-global-v2"))
    fail("expected v2 with explicit support status; old v1 results require review")
  np <- length(x$periods); ns <- length(x$scenario_names); nd <- nrow(x$paired_draw_ids)
  if (!np || !ns || is.null(nd) || nd<1L ||
      !identical(x$metric_names,c("CDD_M0","CDD_M1","delta_CDD"))) fail("period/scenario/draw/metric identity")
  shapes <- list(area_km2_degree_days=c(np,ns,3L,nd),
    mean_CDD_per_area_weighted_valid_supported_day=c(np,ns,3L,nd),
    baseline_all_weather_area_km2_degree_days=c(np,ns,nd),
    baseline_all_weather_mean_per_area_weighted_valid_weather_day=c(np,ns,nd),
    area_km2_valid_supported_days=c(np,ns),area_km2_valid_all_weather_days=c(np,ns),
    eligible_grid_cells=c(np,ns))
  for (nm in names(shapes))
    if (!is.numeric(x[[nm]]) || !identical(as.integer(dim(x[[nm]])),as.integer(shapes[[nm]])))
      fail(paste("shape/type of",nm))
  if (!is.character(x$status) || !identical(as.integer(dim(x$status)),c(np,ns)) ||
      anyNA(x$status) || any(!x$status %in% c("ok","no_eligible_support"))) fail("support status")
  cells <- x$eligible_grid_cells
  den <- x$area_km2_valid_supported_days
  den_all <- x$area_km2_valid_all_weather_days
  if (any(!is.finite(cells) | cells<0 | cells!=floor(cells)) ||
      any(!is.finite(den) | den<0) || any(!is.finite(den_all) | den_all<den))
    fail("counts/denominators")
  max_delta_error <- 0
  for (p in seq_len(np)) for (si in seq_len(ns)) {
    dd <- matrix(x$area_km2_degree_days[p,si,,],3L,nd)
    avg <- matrix(x$mean_CDD_per_area_weighted_valid_supported_day[p,si,,],3L,nd)
    base <- x$baseline_all_weather_area_km2_degree_days[p,si,]
    base_avg <- x$baseline_all_weather_mean_per_area_weighted_valid_weather_day[p,si,]
    values <- c(dd,avg,base,base_avg)
    if (cells[p,si]==0) {
      if (x$status[p,si]!="no_eligible_support" || den[p,si]!=0 || den_all[p,si]!=0 ||
          !all(is.na(values)) || any(is.nan(values))) fail("empty support must have NA draws and zero counts")
      next
    }
    if (x$status[p,si]!="ok" || den[p,si]<=0 || den_all[p,si]<=0 ||
        any(!is.finite(values)) || any(dd[1:2,,drop=FALSE]<0) || any(base<0))
      fail("eligible support must have finite draws and positive denominators")
    err <- max(abs(dd[3L,]-(dd[2L,]-dd[1L,])))
    max_delta_error <- max(max_delta_error,err)
    if (err>1e-7) fail("draw-wise delta identity")
    if (any(abs(avg-dd/den[p,si])>1e-10*pmax(1,abs(avg))) ||
        any(abs(base_avg-base/den_all[p,si])>1e-10*pmax(1,abs(base_avg))) ||
        any(dd[1L,]-base>1e-10*pmax(1,base))) fail("normalization/baseline subset")
  }
  invisible(list(max_delta_error=max_delta_error,
    supported_period_scenarios=sum(cells>0),empty_period_scenarios=sum(cells==0)))
}

# Run the original warm finalizer unchanged, then append the CDD posterior
# summaries. Global CDD intervals are calculated after aggregating SAME-draw
# grid values; marginal per-grid quantiles are NEVER summed.
fast_finalize_global_orig <- fast_finalize_global
fast_finalize_global <- function(config,analysis,bundle,inspected,land,ilon,ilat,
                                 totals,partition_files,out_dir,execution) {
  result <- fast_finalize_global_orig(config,analysis,bundle,inspected,land,ilon,ilat,
                                      totals,partition_files,out_dir,execution)
  scenarios <- fast_stream_validate_scenarios(bundle$g$cdd_scenarios,bundle$ndraw)
  np <- length(bundle$periods); ny <- length(bundle$files)
  nd <- as.integer(bundle$ndraw); ns <- length(scenarios)
  metric_names <- c("CDD_M0", "CDD_M1", "delta_CDD")
  glob <- array(0,c(np,ns,3L,nd),dimnames=list(
    period=bundle$periods,scenario=vapply(scenarios,`[[`,"","name"),
    metric=metric_names,paired_draw_position=as.character(seq_len(nd))))
  base_all <- array(0,c(np,ns,nd))
  denominator <- matrix(0,np,ns)
  denominator_all_t <- matrix(0,np,ns)
  eligible_cells <- matrix(0L,np,ns)
  seen_ids <- integer()
  for (f in partition_files) {
    cpfile <- file.path(dirname(dirname(f)),"checkpoints",basename(f))
    if (!file.exists(cpfile)) stop("Missing committed CDD checkpoint for partition: ",f)
    cp <- readRDS(cpfile)
    if (!identical(cp$fingerprint,execution$fingerprint) || cp$ndraw!=nd ||
        !identical(cp$cdd_scenario_names,vapply(scenarios,`[[`,"","name")))
      stop("CDD checkpoint identity mismatch: ",cpfile)
    cell_ids <- cp$cdd_cell_ids
    if (anyDuplicated(cell_ids) || any(cell_ids %in% seen_ids))
      stop("CDD cell IDs overlap across committed tiles.")
    seen_ids <- c(seen_ids,cell_ids)
    ncells <- length(cell_ids)
    if (!identical(as.integer(dim(cp$cdd_year)),c(ncells,nd,ny,2L,ns)))
      stop("CDD checkpoint shape mismatch.")
    area <- as.numeric(cp$cdd_cell_area_km2)
    if (length(area)!=ncells || any(!is.finite(area)|area<=0)) stop("Invalid physical grid-cell area.")
    for (p in seq_len(np)) {
      years <- if (p<=ny) p else seq_len(ny)
      finite_hist <- rowSums(matrix(cp$cdd_history_finite_days[,years,drop=FALSE],ncells))
      primary <- rowSums(matrix(cp$cdd_primary_days[,years,drop=FALSE],ncells))
      valid_supported_t <- rowSums(matrix(cp$cdd_valid_supported_days[,years,drop=FALSE],ncells))
      valid_all_t <- rowSums(matrix(cp$cdd_valid_weather_days[,years,drop=FALSE],ncells))
      expected <- bundle$planned[p]
      selected <- finite_hist/expected >= bundle$g$min_date_coverage &
        primary >= bundle$g$min_supported_days &
        primary/expected >= bundle$g$min_supported_fraction & valid_supported_t>0
      if (!any(selected)) next
      ix <- which(selected)
      for (si in seq_len(ns)) {
        denominator[p,si] <- denominator[p,si]+sum(area[ix]*valid_supported_t[ix])
        denominator_all_t[p,si] <- denominator_all_t[p,si]+sum(area[ix]*valid_all_t[ix])
        eligible_cells[p,si] <- eligible_cells[p,si]+length(ix)
        for (metric in 1:2) {
          accum <- matrix(0,ncells,nd)
          for (yy in years)
            accum <- accum + matrix(cp$cdd_year[,,yy,metric,si],ncells,nd)
          glob[p,si,metric,] <- glob[p,si,metric,] +
            colSums(accum[ix,,drop=FALSE]*area[ix])
        }
        baseline <- matrix(0,ncells,nd)
        for (yy in years)
          baseline <- baseline + matrix(cp$cdd_all_weather[,,yy,si],ncells,nd)
        base_all[p,si,] <- base_all[p,si,] + colSums(baseline[ix,,drop=FALSE]*area[ix])
      }
    }
  }
  expected_cells <- if (is.null(bundle$g$max_cells)) fast_stream_context()$index$n_cells
    else as.integer(execution$total_land_points)
  if (length(seen_ids)!=expected_cells || anyDuplicated(seen_ids))
    stop("Combined CDD output does not cover every selected run grid cell.")
  glob[,,3L,] <- glob[,,2L,] - glob[,,1L,]
  normalized <- glob
  baseline_normalized <- base_all
  status <- matrix("no_eligible_support",np,ns,
    dimnames=list(period=bundle$periods,scenario=vapply(scenarios,`[[`,"","name")))
  for (p in seq_len(np)) for (si in seq_len(ns)) {
    if (eligible_cells[p,si]>0L) {
      if (denominator[p,si]<=0 || denominator_all_t[p,si]<=0)
        stop("Eligible CDD cells require positive valid-day denominators.")
      status[p,si] <- "ok"
      normalized[p,si,,] <- glob[p,si,,] / denominator[p,si]
      baseline_normalized[p,si,] <- base_all[p,si,] / denominator_all_t[p,si]
    } else {
      # Internal empty sums start at zero; exported no-support results are missing.
      glob[p,si,,] <- normalized[p,si,,] <- NA_real_
      base_all[p,si,] <- baseline_normalized[p,si,] <- NA_real_
    }
  }
  paired <- data.frame(position=seq_len(nd),m1=analysis$draw_ids_m1,m0=analysis$draw_ids_m0)
  if (!identical(paired,result$aggregate_draws$paired_draw_ids))
    stop("Global warm/CDD posterior pairing differs.")
  output <- list(schema="ERA5-stream-cdd-global-v2",status=status,
    periods=bundle$periods,scenario_names=vapply(scenarios,`[[`,"","name"),
    metric_names=metric_names,
    area_km2_degree_days=glob,
    mean_CDD_per_area_weighted_valid_supported_day=normalized,
    baseline_all_weather_area_km2_degree_days=base_all,
    baseline_all_weather_mean_per_area_weighted_valid_weather_day=baseline_normalized,
    area_km2_valid_supported_days=denominator,
    area_km2_valid_all_weather_days=denominator_all_t,
    eligible_grid_cells=eligible_cells,
    paired_draw_ids=paired,
    units="area_km2_degree_days are full-cell area (km2) multiplied by CDD (degree-C days); NOT kWh or energy",
    support="Historical model primary support + finite same-day t2m; no model prediction on unsupported dates",
    note="Joint posterior per-draw aggregation. M0/M1 are separately fitted; draw coupling is the fixed product-posterior convention.",
    scenario_parameters=scenarios, pack_identity=fast_stream_context()$index$identity,
    pilot_not_global=!is.null(bundle$g$max_cells),
    fingerprint=execution$fingerprint)
  fast_stream_validate_global_cdd(output)
  path <- file.path(out_dir,"cdd_global_posterior_draws.rds")
  fast_atomic_save_rds(output,path,replace=TRUE)
  rows <- list(); ir <- 0L
  for (p in seq_len(np)) for (si in seq_len(ns)) for (metric in seq_along(metric_names)) {
    x <- as.numeric(normalized[p,si,metric,]); s <- warm_draw_summary(x)
    ir <- ir+1L
    rows[[ir]] <- data.frame(period=bundle$periods[p],scenario=scenarios[[si]]$name,
      metric=metric_names[metric],status=status[p,si],lower=s["lower"],median=s["median"],upper=s["upper"],
      prob_positive=s["prob_positive"],prob_negative=s["prob_negative"],
      eligible_grid_cells=eligible_cells[p,si],area_km2_valid_supported_days=denominator[p,si])
  }
  utils::write.csv(do.call(rbind,rows),file.path(out_dir,"cdd_global_posterior_summary.csv"),row.names=FALSE)
  # Separately report credible intervals for the UNNORMALIZED global sum.
  # Units: km2 * degree-C * days. This is a climatic indicator, not kWh.
  total_rows <- list(); it <- 0L
  for (p in seq_len(np)) for (si in seq_len(ns)) for (metric in seq_along(metric_names)) {
    tval <- as.numeric(glob[p,si,metric,]); st <- warm_draw_summary(tval)
    it <- it+1L
    total_rows[[it]] <- data.frame(period=bundle$periods[p],scenario=scenarios[[si]]$name,
      metric=metric_names[metric],status=status[p,si],lower=st["lower"],median=st["median"],upper=st["upper"],
      prob_positive=st["prob_positive"],prob_negative=st["prob_negative"],
      units="km2_degreeC_day",eligible_grid_cells=eligible_cells[p,si])
  }
  utils::write.csv(do.call(rbind,total_rows),
    file.path(out_dir,"cdd_global_total_posterior_summary.csv"),row.names=FALSE)
  # Correct warm02's cache-specific metadata. warm02's numerical outputs remain untouched.
  result$metadata$energy_statement <- paste(
    "Predeclared climatic CDD scenarios computed daily then accumulated per fixed paired posterior draw.",
    "No daily posterior cache is written. Stored annual per-cell draws can be reaggregated to regions/years.",
    "Changing nonlinear CDD definitions requires new model predictions; CDD is not electricity.")
  result$metadata$daily_cache_directory <- NA_character_
  result$metadata$same_day_temperature_status <- "prepacked Ta_C, complete source-month audit"
  result$metadata$stream_pack_identity <- fast_stream_context()$index$identity
  result$metadata$stream_cdd_file <- path
  result$metadata$stream_scenario_parameters <- scenarios
  saveRDS(result$metadata,file.path(out_dir,"global_run_metadata.rds"))
  dput(result$metadata[setdiff(names(result$metadata),
     c("input_audit","partition_files","draw_metadata"))],
     file=file.path(out_dir,"global_run_methods.txt"))
  result$cdd <- output
  result
}
