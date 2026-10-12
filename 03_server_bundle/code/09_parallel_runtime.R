# Continuous PSOCK scheduling for the streamed server runner.
# Source AFTER 02_stream_overrides.R. Scientific setup, checkpoint validation,
# fixed draw IDs and final tile-ID reduction follow warm02_reference.R unchanged.
# R parallel's sendCall/recvOneResult are internal APIs, isolated and exercised
# by 10_test_parallel_runtime.R on the installed R version before production.

fast_stream_scheduling <- function(fast) {
  method <- if (is.null(fast$scheduling)) "lpt" else fast$scheduling
  interval <- if (is.null(fast$progress_seconds)) 30 else fast$progress_seconds
  if (!is.character(method) || length(method)!=1L || is.na(method) ||
      !method %in% c("lpt","fifo")) stop("fast$scheduling must be 'lpt' or 'fifo'.")
  if (!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0)
    stop("fast$progress_seconds must be positive.")
  list(method=method,progress_seconds=interval,version="continuous-psock-v1")
}

fast_stream_order_tiles <- function(ids,cell_counts,method="lpt") {
  if (!is.numeric(ids) || anyNA(ids) || any(ids!=as.integer(ids)) ||
      any(ids<1L | ids>length(cell_counts)) || anyDuplicated(ids) ||
      any(!is.finite(cell_counts[ids]) | cell_counts[ids]<=0)) stop("Invalid scheduling IDs/costs.")
  ids <- as.integer(ids)
  if (!method %in% c("lpt","fifo")) stop("Unknown scheduling method.")
  if (method=="lpt") ids[order(-cell_counts[ids],ids)] else ids
}

fast_stream_append_csv <- function(row,path) {
  exists <- file.exists(path) && file.info(path)$size>0
  utils::write.table(row,path,sep=",",row.names=FALSE,col.names=!exists,
    append=exists,quote=TRUE,na="NA")
  invisible(NULL)
}

fast_stream_timing_row <- function(answer,session_id,phase,worker,rank) {
  t <- answer$timing
  if (!is.list(t) || !identical(t$schema,"ERA5-stream-timing-v1") ||
      any(!is.finite(t$stage_seconds) | t$stage_seconds<0) ||
      !is.finite(t$total_seconds) || t$total_seconds<0) stop("Missing/invalid worker timing.")
  stages <- fast_stream_timer()$seconds
  if (any(!names(t$stage_seconds) %in% names(stages))) stop("Unknown timing stage.")
  stages[names(t$stage_seconds)] <- t$stage_seconds
  rows <- t$history_row_counts
  row <- data.frame(session_id=session_id,phase=phase,tile_id=answer$tile_id,
    n_cells=answer$ncells,worker=worker,pid=t$pid,dispatch_rank=rank,
    resumed=isTRUE(answer$resumed),started_at=as.character(t$started_at),
    finished_at=as.character(t$finished_at),total_seconds=t$total_seconds,
    history_rows=if(is.null(rows)) NA_real_ else unname(rows["all"]),
    predicted_rows=if(is.null(rows)) NA_real_ else unname(rows["predicted"]),
    checkpoint_bytes=if(is.null(t$checkpoint_bytes)) NA_real_ else t$checkpoint_bytes)
  for (nm in names(stages)) row[[paste0(nm,"_seconds")]] <- unname(stages[nm])
  row
}

fast_stream_task_envelope <- function(task,fun) {
  tryCatch(list(ok=TRUE,value=fun(task)),error=function(e)
    list(ok=FALSE,tile_id=task$tile_id,message=conditionMessage(e)))
}

fast_stream_queue <- function(cl,tasks,fun,on_result=function(...) NULL) {
  n <- length(tasks); workers <- length(cl)
  if (!n) return(invisible(list()))
  if (!workers) stop("Continuous queue requires at least one worker.")
  ids <- vapply(tasks,function(x) as.integer(x$tile_id),integer(1))
  if (anyNA(ids) || anyDuplicated(ids)) stop("Queue task IDs must be unique.")
  send <- getFromNamespace("sendCall","parallel")
  receive <- getFromNamespace("recvOneResult","parallel")
  active <- rep(NA_integer_,workers); next_job <- 1L
  result <- vector("list",n); failures <- character(); failure_announced <- FALSE
  submit <- function(node) {
    j <- next_job; next_job <<- next_job+1L
    send(cl[[node]],fast_stream_task_envelope,list(tasks[[j]],fun),tag=j)
    active[node] <<- j
  }
  for (node in seq_len(min(workers,n))) submit(node)
  while (any(!is.na(active))) {
    received <- receive(cl)
    node <- received$node; j <- received$tag
    if (length(node)!=1L || is.na(node) || node<1L || node>workers ||
        length(j)!=1L || is.na(j) || is.na(active[node]) || active[node]!=j)
      stop("Continuous queue received an invalid worker/tag.")
    active[node] <- NA_integer_
    value <- received$value
    if (is.list(value) && isTRUE(value$ok)) {
      error <- tryCatch({
        on_result(value$value,tasks[[j]],node,j)
        NULL
      },error=function(e) conditionMessage(e))
      if (is.null(error)) result[j] <- list(value$value) else
        failures <- c(failures,paste0("tile ",ids[j]," callback: ",error))
    } else {
      msg <- if (is.list(value) && is.character(value$message)) value$message else as.character(value)
      failures <- c(failures,paste0("tile ",ids[j],": ",paste(msg,collapse=" ")))
    }
    if (length(failures) && !failure_announced && any(!is.na(active))) {
      message("Queue stopped dispatching after an error; waiting for running tiles to commit.")
      failure_announced <- TRUE
    }
    # On error, drain outstanding work and retain all successful commits. Never
    # silently drop an error, retry within this session, or enqueue fresh tasks.
    if (!length(failures) && next_job<=n) submit(node)
  }
  if (length(failures)) stop("Continuous queue failed; committed tiles can be resumed. ",
    paste(unique(failures),collapse="; "),call.=FALSE)
  invisible(result)
}

run_global_era5 <- function(config, analysis = NULL) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Install ncdf4 before global analysis.")
  g <- global_defaults(config); f <- fast_defaults(config)
  scheduling <- fast_stream_scheduling(config$fast)
  if (is.null(analysis)) analysis <- load_analysis(config)
  if (isTRUE(g$joint_support) && is.null(analysis$joint_reference))
    stop("Global joint support requires the original cached training reference.")
  fixed_set <- analysis$fixed_set["P50"]
  if (length(fixed_set) != 1L || !identical(names(fixed_set), "P50")) stop("This analysis requires P50 only.")
  set_bounds <- as.numeric(stats::quantile(analysis$data$SET_raw, c(.01,.99), na.rm=TRUE, type=7))
  if (length(set_bounds)!=2L || any(!is.finite(set_bounds)) || diff(set_bounds)<=0) stop("Invalid training SET bounds.")
  if (!length(fixed_set) || any(!is.finite(fixed_set))) stop("Invalid fixed SET values.")
  if (is.null(names(fixed_set))) names(fixed_set) <- paste0("SET", seq_along(fixed_set))
  available <- nrow(analysis$draws_m1)
  if (nrow(analysis$draws_m0) != available) stop("M1/M0 posterior row counts differ.")
  requested <- if (is.infinite(g$n_draws)) available else as.integer(g$n_draws)
  ndraw <- min(requested, available)
  if (requested > available) warning("Using all available paired draws: ", ndraw)
  if (ndraw < 50L) stop("At least 50 paired posterior draws are required.")
  positions <- seq_len(ndraw)
  analysis$draws_m1 <- analysis$draws_m1[positions, , drop = FALSE]
  analysis$draws_m0 <- analysis$draws_m0[positions, , drop = FALSE]
  analysis$draw_ids_m1 <- analysis$draw_ids_m1[positions]
  analysis$draw_ids_m0 <- analysis$draw_ids_m0[positions]
  analysis$draw_ids <- analysis$draw_ids_m1
  analysis$draw_metadata$global_draws_requested <- g$n_draws
  analysis$draw_metadata$global_draws_used <- ndraw
  analysis$draw_metadata$paired_draws_before_global_subset <- available
  analysis$draw_metadata$n_draws <- ndraw
  out_dir <- file.path(config$output_dir, "global")
  part_dir <- file.path(out_dir, "grid_partitions")
  checkpoint_dir <- file.path(out_dir, "checkpoints")
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  .set_analysis_stage("inspect_global_inputs", announce = TRUE)
  inspected <- inspect_history_files(g)
  daily_settings <- fast_daily_defaults(config)
  config$daily_cache <- daily_settings
  .set_analysis_stage("inspect_same_day_temperature_inputs", announce = TRUE)
  temperature <- fast_daily_inspect_temperature(config, inspected)
  daily_settings$dir <- file.path(out_dir, "energy_cache")
  dir.create(daily_settings$dir, recursive = TRUE, showWarnings = FALSE)
  first <- inspected$files[[1L]]
  # Enforce the native 0.25-degree spacing and stride=1.
  # ERA5 sample centres are retained; no interpolation or cell-area averaging is done.
  if (g$spatial_stride != 1L ||
      length(first$lon) < 2L || length(first$lat) < 2L ||
      any(!is.finite(diff(first$lon))) || any(!is.finite(diff(first$lat))) ||
      any(abs(abs(diff(first$lon)) - 0.25) > 1e-6) ||
      any(abs(abs(diff(first$lat)) - 0.25) > 1e-6))
    stop("This 0.25-degree script requires stride=1 on regular native 0.25-degree ",
         "longitude and latitude axes; check the ERA5 input coordinates.", call. = FALSE)
  ilon <- seq.int(1L, length(first$lon), by = g$spatial_stride)
  latitude_candidates <- which(first$lat >= g$lat_range[1L] & first$lat <= g$lat_range[2L])
  if (!length(latitude_candidates)) stop("No latitude coordinates in requested domain.")
  ilat <- latitude_candidates[seq.int(1L, length(latitude_candidates), by = g$spatial_stride)]
  .set_analysis_stage("land_mask", announce = TRUE)
  land <- load_land_fraction(g, first$lon, first$lat, ilon, ilat)
  grid <- land$points
  if (!any(grid$land)) stop("No land grid points selected.")
  if (!is.null(g$max_cells)) {
    chosen <- which(grid$land)
    chosen <- chosen[unique(round(seq(1, length(chosen), length.out = min(length(chosen), g$max_cells))))]
    grid$land <- seq_len(nrow(grid)) %in% chosen
    warning("PILOT subset: these are not full-domain global estimates.")
  }
  grid$selected_lon_index <- rep(seq_along(ilon), times = length(ilat))
  grid$selected_lat_index <- rep(seq_along(ilat), each = length(ilon))
  tasks <- fast_make_tasks(grid, ilon, ilat, g$tile_side)
  cell_counts <- vapply(tasks, function(x) nrow(x$tile_points), integer(1))
  planned <- vapply(inspected$files, function(x) length(x$dates), integer(1))
  periods <- c(as.character(g$years), paste0(min(g$years), "-", max(g$years)))
  planned <- c(planned, sum(planned))
  .set_analysis_stage("scan_exact_daily_cache_support", announce = TRUE)
  capacity <- fast_daily_support_preflight(daily_settings, tasks, inspected$files,
    analysis[c("pred", "training_history", "joint_reference")], g, inspected$audit, ndraw, out_dir)
  message("Formal global run: ", length(tasks), " spatial tiles; ", sum(cell_counts),
    " selected land points; ", sum(head(planned, -1L)), " daily dates; ", ndraw,
    " fixed paired posterior draws; lag1-3 retained. Only P50 TSV; M, positive N=1/S and U=M+N computed analytically.")
  .set_analysis_stage("initialize_fast_kernel", announce = TRUE)
  cache <- fast_make_cache(analysis, g$weight_scheme)
  cache_dir <- file.path(config$output_dir, "kernel_cache", paste0("R-", getRversion(), "-", R.version$arch))
  # A resume retains its original backend even if a compiler is installed later.
  manifest_path <- file.path(out_dir, "fast_manifest.rds")
  old_manifest <- if (file.exists(manifest_path)) readRDS(manifest_path) else NULL
  requested_backend <- f$backend
  if (!is.null(old_manifest) && identical(requested_backend, "auto")) requested_backend <- old_manifest$backend
  backend <- fast_init_backend(cache_dir, requested_backend)
  fast_kernel_QC(analysis, cache, config$output_dir)
  fast_boundary_preflight(cache, analysis, set_bounds, out_dir)
  bundle <- list(g = g, files = inspected$files, fixed_set = fixed_set, cache = cache,
    support_analysis = analysis[c("pred", "training_history", "joint_reference")],
    ndraw = ndraw, periods = periods, planned = planned, set_bounds = set_bounds,
    daily_cache = daily_settings, temperature = temperature,
    capacity_supported_rows_by_tile = capacity$supported_rows_by_tile)
  package_names <- c("dlnm", "FNN", "ncdf4", "maps", "Rcpp", "digest")
  packages <- vapply(package_names, function(p) if (requireNamespace(p, quietly = TRUE))
    as.character(utils::packageVersion(p)) else "unavailable", character(1))
  scientific_functions <- c("fast_compute_tile", "fast_make_cache", "fast_predict_ab", ".fast_r_kernel",
    ".fast_shift_probability_r", ".fast_cached_exp", "make_history_design", ".b21_basis", ".tv_basis",
    "history_support", "evaluate_support", "read_rectangular_history", "sum_by_cell", "posterior_rows", "fast_warm_metrics")
  all_symbols <- ls(envir = environment(), all.names = TRUE)
  daily_env <- environment(fast_daily_open)
  daily_symbols <- ls(envir = daily_env, all.names = TRUE)
  scientific_functions <- unique(c(scientific_functions,
    daily_symbols[grepl("^fast_daily_", daily_symbols) &
      vapply(daily_symbols, function(nm) is.function(get(nm, envir = daily_env)), logical(1))]))
  implementation <- lapply(scientific_functions, function(nm) deparse(body(get(nm, mode = "function"))))
  names(implementation) <- scientific_functions
  signature <- list(schema = "ERA5-warm02-run-v1", revision = ANALYSIS_REVISION,
    implementation = implementation, cpp_source = FAST_CPP_CODE,
    global_config = g, daily_cache = fast_daily_identity_settings(daily_settings), temperature = temperature,
    model_input_md5 = analysis$input_md5,
    model_cache = cache, support_reference = bundle$support_analysis, fixed_set = fixed_set,
    set_bounds = set_bounds, input_audit = inspected$audit, files = inspected$files, grid = grid,
    land_method = land$method, backend = backend, packages = packages, R = as.character(getRversion()))
  fingerprint <- fast_object_hash(signature)
  rm(signature)
  if (!is.null(old_manifest) && (!identical(old_manifest$schema, "ERA5-warm02-run-v1") ||
      !identical(old_manifest$fingerprint, fingerprint)))
    stop("Resume fingerprint differs (inputs, scientific settings, software, backend, or draw IDs). ",
         "Set config$fast$resume = 'new' for a separate run; committed results are unchanged.")
  if (is.null(old_manifest)) {
    if (length(list.files(checkpoint_dir, pattern = "^tile_[0-9]+[.]rds$")))
      stop("Checkpoints exist without an immutable fast manifest; cannot establish draw pairing.")
    manifest <- list(schema = "ERA5-warm02-run-v1", revision = ANALYSIS_REVISION,
      fingerprint = fingerprint, backend = backend, created = Sys.time(),
      total_tiles = length(tasks), total_land_points = sum(cell_counts), g = g,
      paired_draw_ids = data.frame(m1 = cache$draw_ids_m1, m0 = cache$draw_ids_m0),
      fixed_set = fixed_set, set_bounds = set_bounds, packages = packages, input_audit = inspected$audit,
      model_input_md5 = analysis$input_md5)
    fast_atomic_save_rds(manifest, manifest_path)
  }
  fast_daily_snapshot(config, analysis, bundle, tasks, fingerprint, out_dir)
  utils::write.csv(inspected$audit, file.path(out_dir, "input_audit.csv"), row.names = FALSE)
  utils::write.csv(grid[grid$land, ], file.path(out_dir, "selected_land_points.csv"), row.names = FALSE)
  checkpoint_paths <- file.path(checkpoint_dir, sprintf("tile_%05d.rds", seq_along(tasks)))
  done <- logical(length(tasks))
  .set_analysis_stage("validate_resume_checkpoints", announce = TRUE)
  for (i in seq_along(tasks))
    fast_daily_recover_warm(tasks[[i]], bundle, fingerprint, checkpoint_paths[i], part_dir)
  for (i in which(file.exists(checkpoint_paths))) {
    cp <- fast_read_checkpoint(checkpoint_paths[i], fingerprint, tasks[[i]], bundle)
    fast_restore_partition(cp, part_dir); done[i] <- TRUE
  }
  initial_done <- sum(done)
  message("Validated resume: ", initial_done, "/", length(tasks), " committed tiles; ",
          sum(cell_counts[done]), "/", sum(cell_counts), " land points.")
  start <- Sys.time(); benchmark <- data.frame(); progress <- list(); chosen_workers <- f$workers
  rate <- NA_real_; elapsed_new <- 0; cells_new <- 0L
  progress_interval <- scheduling$progress_seconds
  last_progress <- -Inf
  session_id <- paste0(format(Sys.time(),"%Y%m%dT%H%M%S"),"_",Sys.getpid(),"_",
    basename(tempfile("session_",tmpdir=out_dir)))
  dispatch_file <- file.path(out_dir,paste0("dispatch_",session_id,".csv"))
  timing_file <- file.path(out_dir,paste0("tile_timings_",session_id,".csv"))
  record_progress <- function(label, workers, seconds = NA_real_, new_cells = 0L) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    remaining <- sum(cell_counts[!done])
    eta <- if (remaining == 0) 0 else if (is.finite(rate) && rate > 0) remaining / rate else NA_real_
    row <- data.frame(time = as.character(Sys.time()), phase = label, backend = backend, workers = workers,
      completed_tiles = sum(done), total_tiles = length(tasks), completed_land_points = sum(cell_counts[done]),
      total_land_points = sum(cell_counts), seconds_this_wave = seconds, land_points_this_wave = new_cells,
      elapsed_this_session_seconds = elapsed, land_points_per_second = rate,
      estimated_remaining_compute_hours = eta / 3600,
      estimated_remaining_recent_hours = if (remaining==0) 0 else if (is.finite(seconds) && seconds>0 && new_cells>0) remaining/(new_cells/seconds)/3600 else NA_real_,
      elapsed_from_launch_hours = as.numeric(difftime(Sys.time(), .warm_session_started, units="hours")),
      projected_total_session_hours_excluding_finalize = as.numeric(difftime(Sys.time(), .warm_session_started, units="hours")) + eta/3600)
    progress[[length(progress) + 1L]] <<- row
    utils::write.csv(do.call(rbind, progress), file.path(out_dir, "progress.csv"), row.names = FALSE)
    last_progress <<- unname(proc.time()[["elapsed"]])
    message("Progress ", sum(done), "/", length(tasks), " tiles, ", sum(cell_counts[done]), "/", sum(cell_counts),
      " land points; ", workers, " workers; remaining compute ~",
      if (is.finite(eta)) sprintf("%.2f h", eta / 3600) else "pending benchmark",
      ". ETA includes measured queue compute/read/write time; excludes final aggregation/plots and may change with geography/I/O.")
    if (is.finite(eta) && remaining>0) {
      forecast <- row$projected_total_session_hours_excluding_finalize
      message(sprintf("Measured forecast for this session: %.2f h including elapsed setup, excluding final aggregation/plots.", forecast))
      if (is.finite(f$target_total_hours) && forecast > f$target_total_hours)
        message("Current measured speed exceeds the ", f$target_total_hours, " h planning target; calculation continues without reducing draws or dates.")
    }
    invisible(row)
  }
  cl <- NULL
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (any(!done)) {
    candidates <- if (isTRUE(f$autotune)) f$worker_candidates else f$workers
    physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (is.finite(physical) && physical > 0) candidates <- pmin(candidates, as.integer(physical))
    candidates <- sort(unique(pmin(candidates, sum(!done))))
    chosen_workers <- candidates[1L]
    .set_analysis_stage("start_PSOCK_workers", workers = max(candidates), announce = TRUE)
    cl <- parallel::makePSOCKcluster(max(candidates), rscript_args = "--vanilla",
      outfile = file.path(out_dir, "worker_log.txt"))
    source_env <- environment(run_global_era5)
    symbols <- ls(envir = source_env, all.names = TRUE)
    fast_functions <- symbols[grepl("^(fast_|[.]fast_)", symbols) &
      vapply(symbols, function(x) is.function(get(x, envir = source_env)), logical(1))]
    exports <- unique(c(fast_functions, "FAST_CPP_CODE", "era5_required_variables", "history_columns",
      ".history_frame", ".b21_basis", ".tv_basis", "make_history_design", ".make_joint_reference",
      "history_support", "evaluate_support", "normal_axis", "read_rectangular_history", "sum_by_cell",
      "posterior_rows", ".set_analysis_stage", ".analysis_error_details"))
    parallel::clusterExport(cl, exports, envir = source_env)
    qc_history <- .fast_qc_histories(analysis)$training
    qc_history <- qc_history[seq_len(min(2L, nrow(qc_history))), , drop = FALSE]
    qc_expected <- fast_predict_ab(cache, qc_history)[c("A1", "B1", "A0", "B0")]
    worker_audit <- vector("list", length(cl))
    # Serial initialization prevents competing sourceCpp cache writes.
    for (i in seq_along(cl)) worker_audit[[i]] <- parallel::clusterCall(cl[i], fast_worker_initialize,
      .libPaths(), cache_dir, backend, bundle, part_dir, checkpoint_dir, fingerprint, qc_history, qc_expected)[[1L]]
    fast_register_workers(config$output_dir, worker_audit)
    saveRDS(worker_audit, file.path(out_dir, "worker_QC.rds"))
    message("All ", length(cl), " workers passed prediction QC; backend=", backend,
            ". Complete production tiles will be used for speed measurement.")
    run_wave <- function(ids, workers, phase) {
      # IDs remain canonical; only the dispatch permutation changes.
      ids <- fast_stream_order_tiles(ids,cell_counts,scheduling$method)
      .set_analysis_stage(phase,tile_ids=ids,workers=workers,announce=TRUE)
      message("Continuous ",scheduling$method," queue: ",length(ids)," tiles (",
        sum(cell_counts[ids])," cells), ",workers," workers; no 2W batch barrier.")
      t0 <- unname(proc.time()[["elapsed"]]); completed_cells <- 0L
      fast_stream_append_csv(data.frame(session_id=session_id,phase=phase,
        dispatch_rank=seq_along(ids),tile_id=ids,n_cells=cell_counts[ids],
        cost_proxy=cell_counts[ids],method=scheduling$method),dispatch_file)
      on_result <- function(answer,task,worker,rank) {
        i <- as.integer(task$tile_id)
        if (!identical(as.integer(answer$tile_id),i) ||
            !identical(as.integer(answer$ncells),as.integer(cell_counts[i])) ||
            !file.exists(checkpoint_paths[i])) stop("Worker returned invalid tile identity/commit: ",i)
        done[i] <<- TRUE
        completed_cells <<- completed_cells+cell_counts[i]
        seconds <- unname(proc.time()[["elapsed"]])-t0
        rate <<- (cells_new+completed_cells)/max(elapsed_new+seconds,.001)
        fast_stream_append_csv(fast_stream_timing_row(answer,session_id,phase,worker,rank),timing_file)
        if (unname(proc.time()[["elapsed"]])-last_progress >= progress_interval || all(done))
          record_progress(phase,workers,seconds,completed_cells)
      }
      fast_stream_queue(cl[seq_len(workers)],tasks[ids],fast_worker_task,on_result)
      seconds <- unname(proc.time()[["elapsed"]])-t0
      list(seconds=seconds,cells=completed_cells,rate=completed_cells/max(seconds,.001))
    }
    benchmark_file <- file.path(out_dir, "worker_benchmark.rds")
    saved_benchmark <- if (file.exists(benchmark_file)) readRDS(benchmark_file) else NULL
    reuse_benchmark <- isTRUE(f$autotune) && !is.null(saved_benchmark) &&
      identical(saved_benchmark$fingerprint, fingerprint) && identical(saved_benchmark$candidates, candidates)
    if (reuse_benchmark) {
      benchmark <- saved_benchmark$table
      chosen_workers <- saved_benchmark$chosen_workers
      rate <- saved_benchmark$chosen_rate
      message("Reusing completed worker selection from this unchanged run: ", chosen_workers, " workers.")
    } else if (isTRUE(f$autotune) && length(candidates) > 1L && sum(!done) >= sum(candidates)) {
      for (workers in candidates) {
        remaining <- which(!done)
        full <- remaining[cell_counts[remaining] == max(cell_counts[remaining])]
        pool <- if (length(full) >= workers) full else remaining
        ids <- pool[unique(as.integer(round(seq(1, length(pool), length.out = workers))))]
        wave <- run_wave(ids, workers, "benchmark_production_tiles")
        row <- data.frame(workers = workers, tiles = length(ids), land_points = wave$cells,
          wall_seconds = wave$seconds, land_points_per_second = wave$rate)
        benchmark <- rbind(benchmark, row); rate <- wave$rate
        utils::write.csv(benchmark, file.path(out_dir, "worker_benchmark.csv"), row.names = FALSE)
        record_progress("benchmark", workers, wave$seconds, wave$cells)
      }
      best <- which.max(benchmark$land_points_per_second)
      chosen_workers <- benchmark$workers[best]; rate <- benchmark$land_points_per_second[best]
      fast_atomic_save_rds(list(fingerprint = fingerprint, candidates = candidates, table = benchmark,
        chosen_workers = chosen_workers, chosen_rate = rate), benchmark_file, replace = TRUE)
      message("Measured worker selection: ", chosen_workers,
        ". Benchmark tiles are retained as completed scientific output; speeds can vary by region.")
    }
    record_progress("production_start", chosen_workers)
    if (any(!done)) {
      wave <- run_wave(which(!done), chosen_workers, "parallel_global_projection")
      elapsed_new <- elapsed_new + wave$seconds; cells_new <- cells_new + wave$cells
      rate <- cells_new / max(elapsed_new, .001)
      record_progress("production", chosen_workers, wave$seconds, wave$cells)
    }
    parallel::stopCluster(cl); cl <- NULL
  }
  .set_analysis_stage("combine_committed_draw_sums", announce = TRUE)
  totals <- NULL; partition_files <- character(length(tasks))
  fields <- c("aggregate_sums", "aggregate_denominator", "aggregate_cell_counts",
              "support_sums", "support_denominator", "support_cell_counts")
  # Fixed tile-ID order makes final sums independent of worker scheduling.
  for (i in seq_along(tasks)) {
    cp <- fast_read_checkpoint(checkpoint_paths[i], fingerprint, tasks[[i]], bundle)
    partition_files[i] <- fast_restore_partition(cp, part_dir)
    if (is.null(totals)) totals <- cp[fields] else
      for (nm in fields) totals[[nm]] <- totals[[nm]] + cp[[nm]]
  }
  cache_completion <- fast_daily_finalize_inventory(bundle, tasks, fingerprint, out_dir)
  execution <- list(backend = backend, selected_workers = chosen_workers, benchmark = benchmark,
    scheduling = scheduling, scheduling_cost_proxy = "selected cells per tile; not measured runtimes",
    dispatch_log = dispatch_file, tile_timing_log = timing_file,
    daily_cache_completion = cache_completion,
    fingerprint = fingerprint, checkpoint_schema = "ERA5-warm02-tile-v1",
    resumed_tiles = initial_done, total_tiles = length(tasks), total_land_points = sum(cell_counts),
    elapsed_projection_session_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
    numerical_QC = "Original official posterior_epred checks plus fast kernel and every worker checked at tolerance <=1e-8.",
    resume = "Only validated per-draw checkpoints of this exact input/configuration/draw fingerprint are combined.",
    reproducibility = "No worker RNG; fixed draw IDs; final tile-ID summation order.")
  fast_finalize_global(config, analysis, bundle, inspected, land, ilon, ilat,
                       totals, partition_files, out_dir, execution)
}
