# Continuous PSOCK scheduling for the streamed server runner.
# Source AFTER 02_stream_overrides.R. Scientific setup, checkpoint validation,
# fixed draw IDs and final tile-ID reduction follow warm02_reference.R unchanged.
# R parallel's sendCall/recvOneResult are internal APIs, isolated and exercised
# by 10_test_parallel_runtime.R on the installed R version before production.
sys.source(file.path(stream_this_dir(),"11_memory_resources.R"),envir=environment())

fast_stream_scheduling <- function(fast) {
  method <- if (is.null(fast$scheduling)) "lpt" else fast$scheduling
  interval <- if (is.null(fast$progress_seconds)) 30 else fast$progress_seconds
  if (!is.character(method) || length(method)!=1L || is.na(method) ||
      !method %in% c("lpt","fifo","memory_lpt")) stop("fast$scheduling must be 'lpt', 'fifo' or 'memory_lpt'.")
  if (!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0)
    stop("fast$progress_seconds must be positive.")
  list(method=method,progress_seconds=interval,version="memory-psock-v2.1",
       memory=fast_stream_memory_settings(fast$memory))
}

fast_stream_order_tiles <- function(ids,cell_counts,method="lpt") {
  if (!is.numeric(ids) || anyNA(ids) || any(ids!=as.integer(ids)) ||
      any(ids<1L | ids>length(cell_counts)) || anyDuplicated(ids) ||
      any(!is.finite(cell_counts[ids]) | cell_counts[ids]<=0)) stop("Invalid scheduling IDs/costs.")
  ids <- as.integer(ids)
  if (!method %in% c("lpt","fifo","memory_lpt")) stop("Unknown scheduling method.")
  if (method!="fifo") ids[order(-cell_counts[ids],ids)] else ids
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

fast_stream_queue <- function(cl,tasks,fun,on_result=function(...) NULL,control=NULL) {
  n <- length(tasks); workers <- length(cl)
  if (!n) return(invisible(list()))
  if (!workers) stop("Continuous queue requires at least one worker.")
  ids <- vapply(tasks,function(x) as.integer(x$tile_id),integer(1))
  if (anyNA(ids) || anyDuplicated(ids)) stop("Queue task IDs must be unique.")
  send <- getFromNamespace("sendCall","parallel")
  receive <- getFromNamespace("recvOneResult","parallel")
  active <- rep(NA_integer_,workers); pending <- seq_len(n); live <- rep(TRUE,workers)
  result <- vector("list",n); failures <- character(); failure_announced <- FALSE
  aware <- !is.null(control); poll <- if(aware) control$poll_seconds else .1
  started <- rep(NA_real_,n); started_at <- vector("list",n); sampled_peak <- rep(NA_real_,n)
  start_rss <- hwm_before <- rep(NA_real_,n); peaks <- numeric(workers)
  bypasses <- integer(n); blocked_at <- rep(NA_real_,n); reserved <- NA_integer_
  reservation_deadline <- NA_real_; reservation_retained <- numeric(workers)
  exhausted <- FALSE
  cleaned_idle <- rep(FALSE,workers); ranks <- integer(n); dispatched <- 0L
  snapshot <- budget <- NULL; retained <- numeric(workers); estimates <- vector("list",n)
  telemetry_valid <- FALSE
  now <- function() unname(proc.time()[["elapsed"]])
  fail <- function(msg) failures <<- unique(c(failures,msg))
  refresh <- function() {
    if(!aware) return(TRUE)
    sample <- tryCatch(control$sample(control$pids),error=function(e) list(ok=FALSE,reason=conditionMessage(e)))
    if(!isTRUE(sample$ok)) {
      telemetry_valid <<- FALSE
      fail(paste0("Memory monitoring unavailable; dispatch paused: ",sample$reason))
      return(FALSE)
    }
    snapshot <<- sample
    telemetry_valid <<- TRUE
    retained <<- pmax(retained,snapshot$worker_rss)
    for(node in which(!is.na(active))) {
      j <- active[node]
      sampled_peak[j] <<- max(sampled_peak[j],snapshot$worker_rss[node],na.rm=TRUE)
      updated <- fast_stream_resource_estimate(control$model,nrow(tasks[[j]]$tile_points),
        control$plan$settings$safety_factor)$peak_bytes
      peaks[node] <<- max(peaks[node],updated,snapshot$worker_rss[node]*control$plan$settings$safety_factor)
    }
    budget <<- fast_stream_memory_budget(snapshot,retained,peaks,control$plan$reserve_bytes,
      control$plan$settings$safety_factor)
    TRUE
  }
  if(aware) {
    if(length(control$pids)!=workers) stop("Memory controller PID count differs from cluster.")
    if(!refresh()) stop(paste(failures,collapse="; "))
    retained <- snapshot$worker_rss*control$plan$settings$safety_factor
    refresh()
  }
  audit <- function(event,j,node,reason,status,answer=NULL) {
    if(!aware || is.null(snapshot) || is.null(budget)) return(invisible(NULL))
    fast_stream_memory_audit(control,event,ranks[j],node,tasks[[j]],estimates[[j]],snapshot,budget,
      sum(!is.na(active)),reason,status,started=started_at[[j]],
      seconds=if(event=="dispatch") NA_real_ else now()-started[j],peak=sampled_peak[j],
      hwm_before=hwm_before[j],checkpoint_bytes=if(is.null(answer$timing$checkpoint_bytes))
        NA_real_ else answer$timing$checkpoint_bytes,telemetry_valid=telemetry_valid)
  }
  submit <- function(j,node,reason) {
    pending <<- setdiff(pending,j); active[node] <<- j
    started[j] <<- now(); started_at[[j]] <<- Sys.time()
    dispatched <<- dispatched+1L; ranks[j] <<- dispatched
    tryCatch({
      if(aware) {
        peaks[node] <<- estimates[[j]]$peak_bytes
        start_rss[j] <<- snapshot$worker_rss[node]; sampled_peak[j] <<- start_rss[j]
        hwm_before[j] <<- snapshot$worker_hwm[node]; cleaned_idle[node] <<- FALSE
        budget <<- fast_stream_memory_budget(snapshot,retained,peaks,control$plan$reserve_bytes,
          control$plan$settings$safety_factor)
        audit("dispatch",j,node,reason,"running")
      }
      send(cl[[node]],fast_stream_task_envelope,list(tasks[[j]],fun),tag=j)
    },error=function(e) {
      active[node] <<- NA_integer_; live[node] <<- FALSE; peaks[node] <<- 0
      fail(paste0("tile ",ids[j]," send/audit failed: ",conditionMessage(e)))
    })
  }
  receive_ready <- function(node) {
    j <- active[node]
    received <- tryCatch(receive(cl[node]),error=function(e)e)
    if(inherits(received,"error")) {
      active[node] <<- NA_integer_; live[node] <<- FALSE; peaks[node] <<- 0
      fail(paste0("tile ",ids[j]," worker disconnected: ",conditionMessage(received)))
      audit("complete",j,node,"worker_disconnected","failed")
      return(invisible(NULL))
    }
    if(length(received$tag)!=1L || is.na(received$tag) || received$tag!=j) {
      fail("Continuous queue received an invalid worker/tag.")
      active[node] <<- NA_integer_; peaks[node] <<- 0
      return(invisible(NULL))
    }
    if(aware) {
      refresh()
      # HWM sharpens this observation only when a NEW lifetime maximum occurred.
      if(is.finite(snapshot$worker_hwm[node]) && snapshot$worker_hwm[node]>hwm_before[j])
        sampled_peak[j] <<- max(sampled_peak[j],snapshot$worker_hwm[node])
    }
    active[node] <<- NA_integer_; peaks[node] <<- 0
    value <- received$value
    if(is.list(value) && isTRUE(value$ok)) {
      err <- tryCatch({on_result(value$value,tasks[[j]],node,j); NULL},error=function(e)conditionMessage(e))
      if(is.null(err)) result[j] <<- list(value$value) else fail(paste0("tile ",ids[j]," callback: ",err))
      if(aware) {
        fast_stream_memory_learn(control$model,nrow(tasks[[j]]$tile_points),now()-started[j],
          sampled_peak[j],start_rss[j],value$value$timing,control$plan$settings$safety_factor)
        refresh()
        audit("complete",j,node,"worker_result",if(is.null(err)) "committed" else "callback_failed",value$value)
      }
    } else {
      msg <- if(is.list(value) && is.character(value$message)) value$message else as.character(value)
      fail(paste0("tile ",ids[j],": ",paste(msg,collapse=" ")))
      audit("complete",j,node,"worker_error","failed")
    }
  }
  step <- function() {
    if(aware) {
      refresh()
      if(!is.null(control$stop_file) && file.exists(control$stop_file))
        fail("Cooperative STOP requested; running tasks drained.")
    }
    while(length(pending) && any(is.na(active) & live) && !length(failures)) {
      idle <- which(is.na(active) & live)
      if(aware) {
        if(!refresh()) break
        for(j in pending) estimates[[j]] <<- fast_stream_resource_estimate(control$model,
          nrow(tasks[[j]]$tile_points),control$plan$settings$safety_factor)
        aged <- pending[!is.na(blocked_at[pending]) &
          (now()-blocked_at[pending]>=control$starvation_seconds | bypasses[pending]>=control$max_bypasses)]
        if(is.na(reserved) && length(aged)) {
          reserved <<- aged[which.min(blocked_at[aged])]
          reservation_retained <<- retained
          reservation_deadline <<- fast_stream_reservation_deadline(reserved,active,started,
            estimates,budget,retained,now())
        }
        choice <- fast_stream_memory_choose(pending,idle,estimates,budget,retained,ids,
          control$plan$settings$backfill,reserved,list(deadline=reservation_deadline,
            now=now(),active=active,started=started,retained=reservation_retained))
        if(is.null(choice)) {
          if(!any(!is.na(active))) {
            # GC is requested only on idle workers. Never assume that it frees RSS.
            clean <- idle[!cleaned_idle[idle]]
            if(length(clean)) {
              for(node in clean) {
                tryCatch(parallel::clusterCall(cl[node],function() {gc(full=TRUE); NULL}),
                  error=function(e) fail(paste0("Idle worker GC failed: ",conditionMessage(e))))
                cleaned_idle[node] <<- TRUE
              }
              if(refresh()) retained[idle] <<- snapshot$worker_rss[idle]*control$plan$settings$safety_factor
              next
            }
            exhausted <<- TRUE
            fail(paste0("No admissible tile with no active work after idle GC (budget/retained RSS). ",
              "Increase free RAM or resume with a smaller worker cap. Checkpoints retained."))
          }
          break
        }
        if(choice$job!=choice$head) {
          j <- choice$head; bypasses[j] <<- bypasses[j]+1L
          if(is.na(blocked_at[j])) blocked_at[j] <<- now()
        }
        if(choice$reason=="bounded_backfill") reservation_retained[choice$node] <<-
          max(reservation_retained[choice$node],estimates[[choice$job]]$peak_bytes)
        if(identical(choice$job,reserved)) {
          reserved <<- NA_integer_; reservation_deadline <<- NA_real_
          reservation_retained <<- numeric(workers)
        }
        submit(choice$job,choice$node,choice$reason)
      } else submit(pending[1L],idle[1L],"legacy_queue")
    }
    if(aware && is.function(control$on_tick)) control$on_tick(sum(!is.na(active)))
    nodes <- which(!is.na(active))
    if(length(nodes)) {
      ready <- tryCatch(socketSelect(lapply(cl[nodes],function(x)x$con),timeout=poll),error=function(e)e)
      if(inherits(ready,"error")) {
        # Identify broken connections without dropping other workers' commits.
        ready <- logical(length(nodes))
        for(k in seq_along(nodes)) {
          probe <- tryCatch(socketSelect(list(cl[[nodes[k]]]$con),timeout=0),error=function(e)e)
          if(inherits(probe,"error")) {
            node <- nodes[k]; fail(paste("Worker socket failed:",conditionMessage(probe)))
            active[node] <<- NA_integer_; peaks[node] <<- 0; live[node] <<- FALSE
          } else ready[k] <- probe[1L]
        }
      }
      # Drain already-ready failures before admitting another task.
      for(node in nodes[ready]) receive_ready(node)
    }
    if(aware && !is.null(snapshot) && !is.null(budget)) {
      if(!is.null(control$monitor_file)) fast_stream_append_csv(data.frame(time=as.character(Sys.time()),
        telemetry_valid=telemetry_valid,sample_time=as.character(snapshot$measured_at),
        active_tasks=sum(!is.na(active)),pending_tasks=length(pending),
        reserved_tile_id=if(is.na(reserved)) NA_integer_ else ids[reserved],
        reservation_seconds_remaining=reservation_deadline-now(),
        available_bytes=snapshot$available_bytes,committed_bytes=budget$committed,
        measured_process_rss_sum_bytes=snapshot$master_rss+sum(snapshot$worker_rss),
        reserved_future_bytes=budget$future,headroom_bytes=budget$headroom,
        memory_scale=control$model$memory_scale,observed_support_fraction=control$model$support_fraction,
        status=if(length(failures)) "draining" else "running"),control$monitor_file)
      if(is.function(control$on_tick)) control$on_tick(sum(!is.na(active)))
    }
  }
  while(length(pending) || any(!is.na(active))) {
    if(length(failures) && !failure_announced) {
      message("Queue stopped dispatching; waiting for running tiles to commit: ",failures[1L])
      failure_announced <- TRUE
    }
    if(length(failures) && !any(!is.na(active))) break
    # First Ctrl+C requests a drain; subsequent interrupts are deferred until drained.
    tryCatch(if(length(failures)) suspendInterrupts(step()) else step(),
      interrupt=function(e) fail("Interrupted; running tasks drained."),error=function(e) fail(conditionMessage(e)))
  }
  if (length(failures)) {
    msg <- paste0("Continuous queue failed; committed tiles can be resumed. ",paste(unique(failures),collapse="; "))
    if(exhausted && length(failures)==1L) stop(structure(list(message=msg,call=NULL),
      class=c("era5_memory_exhausted","error","condition")))
    stop(msg,call.=FALSE)
  }
  invisible(result)
}

# Recover only the explicit all-idle, post-GC capacity condition. STOP, monitoring,
# QC, computation and checkpoint failures must never trigger automatic restarts.
fast_stream_queue_recover <- function(cl,tasks,fun,on_result=function(...) NULL,control=NULL,rebuild=NULL) {
  if(is.null(control) || !is.function(rebuild))
    return(fast_stream_queue(cl,tasks,fun,on_result,control))
  complete <- rep(FALSE,length(tasks)); answers <- vector("list",length(tasks)); attempt <- 0L
  repeat {
    remaining <- which(!complete)
    callback <- function(answer,task,node,j) {
      on_result(answer,task,node,remaining[j])
      complete[remaining[j]] <<- TRUE; answers[remaining[j]] <<- list(answer)
    }
    problem <- tryCatch({
      fast_stream_queue(cl,tasks[remaining],fun,callback,control); NULL
    },era5_memory_exhausted=function(e)e)
    if(is.null(problem)) return(invisible(answers))
    if(attempt>=control$plan$settings$max_rebuilds)
      stop(conditionMessage(problem)," Automatic pool recovery limit reached (",attempt,
        "). Free RAM or lower workers in config_memory.R, then run: bash run_server.sh full",call.=FALSE)
    if(!is.null(control$stop_file) && file.exists(control$stop_file))
      stop("Cooperative STOP requested before idle-pool recovery; checkpoints retained.")
    attempt <- attempt+1L
    message("Idle-pool recovery ",attempt,"/",control$plan$settings$max_rebuilds,
      ": all tasks drained; rebuilding <= ",max(1L,length(cl)%/%2L)," workers; committed tiles retained.")
    replacement <- rebuild(cl,max(1L,length(cl)%/%2L),tasks[!complete],attempt,control)
    cl <- replacement$cluster; control <- replacement$control
  }
}

run_global_era5 <- function(config, analysis = NULL) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Install ncdf4 before global analysis.")
  g <- global_defaults(config); f <- fast_defaults(config)
  scheduling <- fast_stream_scheduling(config$fast)
  memory_mode <- identical(scheduling$method,"memory_lpt")
  if(memory_mode && isTRUE(f$autotune)) stop("Use autotune=FALSE with memory_lpt; admission already adapts continuously.")
  resource_plan <- NULL; memory_control <- NULL; concurrency_sum <- 0; concurrency_samples <- 0
  concurrency_last <- NULL; concurrency_active <- 0L
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
  resume_audit <- data.frame(tile_id=which(done),n_cells=cell_counts[done],status=rep("validated_existing",sum(done)),
    timing_sidecar_present=file.exists(file.path(out_dir,"tile_timings",sprintf("tile_%05d.rds",which(done)))))
  # A validated scientific checkpoint remains authoritative even if a timing CSV/sidecar was lost.
  if(any(done)) utils::write.csv(resume_audit,file.path(out_dir,paste0("resume_",session_id,".csv")),row.names=FALSE)
  record_progress <- function(label, workers, seconds = NA_real_, new_cells = 0L) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    remaining <- sum(cell_counts[!done])
    eta <- if (remaining == 0) 0 else if (is.finite(rate) && rate > 0) remaining / rate else NA_real_
    if(memory_mode && !is.null(resource_plan) && remaining>0 && concurrency_samples>0L) {
      remaining_seconds <- sum(vapply(cell_counts[!done],function(n)
        fast_stream_resource_estimate(resource_plan$model,n,resource_plan$settings$safety_factor)$seconds,numeric(1)))
      eta <- remaining_seconds/max(1,concurrency_sum/concurrency_samples)
    }
    row <- data.frame(time = as.character(Sys.time()), phase = label, backend = backend, workers = workers,
      completed_tiles = sum(done), total_tiles = length(tasks), completed_land_points = sum(cell_counts[done]),
      total_land_points = sum(cell_counts), seconds_this_wave = seconds, land_points_this_wave = new_cells,
      elapsed_this_session_seconds = elapsed, land_points_per_second = rate,
      estimated_remaining_compute_hours = eta / 3600,
      estimated_remaining_recent_hours = if (remaining==0) 0 else if (is.finite(seconds) && seconds>0 && new_cells>0) remaining/(new_cells/seconds)/3600 else NA_real_,
      elapsed_from_launch_hours = as.numeric(difftime(Sys.time(), .warm_session_started, units="hours")),
      projected_total_session_hours_excluding_finalize = as.numeric(difftime(Sys.time(), .warm_session_started, units="hours")) + eta/3600,
      eta_basis=if(memory_mode) "remaining tile-size costs / sampled active concurrency; active tile residuals conservatively uncredited" else "session observed throughput",
      eta_uncertainty=if(sum(done)-initial_done<20L) "high: initial priors/few tiles; geography and support vary" else "approximate: excludes finalization and changing external memory load")
    progress[[length(progress) + 1L]] <<- row
    utils::write.csv(do.call(rbind, progress), file.path(out_dir, "progress.csv"), row.names = FALSE)
    last_progress <<- unname(proc.time()[["elapsed"]])
    message("Progress ", sum(done), "/", length(tasks), " tiles, ", sum(cell_counts[done]), "/", sum(cell_counts),
      " land points; ", workers, " workers; remaining compute ~",
      if (is.finite(eta)) sprintf("%.2f h", eta / 3600) else "pending benchmark",
      ". ETA is approximate, especially early; tile size/support and memory pressure vary; excludes final aggregation/plots.")
    if (is.finite(eta) && remaining>0) {
      forecast <- row$projected_total_session_hours_excluding_finalize
      message(sprintf("Approximate forecast for this session: %.2f h including elapsed setup, excluding final aggregation/plots.", forecast))
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
    if(memory_mode) {
      resource_plan <- fast_stream_memory_plan(f,g,ndraw,cell_counts[!done])
      candidates <- resource_plan$workers
      saveRDS(list(workers=candidates,settings=resource_plan$settings,initial=resource_plan$initial,
        reserve_bytes=resource_plan$reserve_bytes,largest_peak_bytes=resource_plan$largest_peak_bytes,
        prior="User-reported EPYC 9T95 measurements; estimates, NOT guaranteed memory bounds"),
        file.path(out_dir,paste0("resource_plan_",session_id,".rds")))
      message(sprintf("Memory plan: %d CPU slots, %.2f GiB effective RAM, %.2f GiB reserve; pool <= %d, full-tile predicted peak %.2f GiB. Active concurrency is admitted separately.",
        resource_plan$initial$cpu_slots,resource_plan$initial$total_bytes/1024^3,
        resource_plan$reserve_bytes/1024^3,candidates,resource_plan$largest_peak_bytes/1024^3))
    }
    chosen_workers <- candidates[1L]
    .set_analysis_stage("start_PSOCK_workers", workers = max(candidates), announce = TRUE)
    source_env <- environment(run_global_era5)
    symbols <- ls(envir = source_env, all.names = TRUE)
    fast_functions <- symbols[grepl("^(fast_|[.]fast_)", symbols) &
      vapply(symbols, function(x) is.function(get(x, envir = source_env)), logical(1))]
    exports <- unique(c(fast_functions, "FAST_CPP_CODE", "era5_required_variables", "history_columns",
      ".history_frame", ".b21_basis", ".tv_basis", "make_history_design", ".make_joint_reference",
      "history_support", "evaluate_support", "normal_axis", "read_rectangular_history", "sum_by_cell",
      "posterior_rows", ".set_analysis_stage", ".analysis_error_details"))
    qc_history <- .fast_qc_histories(analysis)$training
    qc_history <- qc_history[seq_len(min(2L, nrow(qc_history))), , drop = FALSE]
    qc_expected <- fast_predict_ab(cache, qc_history)[c("A1", "B1", "A0", "B0")]
    # Create/export/initialize one worker at a time. A failed initializer cleans up
    # every already-created connection; registration happens after each successful QC.
    initialize <- function(node,i) {
      parallel::clusterExport(node,exports,envir=source_env)
      parallel::clusterCall(node,fast_worker_initialize,.libPaths(),cache_dir,backend,bundle,
        part_dir,checkpoint_dir,fingerprint,qc_history,qc_expected)[[1L]]
    }
    can_start <- function(audit) {
      if(!memory_mode) return(TRUE)
      if(file.exists(file.path(config$output_dir,"STOP"))) stop("Cooperative STOP requested during worker initialization.")
      pids <- vapply(audit,function(x) as.integer(x$pid),integer(1))
      s <- fast_stream_resources(pids)
      if(!isTRUE(s$ok)) stop("Worker initialization memory monitoring unavailable: ",s$reason)
      b <- fast_stream_memory_budget(s,s$worker_rss*resource_plan$settings$safety_factor,
        numeric(length(pids)),resource_plan$reserve_bytes,resource_plan$settings$safety_factor)
      required <- max(1.5*1024^3*resource_plan$settings$safety_factor,
        resource_plan$idle_prior_bytes+resource_plan$largest_peak_bytes)
      b$headroom>=required
    }
    pool <- fast_stream_start_pool(max(candidates),initialize,
      function(audit) fast_register_workers(config$output_dir,list(audit)),can_start,
      outfile=file.path(out_dir,"worker_log.txt"))
    cl <- pool$cluster; worker_audit <- pool$audit
    candidates <- pmin(candidates,length(cl)); chosen_workers <- candidates[1L]
    if(memory_mode) {
      memory_control <- fast_stream_memory_control(resource_plan,
        vapply(worker_audit,function(x) as.integer(x$pid),integer(1)),out_dir,session_id,
        file.path(config$output_dir,"STOP"))
      memory_control$on_tick <- function(active) {
        tick <- unname(proc.time()[["elapsed"]])
        if(!is.null(concurrency_last)) {
          dt <- tick-concurrency_last
          concurrency_sum <<- concurrency_sum+concurrency_active*dt
          concurrency_samples <<- concurrency_samples+dt
        }
        concurrency_last <<- tick; concurrency_active <<- active
        if(unname(proc.time()[["elapsed"]])-last_progress>=progress_interval)
          record_progress("memory_queue",chosen_workers)
      }
    }
    saveRDS(worker_audit, file.path(out_dir, "worker_QC.rds"))
    message("All ", length(cl), " workers passed prediction QC; backend=", backend,
            ". Complete production tiles will be used for speed measurement.")
    rebuild <- function(old_cluster,cap,pending_tasks,attempt,control) {
      # Called only after the queue has drained every in-flight task and GC failed.
      # Retire the entire idle pool, so native/R heap retention can actually clear.
      parallel::stopCluster(old_cluster); cl <<- NULL
      deadline <- unname(proc.time()[["elapsed"]])+10
      repeat {
        states <- vapply(control$pids,fast_process_state,character(1))
        if(all(states=="dead")) break
        if(unname(proc.time()[["elapsed"]])>=deadline)
          stop("Retired idle workers have not exited within 10 s; recovery stopped. Checkpoints retained.")
        Sys.sleep(.1)
      }
      gc(full=TRUE)
      resource_plan$largest_peak_bytes <<- max(vapply(pending_tasks,function(task)
        fast_stream_resource_estimate(resource_plan$model,nrow(task$tile_points),
          resource_plan$settings$safety_factor)$peak_bytes,numeric(1)))
      replacement <- fast_stream_start_pool(cap,initialize,
        function(audit) fast_register_workers(config$output_dir,list(audit)),can_start,
        outfile=file.path(out_dir,"worker_log.txt"))
      cl <<- replacement$cluster
      chosen_workers <<- length(cl)
      saveRDS(replacement$audit,file.path(out_dir,paste0("worker_QC_",session_id,"_rebuild_",attempt,".rds")))
      control$pids <- vapply(replacement$audit,function(x) as.integer(x$pid),integer(1))
      control$plan <- resource_plan
      fast_stream_append_csv(data.frame(time=as.character(Sys.time()),attempt=attempt,
        old_workers=length(old_cluster),new_workers=length(cl),pending_tiles=length(pending_tasks),
        reason="all_idle_no_fit_after_GC",status="rebuilt_and_QC_passed"),
        file.path(out_dir,paste0("memory_recovery_",session_id,".csv")))
      memory_control <<- control
      list(cluster=cl,control=control)
    }
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
          record_progress(phase,if(memory_mode) chosen_workers else workers,seconds,completed_cells)
      }
      fast_stream_queue_recover(cl[seq_len(workers)],tasks[ids],fast_worker_task,on_result,memory_control,rebuild)
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
    scheduling = scheduling, scheduling_cost_proxy = if(memory_mode)
      "size-bucket seconds/cell updated per completion; conservative full-support matrix memory prior" else "selected cells per tile",
    dispatch_log = dispatch_file, tile_timing_log = timing_file,
    memory_dispatch_log=if(is.null(memory_control)) NULL else memory_control$audit_file,
    memory_monitor_log=if(is.null(memory_control)) NULL else memory_control$monitor_file,
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
