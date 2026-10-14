#!/usr/bin/env Rscript
# Real two-worker PSOCK tests with synthetic weather/posteriors only.
# No fitted model packages, original climate input, or ECS instance required.
cmd <- commandArgs(FALSE)
script <- sub("^--file=","",grep("^--file=",cmd,value=TRUE)[1L])
base <- dirname(normalizePath(script,winslash="/",mustWork=TRUE))

run_parallel_tests <- function() {
  old <- options(era5.stream.keep_fixture=TRUE)
  on.exit(options(old),add=TRUE)
  source(file.path(base,"05_test_stream_tile.R"),local=.GlobalEnv)
  source(file.path(base,"09_parallel_runtime.R"),local=.GlobalEnv)
  fixture <- normalizePath(tmp,winslash="/",mustWork=TRUE)
  stopifnot(startsWith(fixture,paste0(normalizePath(tempdir(),winslash="/"),"/")))
  on.exit(unlink(fixture,recursive=TRUE),add=TRUE)
  checks <- 0L
  check <- function(ok,label) {
    if (!isTRUE(ok)) stop("FAIL: ",label)
    checks <<- checks+1L; message("PASS: ",label)
  }
  rejects <- function(expr,pattern=NULL) {
    e <- tryCatch({force(expr); NULL},error=function(e)e)
    inherits(e,"error") && (is.null(pattern) || grepl(pattern,conditionMessage(e)))
  }
  check(identical(fast_stream_order_tiles(c(5L,2L,4L),c(1,8,2,8,3)),c(2L,4L,5L)),
        "LPT retains canonical IDs and breaks equal costs by tile ID")
  check(identical(fast_stream_order_tiles(c(5L,2L,4L),rep(1,5),"fifo"),c(5L,2L,4L)),
        "FIFO retains input order")
  check(rejects(fast_stream_order_tiles(c(1L,1L),c(2,3))) &&
        rejects(fast_stream_order_tiles(2L,c(2,NA))),"invalid costs and duplicate IDs rejected")
  check(rejects(fast_stream_scheduling(list(scheduling="bad"))) &&
        rejects(fast_stream_scheduling(list(progress_seconds=0))),"invalid scheduling settings rejected")
  check(length(fast_stream_queue(list(),list(),identity))==0L,"empty queue does not need workers")
  check(rejects(fast_stream_queue(list(),list(list(tile_id=1L)),identity)),"nonempty queue requires workers")

  cl <- parallel::makePSOCKcluster(2L,rscript_args="--vanilla",outfile=file.path(fixture,"queue_workers.log"))
  on.exit(if(!is.null(cl)) parallel::stopCluster(cl),add=TRUE)
  # Worker 1 can only finish after job 5: a 2W=4 batch barrier cannot pass this.
  gate <- file.path(fixture,"fifth_job_started")
  jobs <- lapply(1:6,function(i) list(tile_id=as.integer(i),gate=gate))
  gate_fun <- function(task) {
    if (task$tile_id==1L) {
      deadline <- Sys.time()+10
      while (!file.exists(task$gate) && Sys.time()<deadline) Sys.sleep(.01)
      if (!file.exists(task$gate)) stop("2W barrier detected")
    }
    if (task$tile_id==5L) writeLines("started",task$gate)
    list(tile_id=task$tile_id,pid=Sys.getpid())
  }
  order_seen <- integer()
  answers <- fast_stream_queue(cl,jobs,gate_fun,function(answer,task,worker,rank) {
    order_seen <<- c(order_seen,answer$tile_id)
  })
  check(identical(vapply(answers,`[[`,integer(1),"tile_id"),1:6),
        "out-of-order completions return in dispatch order")
  check(match(5L,order_seen)<match(1L,order_seen),"job 5 starts before slow first job finishes: no 2W barrier")
  check(length(unique(vapply(answers,`[[`,integer(1),"pid")))==2L,"tasks execute on two distinct PSOCK workers")

  fail_dir <- file.path(fixture,"queue_failure"); dir.create(fail_dir)
  fail_jobs <- lapply(1:6,function(i) list(tile_id=as.integer(i),root=fail_dir,fail=TRUE))
  fail_fun <- function(task) {
    writeLines("started",file.path(task$root,paste0(task$tile_id,".started")))
    if (task$tile_id==1L && task$fail) stop("deliberate worker failure")
    Sys.sleep(.15)
    writeLines("committed",file.path(task$root,paste0(task$tile_id,".done")))
    list(tile_id=task$tile_id)
  }
  committed <- integer()
  callback <- function(answer,...) committed <<- c(committed,answer$tile_id)
  check(rejects(fast_stream_queue(cl,fail_jobs,fail_fun,callback),"deliberate worker failure"),
        "worker error is propagated after draining active work")
  check(identical(committed,2L) && !file.exists(file.path(fail_dir,"3.started")),
        "active successful tile commits; undispatched jobs remain unstarted")
  old_hash <- tools::md5sum(file.path(fail_dir,"2.done"))
  remaining <- setdiff(1:6,committed)
  retry <- lapply(fail_jobs[remaining],function(x) {x$fail<-FALSE;x})
  fast_stream_queue(cl,retry,fail_fun,callback)
  check(identical(sort(committed),1:6) && identical(old_hash,tools::md5sum(file.path(fail_dir,"2.done"))),
        "resume schedules only missing tiles and retains the previous commit")
  check(rejects(fast_stream_queue(cl,jobs,identity,function(...) stop("callback failure")),"callback failure"),
        "master callback errors also stop dispatch and drain workers")
  check(length(fast_stream_queue(cl,jobs[1],identity))==1L,"cluster remains usable after a drained error")
  parallel::stopCluster(cl); cl <- NULL

  # Exercise the actual global runner, worker exports, callbacks, progress CSV,
  # tile checkpoints and both finalizers. Only model/backend/input inspection
  # are replaced by tiny synthetic equivalents; no scheduler is mocked.
  e <- new.env(parent=.GlobalEnv)
  stream_load_warm02(file.path(base,"warm02_reference.R"),target=e)
  sys.source(file.path(base,"02_stream_overrides.R"),envir=e)
  sys.source(file.path(base,"09_parallel_runtime.R"),envir=e)
  index <- ix; index$n_tiles <- 2L; index$n_cells <- 3L
  index$tiles <- data.frame(n_cells=c(1L,2L))
  index$selected <- data.frame(pack_row=1:3,cell_id=c(1L,3L,4L),
    selected_lon_index=c(1L,3L,4L),selected_lat_index=1L,cell_area_km2=100)
  integration_files <- files
  for (i in 1:2) {
    p <- file.path(fixture,paste0("integration_",i,".nc"))
    dims <- list(ncdf4::ncdim_def("cell","index",1:3),
      ncdf4::ncdim_def("time","days since 2010-01-01",1:2))
    vars <- lapply(c(era5_required_variables,"Ta_C"),function(v)
      ncdf4::ncvar_def(v,"degree_Celsius",dims,missval=-9999,prec="double"))
    nc <- ncdf4::nc_create(p,vars)
    for (v in era5_required_variables) ncdf4::ncvar_put(nc,v,
      if(v=="B21_raw") matrix(c(1,2,1,1,2,1),3,2) else matrix(0,3,2))
    ncdf4::ncvar_put(nc,"Ta_C",matrix(c(28,NA,30,26,29,27),3,2))
    ncdf4::nc_close(nc)
    integration_files[[i]]$path <- p
  }
  saveRDS(index,file.path(fixture,"pack_index.rds"))
  e$fast_stream_set_context(index)
  e$evaluate_support <- env$evaluate_support
  e$fast_predict_ab <- env$fast_predict_ab
  e$fast_make_cache <- function(analysis,...) list(n_draws=nrow(analysis$draws_m1),
    draw_ids_m1=analysis$draw_ids_m1,draw_ids_m0=analysis$draw_ids_m0)
  e$fast_init_backend <- function(...) "r"
  e$fast_kernel_QC <- function(...) invisible(NULL)
  e$fast_boundary_preflight <- function(...) invisible(NULL)
  e$.fast_qc_histories <- function(...) list(training=data.frame(B21_raw=c(1,2),
    TV_lag1_raw=0,TV_lag2_raw=0,TV_lag3_raw=0))
  e$global_defaults <- function(config) config$global
  e$inspect_history_files <- local({
    ff <- integration_files
    for (i in seq_along(ff)) {ff[[i]]$lon<-index$lon_axis;ff[[i]]$lat<-index$lat_axis}
    function(...) list(files=ff,audit=data.frame(year=c(2010L,2011L),source="synthetic"))
  })
  e$fast_worker_initialize <- function(paths,cache_dir,backend,bundle,part_dir,
                                      checkpoint_dir,fingerprint,qc_history,qc_expected) {
    .libPaths(paths)
    fast_stream_set_context(readRDS(file.path(bundle$g$input_dir,"pack_index.rds")))
    assign(".fast_worker_bundle",bundle,envir=.GlobalEnv)
    assign(".fast_worker_part_dir",part_dir,envir=.GlobalEnv)
    assign(".fast_worker_checkpoint_dir",checkpoint_dir,envir=.GlobalEnv)
    assign(".fast_worker_fingerprint",fingerprint,envir=.GlobalEnv)
    list(pid=Sys.getpid(),backend=backend,max_absolute_error=0)
  }
  # Exported functions resolve in the source environment on the master and the
  # worker global environment after parallel::clusterExport.
  for (nm in c("evaluate_support","fast_predict_ab","fast_make_cache","fast_init_backend",
               "fast_kernel_QC","fast_boundary_preflight",".fast_qc_histories",
               "global_defaults","fast_worker_initialize")) environment(e[[nm]]) <- e
  analysis <- list(fixed_set=c(P50=24),data=data.frame(SET_raw=1:50),
    draws_m1=matrix(0,50,1),draws_m0=matrix(0,50,1),draw_ids_m1=1:50,
    draw_ids_m0=101:150,draw_ids=1:50,draw_metadata=list(),input_md5="synthetic",
    pred=list(),training_history=data.frame(),joint_reference=NULL)
  gg <- bundle$g
  gg$input_dir <- fixture; gg$years <- 2010:2011; gg$n_draws <- 50L
  gg$spatial_stride <- 1L; gg$tile_side <- 2L; gg$lat_range <- c(10,10)
  gg$weight_scheme <- "observation"; gg$max_cells <- NULL
  cfg <- list(global=gg,daily_cache=list(enabled=FALSE),fast=list(backend="r",workers=2L,
    autotune=FALSE,scheduling="lpt",progress_seconds=.001))
  cfg$output_dir <- file.path(fixture,"global_lpt")
  dir.create(cfg$output_dir)
  lpt_token <- e$fast_acquire_run_lock(cfg$output_dir)
  lpt_dir <- cfg$output_dir
  on.exit(e$fast_release_run_lock(lpt_dir,lpt_token),add=TRUE)
  lpt <- e$run_global_era5(cfg,analysis)
  cfg$output_dir <- file.path(fixture,"global_fifo"); cfg$fast$scheduling <- "fifo"
  dir.create(cfg$output_dir)
  fifo_token <- e$fast_acquire_run_lock(cfg$output_dir)
  fifo_dir <- cfg$output_dir
  on.exit(e$fast_release_run_lock(fifo_dir,fifo_token),add=TRUE)
  fifo <- e$run_global_era5(cfg,analysis)
  check(identical(lpt$aggregate_draws,fifo$aggregate_draws),"actual runner: FIFO/LPT warm draws and pairing are identical")
  # Run fingerprints also contain run-local cache paths; compare every CDD
  # field except that provenance identity across separate output directories.
  cdd_fields <- setdiff(names(lpt$cdd),"fingerprint")
  check(identical(lpt$cdd[cdd_fields],fifo$cdd[cdd_fields]),
        "actual runner: FIFO/LPT CDD draws and support outputs are identical")
  ex <- lpt$metadata$fast_execution
  timings <- read.csv(ex$tile_timing_log)
  stages <- grep("_seconds$",names(timings),value=TRUE)
  stages <- setdiff(stages,"total_seconds")
  check(nrow(timings)==2L && all(timings$total_seconds>=0) &&
        all(abs(rowSums(timings[,stages])-timings$total_seconds)<1e-8),
        "per-stage elapsed time partitions full tile time, including commit/readback")
  check(all(c("read_history_chunk_seconds","evaluate_support_seconds","predict_supported_ab_seconds",
    "accumulate_supported_cdd_seconds","checkpoint_write_validate_seconds") %in% names(timings)),
    "timing CSV contains I/O, support, prediction, accumulation and checkpoint stages")
  dispatch <- read.csv(ex$dispatch_log)
  check(identical(dispatch$tile_id,2:1) && all(dispatch$method=="lpt") &&
        identical(read.csv(fifo$metadata$fast_execution$dispatch_log)$tile_id,1:2),
        "dispatch audit records LPT reorder and canonical IDs")
  progress <- read.csv(file.path(lpt$output_dir,"progress.csv"))
  check(tail(progress$completed_tiles,1)==2L && tail(progress$completed_land_points,1)==3L,
        "global progress reaches all committed tiles/cells")
  timing_rds <- list.files(file.path(lpt$output_dir,"tile_timings"),full.names=TRUE)
  check(length(timing_rds)==2L && all(vapply(timing_rds,function(p)
    identical(readRDS(p)$fingerprint,ex$fingerprint),logical(1))),"durable per-tile timing sidecars match the run")

  # Actual PSOCK/science pipeline; only host telemetry is deterministic on Windows.
  e$fast_stream_resources <- function(pids=integer()) list(ok=TRUE,reason="synthetic host",
    total_bytes=64*1024^3,available_bytes=60*1024^3,cpu_slots=2L,master_rss=1024^3,
    worker_rss=rep(.5*1024^3,length(pids)),worker_hwm=rep(.6*1024^3,length(pids)),
    pids=pids,limits=character(),measured_at=Sys.time())
  environment(e$fast_stream_resources) <- e
  cfg$output_dir <- file.path(fixture,"global_memory"); cfg$fast$scheduling <- "memory_lpt"
  dir.create(cfg$output_dir)
  memory_dir <- cfg$output_dir
  memory_token <- e$fast_acquire_run_lock(memory_dir)
  on.exit(e$fast_release_run_lock(memory_dir,memory_token),add=TRUE)
  memory_result <- e$run_global_era5(cfg,analysis)
  check(identical(memory_result$aggregate_draws,lpt$aggregate_draws) &&
        identical(memory_result$cdd[cdd_fields],lpt$cdd[cdd_fields]),
        "actual runner: memory-aware/LPT/FIFO scientific draws and pairing are identical")
  memory_ex <- memory_result$metadata$fast_execution
  audit <- read.csv(memory_ex$memory_dispatch_log)
  admitted <- audit[audit$event=="dispatch",]
  check(nrow(admitted)==2L && all(admitted$committed_bytes<=admitted$budget_limit_bytes) &&
        all(admitted$dispatch_headroom_bytes>=0) && sum(audit$status=="committed")==2L,
        "actual runner: memory audit records safe admissions and both commits")
  memory_cp <- file.path(memory_result$output_dir,"checkpoints","tile_00001.rds")
  memory_hash <- tools::md5sum(memory_cp)
  unlink(file.path(memory_result$output_dir,"checkpoints","tile_00002.rds"))
  unlink(file.path(memory_result$output_dir,"tile_timings","tile_00001.rds"))
  unlink(memory_ex$tile_timing_log)
  cfg$fast$workers <- 1L; cfg$fast$memory <- list(reserve_fraction=.2)
  memory_resume <- e$run_global_era5(cfg,analysis)
  check(memory_resume$metadata$fast_execution$resumed_tiles==1L &&
        identical(memory_hash,tools::md5sum(memory_cp)) &&
        identical(memory_resume$cdd,memory_result$cdd) &&
        identical(memory_resume$aggregate_draws,memory_result$aggregate_draws),
        "memory resume: missing timing/CSV and changed execution limits preserve checkpoint/science")
  cfg$fast$workers <- 2L; cfg$fast$memory <- NULL

  # Force all-idle memory exhaustion after the first real scientific checkpoint.
  # A newly initialized one-worker pool must continue the identical computation.
  normal_resources <- e$fast_stream_resources
  e$.normal_estimate <- e$fast_stream_resource_estimate
  e$fast_stream_resource_estimate <- function(model,n_cells,safety_factor) {
    estimate <- .normal_estimate(model,n_cells,safety_factor)
    estimate$peak_bytes <- 32*1024^3
    estimate
  }
  environment(e$fast_stream_resource_estimate) <- e
  cfg$output_dir <- file.path(fixture,"global_recovery")
  dir.create(cfg$output_dir)
  recovery_dir <- cfg$output_dir
  recovery_token <- e$fast_acquire_run_lock(recovery_dir)
  on.exit(e$fast_release_run_lock(recovery_dir,recovery_token),add=TRUE)
  e$.recovery_checkpoint_dir <- file.path(recovery_dir,"global","checkpoints")
  e$.recovery_committed_hash <- NULL
  e$fast_stream_resources <- function(pids=integer()) {
    cp <- list.files(.recovery_checkpoint_dir,pattern="^tile_[0-9]+[.]rds$",full.names=TRUE)
    if(length(cp) && is.null(.recovery_committed_hash)) .recovery_committed_hash <<- tools::md5sum(cp)
    list(ok=TRUE,reason="synthetic retained heaps after first commit",total_bytes=64*1024^3,
      available_bytes=if(length(cp) && length(pids)>1L) 4*1024^3 else 60*1024^3,
      cpu_slots=2L,master_rss=1024^3,worker_rss=rep(.5*1024^3,length(pids)),
      worker_hwm=rep(.6*1024^3,length(pids)),pids=pids,limits=character(),measured_at=Sys.time())
  }
  environment(e$fast_stream_resources) <- e
  recovered <- e$run_global_era5(cfg,analysis)
  check(identical(recovered$aggregate_draws,lpt$aggregate_draws) &&
    identical(recovered$cdd[cdd_fields],lpt$cdd[cdd_fields]),
    "actual runner: idle-pool rebuild preserves every paired scientific result")
  check(length(e$.recovery_committed_hash)==1L && identical(e$.recovery_committed_hash,
    tools::md5sum(names(e$.recovery_committed_hash))),
    "actual runner: checkpoint committed before rebuild is unchanged and never recomputed")
  recovery_log <- read.csv(list.files(file.path(recovery_dir,"global"),pattern="^memory_recovery_",full.names=TRUE))
  recovery_qc <- readRDS(list.files(file.path(recovery_dir,"global"),pattern="^worker_QC_.*_rebuild_1[.]rds$",full.names=TRUE))
  owner <- e$fast_lock_owner(file.path(recovery_dir,".run_lock"))
  check(nrow(recovery_log)==1L && recovery_log$old_workers==2L && recovery_log$new_workers==1L &&
    length(recovery_qc)==1L && recovery_qc[[1]]$pid %in% owner$worker_pids &&
    recovered$metadata$fast_execution$selected_workers==1L,
    "actual runner: rebuilt pool passes initialization QC, registers its PID and logs the reduced size")
  e$fast_stream_resources <- normal_resources
  e$fast_stream_resource_estimate <- e$.normal_estimate

  # Simulate an interrupted run with tile 1 already committed. A resumed runner
  # validates it and sends only tile 2, retaining exactly the same final draws.
  cp1 <- file.path(lpt$output_dir,"checkpoints","tile_00001.rds")
  cp2 <- file.path(lpt$output_dir,"checkpoints","tile_00002.rds")
  cp1_hash <- tools::md5sum(cp1)
  unlink(cp2)
  cfg$output_dir <- dirname(lpt$output_dir); cfg$fast$scheduling <- "lpt"
  resumed <- e$run_global_era5(cfg,analysis)
  rex <- resumed$metadata$fast_execution
  check(rex$resumed_tiles==1L && identical(cp1_hash,tools::md5sum(cp1)),
        "actual runner resume retains validated checkpoint without prediction")
  check(identical(read.csv(rex$dispatch_log)$tile_id,2L),"actual runner resume dispatches only missing tile")
  check(identical(lpt$aggregate_draws,resumed$aggregate_draws) && identical(lpt$cdd,resumed$cdd),
        "resumed warm/CDD final outputs match uninterrupted run exactly")
  message("PARALLEL RUNTIME TESTS PASS: ",checks," checks; real PSOCK, synthetic science only.")
}
run_parallel_tests()
