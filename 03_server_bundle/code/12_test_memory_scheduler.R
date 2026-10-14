#!/usr/bin/env Rscript
# Linux parsers via fixtures + real PSOCK queues using synthetic memory telemetry.
cmd <- commandArgs(FALSE)
script <- sub("^--file=","",grep("^--file=",cmd,value=TRUE)[1L])
base <- dirname(normalizePath(script,winslash="/",mustWork=TRUE))
source(file.path(base,"00_core_loader.R"))
stream_load_warm02(file.path(base,"warm02_reference.R"))
source(file.path(base,"02_stream_overrides.R"))
source(file.path(base,"09_parallel_runtime.R"))

run_memory_tests <- function() {
  fixture <- tempfile("era5_memory_tests_"); dir.create(fixture)
  stopifnot(startsWith(normalizePath(fixture,winslash="/"),paste0(normalizePath(tempdir(),winslash="/"),"/")))
  on.exit(unlink(fixture,recursive=TRUE),add=TRUE)
  checks <- 0L
  check <- function(ok,label) {
    if(!isTRUE(ok)) stop("FAIL: ",label)
    checks <<- checks+1L; message("PASS: ",label)
  }
  rejects <- function(expr,pattern="") {
    err <- tryCatch({force(expr);NULL},error=function(e)e)
    inherits(err,"error") && grepl(pattern,conditionMessage(err))
  }
  gib <- 1024^3
  write <- function(path,lines) {
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE); writeLines(lines,path)
  }
  proc <- file.path(fixture,"proc"); cg <- file.path(fixture,"cgroup")
  write(file.path(proc,"meminfo"),c("MemTotal: 67108864 kB","MemAvailable: 62914560 kB"))
  write(file.path(proc,"self/status"),c("VmRSS: 1048576 kB","Cpus_allowed_list: 0-7"))
  write(file.path(proc,"self/cgroup"),"0::/job")
  write(file.path(proc,"123/status"),c("VmRSS: 524288 kB","VmHWM: 4194304 kB"))
  for(path in c(cg,file.path(cg,"job"))) {
    write(file.path(path,"memory.max"),"max"); write(file.path(path,"memory.high"),"max")
    write(file.path(path,"cpu.max"),"max 100000"); write(file.path(path,"memory.current"),"1073741824")
  }
  write(file.path(cg,"memory.max"),"17179869184")
  write(file.path(cg,"memory.current"),"4294967296")
  write(file.path(cg,"cpu.max"),"400000 100000")
  s <- fast_stream_resources(123L,proc,cg)
  check(s$ok && s$total_bytes==16*gib && s$available_bytes==12*gib && s$cpu_slots<=4L,
    "Linux parser respects ancestor cgroup RAM/current/quota and CPU affinity")
  check(s$worker_rss==.5*gib && s$worker_hwm==4*gib,
    "RSS and lifetime HWM remain distinct measurements")
  write(file.path(cg,"job/memory.high"),"12884901888")
  check(fast_stream_resources(123L,proc,cg)$total_bytes==12*gib,"memory.high is also a conservative ceiling")
  check(!fast_stream_resources(c(123L,999L),proc,cg)$ok,"missing process telemetry is rejected")
  write(file.path(cg,"cpu.max"),"invalid")
  check(!fast_stream_resources(123L,proc,cg)$ok,"malformed cgroup quota is rejected")
  write(file.path(proc,"self/cgroup"),"3:memory:/job")
  check(!fast_stream_resources(123L,proc,cg)$ok,"unsupported cgroup v1 fails closed")
  write(file.path(proc,"self/cgroup"),"0::/")
  unlink(file.path(cg,c("memory.max","memory.high","cpu.max","memory.current")))
  root_snapshot <- fast_stream_resources(123L,proc,cg)
  check(root_snapshot$ok && root_snapshot$total_bytes==64*gib,
    "normal cgroup root without resource-limit files uses host memory and affinity")
  check(fast_stream_cpu_list("0-3,8,10-11")==7 && is.na(fast_stream_cpu_list("4-1")),
    "CPU list ranges and invalid ranges are handled")
  check(rejects(fast_stream_memory_settings(list(safety_factor=.5))) &&
    rejects(fast_stream_memory_settings(list(backfill=NA))) &&
    rejects(fast_stream_memory_settings(list(unknown=1))),"invalid memory settings fail before workers")

  snapshot <- function(pids=integer(),ram=16,available=14,cpu=2) list(ok=TRUE,reason="synthetic",
    total_bytes=ram*gib,available_bytes=available*gib,cpu_slots=as.integer(cpu),master_rss=gib,
    worker_rss=rep(.5*gib,length(pids)),worker_hwm=rep(.6*gib,length(pids)),
    pids=pids,limits=character(),measured_at=Sys.time())
  g <- list(day_chunk=31L,years=2010:2019,cdd_scenarios=list(list()))
  f <- list(workers=64L,memory=list())
  for(cpu in c(32L,64L)) {
    plan <- fast_stream_memory_plan(f,g,2000L,rep(144L,100),snapshot(ram=cpu*4,available=cpu*4-2,cpu=cpu))
    check(plan$workers==cpu && abs(plan$largest_peak_bytes/gib-6.5)<1e-12,
      paste("resource plan for",cpu,"CPU slots and",cpu*4,"GiB: bounded pool, 6.5 GiB full-tile prior"))
  }
  check(rejects(fast_stream_memory_plan(f,g,2000L,144L,snapshot(ram=8,available=3)),"cannot fit"),
    "host too small for one tile fails without starting workers")
  model <- fast_stream_memory_model(g,2000L)
  peak <- fast_stream_resource_estimate(model,144L,1.3)$peak_bytes
  g2 <- g; g2$day_chunk <- 62L
  check(fast_stream_resource_estimate(fast_stream_memory_model(g2,2000L),144L,1.3)$peak_bytes>peak &&
    fast_stream_resource_estimate(fast_stream_memory_model(g,4000L),144L,1.3)$peak_bytes>peak,
    "memory estimate accounts for chunk size and posterior matrix dimensions")
  fast_stream_memory_learn(model,144L,400,8*gib,.5*gib,list(history_row_counts=c(all=100,predicted=40)),1.3)
  raised <- model$memory_scale
  fast_stream_memory_learn(model,144L,300,2*gib,.5*gib,list(),1.3)
  check(raised>1 && model$memory_scale==raised && model$support_fraction==.4,
    "observations increase uncertain memory estimates; low observations cannot lower the floor")
  ss <- snapshot(1:2)
  budget <- fast_stream_memory_budget(ss,c(.5,4)*gib,c(6.5,0)*gib,4*gib,1.3)
  check(budget$committed==12.5*gib && budget$future==10.5*gib && budget$headroom<0,
    "budget counts retained idle RSS and all active future growth, not just MemAvailable")
  est <- list(list(peak_bytes=6.5*gib,seconds=350),list(peak_bytes=gib,seconds=10))
  small_budget <- list(hold=c(.5,.5)*gib,headroom=2*gib)
  pick <- fast_stream_memory_choose(1:2,1:2,est,small_budget,c(.5,.5)*gib,1:2)
  check(pick$job==2L && pick$reason=="backfill","smaller task fills a gap when the longest task cannot fit")
  check(is.null(fast_stream_memory_choose(1:2,1:2,est,small_budget,c(.5,.5)*gib,1:2,reserved=1L)) &&
    is.null(fast_stream_memory_choose(1:2,1:2,est,small_budget,c(.5,.5)*gib,1:2,backfill=FALSE)),
    "reserved oldest large task and disabled backfill both prevent overtaking")

  # Partial initializer failure must close every process connection acquired so far.
  opened <- list()
  init <- function(node,i) {
    opened[[i]] <<- node[[1L]]$con
    if(i==2L) stop("deliberate initialization failure")
    list(pid=parallel::clusterCall(node,Sys.getpid)[[1L]])
  }
  check(rejects(fast_stream_start_pool(2L,init,outfile=file.path(fixture,"init.log")),"initialization failure") &&
    all(vapply(opened,function(con) !isTRUE(tryCatch(isOpen(con),error=function(e)FALSE)),logical(1))),
    "partial pool initialization failure closes all acquired connections")

  cl <- parallel::makePSOCKcluster(2L,rscript_args="--vanilla",outfile=file.path(fixture,"workers.log"))
  on.exit(if(!is.null(cl)) try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  pids <- unlist(parallel::clusterCall(cl,Sys.getpid))
  if(Sys.info()[["sysname"]]=="Linux") {
    live <- fast_stream_resources(pids)
    check(live$ok && all(live$worker_rss>0) && all(live$worker_hwm>0) &&
      live$available_bytes<=live$total_bytes,"live Linux /proc and cgroup telemetry on actual PSOCK worker PIDs")
  }
  make_control <- function(name,counts) {
    plan <- fast_stream_memory_plan(list(workers=2L),g,2000L,counts,snapshot())
    control <- fast_stream_memory_control(plan,pids,fixture,name,file.path(fixture,paste0(name,".STOP")))
    control$sample <- snapshot; control$poll_seconds <- .01
    control$max_bypasses <- 2L
    control
  }
  jobs <- lapply(seq_len(8L),function(i) list(tile_id=as.integer(i),
    tile_points=data.frame(cell_id=seq_len(if(i<=2) 144L else 1L)),root=fixture))
  fun <- function(task) {
    Sys.sleep(if(task$tile_id<=2L) .4 else .01)
    list(tile_id=task$tile_id,pid=Sys.getpid(),timing=list(checkpoint_bytes=100))
  }
  control <- make_control("mixed",c(144,144,rep(1,6)))
  answers <- fast_stream_queue(cl,jobs,fun,control=control)
  audit <- read.csv(control$audit_file); dispatched <- audit[audit$event=="dispatch",]
  check(identical(vapply(answers,function(x)x$tile_id,integer(1)),seq_len(8L)),
    "mixed real PSOCK queue returns every result in original task order")
  check(all(dispatched$committed_bytes<=dispatched$budget_limit_bytes) &&
    all(dispatched$dispatch_headroom_bytes>=0) &&
    all(dispatched$reserved_future_bytes<=dispatched$available_bytes-control$plan$reserve_bytes),
    "every mixed-size dispatch satisfies static budget and live future-growth reservation")
  check(any(dispatched$reason=="backfill") && any(dispatched$reason=="reserved_large") &&
    match(2L,dispatched$tile_id)<=4L,"large tile is reserved after two backfills and does not starve")
  control <- make_control("impossible",144)
  control$sample <- function(pids) snapshot(pids,available=4)
  check(rejects(fast_stream_queue(cl,jobs[1L],fun,control=control),"No admissible tile"),
    "no active tasks and no fitting tile terminates after one idle GC, never deadlocks")

  control <- make_control("stop",rep(1,8))
  completed <- integer()
  small_jobs <- lapply(jobs,function(x) {x$tile_points<-x$tile_points[1,,drop=FALSE];x})
  checkpoint_fun <- function(task) {
    Sys.sleep(.08); writeLines("committed",file.path(task$root,paste0(task$tile_id,".done")))
    list(tile_id=task$tile_id)
  }
  callback <- function(answer,...) {
    completed <<- c(completed,answer$tile_id)
    writeLines("stop",control$stop_file)
  }
  check(rejects(fast_stream_queue(cl,small_jobs,checkpoint_fun,callback,control),"STOP requested") &&
    length(completed)==2L && !file.exists(file.path(fixture,"3.done")),
    "cooperative STOP drains active commits and leaves all pending tasks untouched")
  old_hash <- tools::md5sum(file.path(fixture,"1.done"))
  control <- make_control("resume",rep(1,6))
  rest <- small_jobs[setdiff(1:8,completed)]
  fast_stream_queue(cl,rest,checkpoint_fun,control=control)
  check(identical(old_hash,tools::md5sum(file.path(fixture,"1.done"))) &&
    all(file.exists(file.path(fixture,paste0(1:8,".done")))),
    "restart schedules only uncommitted tasks without modifying retained commits")
  control <- make_control("interrupt",rep(1,8)); interrupted <- FALSE; completed <- integer()
  control$on_tick <- function(active) if(!interrupted && active>0) {
    interrupted <<- TRUE
    signalCondition(structure(list(message="synthetic Ctrl+C"),class=c("interrupt","condition")))
  }
  check(rejects(fast_stream_queue(cl,small_jobs,checkpoint_fun,function(answer,...) {
    completed <<- c(completed,answer$tile_id)
  },control),"Interrupted") && length(completed)==2L,
    "interrupt condition drains active tasks before propagating failure")

  control <- make_control("lost_monitor",rep(1,8)); samples <- 0L; completed <- integer()
  control$sample <- function(pids) {
    samples <<- samples+1L
    if(samples>8L) list(ok=FALSE,reason="synthetic read failure") else snapshot(pids)
  }
  check(rejects(fast_stream_queue(cl,small_jobs,checkpoint_fun,function(answer,...) {
    completed <<- c(completed,answer$tile_id)
  },control),"monitoring unavailable") && length(completed)<=2L,
    "telemetry loss stops new admissions and drains already active work")
  control <- make_control("worker_error",rep(1,8)); completed <- integer()
  error_fun <- function(task) {
    if(task$tile_id==1L) stop("deliberate memory worker error")
    Sys.sleep(.1); list(tile_id=task$tile_id)
  }
  check(rejects(fast_stream_queue(cl,small_jobs,error_fun,function(answer,...) {
    completed <<- c(completed,answer$tile_id)
  },control),"deliberate memory worker error") && identical(completed,2L),
    "memory queue worker error drains the successful sibling")
  control <- make_control("disconnect",rep(1,8)); completed <- integer()
  crash_fun <- function(task) {
    if(task$tile_id==1L) quit(save="no",status=17L,runLast=FALSE)
    Sys.sleep(.15); list(tile_id=task$tile_id)
  }
  check(rejects(fast_stream_queue(cl,small_jobs,crash_fun,function(answer,...) {
    completed <<- c(completed,answer$tile_id)
  },control),"disconnected") && identical(completed,2L),
    "abnormal worker process exit is detected and surviving commits are drained")
  try(parallel::stopCluster(cl),silent=TRUE); cl <- NULL
  message("MEMORY SCHEDULER TESTS PASS: ",checks," checks; real PSOCK, synthetic resource/science inputs.")
  if(Sys.info()[["sysname"]]!="Linux") message("LIVE LINUX /proc/cgroup integration: NOT RUN on this host; parser fixtures passed.")
}
run_memory_tests()
