# Execution-only resource accounting. No scientific arrays or checkpoints change.
# Linux semantics: docs.kernel.org/filesystems/proc.html and admin-guide/cgroup-v2.html.
fast_stream_read_lines <- function(path) {
  tryCatch(suppressWarnings(readLines(path,warn=FALSE)),error=function(e) character())
}

fast_stream_kb <- function(lines,key) {
  value <- grep(paste0("^",key,":\\s+[0-9]+\\s+kB\\s*$"),lines,value=TRUE)
  if (length(value)!=1L) return(NA_real_)
  as.numeric(sub(paste0("^",key,":\\s+([0-9]+).*"),"\\1",value))*1024
}

fast_stream_cpu_list <- function(x) {
  if (length(x)!=1L || is.na(x) || !grepl("^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$",x))
    return(NA_real_)
  pieces <- strsplit(x,",",fixed=TRUE)[[1L]]
  counts <- vapply(strsplit(pieces,"-",fixed=TRUE),function(p) {
    p <- as.numeric(p)
    if(length(p)==1L) 1 else if(p[2]>=p[1]) p[2]-p[1]+1 else NA_real_
  },numeric(1))
  sum(counts)
}

fast_stream_resources <- function(pids=integer(),proc_root="/proc",cgroup_root="/sys/fs/cgroup") {
  mem <- fast_stream_read_lines(file.path(proc_root,"meminfo"))
  total <- fast_stream_kb(mem,"MemTotal"); available <- fast_stream_kb(mem,"MemAvailable")
  self <- fast_stream_read_lines(file.path(proc_root,"self/status"))
  allowed <- sub("^Cpus_allowed_list:\\s*","",grep("^Cpus_allowed_list:",self,value=TRUE))
  cpus <- suppressWarnings(parallel::detectCores(logical=TRUE))
  affinity <- fast_stream_cpu_list(allowed)
  if (is.finite(affinity)) cpus <- min(cpus,affinity,na.rm=TRUE)
  reasons <- character(); limits <- character()
  if(!is.finite(affinity)) reasons <- c(reasons,"CPU affinity unavailable")
  cg <- fast_stream_read_lines(file.path(proc_root,"self/cgroup"))
  v2 <- grep("^0::",cg,value=TRUE)
  if (length(v2)!=1L) reasons <- c(reasons,"cgroup v2 identity unavailable")
  if (length(v2)==1L) {
    root <- normalizePath(cgroup_root,winslash="/",mustWork=FALSE)
    current <- normalizePath(file.path(root,sub("^0::/?","",v2)),winslash="/",mustWork=FALSE)
    if (!identical(current,root) && !startsWith(current,paste0(root,"/"))) {
      reasons <- c(reasons,"cgroup path is outside its mount")
    } else repeat {
      if(!dir.exists(current)) { reasons <- c(reasons,"cgroup directory unavailable"); break }
      for (nm in c("memory.max","memory.high")) {
        # The root cgroup has no limit files. A disabled controller can also omit
        # them at a child; continue to ancestors. Existing but unreadable files fail closed.
        if(!file.exists(file.path(current,nm))) next
        raw <- fast_stream_read_lines(file.path(current,nm))
        if (length(raw)!=1L) { reasons <- c(reasons,paste("unreadable",nm)); next }
        if (raw=="max") next
        cap <- suppressWarnings(as.numeric(raw))
        used <- suppressWarnings(as.numeric(fast_stream_read_lines(file.path(current,"memory.current"))))
        if (length(used)!=1L || !is.finite(used) || used<0 || !is.finite(cap) || cap<=0) {
          reasons <- c(reasons,paste("invalid",nm)); next
        }
        total <- min(total,cap); available <- min(available,max(0,cap-used))
        limits <- c(limits,paste0(nm,"=",format(cap,scientific=FALSE)))
      }
      quota <- fast_stream_read_lines(file.path(current,"cpu.max"))
      if (length(quota)==1L) {
        q <- strsplit(trimws(quota),"\\s+")[[1L]]
        if (length(q)!=2L) reasons <- c(reasons,"invalid cpu.max") else if(q[1L]!="max") {
          q <- suppressWarnings(as.numeric(q))
          if (any(!is.finite(q) | q<=0)) reasons <- c(reasons,"invalid CPU quota") else
            cpus <- min(cpus,max(1,floor(q[1L]/q[2L])))
        }
      } else if(file.exists(file.path(current,"cpu.max"))) reasons <- c(reasons,"unreadable cpu.max")
      if (identical(current,root)) break
      current <- dirname(current)
    }
  }
  workers <- lapply(pids,function(pid) fast_stream_read_lines(file.path(proc_root,as.character(pid),"status")))
  rss <- vapply(workers,fast_stream_kb,numeric(1),key="VmRSS")
  hwm <- vapply(workers,fast_stream_kb,numeric(1),key="VmHWM")
  master <- fast_stream_kb(self,"VmRSS")
  if (length(cpus)!=1L || !is.finite(cpus) || cpus<1) reasons <- c(reasons,"CPU capacity unavailable")
  if (any(!is.finite(c(total,available,master,rss,hwm))) || total<=0 || available<0 ||
      available>total || master<=0 || any(rss<=0)) reasons <- c(reasons,"invalid/missing Linux memory telemetry")
  list(ok=!length(reasons),reason=paste(unique(reasons),collapse="; "),
    total_bytes=total,available_bytes=available,cpu_slots=as.integer(cpus),
    master_rss=master,worker_rss=rss,worker_hwm=hwm,pids=as.integer(pids),
    limits=limits,measured_at=Sys.time())
}

fast_stream_memory_settings <- function(x=NULL) {
  defaults <- list(reserve_fraction=.15,safety_factor=1.3,backfill=TRUE)
  if (is.null(x)) return(defaults)
  if (!is.list(x) || any(!names(x) %in% names(defaults))) stop("Unknown memory scheduling configuration.")
  for(nm in names(x)) defaults[[nm]] <- x[[nm]]
  if (!is.numeric(defaults$reserve_fraction) || length(defaults$reserve_fraction)!=1L ||
      !is.finite(defaults$reserve_fraction) || defaults$reserve_fraction<0 || defaults$reserve_fraction>=.8)
    stop("memory$reserve_fraction must be in [0,0.8).")
  if (!is.numeric(defaults$safety_factor) || length(defaults$safety_factor)!=1L ||
      !is.finite(defaults$safety_factor) || defaults$safety_factor<1)
    stop("memory$safety_factor must be >=1.")
  if (!is.logical(defaults$backfill) || length(defaults$backfill)!=1L || is.na(defaults$backfill))
    stop("memory$backfill must be TRUE or FALSE.")
  defaults
}

fast_stream_memory_model <- function(g,ndraw) {
  # User-reported pilot/full-tile measurements are priors, never hard upper bounds.
  # 0.75 GiB resident + 0.25 annual arrays + 4 GiB transient at 144*31*2000.
  # Transient corresponds to ~60 concurrent double-matrix equivalents.
  # Full support is retained as the memory bound: mean support does not bound a chunk.
  e <- new.env(parent=emptyenv())
  e$day_chunk <- g$day_chunk; e$draws <- ndraw
  e$years <- max(1,length(g$years)); e$scenarios <- max(1,length(g$cdd_scenarios))
  e$memory_scale <- 1
  e$seconds_per_cell <- c(small=2.957,middle=2.957,large=2.437)
  e$observations <- 0L; e$support_fraction <- NA_real_
  e
}

fast_stream_resource_estimate <- function(model,n_cells,safety_factor) {
  size <- n_cells/144; draws <- model$draws/2000
  raw <- (.75 + .25*size*draws*(model$years/10)*model$scenarios +
    4*size*draws*(model$day_chunk/31))*1024^3
  bucket <- if(n_cells<=30) "small" else if(n_cells<=90) "middle" else "large"
  list(peak_bytes=raw*model$memory_scale*safety_factor,raw_peak_bytes=raw,
    seconds=1+n_cells*model$seconds_per_cell[[bucket]],bucket=bucket,
    supported_chunk_rows_upper=n_cells*model$day_chunk)
}

fast_stream_memory_learn <- function(model,n_cells,seconds,observed_peak,start_rss,timing,safety_factor) {
  estimate <- fast_stream_resource_estimate(model,n_cells,safety_factor)
  if(is.finite(seconds) && seconds>0) {
    b <- estimate$bucket
    model$seconds_per_cell[[b]] <- .8*model$seconds_per_cell[[b]]+.2*seconds/n_cells
    model$observations <- model$observations+1L
  }
  # Do not misattribute a large predecessor's retained RSS to a small tile.
  # Low observations never reduce the conservative memory prior.
  if(is.finite(observed_peak) && is.finite(start_rss) && start_rss<=estimate$raw_peak_bytes &&
      observed_peak>estimate$peak_bytes)
    model$memory_scale <- max(model$memory_scale,observed_peak/estimate$raw_peak_bytes)
  rows <- timing$history_row_counts
  if(!is.null(rows) && is.finite(rows[["all"]]) && rows[["all"]]>0 && is.finite(rows[["predicted"]]))
    model$support_fraction <- rows[["predicted"]]/rows[["all"]]
  invisible(NULL)
}

fast_stream_memory_plan <- function(fast,g,ndraw,cell_counts,snapshot=fast_stream_resources()) {
  settings <- fast_stream_memory_settings(fast$memory)
  if(!isTRUE(snapshot$ok)) stop("Memory-aware scheduling requires valid Linux telemetry: ",snapshot$reason)
  model <- fast_stream_memory_model(g,ndraw)
  maximum <- fast_stream_resource_estimate(model,max(cell_counts),settings$safety_factor)$peak_bytes
  reserve <- max(4*1024^3,snapshot$total_bytes*settings$reserve_fraction)
  master <- max(2*1024^3,snapshot$master_rss*settings$safety_factor)
  idle <- .75*1024^3*settings$safety_factor
  # Pool reserves resident memory AND room for at least one largest pending tile.
  capacity <- min(snapshot$total_bytes-reserve-master,
    snapshot$available_bytes-reserve-max(0,master-snapshot$master_rss))
  pool <- min(fast$workers,snapshot$cpu_slots,length(cell_counts),
    floor((capacity-max(0,maximum-idle))/idle))
  if(!is.finite(pool) || pool<1) stop("Memory budget cannot fit one largest tile plus its worker; increase RAM/reserve headroom or reduce the worker cap. No task started.")
  list(workers=as.integer(pool),settings=settings,model=model,reserve_bytes=reserve,
    initial=snapshot,largest_peak_bytes=maximum,idle_prior_bytes=idle)
}

fast_stream_memory_budget <- function(snapshot,retained,active_peaks,reserve_bytes,safety_factor) {
  if(!isTRUE(snapshot$ok) || length(snapshot$worker_rss)!=length(retained) ||
      length(active_peaks)!=length(retained)) stop("Invalid resource snapshot/budget dimensions.")
  current <- snapshot$worker_rss
  hold <- pmax(current,retained,active_peaks)
  master <- max(2*1024^3,snapshot$master_rss*safety_factor)
  future <- sum(pmax(0,hold-current))+max(0,master-snapshot$master_rss)
  committed <- master+sum(hold)
  list(hold=hold,committed=committed,future=future,
    limit=snapshot$total_bytes-reserve_bytes,
    headroom=min(snapshot$total_bytes-reserve_bytes-committed,
      snapshot$available_bytes-reserve_bytes-future))
}

fast_stream_memory_choose <- function(pending,idle,estimates,budget,retained,ids,
                                      backfill=TRUE,reserved=NA_integer_) {
  if(!length(pending) || !length(idle)) return(NULL)
  order <- pending[order(-vapply(estimates[pending],`[[`,numeric(1),"seconds"),ids[pending])]
  if(!is.na(reserved) && reserved %in% pending) order <- c(reserved,setdiff(order,reserved))
  # Prefer reusing a resident worker, avoiding needless expansion of idle heaps.
  idle <- idle[order(-retained[idle],idle)]
  candidates <- if(!is.na(reserved) || !backfill) head(order,1L) else order
  for(j in candidates) for(node in idle) {
    extra <- max(0,estimates[[j]]$peak_bytes-budget$hold[node])
    if(is.finite(budget$headroom) && extra<=budget$headroom) return(list(job=j,node=node,
      extra=extra,head=order[1L],reason=if(!is.na(reserved)) "reserved_large" else
        if(j==order[1L]) "longest_fits" else "backfill"))
  }
  NULL
}

fast_stream_memory_control <- function(plan,pids,out_dir,session_id,stop_file) {
  list(plan=plan,pids=as.integer(pids),sample=fast_stream_resources,
    audit_file=file.path(out_dir,paste0("memory_dispatch_",session_id,".csv")),
    monitor_file=file.path(out_dir,paste0("memory_monitor_",session_id,".csv")),
    session_id=session_id,stop_file=stop_file,poll_seconds=1,
    starvation_seconds=60,max_bypasses=8L,model=plan$model)
}

fast_stream_start_pool <- function(n,initialize,on_ready=function(...) NULL,
                                  can_start=function(...) TRUE,outfile="") {
  cl <- NULL; audit <- list(); success <- FALSE
  on.exit(if(!success && !is.null(cl)) try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  for(i in seq_len(n)) {
    if(!can_start(audit)) break
    node <- parallel::makePSOCKcluster(1L,rscript_args="--vanilla",outfile=outfile)
    cl <- structure(c(unclass(cl),unclass(node)),class=class(node))
    audit[[i]] <- initialize(cl[i],i)
    on_ready(audit[[i]])
  }
  if(!length(cl)) stop("No worker can be initialized within the current memory budget.")
  success <- TRUE
  list(cluster=cl,audit=audit)
}

fast_stream_memory_audit <- function(control,event,job,node,task,estimate,snapshot,budget,
                                    active_count,reason,status,started=NA,seconds=NA_real_,
                                    peak=NA_real_,hwm_before=NA_real_,checkpoint_bytes=NA_real_,telemetry_valid=TRUE) {
  if(is.null(control$audit_file)) return(invisible(NULL))
  fast_stream_append_csv(data.frame(session_id=control$session_id,event=event,time=as.character(Sys.time()),
    telemetry_valid=telemetry_valid,sample_time=as.character(snapshot$measured_at),
    tile_id=task$tile_id,n_cells=nrow(task$tile_points),worker=node,pid=control$pids[node],
    dispatch_rank=job,started_at=as.character(started),actual_seconds=seconds,
    predicted_seconds=estimate$seconds,predicted_peak_bytes=estimate$peak_bytes,
    sampled_task_peak_rss_bytes=peak,worker_rss_bytes=snapshot$worker_rss[node],
    process_lifetime_hwm_before_bytes=hwm_before,process_lifetime_hwm_after_bytes=snapshot$worker_hwm[node],
    available_bytes=snapshot$available_bytes,budget_limit_bytes=budget$limit,
    committed_bytes=budget$committed,reserved_future_bytes=budget$future,
    dispatch_headroom_bytes=budget$headroom,active_tasks=active_count,reason=reason,
    checkpoint_bytes=checkpoint_bytes,status=status),control$audit_file)
}
