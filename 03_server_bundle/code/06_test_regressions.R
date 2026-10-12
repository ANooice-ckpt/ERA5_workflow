#!/usr/bin/env Rscript
# Model-free regression tests. Synthetic support/prediction only; no RStan or ECS.
cmd <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",grep("^--file=",cmd,value=TRUE)[1L])
base <- dirname(normalizePath(script,mustWork=TRUE))

run_regressions <- function() {
  old <- options(era5.stream.keep_fixture=TRUE)
  on.exit(options(old),add=TRUE)
  source(file.path(base,"05_test_stream_tile.R"),local=.GlobalEnv)
  fixture_dir <- normalizePath(tmp,winslash="/",mustWork=TRUE)
  temp_root <- paste0(normalizePath(tempdir(),winslash="/",mustWork=TRUE),"/")
  stopifnot(startsWith(fixture_dir,temp_root))
  on.exit(unlink(fixture_dir,recursive=TRUE),add=TRUE)
  nchecks <- 0L
  check <- function(ok,label) {
    if (!isTRUE(ok)) stop("FAIL: ",label)
    nchecks <<- nchecks+1L
    message("PASS: ",label)
  }
  rejects <- function(expr,pattern=NULL) {
    err <- tryCatch({force(expr); NULL},error=function(e)e)
    inherits(err,"error") && (is.null(pattern) || grepl(pattern,conditionMessage(err)))
  }

  # Golden comparison with the unmodified warm02 implementation on raw-shaped NC.
  core <- new.env(parent=.GlobalEnv)
  stream_load_warm02(file.path(base,"warm02_reference.R"),target=core)
  original_env <- new.env(parent=core)
  original_env$evaluate_support <- env$evaluate_support
  original_env$fast_predict_ab <- env$fast_predict_ab
  original_compute <- core$fast_compute_tile
  environment(original_compute) <- original_env
  raw_bundle <- bundle
  for (i in seq_along(files)) {
    pc <- ncdf4::nc_open(files[[i]]$path)
    weather <- lapply(era5_required_variables,function(v)ncdf4::ncvar_get(pc,v))
    ncdf4::nc_close(pc)
    raw_path <- file.path(tmp,paste0("raw_",i,".nc"))
    dims <- list(ncdf4::ncdim_def("longitude","degrees_east",c(0,.25)),
      ncdf4::ncdim_def("latitude","degrees_north",10),
      ncdf4::ncdim_def("time","days since 2010-01-01",0:1))
    vars <- lapply(era5_required_variables,function(v)
      ncdf4::ncvar_def(v,"degree_Celsius",dims,missval=-9999,prec="double"))
    nc <- ncdf4::nc_create(raw_path,vars)
    for (j in seq_along(vars)) ncdf4::ncvar_put(nc,era5_required_variables[j],array(weather[[j]],c(2,1,2)))
    ncdf4::nc_close(nc)
    raw_bundle$files[[i]]$path <- raw_path
  }
  original_result <- original_compute(task,raw_bundle,file.path(tmp,"original_parts"),
    file.path(tmp,"original_checkpoints"),"synthetic-original")
  original_cp <- readRDS(original_result$checkpoint)
  fields <- c("aggregate_sums","aggregate_denominator","aggregate_cell_counts",
    "support_sums","support_denominator","support_cell_counts","summary_table",
    "joint_draw_means","joint_draw_denominator","joint_draw_cell_ids")
  for (nm in fields) check(isTRUE(all.equal(cp[[nm]],original_cp[[nm]],tolerance=1e-12)),
                          paste("original warm field unchanged:",nm))
  scalar_cdd <- array(0,dim(cp$cdd_year))
  temps <- list(c(28,30),c(27,24))
  for (yy in 1:2) for (d in 1:3) for (t in temps[[yy]]) {
    u1 <- (.4+.03*d+1)/(.12+.01*d)
    u0 <- (.3+.004*d+1)/(.1+.002*d)
    scalar_cdd[1,d,yy,1,1] <- scalar_cdd[1,d,yy,1,1]+max(t-25,0)
    scalar_cdd[1,d,yy,2,1] <- scalar_cdd[1,d,yy,2,1]+max(t-25-(u1-u0),0)
  }
  check(max(abs(cp$cdd_year-scalar_cdd))<1e-12,"CDD agrees with independent daily scalar formula")

  corruptions <- list(
    negative_cdd=function(z) {z$cdd_year[1,1,1,1,1]<- -1; z},
    negative_baseline=function(z) {z$cdd_all_weather[1,1,1,1]<- -1; z},
    infinite_cdd=function(z) {z$cdd_year[1,1,1,1,1]<-Inf; z},
    missing_cdd=function(z) {z$cdd_year[1,1,1,1,1]<-NA_real_; z},
    excess_weather_days=function(z) {z$cdd_valid_weather_days[1,1]<-999; z},
    fractional_days=function(z) {z$cdd_valid_weather_days[2,1]<-1.5; z},
    supported_exceeds_weather=function(z) {z$cdd_valid_weather_days[1,1]<-1; z},
    impossible_intersection=function(z) {z$cdd_valid_supported_days[1,1]<-1; z},
    nonzero_unsupported_cdd=function(z) {z$cdd_year[2,1,1,2,1]<-1; z},
    nonzero_missing_weather_baseline=function(z) {z$cdd_valid_weather_days[2,1]<-0; z},
    baseline_subset=function(z) {z$cdd_year[1,1,1,1,1]<-100; z},
    equal_support_baseline=function(z) {z$cdd_all_weather[1,1,1,1]<-100; z},
    mismatched_warm_support=function(z) {
      z$cdd_history_finite_days[2,1]<-1; z
    },
    changed_cell_area=function(z) {z$cdd_cell_area_km2[1]<-101; z})
  for (nm in names(corruptions)) {
    bad_path <- file.path(tmp,paste0(nm,".rds"))
    saveRDS(corruptions[[nm]](cp),bad_path)
    check(rejects(fast_read_checkpoint(bad_path,"synthetic-stream",task,bundle)),
          paste("reject corrupted checkpoint:",nm))
  }

  # Resume a committed tile without invoking a predictor; restore missing export.
  resume_env <- new.env(parent=env)
  resume_env$fast_predict_ab <- function(...) stop("Unexpected repeated prediction")
  resume_compute <- fast_compute_tile
  environment(resume_compute) <- resume_env
  unlink(out$partition)
  resumed <- resume_compute(task,bundle,file.path(tmp,"parts"),file.path(tmp,"checkpoint"),"synthetic-stream")
  check(isTRUE(resumed$resumed) && identical(readRDS(resumed$partition),cp$summary_table),
        "checkpoint resume restores partition without prediction")

  # Sparse edge tiles can contain one cell; production omits annual warm tables.
  single_task <- task
  single_task$lon_indices <- 1L; single_task$keep <- 1L
  single_task$tile_points <- task$tile_points[1L,,drop=FALSE]
  multi_bundle <- bundle
  multi_bundle$g$save_annual_grid <- FALSE
  multi_bundle$g$cdd_scenarios[[2L]] <- list(name="draw_thresholds",theta0=c(24,25,26),slope=2,kappa=0)
  multi_out <- compute(single_task,multi_bundle,file.path(tmp,"single_parts"),
    file.path(tmp,"single_cp"),"synthetic-stream")
  multi_cp <- fast_read_checkpoint(multi_out$checkpoint,"synthetic-stream",single_task,multi_bundle)
  check(identical(dim(multi_cp$cdd_year),c(1L,3L,2L,2L,2L)) && nrow(multi_cp$summary_table)==1L,
        "single-cell multi-scenario checkpoint with full-period-only warm table")
  check(identical(as.numeric(multi_cp$cdd_year[1,,1,1,2]),c(10,8,6)) &&
        identical(multi_cp$cdd_year[,,,1,2],multi_cp$cdd_year[,,,2,2]),
        "draw-specific thresholds and zero-transfer scenario stay aligned")

  analysis <- list(draw_ids_m1=1:3,draw_ids_m0=11:13)
  final_env <- new.env(parent=environment(fast_finalize_global))
  final_env$fast_finalize_global_orig <- function(...,out_dir) list(
    aggregate_draws=list(paired_draw_ids=data.frame(position=1:3,m1=1:3,m0=11:13)),
    metadata=list(),output_dir=out_dir)
  # Match the original argument names because the wrapper calls it positionally.
  formals(final_env$fast_finalize_global_orig) <- formals(fast_finalize_global)
  final <- fast_finalize_global
  environment(final) <- final_env
  finalize_case <- function(name,b=bundle,z=cp) {
    out_dir <- file.path(tmp,name,"global")
    dir.create(file.path(out_dir,"grid_partitions"),recursive=TRUE)
    dir.create(file.path(out_dir,"checkpoints"))
    part <- file.path(out_dir,"grid_partitions","tile_00001.rds")
    saveRDS(z$summary_table,part)
    saveRDS(z,file.path(out_dir,"checkpoints","tile_00001.rds"))
    final(list(),analysis,b,list(),list(),1:2,1L,list(),part,out_dir,
          list(fingerprint="synthetic-stream",total_land_points=2L))
  }
  normal <- finalize_case("normal")
  check(all(normal$cdd$status=="ok") &&
        abs(normal$cdd$area_km2_degree_days[1,1,1,1]-800)<1e-12,
        "eligible result retains exact area-weighted totals")
  empty_bundle <- bundle
  empty_bundle$g$min_supported_days <- 99L
  empty <- finalize_case("empty",empty_bundle)
  check(all(empty$cdd$status=="no_eligible_support") &&
        all(empty$cdd$eligible_grid_cells==0) &&
        all(is.na(empty$cdd$area_km2_degree_days)) &&
        all(is.na(empty$cdd$baseline_all_weather_area_km2_degree_days)),
        "empty totals and baseline are NA with explicit status")
  for (name in c("cdd_global_posterior_summary.csv","cdd_global_total_posterior_summary.csv")) {
    tab <- read.csv(file.path(empty$output_dir,name))
    check(all(tab$status=="no_eligible_support") &&
          all(is.na(as.matrix(tab[,c("lower","median","upper","prob_positive","prob_negative")]))),
          paste("empty CSV intervals are NA:",name))
  }
  mixed_bundle <- bundle
  mixed_bundle$g$min_supported_days <- 3L
  mixed <- finalize_case("mixed",mixed_bundle)
  check(identical(as.vector(mixed$cdd$status),c("no_eligible_support","no_eligible_support","ok")),
        "annual empty periods coexist with eligible full period")
  cold_bundle <- bundle
  cold_bundle$g$cdd_scenarios[[1]]$theta0 <- 100
  cold_out <- compute(task,cold_bundle,file.path(tmp,"cold_parts"),file.path(tmp,"cold_cp"),"synthetic-stream")
  cold <- finalize_case("cold",cold_bundle,readRDS(cold_out$checkpoint))
  check(all(cold$cdd$status=="ok") && all(cold$cdd$area_km2_degree_days==0),
        "genuine zero CDD remains zero on eligible cold days")
  false_zero <- empty$cdd
  false_zero$area_km2_degree_days[] <- 0
  check(rejects(fast_stream_validate_global_cdd(false_zero)),"verifier rejects false zero for empty support")
  false_na <- normal$cdd
  false_na$area_km2_degree_days[1,1,1,1] <- NA_real_
  check(rejects(fast_stream_validate_global_cdd(false_na)),"verifier rejects missing draws on eligible support")

  # Exercise the actual verifier entrypoint with 2000 synthetic paired positions,
  # including the all-empty case that used to fail max(...,na.rm=TRUE).
  verify_case <- function(result,label) {
    cdd <- result$cdd
    cdd$paired_draw_ids <- data.frame(position=1:2000,m1=1:2000,m0=2001:4000)
    for (nm in c("area_km2_degree_days","mean_CDD_per_area_weighted_valid_supported_day",
                 "baseline_all_weather_area_km2_degree_days",
                 "baseline_all_weather_mean_per_area_weighted_valid_weather_day")) {
      a <- cdd[[nm]]; dims <- dim(a); dims[length(dims)] <- 2000L
      cdd[[nm]] <- array(rep(as.vector(a),length.out=prod(dims)),dims)
    }
    # Only empty and cold cases are expanded; their NA/zero draws are constant.
    saveRDS(cdd,file.path(result$output_dir,"cdd_global_posterior_draws.rds"))
    saveRDS(list(fingerprint="synthetic-stream",total_tiles=1L,total_land_points=2L),
      file.path(result$output_dir,"fast_manifest.rds"))
    saveRDS(list(status="COMPLETED",fingerprint="synthetic-stream"),
      file.path(result$output_dir,"numerical_completion.rds"))
    saveRDS(list(paired_draw_ids=cdd$paired_draw_ids,estimates=array(0,c(3L,8L,2000L))),
      file.path(result$output_dir,"warm_boundary_area_day_posterior_draws.rds"))
    root <- file.path(tmp,paste0("verify_",label)); dir.create(root)
    writeLines(dirname(result$output_dir),file.path(root,"latest_run.txt"))
    pack <- file.path(root,"pack"); dir.create(pack)
    index <- ix; index$n_tiles <- 1L
    saveRDS(index,file.path(pack,"pack_index.rds"))
    config <- file.path(root,"config.R")
    writeLines(c("ERA5_CONFIG <-",capture.output(dput(list(
      remote=list(output_dir=root,pack_dir=pack),local=list(pack_dir=pack),
      pack=list(expected_selected_cells=2L))))),config)
    log <- system2(file.path(R.home("bin"),"Rscript"),
      c("--vanilla",shQuote(file.path(base,"04_verify_outputs.R")),shQuote(config),"run"),
      stdout=TRUE,stderr=TRUE)
    check(is.null(attr(log,"status")) && any(grepl("RUN VALIDATION PASS",log)),
          paste("04 verifier entrypoint:",label))
  }
  verify_case(empty,"all_empty")
  verify_case(cold,"valid_zero")

  # New, incomplete, failed, completed, and publication-crash-window states.
  root <- file.path(tmp,"run_states")
  dest <- fast_prepare_run(root,list(resume="new"))
  fast_release_run_lock(dest$path,dest$locktoken)
  check(isTRUE(fast_select_run_dir(root,list(resume="auto"))$resumed),"auto resumes incomplete run")
  fast_stream_write_run_status(dest$path,"RUNNING")
  check(fast_existing_run(dest$path,dest$root)$status=="RUNNING","running status published")
  fast_stream_write_run_status(dest$path,"FAILED")
  check(fast_select_run_dir(root,list(resume="auto"))$previous_status=="FAILED","auto resumes failed run")
  global <- file.path(dest$path,"global"); dir.create(global)
  for (nm in c("warm_boundary_area_day_posterior_draws.rds","cdd_global_posterior_draws.rds",
               "global_run_metadata.rds")) saveRDS(list(),file.path(global,nm))
  saveRDS(list(fingerprint="test-complete"),file.path(global,"fast_manifest.rds"))
  receipt <- file.path(global,"numerical_completion.rds")
  saveRDS(list(status="COMPLETED",fingerprint="test-complete"),receipt)
  check(fast_existing_run(dest$path,dest$root)$status=="COMPLETED",
        "numeric receipt overrides stale failed/running status")
  unlink(file.path(dest$path,"standalone_run_status.rds"))
  check(rejects(fast_prepare_run(root,list(resume="auto")),"already complete"),
        "legacy v1 completion receipt prevents automatic rerun")
  check(isTRUE(fast_select_run_dir(root,list(resume=dest$path))$resumed),
        "explicit resume remains available for integrity review")
  saveRDS(list(status="COMPLETED",fingerprint="wrong"),receipt)
  check(rejects(fast_existing_run(dest$path,dest$root),"fingerprint"),"reject mismatched completion receipt")
  unlink(receipt)
  race_env <- new.env(parent=environment(fast_prepare_run))
  race_env$fast_select_run_dir <- function(...) {
    chosen <- fast_select_run_dir(...)
    saveRDS(list(status="COMPLETED",fingerprint="test-complete"),receipt)
    chosen
  }
  prepare_race <- fast_prepare_run; environment(prepare_race) <- race_env
  check(rejects(prepare_race(root,list(resume="auto")),"finished before lock acquisition"),
        "completion is rechecked after lock acquisition")
  token <- fast_acquire_run_lock(dest$path)
  fast_release_run_lock(dest$path,token)
  check(TRUE,"completion race releases the run lock")
  fast_stream_write_run_status(dest$path,"COMPLETED","test-complete")
  check(readRDS(file.path(dest$path,"standalone_run_status.rds"))$status=="COMPLETED",
        "completed status published in inherited selector format")
  message("REGRESSION TESTS PASS: ",nchecks," checks. Synthetic fixtures only; fitted-model/ECS validation remains separate.")
}
run_regressions()
