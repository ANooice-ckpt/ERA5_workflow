#!/usr/bin/env Rscript
# Tiny fully synthetic streamed tile fixture: NO original weather or fitted models.
# Checks yearly 3-draw CDD arrays, support and missing-temperature counts,
# draw alignment, and warm02's original warm-boundary checkpoint validator.
cmd <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",grep("^--file=",cmd,value=TRUE)[1L])
base <- dirname(normalizePath(script,mustWork=TRUE))
source(file.path(base,"00_core_loader.R"),local=FALSE)
stream_load_warm02(file.path(base,"warm02_reference.R"))
if (!requireNamespace("ncdf4",quietly=TRUE)) stop("Install ncdf4 for this synthetic test.")
fast_daily_cdd_selftest()
source(file.path(base,"02_stream_overrides.R"),local=.GlobalEnv)

tmp <- tempfile("era5_stream_test_")
dir.create(tmp)

ix <- list(schema="ERA5-indexed-weather-v1",n_cells=2L,
  lon_axis=seq(0,by=.25,length.out=1440L),lat_axis=c(10,9.75),
  lat_indices=1L,
  tiles=data.frame(n_cells=2L),
  selected=data.frame(pack_row=1:2,cell_id=1:2,selected_lon_index=1:2,
    selected_lat_index=c(1L,1L)))
fast_stream_set_context(ix)

source_nc <- function(year,first=TRUE) {
  path<-file.path(tmp,sprintf("packed_%d.nc",year))
  dcell<-ncdf4::ncdim_def("cell","index",1:2)
  dtime<-ncdf4::ncdim_def("time","days since 2010-01-01",1:2)
  vars<-lapply(c(era5_required_variables,"Ta_C"),function(v)
    ncdf4::ncvar_def(v,"degree_Celsius",list(dcell,dtime),missval=-9999,prec="double"))
  nc<-ncdf4::nc_create(path,vars,force_v4=TRUE)
  tryCatch({
    # cell 1 always primary-supported; cell 2 finite but outside mock support.
    for(v in era5_required_variables) {
      mat <- if (v=="B21_raw") matrix(c(1,2,1,2),2,2) else matrix(0,2,2)
      ncdf4::ncvar_put(nc,v,mat)
    }
    t <- if(first) matrix(c(28,30,30,NA),2,2) else matrix(c(27,29,24,NA),2,2)
    ncdf4::ncvar_put(nc,"Ta_C",t)
  },finally=ncdf4::nc_close(nc))
  path
}
files <- list()
for (i in 1:2) {
  yy<-2009L+i
  path<-source_nc(yy,first=i==1L)
  files[[i]] <- list(path=path,year=yy,
    dates=as.Date(sprintf("%d-01-01",yy))+0:1,
    offsets=stats::setNames(rep(0,4),era5_required_variables))
}

task <- list(tile_id=1L,lon_indices=1:2,lat_indices=1L,keep=1:2,
  tile_points=data.frame(cell_id=1:2,lon_plot=c(0,.25),lat=c(10,10),
    land_fraction=c(1,1),area_weight=c(.9,.9),cell_area_km2=c(100,100)))
sc <- list(list(name="synthetic",theta0=25,slope=1,kappa=1))
g <- list(day_chunk=1L,joint_support=FALSE,
  min_date_coverage=.5,min_supported_days=1L,min_supported_fraction=0,
  support_fraction_sensitivity=c(0,.5,.8),save_annual_grid=TRUE,
  max_cells=NULL,cdd_scenarios=sc)
bundle <- list(g=g,files=files,fixed_set=c(P50=24),set_bounds=c(-100,100),
  cache=list(n_draws=3L),support_analysis=list(),ndraw=3L,
  periods=c("2010","2011","2010-2011"),planned=c(2L,2L,4L),
  daily_cache=list(enabled=FALSE))

# Inject ONLY synthetic posterior prediction and support functions.
env <- new.env(parent=environment(fast_compute_tile))
env$evaluate_support <- function(analysis,history,joint) {
  finite <- rowSums(!is.finite(as.matrix(history)))==0L
  supported <- finite & history$B21_raw<1.5
  data.frame(finite=finite,marginal_inside=supported,joint_inside=supported,
    B21_outside=!supported,TV_lag1_outside=FALSE,TV_lag2_outside=FALSE,TV_lag3_outside=FALSE)
}
env$fast_predict_ab <- function(cache,history) {
  d<-1:cache$n_draws
  list(A1=outer(history$B21_raw,d,function(h,z).12+.01*h*z),
    B1=outer(history$B21_raw,d,function(h,z).4+.03*h*z),
    A0=.1+.002*d,B0=.3+.004*d)
}
compute<-fast_compute_tile
environment(compute)<-env
out<-compute(task,bundle,file.path(tmp,"parts"),file.path(tmp,"checkpoint"),"synthetic-stream")
cp<-readRDS(out$checkpoint)
stopifnot(identical(dim(cp$cdd_year),c(2L,3L,2L,2L,1L)),
  all(cp$cdd_primary_days[1L,]==2L),all(cp$cdd_primary_days[2L,]==0L),
  cp$cdd_valid_weather_days[1L,1L]==2L,
  cp$cdd_valid_weather_days[2L,1L]==1L,
  all(cp$cdd_valid_supported_days[1L,]==2L),
  all(cp$cdd_valid_supported_days[2L,]==0L))
# Primary supported daily t2m for cell 1: 28,30 in 2010 => 3+5 CDD0.
stopifnot(max(abs(cp$cdd_year[1L,,1L,1L,1L]-8)) < 1e-12)
# Non-supported cell 2 still contributes to all-weather M0 baseline, where Ta is finite.
stopifnot(all(cp$cdd_year[2L,,,1L,1L]==0),
  all(cp$cdd_all_weather[2L,,1L,1L]==5))
# Draws remain separate: M1 CDD differs across the three draw positions.
if (length(unique(round(cp$cdd_year[1L,,1L,2L,1L],12)))<2L)
  stop("Synthetic model posterior draws were collapsed.")
# Stream wrapper's checkpoint verifier must accept the freshly committed tile.
validated<-fast_read_checkpoint(out$checkpoint,"synthetic-stream",task,bundle)
stopifnot(identical(validated$joint_draw_cell_ids,1:2))
message("SYNTHETIC STREAM TEST PASS: exact supported CDD, all-weather baseline, paired draws, immutable checkpoint")
if (!isTRUE(getOption("era5.stream.keep_fixture",FALSE))) unlink(tmp,recursive=TRUE)
