options(warn=1)
root <- "E:/ERA5/ERA5_workflow"
code <- file.path(root,"03_server_bundle/code")
source(file.path(code,"00_core_loader.R"))
stream_load_warm02(file.path(code,"warm02_reference.R"))
mask <- readRDS(file.path(root,"02_local/mask/mask_0p25deg_pop_gt1.rds"))
sel <- mask$grid[mask$grid$mask_code==1L,,drop=FALSE]
path <- file.path(root,"01_source_data/02_Thermal history/thermal_history_2010.nc")
nc <- ncdf4::nc_open(path)
lon <- nc_axis(nc,"lon")$values; lat <- nc_axis(nc,"lat")$values
ncdf4::nc_close(nc)
dates <- seq(as.Date("2010-01-01"),as.Date("2010-01-31"),by="day")
history <- list(files=list(list(lon=lon,lat=lat,dates=dates)))
temp <- fast_daily_inspect_temperature(list(daily_cache=list(temperature_mode="required",
  temperature_dir=file.path(root,"01_source_data/01_raw_ERA5"),temperature_checksum_inputs=FALSE)),history)
start <- proc.time()[[3L]]
whole <- fast_daily_temperature_chunk(temp,list(tile_points=sel),dates)
elapsed <- proc.time()[[3L]]-start
# Compare all 104344 cells x 31 dates using the original narrow-band calls.
start <- proc.time()[[3L]]
for (band in unique((sel$selected_lat_index-1L)%/%12L)) {
  rows <- which((sel$selected_lat_index-1L)%/%12L==band)
  old <- fast_daily_temperature_chunk(temp,list(tile_points=sel[rows,,drop=FALSE]),dates)
  stopifnot(identical(as.numeric(old),as.numeric(whole[rows,,drop=FALSE])))
}
old_elapsed <- proc.time()[[3L]]-start
cat("MONTH READER EXACT PASS:",nrow(sel)*length(dates),"values; whole-month seconds:",
    elapsed,"; narrow-band seconds:",old_elapsed,"\n")
# Expanded thermal-history bands retain exactly the same original source values.
nc <- ncdf4::nc_open(path)
ilat <- which(lat >= -60 & lat <= 80)
ii <- unique(round(seq(1,nrow(sel),length.out=40)))
for (v in era5_required_variables) {
  wide <- read_rectangular_history(nc,v,1:1440,ilat[1:90],1L,31L)
  for (first in seq(1L,85L,by=12L)) {
    b <- first:min(first+11L,90L)
    narrow <- read_rectangular_history(nc,v,1:1440,ilat[b],1L,31L)
    rows <- as.vector(outer(1:1440,(b-1L)*1440L,"+"))
    stopifnot(identical(as.numeric(narrow),as.numeric(wide[rows,,drop=FALSE])))
  }
}
ncdf4::nc_close(nc)
cat("HISTORY BAND EXACT PASS: all four variables, 1440 x 90 x 31 values each\n")
