#!/usr/bin/env Rscript
# One-time Windows/local preparation: mask + source ERA5 -> sparse indexed cell x day NetCDF.
# NO model loading and NO posterior draws; no interpolation/resampling of weather.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript --vanilla 01_pack_local.R config.R")
cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE)[1L])
base <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
source(file.path(base, "00_core_loader.R"), local = FALSE)
cfg <- stream_read_config(normalizePath(args[[1L]], mustWork = TRUE))
stream_load_warm02(file.path(base, "warm02_reference.R"))

for (p in c("ncdf4", "digest")) if (!requireNamespace(p, quietly = TRUE)) stop("Install R package: ", p)
p <- cfg$pack; l <- cfg$local
if (!identical(as.integer(p$years), 2010:2019)) stop("This package contract is fixed to 2010-2019.")
if (!identical(as.numeric(p$latitude_range), c(-60,80))) stop("Mismatch: mask latitude range.")
if (p$tile_side < 1 || p$tile_side != as.integer(p$tile_side)) stop("Invalid tile_side.")
if (p$latitude_band < 1 || p$latitude_band != as.integer(p$latitude_band)) stop("Invalid latitude_band.")
if (!file.exists(l$mask_rds)) stop("Missing completed population mask RDS: ", l$mask_rds)
if (!dir.exists(paste0(l$history_dir))) stop("Missing history directory.")
dir.create(l$pack_dir, recursive = TRUE, showWarnings = FALSE)
mask_md5 <- unname(tools::md5sum(l$mask_rds))
mask <- readRDS(l$mask_rds)
if (!all(c("grid", "land_population", "fingerprint") %in% names(mask)))
  stop("Unexpected mask RDS schema: must be output from test_population_mask.R.")
mg <- mask$grid
required <- c("cell_id", "lon", "lat", "lon_plot", "mask_code", "selected_lon_index", "selected_lat_index")
if (!all(required %in% names(mg))) stop("Mask grid lacks: ", paste(setdiff(required, names(mg)), collapse=", "))
if (anyNA(mg$mask_code) || anyDuplicated(mg$cell_id)) stop("Mask cell IDs/status are not unique or finite.")
sel <- mg[mg$mask_code == 1L, , drop = FALSE]
if (!nrow(sel) || nrow(sel) != p$expected_selected_cells)
  stop("Selected mask cell count differs from configuration: ", nrow(sel))
if (anyNA(sel$cell_id) || anyNA(sel$lon) || anyNA(sel$lat)) stop("Bad selected coordinates.")
if (!all(sel$mask_code == 1L)) stop("Unexpected selected mask status.")
if (!all(sel$selected_lon_index >= 1 & sel$selected_lat_index >= 1)) stop("Bad mask indices.")

# Read original source's CF coordinates/dates/units with its validated reader.
g <- global_defaults(list(global=list(input_dir=l$history_dir, years=p$years, spatial_stride=1L,
  lat_range=p$latitude_range, checksum_inputs=isTRUE(p$checksum_original_inputs),
  preprocessing_acknowledged=isTRUE(cfg$run$preprocessing_acknowledged))))
message("Auditing original thermal history NetCDF inventory...")
inspected <- inspect_history_files(g)
first <- inspected$files[[1L]]
lon <- first$lon; lat <- first$lat
ilat <- which(lat >= p$latitude_range[1L] & lat <= p$latitude_range[2L])
stopifnot(length(lon) == 1440L, length(ilat) == 561L)
if (nrow(mg) != length(lon)*length(ilat)) stop("Mask grid dimension differs from source ERA5 axes.")
if (!same_coordinate(unique(mg$lon)[seq_along(lon)], lon)) {
  # unique preserves first-seen order; fail loudly, never rely on nearest neighbour.
  stop("Mask longitude values differ from source ERA5 coordinate order.")
}
if (!same_coordinate(mg$lat[seq(1L, nrow(mg), by=length(lon))], lat[ilat]))
  stop("Mask latitude values differ from source ERA5 coordinate order.")
expected_cell <- sel$selected_lon_index + (sel$selected_lat_index-1L)*length(lon)
if (!identical(as.integer(sel$cell_id), as.integer(expected_cell)))
  stop("Mask cell_id does not agree with source longitude/latitude indices.")

# Preserve the mask's target population/area metadata; area is whole WGS84 cell.
pop <- mask$land_population
if (!all(c("cell_id", "cell_area_km2", "population_lower_bound", "population_density_lower_bound") %in% names(pop)))
  stop("Population/whole-cell area metadata missing from mask.")
lookup <- match(sel$cell_id, pop$cell_id)
if (anyNA(lookup)) stop("Selected mask cells absent from population table.")
sel$cell_area_km2 <- pop$cell_area_km2[lookup]
sel$population_lower_bound <- pop$population_lower_bound[lookup]
sel$population_density_lower_bound <- pop$population_density_lower_bound[lookup]
if (any(!is.finite(sel$cell_area_km2) | sel$cell_area_km2 <= 0)) stop("Bad cell areas.")
if (any(!is.finite(sel$population_density_lower_bound) | sel$population_density_lower_bound <= 1))
  stop("Selected cell no longer satisfies density >1 lower-bound rule.")

# Stable ordering: tile row-major (y before x), then original geographic cell_id.
sel$tile_x <- (as.integer(sel$selected_lon_index) - 1L) %/% as.integer(p$tile_side)
sel$tile_y <- (as.integer(sel$selected_lat_index) - 1L) %/% as.integer(p$tile_side)
sel <- sel[order(sel$tile_y, sel$tile_x, sel$cell_id), , drop=FALSE]
rownames(sel) <- NULL
sel$tile_key <- paste(sel$tile_y, sel$tile_x, sep=":")
key_run <- rle(sel$tile_key)
sel$tile_id <- rep(seq_along(key_run$lengths), key_run$lengths)
sel$pack_row <- seq_len(nrow(sel))
tile_end <- cumsum(key_run$lengths)
tiles <- data.frame(tile_id=seq_along(key_run$lengths), tile_key=key_run$values,
  start_row=as.integer(tile_end-key_run$lengths+1L), n_cells=as.integer(key_run$lengths),
  stringsAsFactors=FALSE)
if (any(tiles$n_cells > p$tile_side^2)) stop("Packed tile exceeds capacity.")
message("Selected ", nrow(sel), " cells in ", nrow(tiles), " nonempty tiles (", p$tile_side, "x", p$tile_side, ").")

message("Auditing original t2m inventory (daily UTC, all 120 months)...")
temp <- fast_daily_inspect_temperature(list(daily_cache=list(
  temperature_mode="required", temperature_dir=l$temperature_dir,
  temperature_checksum_inputs=isTRUE(p$checksum_original_inputs))), inspected)
if (!identical(temp$status, "embedded")) stop("The 120 monthly t2m files are required for cooling scenarios.")

scientific_identity <- list(schema="ERA5-indexed-weather-v1", mask_md5=mask_md5,
  input_history_audit=inspected$audit, input_temperature_audit=temp$audit,
  source_warm02_md5=unname(tools::md5sum(file.path(base, "warm02_reference.R"))),
  selected_cell_ids=sel$cell_id, selected_lons=sel$lon, selected_lats=sel$lat,
  source_lon=lon, source_lat=lat, tile_side=as.integer(p$tile_side),
  longitude_latitude_rule="Native 0.25 deg centres unchanged; land any overlap; density lower bound >1",
  temperature_rule="Daily Ta_C is source daily t2m converted to Celsius; no interpolation")
identity <- fast_object_hash(scientific_identity)

month_batches <- function(dates) {
  key <- format(dates, "%Y-%m")
  split(seq_along(dates), factor(key, levels=unique(key)))
}

check_year_samples <- function(packed_path, source_item, temp, selected, ntry=12L) {
  pc <- ncdf4::nc_open(packed_path)
  oc <- ncdf4::nc_open(source_item$path)
  on.exit({ncdf4::nc_close(pc); ncdf4::nc_close(oc)}, add=TRUE)
  ii <- unique(as.integer(round(seq(1,nrow(selected),length.out=min(ntry,nrow(selected))))))
  dd <- unique(as.integer(round(seq(1,length(source_item$dates),length.out=3L))))
  for (day in dd) {
    coords <- selected[ii, , drop=FALSE]
    temps <- fast_daily_temperature_chunk(temp, list(tile_points=coords), source_item$dates[day])[,1L]
    packed_temp <- ncdf4::ncvar_get(pc,"Ta_C",start=c(1L,day),
      count=c(nrow(selected),1L),collapse_degen=FALSE)[ii,1L]
    if (!isTRUE(all.equal(as.numeric(temps), as.numeric(packed_temp), tolerance=0, check.attributes=FALSE)))
      stop("Packed temperature sample differed at day ",day, " (no rounding tolerated).")
    for (v in era5_required_variables) {
      exact <- vapply(seq_along(ii),function(j) {
        lon_i <- as.integer(coords$selected_lon_index[j]); lat_i <- as.integer(ilat[coords$selected_lat_index[j]])
        read_rectangular_history(oc,v,lon_i,lat_i,day,1L,source_item$offsets[[v]])[1L,1L]
      },numeric(1))
      actual <- ncdf4::ncvar_get(pc,v,start=c(1L,day),
        count=c(nrow(selected),1L),collapse_degen=FALSE)[ii,1L]
      if (!isTRUE(all.equal(exact, as.numeric(actual), tolerance=0, check.attributes=FALSE)))
        stop("Packed data sample differed: ",v," / day ",day)
    }
  }
  TRUE
}

results <- vector("list",length(inspected$files))
for (yy in seq_along(inspected$files)) {
  source_item <- inspected$files[[yy]]
  year <- source_item$year; final <- file.path(l$pack_dir,sprintf("packed_%04d.nc",year))
  receipt <- file.path(l$pack_dir,sprintf("packed_%04d.done.rds",year))
  if (file.exists(final) || file.exists(receipt)) {
    if (!file.exists(final) || !file.exists(receipt)) stop("Orphan packed file/receipt: ",year)
    checked <- readRDS(receipt)
    if (!identical(checked$identity,identity) || !identical(checked$md5,unname(tools::md5sum(final))))
      stop("Existing year is not validated for current mask/inputs. Use a new pack directory: ",year)
    message("Skipping previously committed, MD5-verified year: ",year)
    results[[yy]] <- checked
    next
  }
  message("Packing year ",year," (",yy,"/",length(inspected$files),")...")
  tmp <- paste0(final,".partial")
  if (file.exists(tmp)) {
    warning("Removing incomplete uncommitted pack (year ",year,"): ",tmp)
    unlink(tmp)
  }
  src_nc <- ncdf4::nc_open(source_item$path)
  sc <- tryCatch({
    time_axis <- nc_axis(src_nc,"time")
    dcell <- ncdf4::ncdim_def("cell","index",as.numeric(sel$pack_row),create_dimvar=TRUE)
    dtime <- ncdf4::ncdim_def("time",time_axis$units,time_axis$values,create_dimvar=TRUE)
    dn <- list(dcell,dtime)
    chunks <- as.integer(c(min(256L,nrow(sel)),min(31L,length(source_item$dates))))
    varlist <- lapply(c(era5_required_variables,"Ta_C"),function(v) ncdf4::ncvar_def(
      name=v, units=if (v=="Ta_C") "degree_Celsius" else "degree_Celsius",
      dim=dn, missval=-9.969209968386869e36, prec="double",
      compression=as.integer(p$netcdf_compression_level), chunksizes=chunks,
      longname=if(v=="Ta_C") "daily t2m at original ERA5 grid centres" else v))
    varlist <- c(varlist, list(ncdf4::ncvar_def("cell_id","1",list(dcell),prec="integer",
      missval=-2147483647L)))
    nc <- ncdf4::nc_create(tmp, varlist, force_v4=TRUE)
    tryCatch({
      ncdf4::ncatt_put(nc,0,"pack_schema","ERA5-indexed-weather-v1")
      ncdf4::ncatt_put(nc,0,"pack_identity",identity)
      ncdf4::ncatt_put(nc,0,"land_rule","any Natural Earth 10m land/minor-island overlap")
      ncdf4::ncatt_put(nc,0,"population_rule","2015 density lower bound >1 per full WGS84 cell km2")
      ncdf4::ncatt_put(nc,0,"source_history_file",basename(source_item$path))
      ncdf4::ncvar_put(nc,"cell_id",as.integer(sel$cell_id))
      # Work on one month at a time: ~130 MB for five double buffers, not ten years in memory.
      for (days in month_batches(source_item$dates)) {
        nd <- length(days); ds <- days[1L]
        dates <- source_item$dates[days]
        matrices <- setNames(lapply(c(era5_required_variables,"Ta_C"),function(x)
          matrix(NA_real_,nrow(sel),nd)),c(era5_required_variables,"Ta_C"))
        # Read this month's daily t2m once. Source chunks cover 720x361x11;
        # reopening the same file for every narrow latitude band repeatedly
        # decompresses those large chunks. The exact same coordinate/date reader
        # is used; only I/O batching changes (transient memory below a few GiB).
        matrices$Ta_C <- fast_daily_temperature_chunk(temp,list(tile_points=sel),dates)
        # Geographic read bands prevent expensive per-cell NetCDF random I/O.
        bands <- split(seq_along(ilat),ceiling(seq_along(ilat)/p$latitude_band))
        for (b in bands) {
          rows <- which(sel$selected_lat_index %in% b)
          if (!length(rows)) next
          lon_src <- seq.int(min(sel$selected_lon_index[rows]),max(sel$selected_lon_index[rows]))
          lat_src <- ilat[seq.int(min(sel$selected_lat_index[rows]),max(sel$selected_lat_index[rows]))]
          local_row <- match(sel$selected_lon_index[rows],lon_src) +
            (match(ilat[sel$selected_lat_index[rows]],lat_src)-1L)*length(lon_src)
          if (anyNA(local_row)) stop("Band index alignment failed.")
          for (v in era5_required_variables) {
            rect <- read_rectangular_history(src_nc,v,lon_src,lat_src,ds,nd,source_item$offsets[[v]])
            matrices[[v]][rows,] <- rect[local_row,,drop=FALSE]
          }
        }
        for (v in names(matrices)) ncdf4::ncvar_put(nc,v,matrices[[v]],
          start=c(1L,ds),count=c(nrow(sel),nd))
        rm(matrices); gc(verbose=FALSE)
        message("  ",format(dates[1L])," to ",format(tail(dates,1L))," committed to open NetCDF")
      }
    }, finally=ncdf4::nc_close(nc))
    TRUE
  },finally=ncdf4::nc_close(src_nc))
  if (!isTRUE(sc)) stop("Packing failed for year: ",year)
  check_year_samples(tmp, source_item, temp, sel)
  if (!file.rename(tmp,final)) stop("Cannot publish completed annual pack: ",final)
  receipt_obj <- list(schema="ERA5-indexed-weather-year-v1",year=year,identity=identity,
    md5=unname(tools::md5sum(final)),bytes=unname(file.info(final)$size),
    n_cells=nrow(sel),n_days=length(source_item$dates),timestamp=Sys.time())
  saveRDS(receipt_obj,receipt)
  results[[yy]] <- receipt_obj
  message("Year ",year," verified and committed: ",round(receipt_obj$bytes/1e9,2)," GB")
}

manifest <- list(schema="ERA5-indexed-weather-v1",identity=identity,mask_md5=mask_md5,
  mask_source=normalizePath(l$mask_rds,winslash="/",mustWork=TRUE),
  n_cells=nrow(sel),n_tiles=nrow(tiles),tile_side=as.integer(p$tile_side),
  latitude_range=as.numeric(p$latitude_range),selected=sel,tiles=tiles,
  lon_axis=lon,lat_axis=lat,lat_indices=ilat,years=p$years,
  dates=lapply(inspected$files,function(x)as.character(x$dates)),
  year_files=results,history_audit=inspected$audit,temperature_audit=temp$audit,
  land_method=mask$fingerprint$land_method,
  area_weight_statement="Within selected mask, original warm02 cos(latitude) weights; binary overlap does NOT measure actual land fraction.",
  population_area_statement="cell_area_km2 is the entire 0.25-degree cell area, not fractional land area.",
  output_values="Source ERA5 thermal-history values and same-day Celsius Ta stored as float64; no grid/time interpolation.")
# Always publish index last, after all ten year receipts are verified.
saveRDS(manifest,file.path(l$pack_dir,"pack_index.rds"),compress="gzip")
utils::write.csv(sel,file.path(l$pack_dir,"selected_cells.csv"),row.names=FALSE)
utils::write.csv(tiles,file.path(l$pack_dir,"tiles.csv"),row.names=FALSE)
message("PREPACK COMPLETE: ",l$pack_dir,"; years=",length(results),
        "; cells=",nrow(sel),"; nonempty tiles=",nrow(tiles))
