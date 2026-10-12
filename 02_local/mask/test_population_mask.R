# Independent ERA5 spatial-mask test; does not fit models or evaluate weather.
# Run with Rscript --vanilla test_population_mask.R [ERA5_ROOT] [OUTPUT_DIR]
# Uses only selected definitions from the existing scripts, never their runners.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else "E:/ERA5/ERA5_workflow"
out <- if (length(args) >= 2L) args[2L] else
  file.path(root, "02_local/mask")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
local_lib <- file.path(root, "02_local/runtime/Rlibrary")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", GDAL_CACHEMAX = "128")
needed <- c("ncdf4", "maps", "terra", "sf", "ggplot2", "ragg", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
terra::terraOptions(tempdir = file.path(out, "_tmp"), memfrac = 0.2, progress = 0)
dir.create(file.path(out, "_tmp"), showWarnings = FALSE)
start_time <- Sys.time()

warm_file <- file.path(root, "03_server_bundle/code/warm02_reference.R")
population_file <- file.path(root, "90_archive/scripts/0.5°_ERA5-Figure1-population.R")
pop_tif <- file.path(root, "01_source_data/07_Population_2015",
                    "GHS_POP_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif")

# Evaluate explicitly named assignments only. In particular, source() is not used:
# the warm script's quoted program and both automatic runners must never execute.
assignment_name <- function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]))
    as.character(x[[2L]]) else ""
}
select_assignments <- function(expressions, names, envir) {
  keys <- vapply(expressions, assignment_name, character(1))
  for (nm in names) {
    hit <- which(keys == nm)
    if (length(hit) != 1L) stop("Expected one source definition: ", nm)
    eval(expressions[[hit]], envir = envir)
  }
}
warm_top <- as.list(parse(warm_file, encoding = "UTF-8"))
warm_key <- vapply(warm_top, assignment_name, character(1))
warm_expr <- warm_top[[which(warm_key == "WARM02_PROGRAM")]][[3L]]
stopifnot(identical(warm_expr[[1L]], as.name("quote")),
          identical(warm_expr[[2L]][[1L]], as.name("{")))
warm_definitions <- as.list(warm_expr[[2L]])[-1L]
e <- new.env(parent = globalenv())
select_assignments(warm_definitions,
  c("config", "normal_axis", "same_coordinate", "nc_axis", "load_land_fraction"), e)
population_definitions <- as.list(parse(population_file, encoding = "UTF-8"))
select_assignments(population_definitions,
  c("p1_abort", "p1_wgs84_rectangle_km2", "p1_longitude_overlap",
    "p1_latitude_overlap", "p1_weighted_band_statistics",
    "p1_collect_weighted_targets", "p1_aggregate_population"), e)

g <- e$config$global
g$input_dir <- file.path(root,"01_source_data/02_Thermal history")
stopifnot(g$stride == 1L, identical(as.numeric(g$lat_range), c(-60, 80)))
history_nc <- file.path(g$input_dir, sprintf("thermal_history_%d.nc", min(g$years)))
nc <- ncdf4::nc_open(history_nc)
lon <- e$nc_axis(nc, "lon")$values
lat <- e$nc_axis(nc, "lat")$values
ncdf4::nc_close(nc)
stopifnot(all(abs(abs(diff(lon)) - 0.25) < 1e-8),
          all(abs(abs(diff(lat)) - 0.25) < 1e-8))
ilon <- seq_along(lon)
ilat <- which(lat >= g$lat_range[1L] & lat <= g$lat_range[2L])
message("Constructing exactly the production land grid: ", length(ilon), " x ", length(ilat))
land <- e$load_land_fraction(g, lon, lat, ilon, ilat)
grid <- land$points
grid$selected_lon_index <- rep(seq_along(ilon), times = length(ilat))
grid$selected_lat_index <- rep(seq_along(ilat), each = length(ilon))
grid$original_center_land <- grid$land

# The requested land rule: any overlap with Natural Earth 1:10m land + minor
# islands. Build cells around the ORIGINAL climate sample centres, never shift
# the ERA5 grid by half a cell. Duplicate the dateline edge for mask construction.
coast_dir <- file.path(out, "coastline")
land_shp <- file.path(coast_dir, "ne_10m_land", "ne_10m_land.shp")
island_shp <- file.path(coast_dir, "ne_10m_minor_islands", "ne_10m_minor_islands.shp")
if (!all(file.exists(c(land_shp, island_shp))))
  stop("Download Natural Earth 10m land and minor-island shapefiles into coastline/ first.")
land_sf <- sf::st_read(land_shp, quiet = TRUE)
island_sf <- sf::st_read(island_shp, quiet = TRUE)
sf::sf_use_s2(FALSE)
land_geometry <- sf::st_sfc(c(sf::st_geometry(land_sf), sf::st_geometry(island_sf)), crs = 4326)
land_sf <- sf::st_sf(geometry = sf::st_make_valid(land_geometry))
land_vector <- terra::vect(land_sf)
template <- terra::rast(ncols = length(ilon) + 1L, nrows = length(ilat),
  xmin = -180.125, xmax = 180.125, ymin = min(lat[ilat]) - 0.125,
  ymax = max(lat[ilat]) + 0.125, crs = "EPSG:4326")
message("Rasterizing Natural Earth 1:10m land + minor islands (any overlap).")
touch_raster <- terra::rasterize(land_vector, template, field = 1,
                                background = 0, touches = TRUE)
center_raster <- terra::rasterize(land_vector, template, field = 1,
                                 background = 0, touches = FALSE)
xy <- cbind(grid$lon_plot, grid$lat)
grid$land <- terra::extract(touch_raster, xy)[, 1L] > 0
grid$ne_center_land <- terra::extract(center_raster, xy)[, 1L] > 0
dateline <- which(grid$lon_plot == -180)
if (length(dateline)) {
  opposite <- cbind(180, grid$lat[dateline])
  grid$land[dateline] <- grid$land[dateline] |
    terra::extract(touch_raster, opposite)[, 1L] > 0
  grid$ne_center_land[dateline] <- grid$ne_center_land[dateline] |
    terra::extract(center_raster, opposite)[, 1L] > 0
}
stopifnot(all(grid$land[grid$ne_center_land]))
land$method <- "Natural Earth 1:10m land v5.1.1 + minor islands v4.1.0; terra rasterize touches=TRUE; dateline halves united"
# This is a binary candidate mask, not a measured fractional land cover.
grid$land_fraction <- NULL
grid$area_weight <- NULL
land_points <- grid[grid$land, , drop = FALSE]
stopifnot(nrow(land_points) > 0L, !anyDuplicated(grid$cell_id))
message("Land-overlap cells: ", nrow(land_points),
        "; same-coastline centre cells: ", sum(grid$ne_center_land),
        "; original-code centre cells: ", sum(grid$original_center_land))

source_paths <- c(warm_script = warm_file, population_script = population_file,
  population_raster = pop_tif, coastline = land_shp, minor_islands = island_shp)
fingerprint <- list(md5 = as.list(tools::md5sum(source_paths)),
  history_coordinate_file = normalizePath(history_nc, winslash = "/"),
  lon = lon, lat = lat[ilat], land_cell_ids = land_points$cell_id,
  density_threshold = 1, comparison = ">", dx = 0.25, dy = 0.25,
  latitude_range = g$lat_range, land_method = land$method)
cache_file <- file.path(out, "population_matching.rds")
if (file.exists(cache_file)) {
  cached <- readRDS(cache_file)
  if (!identical(cached$fingerprint, fingerprint))
    stop("Existing population cache does not match inputs. Use a new output directory.")
  matching <- cached$matching
  message("Reusing verified population aggregation.")
} else {
  message("Aggregating 2015 population counts to native 0.25-degree land cells.")
  matching <- e$p1_aggregate_population(pop_tif, land_points,
    dx = 0.25, dy = 0.25, out_dir = file.path(out, "population_audit"),
    config = list(progress_every = 20L, max_band_mb = 256))
  saveRDS(list(fingerprint = fingerprint, matching = matching), cache_file)
}
p <- matching$table
p$original_center_land <- land_points$original_center_land
p$ne_center_land <- land_points$ne_center_land
stopifnot(identical(p$cell_id, land_points$cell_id),
  max(abs(p$population_density_lower_bound - p$population_lower_bound / p$cell_area_km2)) < 1e-9,
  all(is.finite(p$population_density_lower_bound)), all(p$population_lower_bound >= 0))

# NoData is unknown, not zero. A known-contribution lower bound > 1 suffices
# to pass even when some source pixels are unknown; otherwise keep NA if incomplete.
p$selected <- ifelse(p$population_density_lower_bound > 1, TRUE,
                     ifelse(p$complete, FALSE, NA))
p$mask_status <- ifelse(p$selected %in% TRUE, "retained_gt1",
                       ifelse(p$selected %in% FALSE, "excluded_le1", "unknown_population"))
grid$mask_code <- 0L
grid$mask_code[grid$land] <- match(p$mask_status,
  c("retained_gt1", "excluded_le1", "unknown_population"))
grid$selected <- FALSE
grid$selected[grid$land] <- p$selected
grid$population_density <- NA_real_
grid$population_density[grid$land] <- p$population_density
grid$population_density_lower_bound <- NA_real_
grid$population_density_lower_bound[grid$land] <- p$population_density_lower_bound
grid$population_complete <- NA
grid$population_complete[grid$land] <- p$complete
stopifnot(sum(grid$mask_code == 1L) == sum(p$selected %in% TRUE),
          all(!grid$selected[!grid$land]),
          sum(table(factor(grid$mask_code, levels = 0:3))) == nrow(grid))

write_gzip_csv <- function(x, path) {
  con <- gzfile(path, "wt", encoding = "UTF-8")
  on.exit(close(con))
  write.csv(x, con, row.names = FALSE, na = "")
}
write_gzip_csv(p, file.path(out, "mask_land_cells.csv.gz"))
retained <- p[p$selected %in% TRUE, , drop = FALSE]
write_gzip_csv(retained, file.path(out, "mask_selected_cells.csv.gz"))
saveRDS(list(grid = grid, land_population = p, fingerprint = fingerprint,
             mask_codes = c(no_mapped_land_overlap = 0L, retained_gt1 = 1L,
                            excluded_le1 = 2L, unknown_population = 3L)),
        file.path(out, "mask_0p25deg_pop_gt1.rds"))

# Preserve original coordinate ordering and cell_id = lon-index + (lat-index-1)*nlon.
nc_path <- file.path(out, "mask_0p25deg_pop_gt1.nc")
d_lon <- ncdf4::ncdim_def("longitude", "degrees_east", lon[ilon])
d_lat <- ncdf4::ncdim_def("latitude", "degrees_north", lat[ilat])
dims <- list(d_lon, d_lat)
vars <- list(
  ncdf4::ncvar_def("mask_status", "1", dims, missval = -32767L, prec = "short",
                  longname = "0 no mapped land overlap; 1 retained; 2 excluded <=1; 3 unknown", compression = 4),
  ncdf4::ncvar_def("selected", "1", dims, missval = -32767L, prec = "short",
                  longname = "1 retained, 0 excluded/ocean, missing unknown population", compression = 4),
  ncdf4::ncvar_def("population_density_lower_bound", "persons km-2", dims,
                  missval = -9999, prec = "double", compression = 4))
nco <- ncdf4::nc_create(nc_path, vars, force_v4 = TRUE)
shape <- c(length(ilon), length(ilat))
ncdf4::ncvar_put(nco, "mask_status", array(grid$mask_code, shape))
ncdf4::ncvar_put(nco, "selected", array(as.integer(grid$selected), shape))
ncdf4::ncvar_put(nco, "population_density_lower_bound",
                array(grid$population_density_lower_bound, shape))
ncdf4::ncatt_put(nco, 0, "population_source", basename(pop_tif))
ncdf4::ncatt_put(nco, 0, "land_mask_method", land$method)
ncdf4::ncatt_put(nco, 0, "density_definition",
  "Conservatively allocated population / whole 0.25-degree WGS84 ellipsoid cell area; NOT land-only area.")
ncdf4::ncatt_put(nco, 0, "selection_rule",
  "Native centre latitude -60..80 inclusive; any Natural Earth 10m land/minor-island overlap; known population density lower bound strictly >1.")
ncdf4::ncatt_put(nco, 0, "missing_policy",
  "Unknown source pixels never imputed as zero; incomplete cells not demonstrably >1 remain unknown.")
ncdf4::ncatt_put(nco, 0, "scope", "Spatial mask only; no daily weather/model support filtering.")
ncdf4::nc_close(nco)
check <- ncdf4::nc_open(nc_path)
stopifnot(identical(as.integer(ncdf4::ncvar_get(check, "mask_status")), grid$mask_code),
          identical(as.integer(ncdf4::ncvar_get(check, "selected")), as.integer(grid$selected)))
ncdf4::nc_close(check)

n_keep <- nrow(retained)
n_unknown <- sum(is.na(p$selected))
counts <- data.frame(category = c("all_grid_centres", "no_land_overlap",
  "land_overlap_cells", "land_retained_gt1", "land_excluded_le1", "land_unknown_population",
  "retained_with_complete_population", "retained_using_lower_bound"),
  n_cells = c(nrow(grid), sum(!grid$land), nrow(p), n_keep,
    sum(p$selected %in% FALSE), n_unknown, sum(retained$complete), sum(!retained$complete)))
write.csv(counts, file.path(out, "mask_counts.csv"), row.names = FALSE)
tile_id <- function(x) paste((x$selected_lon_index - 1L) %/% g$tile_side,
                             (x$selected_lat_index - 1L) %/% g$tile_side, sep = ":")
summary <- list(latitude_range = g$lat_range, resolution_degrees = 0.25,
  population_epoch = 2015L, density_threshold_persons_per_km2 = 1,
  threshold_comparison = ">", full_grid_cells = nrow(grid), land_cells = nrow(p),
  retained_cells = n_keep, excluded_low_density_cells = sum(p$selected %in% FALSE),
  unresolved_population_cells = n_unknown,
  retained_pct_of_land_cells = 100 * n_keep / nrow(p),
  land_cell_reduction_pct = 100 * (1 - n_keep / nrow(p)),
  known_population_retained_pct = 100 * sum(retained$population_lower_bound) / sum(p$population_lower_bound),
  original_nonempty_tiles = length(unique(tile_id(land_points))),
  retained_nonempty_tiles = length(unique(tile_id(grid[grid$mask_code == 1L, ]))),
  retained_with_complete_population = sum(retained$complete),
  retained_using_lower_bound = sum(!retained$complete),
  original_code_center_land_cells = sum(grid$original_center_land),
  same_coastline_center_land_cells = sum(grid$ne_center_land),
  same_coastline_center_retained_gt1 = sum(retained$ne_center_land),
  retained_added_by_any_overlap = sum(!retained$ne_center_land),
  original_center_retained_within_any_overlap = sum(retained$original_center_land),
  original_center_outside_new_coastline = sum(grid$original_center_land & !grid$land),
  land_method = land$method,
  max_population_conservation_difference = matching$audit$max_abs_population_conservation_difference,
  spatial_only = TRUE, elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
jsonlite::write_json(summary, file.path(out, "mask_summary.json"), pretty = TRUE, auto_unbox = TRUE)

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
message("Building exact coastline clipping for map display.")
land_union <- sf::st_union(sf::st_geometry(land_sf))
window <- sf::st_as_sfc(sf::st_bbox(c(xmin = -180, ymin = -60.125,
                                    xmax = 180, ymax = 80.125), crs = sf::st_crs(4326)))
ocean_sf <- sf::st_sf(geometry = sf::st_difference(window, land_union))
display <- data.frame(lon = grid$lon_plot, lat = grid$lat,
                      status = factor(grid$mask_code, levels = 0:3))
# Display the eastern half of the centre at -180 also at +180; counts stay unique.
edge <- display[display$lon == -180, ]
edge$lon <- 180
display <- rbind(display, edge)
map <- ggplot2::ggplot(display, ggplot2::aes(lon, lat)) +
  ggplot2::geom_raster(ggplot2::aes(fill = status)) +
  ggplot2::geom_sf(data = ocean_sf, inherit.aes = FALSE, fill = "#EDF3F7", colour = NA) +
  ggplot2::geom_sf(data = land_sf, inherit.aes = FALSE, fill = NA,
                   linewidth = 0.10, colour = "#68747A") +
  ggplot2::scale_fill_manual(values = c("0" = "#EDF3F7", "1" = "#087F8C",
                                       "2" = "#D9DEE1", "3" = "#D79A53"),
    breaks = c("1", "2", "3"), drop = FALSE, name = NULL,
    labels = c(sprintf("保留：人口密度 >1（%s 格）", fmt(n_keep)),
               sprintf("排除：人口密度 ≤1（%s 格）", fmt(sum(p$selected %in% FALSE))),
               sprintf("人口数据不足（%s 格）", fmt(n_unknown)))) +
  ggplot2::scale_x_continuous(breaks = seq(-180, 180, 60)) +
  ggplot2::scale_y_continuous(breaks = seq(-60, 80, 20)) +
  ggplot2::coord_sf(crs = 4326, xlim = c(-180, 180), ylim = c(-60.125, 80.125), expand = FALSE) +
  ggplot2::labs(title = "ERA5 0.25° 陆地人口筛选",
    subtitle = sprintf("2015 年人口密度 >1 人/km²  |  保留 %s / %s 个陆地网格（%.1f%%）",
                       fmt(n_keep), fmt(nrow(p)), 100 * n_keep / nrow(p)),
    x = "经度", y = "纬度",
    caption = paste("范围：60°S–80°N；网格内有陆地即纳入候选，绘图沿海岸线裁切。人口密度分母为完整网格面积。",
      "数据：GHS-POP R2023A（2015）；缺失人口保留未知状态。此图尚未筛选逐日气象或模型支持范围。", sep = "\n")) +
  ggplot2::theme_minimal(base_size = 12, base_family = "Microsoft YaHei") +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    legend.position = "bottom", legend.text = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(face = "bold", size = 19),
    plot.subtitle = ggplot2::element_text(size = 12, margin = ggplot2::margin(b = 12)),
    plot.caption = ggplot2::element_text(hjust = 0, size = 9, lineheight = 1.4),
    plot.margin = ggplot2::margin(15, 18, 15, 15))
ggplot2::ggsave(file.path(out, "mask_global.png"), map, width = 15, height = 8,
                dpi = 220, device = ragg::agg_png, bg = "white")
ggplot2::ggsave(file.path(out, "mask_global.pdf"), map, width = 15, height = 8,
                device = grDevices::cairo_pdf, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(c(
  "独立空间 mask 测试（不运行原脚本主程序，不计算模型预测）",
  sprintf("网格中心纬度范围：%s；网格间隔：0.25°；人口年份：2015。", paste(g$lat_range, collapse = " 至 ")),
  sprintf("全部格点 %s；陆地格点 %s；确认 >1 人/km² 的陆地格点 %s。", fmt(nrow(grid)), fmt(nrow(p)), fmt(n_keep)),
  sprintf("确认 ≤1 的陆地格点 %s；人口数据不足 %s。", fmt(sum(p$selected %in% FALSE)), fmt(n_unknown)),
  sprintf("保留陆地格点的 %.4f%%；格点数量减少 %.4f%%（不代表运行时间同比减少）。", summary$retained_pct_of_land_cells, summary$land_cell_reduction_pct),
  "方法：读取新版 warm02 的 config 和原始坐标。原中心陆地规则仅作对照；主 mask 改用 Natural Earth 1:10m 陆地及小岛多边形，与网格有交集即入选陆地候选。",
  "仅复用人口脚本中可接受任意目标分辨率的保守面积分配函数，不运行其半度主程序。",
  "人口原始栅格值是每个像元的人数，不是人口密度。依据真实像元边界和 WGS84 椭球面积分配人数，除以完整 0.25° 网格面积。",
  "没有按网格内陆地面积校正分母。海岸网格也使用完整网格面积。",
  "陆地使用 terra::rasterize(touches=TRUE)，跨日期变更线网格合并两半；不是 ERA5 官方陆海比例 mask。结果受边界数据分辨率影响，并非保证包含世界上每一个微小岛礁。",
  "绘图用真实矢量海岸线遮去海面，保持原始 0.25° 计算网格；没有对科学数据插值或平滑。",
  "NoData 和源数据范围外区域不视为零人口。已知人口贡献足以令密度 >1 的格点可确认保留；否则人口覆盖不全时标为未知。",
  "严格使用 >1，不是旧人口脚本的 >=1。气候/模型逐日支持、有效天数等条件尚未检查。",
  "NetCDF 保留源经纬度排序。mask_status: 0 无陆地交集；1 保留；2 低密度排除；3 人口未知。",
  "selected: 1 保留；0 排除/海洋；缺失为人口未知。RDS 和 CSV 保留原生产 cell_id。",
  "Windows 重跑：在 PowerShell 中执行本目录的 run_mask.ps1，它配置 UTF-8 和临时目录后运行测试。",
  "Windows 若 R 无法创建临时目录，请先将 TMPDIR、TMP、TEMP 设为本测试目录下的 _tmp。",
  "依赖 ncdf4、maps、terra、sf、ggplot2、ragg、jsonlite；本次新增包安装在 _Rlibrary，不修改系统库。",
  "人口数据说明：https://human-settlement.emergency.copernicus.eu/documents/GHSL_Data_Package_2023.pdf",
  "海岸线来源：https://www.naturalearthdata.com/downloads/10m-physical-vectors/"
), file.path(out, "README.txt"), useBytes = TRUE)
print(counts, row.names = FALSE)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
message("Completed: ", normalizePath(out, winslash = "/"))
