# Independent geometry and artifact checks for the completed spatial mask.
if (.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
root <- "E:/ERA5/ERA5_workflow"
out <- file.path(root,"02_local/mask")
.libPaths(c(file.path(root, "02_local/runtime/Rlibrary"), .libPaths()))
x <- readRDS(file.path(out, "mask_0p25deg_pop_gt1.rds"))
g <- x$grid
p <- x$land_population
sf::sf_use_s2(FALSE)
land <- sf::st_read(file.path(out, "coastline/ne_10m_land/ne_10m_land.shp"), quiet = TRUE)
islands <- sf::st_read(file.path(out, "coastline/ne_10m_minor_islands/ne_10m_minor_islands.shp"), quiet = TRUE)
coast <- sf::st_make_valid(sf::st_sfc(c(sf::st_geometry(land), sf::st_geometry(islands)), crs = 4326))
set.seed(20261012)
ix <- unique(c(sample(which(g$land & !g$ne_center_land), 100),
  sample(which(g$ne_center_land), 50), sample(which(!g$land), 50), which(g$lon_plot == -180)))
rectangle <- function(w, e, s, n) matrix(c(w,s,e,s,e,n,w,n,w,s), ncol = 2, byrow = TRUE)
cells <- sf::st_sfc(lapply(ix, function(i) {
  w <- g$lon_plot[i] - 0.125; e <- g$lon_plot[i] + 0.125
  s <- g$lat[i] - 0.125; n <- g$lat[i] + 0.125
  if (w < -180) sf::st_multipolygon(list(list(rectangle(w + 360,180,s,n)),
                                        list(rectangle(-180,e,s,n)))) else
    sf::st_polygon(list(rectangle(w,e,s,n)))
}), crs = 4326)
reference <- lengths(sf::st_intersects(cells, coast)) > 0L
stopifnot(identical(reference, g$land[ix]))
nc <- ncdf4::nc_open(file.path(out, "mask_0p25deg_pop_gt1.nc"))
stopifnot(identical(as.integer(ncdf4::ncvar_get(nc, "mask_status")), g$mask_code),
          identical(as.integer(ncdf4::ncvar_get(nc, "selected")), as.integer(g$selected)),
          all(p$selected == (p$population_density > 1)),
          sum(g$selected) == 104344L, all(p$complete))
ncdf4::nc_close(nc)
csv <- read.csv(gzfile(file.path(out, "mask_selected_cells.csv.gz")))
stopifnot(nrow(csv) == sum(g$selected), identical(csv$cell_id, g$cell_id[g$selected]))
fingerprints <- unlist(x$fingerprint$md5)
# Preserve the frozen mask and its original provenance; resolve only file locations.
migration <- read.csv(file.path(root,"02_local/audit/path_migration.csv"),fileEncoding="UTF-8")
slash <- function(x) gsub("\\", "/", x, fixed=TRUE)
old <- slash(migration$original_path); new <- slash(migration$current_path)
relocated <- vapply(slash(names(fingerprints)),function(path) {
  hits <- which(path==old | startsWith(path,paste0(old,"/")))
  if (!length(hits)) return(path)
  hit <- hits[which.max(nchar(old[hits]))]
  paste0(new[hit],substring(path,nchar(old[hit])+1L))
},character(1))
stopifnot(all(file.exists(relocated)),
          identical(unname(tools::md5sum(relocated)), unname(fingerprints)))
audit <- readRDS(file.path(out, "population_audit/population_source_audit.rds"))
stopifnot(audit$status == "COMPLETED", audit$max_abs_population_conservation_difference < 1e-6)
result <- list(status = "PASS", independent_polygon_checks = length(ix),
  dateline_cells_checked = sum(g$lon_plot[ix] == -180),
  raster_netcdf_rds_csv_counts_agree = TRUE, source_file_hashes_unchanged = TRUE,
  population_conservation_max_abs_persons = audit$max_abs_population_conservation_difference)
jsonlite::write_json(result, file.path(out, "validation.json"), pretty = TRUE, auto_unbox = TRUE)
print(result)
