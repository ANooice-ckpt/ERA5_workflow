# 0.5° ERA5-Figure1-population.R — completed warm02 run -> new GHS-POP matching -> Figure 4.
# 先把 HALF_DEGREE_RUN_DIR 指向新生成、已完成的 stride=2 warm02 run_... 目录，再 Source。
# 默认首次运行会重新匹配 2015 年 GHS-POP，并绘制人口限制后的 Figure 4。
# 如果已完成匹配、只想重画三幅图，把 HALF_DEGREE_MODE 改为 "three_maps"。
# 0.5° 是原生 0.25° ERA5 点每 2 个坐标抽样，并非 0.5° 区域的气温平均。
# 上左a: delta_M；上右b: delta_N；下方居中c: delta_U。三幅地图完全等大。
# 论文图号为Figure 4；输出文件名沿用Figure1_thermal_response以兼容已有流程。
# 保留a/b/c图号，不添加顶部说明或英文子图标题；共同色带与灰色说明上下居中排列。
# 本脚本不拟合热史模型、不重复 warm02 全球 ERA5 投影，也不读取逐日能耗缓存。
# 正式分析严格校验输入 run 的分辨率、stride、指纹和格点清单。

HALF_DEGREE_RUN_DIR <- NULL  # 必填：新 0.5° warm02 的 run_... 绝对路径；绝不可填写旧 1° run
HALF_DEGREE_MODE <- "population_then_maps" # 可改为 "population" 或 "three_maps"

THREE_MAP_REVISION <- "ERA5-Figure1-population-0.5deg-v1-components-boundary"
THREE_MAP_CONFIG <- list(
  RUN_DIR = HALF_DEGREE_RUN_DIR,
  FIGURE_DIR = NULL, # 自动选择此run内最新成功的Figure1_population_...目录
  main_domain = "pop_ge1",
  install_missing = TRUE,
  verify_saved_grid_hash = TRUE,
  draw_map_frame = TRUE,
  frame_linewidth_mm = 0.14, # 约0.4 pt；仅地图面板框，不套整个复合图
  frame_colour = "#666666",
  show_metric_titles = FALSE, # 若设TRUE，按ΔM/ΔN/ΔU顺序显示标题；默认全部取消
  show_panel_letters = TRUE,
  panel_letter_size_pt = 8,
  panel_letter_space_mm = 3.5,
  figure_width_mm = 180,
  outer_margin_mm = 1.5,
  column_gap_mm = 4,
  row_gap_mm = 4,
  map_margin_mm = 0.7,
  legend_height_mm = 19,
  colourbar_width_mm = 78,
  colourbar_height_mm = 2.6,
  base_font_size = 7,
  font_family = "sans",
  dpi = 600L,
  save_tiff = TRUE,
  save_individual_maps = FALSE
)

# Population matching, posterior aggregation, and map functions are self-contained below.

POPULATION_REVISION <- "ERA5-Figure1-population-0.5deg-v1-real-extent"
POPULATION_CONFIG <- list(
  RUN_DIR = HALF_DEGREE_RUN_DIR,
  population_dir = "I:/2026-9 ERA5/07_Population_2015",
  population_tif = NULL, # 自动查找唯一的2015/4326/30ss人口TIF；不读.ovr或压缩包
  matching_rds = NULL,  # 新0.5°首次运行必须重算；仅重用同一run指纹与0.5°分辨率匹配RDS
  main_domain = "pop_ge1",
  install_missing = TRUE,
  save_tiff = TRUE,
  dpi = 600L,
  save_all_domain_figures = TRUE,
  progress_every = 20L,
  max_band_mb = 512,
  gdal_cache_mb = 128L,
  strict_product = TRUE # 正式0.5度分析始终为TRUE；小栅格函数测试应另行调用底层函数
)


FIGURE1_REVISION <- "Figure1-0.5deg-v1-population-real-extent"
FIGURE1_CONFIG <- list(
  RUN_DIR = HALF_DEGREE_RUN_DIR,
  install_missing = TRUE,
  verify_checkpoints = TRUE, # reads saved checkpoint RDS only; never daily .f32 payloads
  save_tiff = TRUE,
  dpi = 600L,
  figure_width_mm = 180,
  figure_height_mm = 125,
  base_font_size = 7,
  font_family = "sans",
  map_scale_limit = NULL, # shared full symmetric range; no silent clipping
  save_individual_panels = FALSE,
  main_domain = "original",
  root_range_thresholds = c(0.50, 0.80, 0.95),
  history_fraction_thresholds = c(0.25, 0.50, 0.75),
  save_all_domain_figures = TRUE,
  external_domain_csv = NULL,
  external_domain_basis = NULL,
  elevation_csv = NULL,
  elevation_thresholds_m = numeric(0), # 如 c(2000,3000)，仅在填入海拔CSV和依据后运行
  elevation_basis = NULL
)

f1_stop <- function(...) stop(..., call. = FALSE)
f1_assert <- function(ok, ...) if (!isTRUE(ok)) f1_stop(...)
f1_close <- function(a, b, label, tol = 1e-8) {
  a <- as.numeric(a); b <- as.numeric(b)
  f1_assert(length(a) == length(b), "Length mismatch: ", label)
  f1_assert(identical(is.na(a), is.na(b)), "Missingness mismatch: ", label)
  good <- !is.na(a)
  f1_assert(all(is.finite(a[good])) && all(is.finite(b[good])), "Nonfinite values: ", label)
  err <- if (any(good)) max(abs(a[good] - b[good]) / pmax(1, abs(a[good]), abs(b[good]))) else 0
  f1_assert(err <= tol, "Numerical mismatch: ", label, "; relative/absolute scaled error = ", signif(err, 6))
  invisible(err)
}
f1_read <- function(path) {
  f1_assert(file.exists(path), "Required saved output is missing: ", path)
  tryCatch(readRDS(path), error = function(e) f1_stop("Cannot read ", path, ": ", conditionMessage(e)))
}
f1_packages <- function(config) {
  required <- c("ggplot2", "maps", "patchwork", "scales")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && isTRUE(config$install_missing)) {
    message("Installing missing plotting packages: ", paste(missing, collapse = ", "))
    utils::install.packages(missing, repos = "https://cloud.r-project.org",
      type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType"))
  }
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  f1_assert(!length(missing), "Install these packages first: ", paste(missing, collapse = ", "))
  f1_assert(utils::packageVersion("ggplot2") >= "3.4.0", "ggplot2 >= 3.4.0 is required. Update it in a clean R session.")
  invisible(required)
}
f1_resolve_run <- function(config) {
  run <- config$RUN_DIR
  f1_assert(is.character(run) && length(run) == 1L && !is.na(run) && nzchar(run),
    "Set HALF_DEGREE_RUN_DIR to the completed NEW stride=2 warm02 run_... directory. Automatic latest_run selection is disabled.")
  f1_assert(dir.exists(run), "Run directory not found: ", run)
  run <- normalizePath(run, winslash = "/", mustWork = TRUE)
  if (basename(run) == "global" && file.exists(file.path(run, "global_run_metadata.rds"))) run <- dirname(run)
  f1_assert(dir.exists(file.path(run, "global")), "Selected folder must contain global/: ", run)
  run
}

f1_verify_half_degree_run <- function(b) {
  resolution <- b$meta$grid_resolution[c("lon", "lat")]
  f1_assert(length(resolution) == 2L && all(is.finite(resolution)), "Invalid saved map resolution.")
  f1_close(resolution, c(0.5, 0.5), "0.5-degree representative grid spacing")
  stride <- b$meta$global_config$spatial_stride
  if (is.null(stride)) stride <- b$meta$global_config$stride
  f1_assert(length(stride) == 1L && !is.na(stride) && stride == 2L,
    "This script requires the native 0.25-degree ERA5 grid sampled with stride=2. Do not use the old stride=4 run.")
  f1_assert(identical(as.integer(stride), as.integer(b$manifest$g$spatial_stride)),
    "Saved metadata and run manifest disagree on the ERA5 spatial stride.")
  f1_assert(length(b$land$cell_id) == b$manifest$total_land_points &&
    length(b$land$cell_id) == b$meta$fast_execution$total_land_points,
    "Saved half-degree grid-point count differs among selected-land inventory, manifest, and metadata.")
  invisible(resolution)
}
f1_draw_ids <- function(x, label) {
  f1_assert(is.data.frame(x) && all(c("m1", "m0") %in% names(x)), "Missing paired draw IDs: ", label)
  a <- x[c("m1", "m0")]
  for (nm in names(a)) {
    f1_assert(is.numeric(a[[nm]]) && all(is.finite(a[[nm]])) &&
      all(a[[nm]] > 0 & a[[nm]] == floor(a[[nm]])) && !anyDuplicated(a[[nm]]), "Invalid draw IDs: ", label, " / ", nm)
    a[[nm]] <- as.integer(a[[nm]])
  }
  rownames(a) <- NULL
  a
}
f1_quantiles <- function(x) {
  if (!length(x) || all(is.na(x))) return(c(lower = NA_real_, median = NA_real_, upper = NA_real_, prob_positive = NA_real_, prob_negative = NA_real_))
  f1_assert(all(is.finite(x)), "Partially missing posterior draws are not allowed.")
  c(stats::setNames(stats::quantile(x, c(.025, .5, .975), type = 8, names = FALSE), c("lower", "median", "upper")),
    prob_positive = mean(x > 0), prob_negative = mean(x < 0))
}
f1_load_identity <- function(run) {
  gdir <- file.path(run, "global")
  paths <- file.path(gdir, c("global_run_metadata.rds", "numerical_completion.rds", "fast_manifest.rds",
    "warm_boundary_area_day_posterior_draws.rds", "warm_boundary_area_day_posterior_summary.rds", "selected_land_points.csv"))
  f1_assert(all(file.exists(paths)), "Incomplete numerical outputs. Missing: ", paste(paths[!file.exists(paths)], collapse = ", "))
  meta <- f1_read(paths[1]); complete <- f1_read(paths[2]); manifest <- f1_read(paths[3])
  f1_assert(identical(complete$status, "COMPLETED"), "Numerical finalization has not completed; do not run this script yet.")
  f1_assert(identical(manifest$schema, "ERA5-warm02-run-v1") && identical(meta$checkpoint_schema, "ERA5-warm02-tile-v1"), "This script requires completed ERA5-warm02 output, not warm01 or the older global01 run.")
  fp <- manifest$fingerprint
  f1_assert(is.character(fp) && length(fp) == 1L && nzchar(fp) &&
    identical(fp, complete$fingerprint) && identical(fp, meta$fast_execution$fingerprint), "Run fingerprint mismatch.")
  f1_assert(isFALSE(meta$pilot_not_global) && is.null(meta$global_config$max_cells), "Pilot outputs cannot be used as the formal global Figure 1.")
  f1_assert(length(meta$fixed_set) == 1L && identical(names(meta$fixed_set), "P50"), "P50-only metadata contract failed.")
  f1_assert(is.numeric(meta$set_bounds) && length(meta$set_bounds) == 2L &&
    all(is.finite(meta$set_bounds)) && diff(meta$set_bounds) > 0,
    "Invalid saved marginal training SET P1-P99 bounds.")
  f1_close(meta$set_bounds, manifest$set_bounds, "Metadata/manifest SET bounds")
  f1_assert(length(meta$periods) > 1L && length(meta$years) > 0L, "Missing time domain.")
  period <- tail(as.character(meta$periods), 1L)
  expected <- paste0(min(meta$years), "-", max(meta$years))
  f1_assert(identical(period, expected), "The final period must be the complete study interval.")
  draw <- f1_read(paths[4]); saved_summary <- f1_read(paths[5])
  ids <- f1_draw_ids(draw$paired_draw_ids, "posterior draws")
  f1_assert(identical(ids, f1_draw_ids(meta$paired_draw_ids, "metadata")) &&
    identical(ids, f1_draw_ids(manifest$paired_draw_ids, "manifest")), "M1/M0 draw pairing differs between saved outputs.")
  arr <- draw$estimates
  metrics <- c("delta_M", "delta_N", "delta_U", "delta_TSV_P50", "U_M1", "U_M0", "U_pair_in_SET_range", "absolute_delta_U")
  f1_assert(is.numeric(arr) && identical(as.integer(dim(arr)), c(as.integer(length(meta$periods)), 8L, as.integer(nrow(ids)))) &&
    identical(dimnames(arr)[[1]], as.character(meta$periods)) && identical(dimnames(arr)[[2]], metrics), "Unexpected global posterior array schema.")
  f1_assert(nrow(ids) == meta$global_draws_actual, "Draw count mismatch.")
  f1_assert(all(is.finite(arr) | is.na(arr)) && !any(is.nan(arr)), "Nonfinite posterior array.")
  f1_close(arr[, "delta_U", ], arr[, "delta_M", ] + arr[, "delta_N", ], "delta_U = delta_M + delta_N")
  f1_close(arr[, "delta_U", ], arr[, "U_M1", ] - arr[, "U_M0", ], "delta_U = U_M1 - U_M0")
  f1_assert(all(arr[, "absolute_delta_U", ] >= abs(arr[, "delta_U", ]) - 1e-8, na.rm = TRUE), "Daily absolute-shift identity failed.")
  rows <- list(); k <- 0L
  for (p in seq_along(meta$periods)) for (m in seq_along(metrics)) {
    k <- k + 1L; q <- f1_quantiles(arr[p, m, ])
    rows[[k]] <- data.frame(period = meta$periods[p], metric = metrics[m],
      lower = q["lower"], median = q["median"], upper = q["upper"],
      prob_positive = q["prob_positive"], prob_negative = q["prob_negative"],
      eligible_grid_points = draw$aggregate_cell_counts[p],
      weighted_supported_cell_day_denominator = draw$aggregate_denominator[p], row.names = NULL)
  }
  summary <- do.call(rbind, rows)
  key <- function(x) paste(x$period, x$metric, sep = "|")
  f1_assert(is.data.frame(saved_summary) && !anyDuplicated(key(saved_summary)) && setequal(key(summary), key(saved_summary)), "Saved posterior summary rows differ.")
  jj <- match(key(summary), key(saved_summary))
  for (nm in c("lower", "median", "upper", "prob_positive", "prob_negative", "eligible_grid_points", "weighted_supported_cell_day_denominator"))
    f1_close(summary[[nm]], saved_summary[[nm]][jj], paste("Global summary", nm))
  land <- utils::read.csv(paths[6], check.names = FALSE)
  f1_assert(all(c("cell_id", "lon_plot", "lat", "area_weight", "land_fraction") %in% names(land)) && !anyDuplicated(land$cell_id), "Invalid selected land-point inventory.")
  f1_assert(nrow(land) == manifest$total_land_points && nrow(land) == meta$fast_execution$total_land_points, "Land-point count mismatch.")
  list(run = run, directory = gdir, meta = meta, manifest = manifest, draws = draw,
    summary = summary, period = period, land = land, input_files = paths, fingerprint = fp)
}

f1_validate_grid <- function(d, b) {
  metrics <- dimnames(b$draws$estimates)[[2]]
  qcols <- as.vector(outer(metrics, c("lower", "median", "upper", "prob_positive", "prob_negative"), paste, sep = "_"))
  req <- c("cell_id", "lon", "lat", "area_weight", "land_fraction", "period", "fixed_set_name", "fixed_set",
    "n_expected_days", "n_valid_history_days", "n_common_support_days", "n_estimation_days", "date_coverage",
    "common_support_fraction", "eligible_for_area_day_summary", "eligible_for_supported_day_map",
    "SET_range_lower", "SET_range_upper", qcols)
  f1_assert(is.data.frame(d) && nrow(d) > 0 && all(req %in% names(d)), "Missing grid-summary fields.")
  f1_assert(!anyDuplicated(d$cell_id) && setequal(d$cell_id, b$land$cell_id), "Duplicated/missing/mixed grid cell IDs.")
  f1_assert(all(d$period == b$period) && all(d$fixed_set_name == "P50") &&
    all(d$fixed_set == unname(b$meta$fixed_set)), "Wrong period or fixed-SET reference in grid output.")
  f1_close(d$SET_range_lower, rep(b$meta$set_bounds[1L], nrow(d)), "Grid lower SET bound")
  f1_close(d$SET_range_upper, rep(b$meta$set_bounds[2L], nrow(d)), "Grid upper SET bound")
  j <- match(d$cell_id, b$land$cell_id)
  for (nm in c("lat", "area_weight", "land_fraction")) f1_close(d[[nm]], b$land[[nm]][j], paste("Land inventory", nm))
  f1_close(d$lon, b$land$lon_plot[j], "Land inventory longitude")
  f1_assert(all(is.finite(d$area_weight) & d$area_weight > 0), "Invalid grid area weights.")
  for (nm in c("n_expected_days", "n_valid_history_days", "n_common_support_days", "n_estimation_days"))
    f1_assert(is.numeric(d[[nm]]) && all(is.finite(d[[nm]]) & d[[nm]] >= 0 & d[[nm]] == floor(d[[nm]])), "Invalid day counts: ", nm)
  f1_assert(all(d$n_expected_days == b$meta$n_expected_total_days) &&
    all(d$n_estimation_days == d$n_common_support_days) &&
    all(d$n_common_support_days <= d$n_valid_history_days & d$n_valid_history_days <= d$n_expected_days), "Inconsistent day counts.")
  f1_close(d$date_coverage, d$n_valid_history_days / d$n_expected_days, "Finite-history coverage")
  f1_close(d$common_support_fraction, d$n_common_support_days / d$n_expected_days, "History support fraction")
  g <- b$meta$global_config
  f1_assert(all(c("min_date_coverage", "min_supported_days", "min_supported_fraction") %in% names(g)), "Saved eligibility thresholds are missing.")
  eligible <- d$date_coverage >= g$min_date_coverage & d$n_common_support_days >= g$min_supported_days &
    d$common_support_fraction >= g$min_supported_fraction
  f1_assert(!anyNA(d$eligible_for_area_day_summary) &&
    identical(as.logical(d$eligible_for_area_day_summary), eligible) &&
    identical(as.logical(d$eligible_for_supported_day_map), eligible), "Saved map/summary eligibility differs from original rules.")
  f1_assert(any(eligible), "No eligible cells. Inspect support diagnostics in the completed run.")
  for (m in metrics) {
    lower <- d[[paste0(m, "_lower")]]; middle <- d[[paste0(m, "_median")]]; upper <- d[[paste0(m, "_upper")]]
    supported <- d$n_common_support_days > 0
    f1_assert(all(is.finite(lower[supported])) && all(is.finite(middle[supported])) && all(is.finite(upper[supported])), "Missing supported-cell summaries: ", m)
    f1_assert(all(is.na(lower[!supported])) && all(is.na(middle[!supported])) && all(is.na(upper[!supported])), "Unsupported cells must have NA estimates: ", m)
    f1_assert(all(lower[supported] <= middle[supported] + 1e-10 & middle[supported] <= upper[supported] + 1e-10), "Invalid credible interval ordering: ", m)
    p <- d[[paste0(m, "_prob_positive")]]; n <- d[[paste0(m, "_prob_negative")]]
    f1_assert(all(is.finite(p[supported]) & p[supported] >= 0 & p[supported] <= 1) &&
      all(is.finite(n[supported]) & n[supported] >= 0 & n[supported] <= 1) &&
      all(p[supported] + n[supported] <= 1 + 1e-10), "Invalid posterior sign probabilities: ", m)
  }
  f1_assert(all(d$U_pair_in_SET_range_median[d$n_common_support_days > 0] >= 0 &
    d$U_pair_in_SET_range_median[d$n_common_support_days > 0] <= 1), "Invalid root-range diagnostic.")
  p <- match(b$period, b$meta$periods)
  f1_close(sum(d$area_weight[eligible] * d$n_common_support_days[eligible]), b$draws$aggregate_denominator[p], "Global area-day denominator")
  f1_close(sum(eligible), b$draws$aggregate_cell_counts[p], "Global eligible cell count")
  d[order(d$cell_id), , drop = FALSE]
}

f1_load_grid <- function(b, config) {
  saved <- b$meta$partition_files
  n <- as.integer(b$manifest$total_tiles)
  expected <- sprintf("tile_%05d.rds", seq_len(n))
  f1_assert(length(saved) == n && !anyDuplicated(basename(saved)) && setequal(basename(saved), expected), "Partition inventory differs from the run manifest.")
  # Prefer colocated files to allow moving a complete run to a different drive.
  paths <- file.path(b$directory, "grid_partitions", expected)
  f1_assert(all(file.exists(paths)), "Missing colocated grid partitions: ", paste(paths[!file.exists(paths)], collapse = ", "))
  pieces <- vector("list", n); checks <- vector("list", n)
  aggregate_sums <- NULL; aggregate_den <- NULL
  nd <- nrow(b$draws$paired_draw_ids)
  domain_state <- f2_init_domain_state(config$domain_spec, nd)
  area_draw_sums <- matrix(0, nrow = nd, ncol = 3L, dimnames = list(NULL, c("positive", "negative", "zero")))
  eligible_area <- 0
  for (i in seq_len(n)) {
    tab <- f1_read(paths[i]); f1_assert(is.data.frame(tab) && "period" %in% names(tab), "Invalid partition: ", paths[i])
    x <- tab[tab$period == b$period, , drop = FALSE]
    f1_assert(nrow(x) > 0L && all(c("cell_id", "area_weight", "eligible_for_area_day_summary") %in% names(x)), "Missing full-period tile rows: ", i)
    cp_path <- file.path(b$directory, "checkpoints", expected[i])
    if (isTRUE(config$verify_checkpoints)) {
      cp <- f1_read(cp_path)
      f1_assert(identical(cp$schema, "ERA5-warm02-tile-v1") && identical(cp$fingerprint, b$fingerprint) &&
        cp$tile_id == i && cp$ndraw == nd && identical(as.character(cp$periods), as.character(b$meta$periods)), "Checkpoint identity mismatch: ", cp_path)
      f1_assert(identical(cp$summary_table, tab), "Partition differs from its committed checkpoint: ", paths[i])
      f1_assert(identical(dim(cp$aggregate_sums), dim(b$draws$estimates)), "Checkpoint aggregate dimensions differ: ", i)
      f1_assert(all(is.finite(cp$aggregate_sums)) && all(is.finite(cp$aggregate_denominator)), "Nonfinite checkpoint aggregates: ", i)
      if (is.null(aggregate_sums)) {
        aggregate_sums <- cp$aggregate_sums; aggregate_den <- cp$aggregate_denominator
      } else {
        aggregate_sums <- aggregate_sums + cp$aggregate_sums
        aggregate_den <- aggregate_den + cp$aggregate_denominator
      }
      joint <- cp$joint_draw_means
      f1_assert(identical(as.integer(dim(joint)), c(as.integer(nrow(x)), as.integer(nd), 4L)) &&
        identical(dimnames(joint)[[3]], c("delta_M", "delta_N", "delta_U", "delta_TSV_P50")) &&
        identical(as.integer(cp$joint_draw_cell_ids), as.integer(x$cell_id)), "Checkpoint joint-grid posterior dimensions/order differ: ", i)
      f1_close(cp$joint_draw_denominator, x$n_common_support_days, "Checkpoint supported-day counts")
      take <- which(as.logical(x$eligible_for_area_day_summary))
      if (length(take)) {
        ju <- matrix(joint[take, , 3L, drop = FALSE], nrow = length(take), ncol = nd)
        jm <- matrix(joint[take, , 1L, drop = FALSE], nrow = length(take), ncol = nd)
        jn <- matrix(joint[take, , 2L, drop = FALSE], nrow = length(take), ncol = nd)
        f1_assert(all(is.finite(c(ju, jm, jn))), "Missing eligible-cell joint draws: ", i)
        f1_close(ju, jm + jn, "Cell-draw decomposition")
        for (metric in c("delta_M", "delta_N", "delta_U")) {
          mat <- switch(metric, delta_M = jm, delta_N = jn, delta_U = ju)
          q <- t(apply(mat, 1L, f1_quantiles))
          for (field in colnames(q)) f1_close(q[, field], x[[paste0(metric, "_", field)]][take], paste("Grid posterior summary", metric, field))
        }
        w <- x$area_weight[take]
        area_draw_sums[, "positive"] <- area_draw_sums[, "positive"] + colSums((ju > 0) * w)
        area_draw_sums[, "negative"] <- area_draw_sums[, "negative"] + colSums((ju < 0) * w)
        area_draw_sums[, "zero"] <- area_draw_sums[, "zero"] + colSums((ju == 0) * w)
        eligible_area <- eligible_area + sum(w)
      }
      domain_state <- f2_accumulate_domains(domain_state, x, joint, config$domain_spec)
      rm(cp, joint)
    }
    pieces[[i]] <- x
    checks[[i]] <- data.frame(tile_id = i, cells = nrow(x), checkpoint_verified = isTRUE(config$verify_checkpoints),
      partition = paths[i], checkpoint = if (isTRUE(config$verify_checkpoints)) cp_path else NA_character_)
    if (i %% 50L == 0L || i == n) message("Read/validated saved tiles: ", i, "/", n)
  }
  grid <- f1_validate_grid(do.call(rbind, pieces), b)
  area_draws <- NULL
  if (isTRUE(config$verify_checkpoints)) {
    f1_close(aggregate_den, b$draws$aggregate_denominator, "Checkpoint/global denominators")
    for (p in seq_along(b$meta$periods)) {
      a <- if (aggregate_den[p] > 0) aggregate_sums[p, , ] / aggregate_den[p] else
        matrix(NA_real_, nrow = 8L, ncol = nd)
      f1_close(a, b$draws$estimates[p, , ], paste("Checkpoint/global posterior", b$meta$periods[p]))
    }
    f1_assert(eligible_area > 0, "No eligible area for posterior fractions.")
    area_draws <- area_draw_sums / eligible_area
    f1_close(rowSums(area_draws), rep(1, nd), "Posterior area fractions sum to one")
  }
  list(data = grid, checks = do.call(rbind, checks), partition_files = paths, area_draws = area_draws,
    domain_state = domain_state)
}



# Range-screened scenarios are not a validated area of applicability (AOA).
# Their posterior intervals are conditional on the fixed selected locations.
f2_thresholds <- function(x, label, bounded = TRUE) {
  f1_assert(is.numeric(x) && all(is.finite(x)) && !anyDuplicated(x), "Invalid thresholds: ", label)
  if (bounded) f1_assert(all(x > 0 & x <= 1), label, " must lie in (0,1].") else
    f1_assert(all(x > 0 & x == floor(x)), label, " must be positive whole metres.")
  sort(x)
}
f2_nonempty <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
f2_read_spatial_input <- function(path, b, kind, basis) {
  if (is.null(path)) return(NULL)
  f1_assert(f2_nonempty(path) && file.exists(path), "Cannot find ", kind, " CSV: ", path)
  f1_assert(f2_nonempty(basis), "Supply a nonempty provenance/selection basis for ", kind, ".")
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  req <- if (kind == "external_domain") c("cell_id", "include", "reason") else c("cell_id", "elevation_m")
  f1_assert(all(req %in% names(d)), kind, " CSV requires: ", paste(req, collapse = ", "))
  f1_assert(is.numeric(d$cell_id) && all(is.finite(d$cell_id) & d$cell_id == floor(d$cell_id)) &&
    !anyDuplicated(d$cell_id) && setequal(d$cell_id, b$land$cell_id),
    kind, " must contain every selected cell_id exactly once; no absent, extra or duplicate IDs.")
  d <- d[match(b$land$cell_id, d$cell_id), , drop = FALSE]
  if ("lon" %in% names(d)) f1_close(d$lon, b$land$lon_plot, paste(kind, "longitude"))
  if ("lat" %in% names(d)) f1_close(d$lat, b$land$lat, paste(kind, "latitude"))
  if ("run_fingerprint" %in% names(d)) f1_assert(!anyNA(d$run_fingerprint) && all(d$run_fingerprint == b$fingerprint), "External CSV run fingerprint differs: ", kind)
  if (kind == "external_domain") {
    raw <- toupper(trimws(as.character(d$include)))
    f1_assert(!anyNA(raw) && all(raw %in% c("TRUE", "FALSE", "0", "1")), "include must be explicit TRUE/FALSE or 1/0.")
    d$include <- raw %in% c("TRUE", "1")
    f1_assert(!anyNA(d$reason) && all(nzchar(trimws(as.character(d$reason)))), "Every external-mask row requires a reason.")
    d$reason <- as.character(d$reason)
  } else f1_assert(is.numeric(d$elevation_m) && all(is.finite(d$elevation_m)), "elevation_m must be finite metres for every selected cell; missing elevations cannot be treated as low altitude.")
  attr(d, "source_path") <- normalizePath(path, winslash = "/", mustWork = TRUE)
  attr(d, "basis") <- basis
  d
}

f2_prepare_domains <- function(b, config) {
  if (!is.null(config$population_table)) return(p1_prepare_domains(b, config))
  roots <- f2_thresholds(config$root_range_thresholds, "root_range_thresholds")
  histories <- f2_thresholds(config$history_fraction_thresholds, "history_fraction_thresholds")
  elevations <- f2_thresholds(config$elevation_thresholds_m, "elevation_thresholds_m", FALSE)
  ext <- f2_read_spatial_input(config$external_domain_csv, b, "external_domain", config$external_domain_basis)
  elev <- f2_read_spatial_input(config$elevation_csv, b, "elevation", config$elevation_basis)
  if (length(elevations)) f1_assert(!is.null(elev), "Elevation scenarios require elevation_csv and elevation_basis. No elevation data are present in warm02 summaries.")
  defs <- data.frame(domain_id = "original", domain_label = "Original history-supported domain", kind = "original", threshold = NA_real_,
    basis = "Original saved eligibility and supported dates; not a claim of fully validated global transferability", stringsAsFactors = FALSE)
  add <- function(id, label, kind, cut, basis) {
    defs <<- rbind(defs, data.frame(domain_id = id, domain_label = label, kind = kind, threshold = cut, basis = basis))
  }
  for (z in histories) add(sprintf("history_%03d", round(z * 100)), sprintf("History-supported dates >= %g%%", 100 * z), "history", z,
    "Post hoc location-coverage sensitivity; preserves all original supported dates at retained locations")
  for (z in roots) add(sprintf("root_%03d", round(z * 100)), sprintf("Median root-range fraction >= %g%%", 100 * z), "root", z,
    "Post hoc model-derived SET-range diagnostic; NOT a validated SET-by-history AOA; fixed location selection")
  if (!is.null(ext)) add("external", "User-defined external domain", "external", NA_real_, config$external_domain_basis)
  for (z in elevations) add(paste0("elev_", z), sprintf("Elevation <= %g m (sensitivity)", z), "elevation", z, config$elevation_basis)
  f1_assert(!anyDuplicated(defs$domain_id), "Thresholds create duplicate domain IDs; use distinct percentage points.")
  f1_assert(config$main_domain %in% defs$domain_id, "Unknown main_domain: ", config$main_domain,
    ". Available: ", paste(defs$domain_id, collapse = ", "))
  list(definitions = defs, external = ext, elevation = elev)
}

f2_domain_membership <- function(d, spec) {
  if (!is.null(spec$population)) return(p1_domain_membership(d, spec))
  base <- as.logical(d$eligible_for_area_day_summary)
  f1_assert(!anyNA(base), "Missing original eligibility.")
  ans <- data.frame(cell_id = d$cell_id)
  for (j in seq_len(nrow(spec$definitions))) {
    def <- spec$definitions[j, ]; cut <- def$threshold
    take <- switch(def$kind,
      original = base,
      history = base & is.finite(d$common_support_fraction) & d$common_support_fraction >= cut,
      root = base & is.finite(d$U_pair_in_SET_range_median) & d$U_pair_in_SET_range_median >= cut,
      external = base & spec$external$include[match(d$cell_id, spec$external$cell_id)],
      elevation = base & spec$elevation$elevation_m[match(d$cell_id, spec$elevation$cell_id)] <= cut,
      f1_stop("Unknown domain kind: ", def$kind))
    f1_assert(length(take) == nrow(d) && !anyNA(take), "Invalid membership: ", def$domain_id)
    ans[[def$domain_id]] <- take
  }
  ans
}

f2_init_domain_state <- function(spec, nd) {
  ids <- spec$definitions$domain_id; metrics <- c("delta_M", "delta_N", "delta_U", "delta_TSV_P50")
  list(sums = array(0, c(length(ids), 4L, nd), dimnames = list(domain_id = ids, metric = metrics, paired_draw_position = as.character(seq_len(nd)))),
    area_sign_sums = array(0, c(length(ids), 3L, nd), dimnames = list(domain_id = ids, direction = c("positive", "negative", "zero"), paired_draw_position = as.character(seq_len(nd)))),
    denominator = setNames(numeric(length(ids)), ids), area = setNames(numeric(length(ids)), ids),
    cells = setNames(integer(length(ids)), ids), expected_area_days = setNames(numeric(length(ids)), ids))
}

f2_accumulate_domains <- function(state, d, joint, spec) {
  membership <- f2_domain_membership(d, spec)
  nd <- dim(joint)[2L]
  for (id in spec$definitions$domain_id) {
    ii <- which(membership[[id]])
    if (!length(ii)) next
    w <- d$area_weight[ii]; wd <- w * d$n_common_support_days[ii]
    for (m in seq_len(4L)) {
      mat <- matrix(joint[ii, , m, drop = FALSE], nrow = length(ii), ncol = nd)
      f1_assert(all(is.finite(mat)), "Nonfinite selected posterior draws: ", id)
      state$sums[id, m, ] <- state$sums[id, m, ] + colSums(mat * wd)
    }
    u <- matrix(joint[ii, , 3L, drop = FALSE], nrow = length(ii), ncol = nd)
    state$area_sign_sums[id, "positive", ] <- state$area_sign_sums[id, "positive", ] + colSums((u > 0) * w)
    state$area_sign_sums[id, "negative", ] <- state$area_sign_sums[id, "negative", ] + colSums((u < 0) * w)
    state$area_sign_sums[id, "zero", ] <- state$area_sign_sums[id, "zero", ] + colSums((u == 0) * w)
    state$denominator[id] <- state$denominator[id] + sum(wd)
    state$area[id] <- state$area[id] + sum(w)
    state$cells[id] <- state$cells[id] + length(ii)
    state$expected_area_days[id] <- state$expected_area_days[id] + sum(w * d$n_expected_days[ii])
  }
  state
}

f2_finalize_domains <- function(b, loaded, spec) {
  state <- loaded$domain_state; ids <- spec$definitions$domain_id
  arr <- state$sums; fractions <- state$area_sign_sums
  for (id in ids) {
    arr[id, , ] <- if (state$denominator[id] > 0) state$sums[id, , ] / state$denominator[id] else NA_real_
    fractions[id, , ] <- if (state$area[id] > 0) state$area_sign_sums[id, , ] / state$area[id] else NA_real_
  }
  f1_close(arr[, "delta_U", ], arr[, "delta_M", ] + arr[, "delta_N", ], "Domain-wise decomposition")
  p <- match(b$period, b$meta$periods)
  f1_close(arr["original", , ], b$draws$estimates[p, seq_len(4L), ], "Original domain unchanged")
  f1_close(t(fractions["original", , ]), loaded$area_draws, "Original posterior area fractions unchanged")
  membership <- f2_domain_membership(loaded$data, spec)
  coverage <- spec$definitions
  coverage$n_cells <- unname(state$cells[ids])
  coverage$status <- ifelse(coverage$n_cells > 0, "AVAILABLE", "EMPTY")
  coverage$area_weight_sum <- unname(state$area[ids])
  coverage$weighted_supported_cell_day_denominator <- unname(state$denominator[ids])
  coverage$retained_area_fraction <- unname(state$area[ids] / state$area["original"])
  coverage$retained_area_day_fraction <- unname(state$denominator[ids] / state$denominator["original"])
  coverage$supported_fraction_within_selected_locations <- ifelse(state$expected_area_days[ids] > 0,
    unname(state$denominator[ids] / state$expected_area_days[ids]), NA_real_)
  coverage$period <- b$period
  for (id in ids) {
    ii <- membership[[id]]
    f1_close(sum(ii), state$cells[id], paste("Domain selected cells", id))
    f1_close(sum(loaded$data$area_weight[ii] * loaded$data$n_common_support_days[ii]), state$denominator[id], paste("Domain denominator", id))
    f1_close(sum(loaded$data$area_weight[ii]), state$area[id], paste("Domain selected area", id))
    if (state$area[id] > 0) f1_close(colSums(fractions[id, , ]), rep(1, dim(arr)[3L]), paste("Domain sign fractions sum", id))
  }
  rows <- list(); comparisons <- list(); signs <- list(); k <- 0L
  for (id in ids) for (m in dimnames(arr)[[2L]]) {
    k <- k + 1L; q <- f1_quantiles(arr[id, m, ]); cov <- coverage[coverage$domain_id == id, ]
    rows[[k]] <- data.frame(domain_id = id, domain_label = cov$domain_label, kind = cov$kind, period = b$period,
      metric = m, lower = q["lower"], median = q["median"], upper = q["upper"],
      prob_positive = q["prob_positive"], prob_negative = q["prob_negative"], n_cells = cov$n_cells,
      retained_area_fraction = cov$retained_area_fraction, retained_area_day_fraction = cov$retained_area_day_fraction,
      weighted_supported_cell_day_denominator = cov$weighted_supported_cell_day_denominator,
      status = cov$status, unit = if (m == "delta_TSV_P50") "TSV" else "degC_SET", row.names = NULL)
    dq <- f1_quantiles(arr[id, m, ] - arr["original", m, ])
    comparisons[[k]] <- data.frame(domain_id = id, metric = m, lower = dq["lower"], median = dq["median"], upper = dq["upper"],
      interpretation = "Change in aggregate estimand after selecting locations; same posterior draw pairing", row.names = NULL)
  }
  for (id in ids) for (direction in dimnames(fractions)[[2L]]) {
    q <- f1_quantiles(fractions[id, direction, ])
    signs[[length(signs) + 1L]] <- data.frame(domain_id = id, direction = direction,
      posterior_mean = if (state$area[id] > 0) mean(fractions[id, direction, ]) else NA_real_,
      lower = q["lower"], median = q["median"], upper = q["upper"],
      interpretation = "Area fraction for sign of location supported-date mean, not daily sign frequency or population fraction", row.names = NULL)
  }
  list(draws = arr, sign_fraction_draws = fractions, coverage = coverage,
    membership = membership, summary = do.call(rbind, rows), contrasts = do.call(rbind, comparisons), sign_summary = do.call(rbind, signs))
}


# Figure 1 v2 plotting module. This file performs no I/O and no model prediction.
# API:
#   grid: validated one-row-per-cell table, with effective_eligible for the main
#         selected domain; falls back to eligible_for_area_day_summary.
#   summary: one primary-period summary per metric in the SAME main domain.
#   domains: optional wide data.frame(cell_id, original, root_050, ...), a
#            named list of logical vectors, or long cell_id/domain_id/eligible.
#   sensitivity_summary: optional data.frame(domain_id[,domain_label], metric,
#            median, lower, upper[,retained_area_fraction]); all values must be
#            computed using the original paired posterior draws, not cell medians.
#            retained_area_fraction is relative to original eligible land area.
#   config: optional main_domain, main_domain_label, base_font_size, font_family,
#           map_scale_limit, period_label, annual_summary and domain_labels.
# Returns ggplot/patchwork objects, the plotted tables, limits, and full captions.
f2_make_plots <- function(grid, summary, meta, config = list(), domains = NULL,
                          sensitivity_summary = NULL) {
  need <- c("ggplot2", "maps", "patchwork", "scales")
  missing <- need[!vapply(need, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) stop("Missing plotting packages: ", paste(missing, collapse = ", "), call. = FALSE)
  cfg <- function(name, fallback) if (is.null(config[[name]])) fallback else config[[name]]
  base_size <- cfg("base_font_size", 7)
  font_family <- cfg("font_family", "sans")
  if (!is.numeric(base_size) || length(base_size) != 1L || !is.finite(base_size) || base_size < 6) stop("base_font_size must be finite and >= 6.")
  grid <- as.data.frame(grid)
  metrics <- c("delta_U", "delta_M", "delta_N")
  required <- c("cell_id", "lon", "lat", "eligible_for_area_day_summary", "common_support_fraction",
                "U_pair_in_SET_range_median", "delta_U_lower", "delta_U_upper", paste0(metrics, "_median"))
  if (!all(required %in% names(grid))) stop("Plot data lack: ", paste(setdiff(required, names(grid)), collapse = ", "))
  eligibility_field <- if ("effective_eligible" %in% names(grid)) "effective_eligible" else "eligible_for_area_day_summary"
  if (!is.logical(grid[[eligibility_field]]) || anyNA(grid[[eligibility_field]])) stop("Main-domain eligibility must be complete logical values.")
  eligible <- grid[[eligibility_field]]
  original_eligible <- !is.na(grid$eligible_for_area_day_summary) & grid$eligible_for_area_day_summary
  if (any(eligible & !original_eligible)) stop("The main domain cannot add previously ineligible cells.")
  if (!any(eligible)) stop("The selected main domain contains no eligible cells.")
  main_domain <- cfg("main_domain", "original")
  main_domain_label <- cfg("main_domain_label", if (main_domain == "original") "Original eligible domain" else main_domain)
  res <- meta$grid_resolution
  if (is.null(names(res)) || !all(c("lon", "lat") %in% names(res)) || any(!is.finite(res[c("lon", "lat")]) | res[c("lon", "lat")] <= 0)) stop("Invalid named grid resolution.")
  lat_range <- meta$global_config$lat_range
  if (is.null(lat_range)) lat_range <- range(grid$lat, finite = TRUE) + c(-1, 1) * res[["lat"]] / 2
  if (length(lat_range) != 2L || any(!is.finite(lat_range)) || lat_range[1] >= lat_range[2]) stop("Invalid latitude range.")
  if (any(!is.finite(grid$lon) | !is.finite(grid$lat)) || any(grid$lon < -180 | grid$lon > 180)) stop("Invalid map coordinates.")
  validate_summary <- function(s, label) {
    if (!all(c("metric", "median", "lower", "upper") %in% names(s))) stop(label, " requires metric, median, lower and upper.")
    s <- s[s$metric %in% metrics, , drop = FALSE]
    if (nrow(s) != 3L || anyDuplicated(s$metric)) stop(label, " requires exactly one row for each delta_U, delta_M and delta_N.")
    s <- s[match(metrics, s$metric), , drop = FALSE]
    if (any(!is.finite(as.matrix(s[c("lower", "median", "upper")])))) stop(label, " intervals must be finite.")
    if (any(s$lower > s$median | s$median > s$upper)) stop(label, " contains invalid quantile ordering.")
    s
  }
  forest <- validate_summary(summary, "Main-domain summary")
  vals <- unlist(lapply(metrics, function(m) grid[[paste0(m, "_median")]][eligible]), use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) stop("There are no finite, eligible map estimates.")
  observed_limit <- max(abs(vals))
  requested_limit <- cfg("map_scale_limit", NULL)
  if (!is.null(requested_limit) && (length(requested_limit) != 1L || !is.numeric(requested_limit) || !is.finite(requested_limit) || requested_limit <= 0)) stop("map_scale_limit must be NULL or one positive finite number.")
  if (!is.null(requested_limit) && requested_limit < observed_limit) warning("map_scale_limit was expanded to include every displayed value; no value was clipped.", call. = FALSE)
  lim <- max(observed_limit, if (is.null(requested_limit)) 0 else requested_limit)
  if (lim == 0) lim <- 0.01
  period <- cfg("period_label", NULL)
  if (is.null(period)) {
    yrs <- meta$global_config$years
    if (is.null(yrs)) yrs <- meta$years
    period <- if (length(yrs)) paste(range(yrs), collapse = "-") else "Saved analysis period"
  }
  blue <- "#2166AC"; orange <- "#D55E00"; grey_na <- "#D9D9D9"
  world <- ggplot2::map_data("world")
  theme_map <- ggplot2::theme_minimal(base_size = base_size, base_family = font_family) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + .5, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - .5, hjust = 0),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"),
      plot.margin = ggplot2::margin(2, 2, 2, 2), legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size - .5),
      legend.text = ggplot2::element_text(size = base_size - .5),
      legend.key.height = grid::unit(2.4, "mm"), legend.key.width = grid::unit(4, "mm"),
      legend.spacing.x = grid::unit(1.5, "mm"), legend.box.spacing = grid::unit(1, "mm"))
  map_base <- function(d) ggplot2::ggplot() +
    ggplot2::geom_polygon(data = world, ggplot2::aes(x = long, y = lat, group = group),
      inherit.aes = FALSE, fill = grey_na, colour = NA, show.legend = FALSE) +
    ggplot2::geom_tile(data = d, ggplot2::aes(x = lon, y = lat, fill = value),
      width = res[["lon"]], height = res[["lat"]], inherit.aes = FALSE) +
    ggplot2::geom_polygon(data = world, ggplot2::aes(x = long, y = lat, group = group),
      inherit.aes = FALSE, fill = NA, colour = "#777777", linewidth = .09, show.legend = FALSE) +
    ggplot2::coord_quickmap(xlim = c(-180, 180), ylim = lat_range, expand = FALSE, clip = "on") + theme_map
  # The invisible dummy point creates a genuine, separate categorical legend.
  # Its alpha is overridden in the legend only; it adds no visible map mark.
  missing_legend <- function(p) p +
    ggplot2::geom_point(data = data.frame(lon = 0, lat = 0, missingness = "No estimate / outside selected domain"),
      ggplot2::aes(x = lon, y = lat, shape = missingness), inherit.aes = FALSE,
      colour = grey_na, size = 2, alpha = 0, show.legend = TRUE) +
    ggplot2::scale_shape_manual(values = 15, name = NULL,
      guide = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1, size = 2.3)))
  map_data <- list()
  map_breaks <- unique(c(-lim, pretty(c(-lim, lim), n = 4), 0, lim))
  map_breaks <- map_breaks[map_breaks >= -lim & map_breaks <= lim]
  if (length(map_breaks) > 7L) map_breaks <- seq(-lim, lim, length.out = 5)
  main_map <- function(metric, title) {
    d <- grid[c("cell_id", "lon", "lat")]
    d$value <- grid[[paste0(metric, "_median")]]
    d$value[!eligible | !is.finite(d$value)] <- NA_real_
    d$metric <- metric
    map_data[[metric]] <<- d
    missing_legend(map_base(d) + ggplot2::scale_fill_gradient2(low = blue, mid = "white", high = orange,
      midpoint = 0, limits = c(-lim, lim), na.value = grey_na, breaks = map_breaks,
      labels = scales::label_number(accuracy = if (lim >= 1) .1 else NULL), name = "Posterior median shift (\u00b0C SET)",
      guide = ggplot2::guide_colourbar(order = 1, title.position = "top", title.hjust = .5,
        barwidth = grid::unit(62, "mm"), barheight = grid::unit(2.4, "mm"))) + ggplot2::labs(title = title))
  }
  p_u <- main_map("delta_U", expression("Warm boundary, " * Delta * U))
  p_m <- main_map("delta_M", expression("Neutral position, " * Delta * M))
  p_n <- main_map("delta_N", expression("Positive half-width, " * Delta * N))
  numeric_label <- function(x) sprintf("%.3f", ifelse(abs(x) < .0005, 0, x))
  forest$y <- c(3, 2, 1)
  forest$label <- paste0(numeric_label(forest$median), " [", numeric_label(forest$lower), ", ", numeric_label(forest$upper), "]")
  xr <- range(c(forest$lower, forest$upper, 0)); span <- diff(xr)
  if (span == 0) span <- 1
  xr <- xr + c(-.07, .07) * span
  breaks <- if (diff(xr) <= 5 && diff(xr) >= 1) seq(floor(xr[1] * 2) / 2, ceiling(xr[2] * 2) / 2, by = .5) else pretty(xr, n = 5)
  breaks <- breaks[breaks >= xr[1] & breaks <= xr[2]]
  forest$annotation_x <- xr[2] - .015 * span
  theme_forest <- ggplot2::theme_classic(base_size = base_size, base_family = font_family) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = base_size + .5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size - .5),
      axis.title = ggplot2::element_text(size = base_size), axis.text = ggplot2::element_text(size = base_size - .5, colour = "#333333"),
      axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(linewidth = .25),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"), plot.margin = ggplot2::margin(2, 3, 2, 5))
  p_forest <- ggplot2::ggplot(forest) +
    ggplot2::geom_vline(xintercept = 0, linewidth = .25, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, y = y, yend = y), linewidth = .6, colour = "#333333") +
    ggplot2::geom_point(ggplot2::aes(x = median, y = y), size = 1.9, shape = 21, stroke = .3, colour = "#222222", fill = "white") +
    ggplot2::geom_text(ggplot2::aes(x = annotation_x, y = y + .24, label = label), hjust = 1,
      size = (base_size - .5) / 2.845, family = font_family, colour = "#333333") +
    ggplot2::scale_y_continuous(breaks = c(1, 2, 3), labels = expression(Delta*N, Delta*M, Delta*U), limits = c(.5, 3.65)) +
    ggplot2::scale_x_continuous(limits = xr, breaks = breaks, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::labs(title = "Area-day weighted response", subtitle = "Median [95% credible interval]",
      x = "Mean shift (\u00b0C SET)", y = NULL) + theme_forest + ggplot2::theme(legend.position = "none")
  main <- patchwork::wrap_plots(list(p_u, p_m, p_n, p_forest), ncol = 2, guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a", theme = ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 1, 1)))
  main <- main & ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
  if (main_domain != "original") main <- main + patchwork::plot_annotation(
    subtitle = paste0("Analysis domain: ", main_domain_label),
    theme = ggplot2::theme(plot.subtitle = ggplot2::element_text(size = base_size, family = font_family)))
  main_caption <- paste0("Thermal-history effects on the warm-sensation boundary, ", period, ". Main analysis domain: ", main_domain_label, ". ",
    "M1 minus M0, standardized to study-sample covariates. U = M + N is the SET root at mean TSV = +1; N = 1/A > 0, and delta N has no additional minus sign. ",
    "Maps show posterior medians of each grid point's mean over its original supported dates. Tiles represent regularly sampled points, without spatial interpolation. ",
    "Grey land denotes absent estimates or cells excluded from the selected domain; white denotes approximately zero shift or ocean. All three maps share one linear symmetric colour scale spanning every displayed value. ",
    "Panel d uses the same grid domain as the maps and an area-times-supported-day weighted aggregation of paired posterior draws. It is not an average of posterior medians. ",
    "Intervals are 95% equal-tailed credible intervals. Summaries concern supported dates and are not necessarily calendar-year averages or population-weighted exposures. ",
    "Marginal medians need not add even though delta U = delta M + delta N holds draw by draw. Domain exclusions are sensitivity definitions, not proof of transportability within the remaining domain.")
  support_data <- list()
  fraction_map <- function(field, title, legend_title, mask = FALSE) {
    d <- grid[c("cell_id", "lon", "lat")]; d$value <- grid[[field]]
    if (any(is.finite(d$value) & (d$value < -1e-8 | d$value > 1 + 1e-8))) stop("Invalid fraction in ", field)
    d$value[!is.finite(d$value)] <- NA_real_
    if (mask) d$value[!original_eligible] <- NA_real_
    d$field <- field; support_data[[field]] <<- d
    map_base(d) + ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1),
      breaks = c(0, .25, .5, .75, 1), labels = scales::label_percent(accuracy = 1), na.value = grey_na, name = legend_title,
      guide = ggplot2::guide_colourbar(title.position = "top", title.hjust = .5,
        barwidth = grid::unit(55, "mm"), barheight = grid::unit(2.4, "mm"))) + ggplot2::labs(title = title)
  }
  s_support <- fraction_map("common_support_fraction", "Thermal-history support", "Supported / expected dates")
  s_root <- fraction_map("U_pair_in_SET_range_median", "Warm-root range diagnostic", "Median fraction of supported dates", mask = TRUE)
  d <- grid[c("cell_id", "lon", "lat")]
  levels_dir <- c("Below zero", "Includes zero", "Above zero")
  valid <- original_eligible & is.finite(grid$delta_U_lower) & is.finite(grid$delta_U_upper)
  cl <- ifelse(grid$delta_U_upper < 0, levels_dir[1], ifelse(grid$delta_U_lower > 0, levels_dir[3], levels_dir[2])); cl[!valid] <- NA_character_
  d$value <- factor(cl, levels = levels_dir); support_data$delta_U_interval_direction <- d
  s_dir <- map_base(d) + ggplot2::scale_fill_manual(values = stats::setNames(c(blue, "#F6F6F6", orange), levels_dir),
    breaks = levels_dir, limits = levels_dir, drop = FALSE, na.value = grey_na, name = "Pointwise 95% credible interval", na.translate = FALSE,
    guide = ggplot2::guide_legend(title.position = "top", title.hjust = .5, nrow = 1)) +
    ggplot2::labs(title = expression("Warm-boundary direction, " * Delta * U))
  support <- patchwork::wrap_plots(list(s_support, s_root, s_dir), ncol = 1, guides = "keep") +
    patchwork::plot_annotation(tag_levels = "a")
  support_caption <- paste0("Diagnostics on the original saved grid, before optional domain exclusions. Thermal-history support is the fraction of expected dates satisfying the saved marginal and joint history checks. Root-range and direction maps use the original eligible domain. ",
    "The root-range diagnostic is the posterior median of the fraction of supported dates for which both warm roots lie within the saved marginal training SET P1-P99 bounds. ",
    "It neither validates conditional SET-by-history support nor establishes suitability for high-altitude populations. Direction categories use pointwise posterior intervals, with no simultaneous field-significance claim. Grey denotes unavailable estimates.")
  sensitivity <- NULL; domain_map <- NULL; sensitivity_data <- NULL; domain_data <- NULL
  if (!is.null(sensitivity_summary) && nrow(sensitivity_summary)) {
    ss <- as.data.frame(sensitivity_summary)
    if (!"domain_id" %in% names(ss)) stop("sensitivity_summary requires domain_id.")
    if (!"domain_label" %in% names(ss)) ss$domain_label <- ss$domain_id
    ss <- ss[ss$metric %in% metrics, , drop = FALSE]
    ids <- unique(ss$domain_id)
    ss <- do.call(rbind, lapply(ids, function(id) {
      one <- ss[ss$domain_id == id, , drop = FALSE]
      if (all(is.na(as.matrix(one[c("median", "lower", "upper")])))) {
        if (nrow(one) != 3L || anyDuplicated(one$metric)) stop("Empty domain requires one row per metric: ", id)
        one$plot_empty <- TRUE
        one[match(metrics, one$metric), , drop = FALSE]
      } else {
        one <- validate_summary(one, paste("Domain", id)); one$plot_empty <- FALSE; one
      }
    }))
    ss$domain_display <- ss$domain_label
    if ("retained_area_fraction" %in% names(ss)) {
      ok <- is.finite(ss$retained_area_fraction)
      ss$domain_display[ok] <- paste0(ss$domain_label[ok], " (", sprintf("%.1f", 100 * ss$retained_area_fraction[ok]), "% area)")
    }
    ss$domain_display <- factor(ss$domain_display, levels = rev(unique(ss$domain_display)))
    ss$component <- factor(ss$metric, levels = metrics, labels = c("Warm boundary: delta U", "Neutral position: delta M", "Positive half-width: delta N"))
    sensitivity <- ggplot2::ggplot(ss) +
      ggplot2::geom_vline(xintercept = 0, linewidth = .25, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, y = domain_display, yend = domain_display), linewidth = .5, colour = "#444444", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(x = median, y = domain_display), size = 1.8, shape = 21, fill = "white", stroke = .3, na.rm = TRUE) +
      ggplot2::geom_text(data = ss[ss$plot_empty, , drop = FALSE], ggplot2::aes(x = 0, y = domain_display),
        label = "No eligible cells", inherit.aes = FALSE, size = (base_size - .5) / 2.845,
        hjust = -.1, colour = "#777777", family = font_family) +
      ggplot2::facet_wrap(~component, ncol = 1, scales = "fixed") +
      ggplot2::labs(x = "Area-day weighted mean shift (\u00b0C SET)", y = NULL) + theme_forest +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", hjust = 0),
        panel.spacing.y = grid::unit(4, "mm"))
    sensitivity_data <- ss
  }
  if (!is.null(domains)) {
    if (is.list(domains) && !is.data.frame(domains)) {
      if (is.null(names(domains)) || any(!nzchar(names(domains)))) stop("domains must be a named list.")
      dd <- lapply(names(domains), function(id) {
        e <- domains[[id]]
        if (!is.logical(e) || length(e) != nrow(grid) || anyNA(e)) stop("Invalid logical domain mask: ", id)
        data.frame(cell_id = grid$cell_id, domain_id = id, eligible = e, stringsAsFactors = FALSE)
      }); dd <- do.call(rbind, dd)
    } else {
      dd <- as.data.frame(domains)
      if ("cell_id" %in% names(dd) && !"domain_id" %in% names(dd)) {
        mask_names <- setdiff(names(dd), "cell_id")
        dd <- do.call(rbind, lapply(mask_names, function(id) {
          if (!is.logical(domains[[id]]) || anyNA(domains[[id]])) stop("Invalid wide-domain logical mask: ", id)
          data.frame(cell_id = domains$cell_id, domain_id = id, eligible = domains[[id]], stringsAsFactors = FALSE)
        }))
      }
      if (!all(c("cell_id", "domain_id", "eligible") %in% names(dd))) stop("Long domains table requires cell_id, domain_id and eligible.")
      if (!is.logical(dd$eligible) || anyNA(dd$eligible) || anyDuplicated(paste(dd$domain_id, dd$cell_id))) stop("Invalid or duplicated long-domain mask rows.")
    }
    dd <- dd[dd$domain_id == main_domain, , drop = FALSE]
    if (!nrow(dd)) stop("domains does not contain the selected main_domain: ", main_domain)
    ix <- match(dd$cell_id, grid$cell_id)
    if (anyNA(ix)) stop("Domain masks reference unknown cell IDs.")
    if (any(dd$eligible & !original_eligible[ix])) stop("Domain masks cannot add originally ineligible cells.")
    dd$lon <- grid$lon[ix]; dd$lat <- grid$lat[ix]
    dd$value <- ifelse(!original_eligible[ix], "Originally unavailable", ifelse(dd$eligible, "Retained", "Excluded by domain rule"))
    if ("population_density_lower_bound" %in% names(grid) && grepl("^pop_", main_domain)) {
      population_cut <- if (main_domain == "pop_ge5") 5 else if (main_domain == "pop_gt0") 0 else 1
      population_pass <- p1_include(grid, population_cut, strict = main_domain == "pop_gt0")
      population_unknown <- !grid$complete & !population_pass
      dd$value[original_eligible[ix] & !dd$eligible & population_unknown[ix]] <- "Population data insufficient"
    }
    dd$value <- factor(dd$value, levels = c("Retained", "Excluded by domain rule", "Originally unavailable", "Population data insufficient"))
    if (!"domain_label" %in% names(dd)) dd$domain_label <- dd$domain_id
    labels <- cfg("domain_labels", NULL)
    if (!is.null(labels)) {
      matched <- match(dd$domain_id, names(labels)); use <- !is.na(matched)
      dd$domain_label[use] <- unname(labels[matched[use]])
    }
    dd$domain_label <- factor(dd$domain_label, levels = unique(dd$domain_label))
    domain_map <- map_base(dd) + ggplot2::scale_fill_manual(
      values = c("Retained" = "#3D7E9B", "Excluded by domain rule" = "#E9A55A", "Originally unavailable" = grey_na, "Population data insufficient" = "#A985B0"),
      drop = FALSE, na.value = grey_na, name = "Grid membership", guide = ggplot2::guide_legend(nrow = 1)) +
      ggplot2::labs(title = main_domain_label)
    domain_data <- dd
  }
  annual <- NULL; annual_data <- cfg("annual_summary", NULL)
  if (!is.null(annual_data) && nrow(annual_data)) {
    if (!all(c("period", "metric", "median", "lower", "upper") %in% names(annual_data))) stop("annual_summary is incomplete.")
    annual_data <- annual_data[annual_data$metric %in% metrics & grepl("^[0-9]{4}$", annual_data$period), , drop = FALSE]
    if (nrow(annual_data)) {
      annual_data$year <- as.integer(annual_data$period)
      annual_data$component <- factor(annual_data$metric, levels = metrics, labels = c("Warm boundary: delta U", "Neutral position: delta M", "Positive half-width: delta N"))
      annual <- ggplot2::ggplot(annual_data, ggplot2::aes(x = year, y = median)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = .25, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = blue, alpha = .18) +
        ggplot2::geom_line(colour = blue, linewidth = .5) + ggplot2::geom_point(colour = blue, size = 1.2) +
        ggplot2::facet_wrap(~component, ncol = 1, scales = "fixed") +
        ggplot2::scale_x_continuous(breaks = pretty(range(annual_data$year), n = 5)) +
        ggplot2::labs(title = "Original supported domain", x = "Year", y = "Mean shift (\u00b0C SET)") + theme_forest +
        ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", hjust = 0))
    }
  }
  list(main = main, support = support, sensitivity = sensitivity, domain_map = domain_map, annual = annual,
    panels = list(delta_U = p_u, delta_M = p_m, delta_N = p_n, posterior_summary = p_forest,
      support_fraction = s_support, warm_root_range = s_root, direction = s_dir),
    plot_data = list(maps = do.call(rbind, map_data), forest = forest, support = support_data,
      sensitivity = sensitivity_data, domains = domain_data, annual = annual_data), scale_limit = lim,
    captions = list(main = main_caption, support = support_caption,
      sensitivity = "Sensitivity to fixed grid-domain restrictions. Points and lines show area-times-original-supported-day weighted posterior medians and 95% credible intervals, recomputed from common saved draws for each domain. The area percentages refer to retained area relative to the original eligible grid domain, not population or supported-day fractions. Restriction changes the target domain and does not establish joint SET-by-history support. All original supported dates within retained cells remain in the calculation.",
      domain_map = "Fixed grid-domain membership used for the sensitivity analysis. Grey indicates originally unavailable estimates; orange indicates originally eligible cells excluded by a stated domain rule. Purple, where present, indicates insufficient population evidence for the selected threshold; it does not indicate zero population. These maps document exclusions and do not label the retained domain as externally validated.",
      annual = "Annual area-day weighted posterior medians and pointwise 95% credible intervals on the original eligible domain. No temporal trend model was fitted; this annual panel is not a recomputation under optional spatial exclusions."))
}



# Reporting uses supported-date estimands and the original saved draw pairing.
f1_metric_dictionary <- function() data.frame(
  metric = c("delta_M", "delta_N", "delta_U", "delta_TSV_P50", "U_M1", "U_M0", "U_pair_in_SET_range", "absolute_delta_U"),
  label = c("Neutral SET shift", "Positive half-width change", "Warm-boundary SET shift", "TSV change at sample P50 SET",
    "Warm-boundary SET with history", "Warm-boundary SET without history", "Fraction of supported dates with both roots in training SET bounds", "Mean absolute daily warm-boundary shift"),
  unit = c("degC_SET", "degC_SET", "degC_SET", "TSV", "degC_SET", "degC_SET", "fraction", "degC_SET"),
  stringsAsFactors = FALSE)

f1_write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
f1_write_text <- function(x, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(x, con, useBytes = FALSE)
}
f1_file_inventory <- function(paths, hash = TRUE) {
  z <- file.info(paths)
  data.frame(path = normalizePath(paths, winslash = "/", mustWork = TRUE), bytes = z$size,
    modified = format(z$mtime, "%Y-%m-%d %H:%M:%S %z"),
    md5 = if (hash) unname(tools::md5sum(paths)) else NA_character_, row.names = NULL)
}

f1_export_tables <- function(b, loaded, out) {
  d <- loaded$data; e <- d$eligible_for_area_day_summary
  dictionary <- f1_metric_dictionary()
  summary <- b$summary
  j <- match(summary$metric, dictionary$metric)
  summary$label <- dictionary$label[j]; summary$unit <- dictionary$unit[j]
  summary$estimand <- "Area-weighted supported cell-day mean; saved M1/M0 draw pairing"
  main <- summary[summary$period == b$period, , drop = FALSE]
  annual <- summary[summary$period != b$period, , drop = FALSE]
  f1_write_csv(main, file.path(out, "Figure1_global_summary.csv"))
  f1_write_csv(annual, file.path(out, "Figure1_annual_summary.csv"))
  f1_write_csv(d, file.path(out, "Figure1_grid_data.csv"))
  f1_write_csv(loaded$checks, file.path(out, "Figure1_tile_validation.csv"))
  f1_write_csv(dictionary, file.path(out, "Figure1_metric_dictionary.csv"))
  saveRDS(list(revision = FIGURE1_REVISION, fingerprint = b$fingerprint,
    estimates = b$draws$estimates, paired_draw_ids = b$draws$paired_draw_ids,
    aggregate_denominator = b$draws$aggregate_denominator,
    aggregate_cell_counts = b$draws$aggregate_cell_counts,
    note = "Unchanged saved posterior draws. Cross-model pairing is the original fixed working product-posterior coupling."),
    file.path(out, "Figure1_global_posterior_draws.rds"), compress = FALSE)

  w <- d$area_weight; we <- w[e]; support <- d$n_common_support_days
  domain <- data.frame(
    quantity = c("selected_land_points", "eligible_land_points", "posterior_draws", "expected_days_per_point",
      "eligible_area_fraction_of_selected_domain", "finite_history_area_day_fraction_of_selected_domain",
      "supported_area_day_fraction_of_selected_domain", "supported_date_fraction_within_eligible_area",
      "selected_domain_area_weight_sum", "eligible_area_weight_sum", "eligible_supported_area_day_weight_sum"),
    value = c(nrow(d), sum(e), b$meta$global_draws_actual, b$meta$n_expected_total_days,
      sum(we) / sum(w), sum(w * d$n_valid_history_days) / sum(w * d$n_expected_days),
      sum(w * support) / sum(w * d$n_expected_days), sum(we * support[e]) / sum(we * d$n_expected_days[e]),
      sum(w), sum(we), sum(we * support[e])),
    unit = c(rep("count", 4), rep("fraction", 4), "saved_area_weight", "saved_area_weight", "saved_area_weight_day"),
    domain = "Selected land-point domain; no population weighting and no unsampled-land inference")
  f1_write_csv(domain, file.path(out, "Figure1_domain_coverage.csv"))
  area_mean <- function(x) sum(we * x) / sum(we)
  spatial <- data.frame(
    quantity = c("area_fraction_with_positive_posterior_median_delta_U", "area_fraction_with_negative_posterior_median_delta_U",
      "area_fraction_with_zero_posterior_median_delta_U", "area_fraction_with_pointwise_CrI_above_zero",
      "area_fraction_with_pointwise_CrI_below_zero", "area_fraction_with_pointwise_CrI_including_zero",
      "posterior_expected_area_fraction_positive_delta_U", "posterior_expected_area_fraction_negative_delta_U"),
    value = c(area_mean(d$delta_U_median[e] > 0), area_mean(d$delta_U_median[e] < 0), area_mean(d$delta_U_median[e] == 0),
      area_mean(d$delta_U_lower[e] > 0), area_mean(d$delta_U_upper[e] < 0),
      area_mean(d$delta_U_lower[e] <= 0 & d$delta_U_upper[e] >= 0),
      area_mean(d$delta_U_prob_positive[e]), area_mean(d$delta_U_prob_negative[e])),
    denominator = "Eligible land area weights; each location's supported-date mean delta_U",
    interpretation = c(rep("Descriptive area fraction, without a posterior interval", 6),
      rep("Posterior expectation from weighted marginal sign probabilities", 2)))
  f1_write_csv(spatial, file.path(out, "Figure1_spatial_descriptives.csv"))
  fractions <- NULL
  if (!is.null(loaded$area_draws)) {
    f1_close(mean(loaded$area_draws[, "positive"]), spatial$value[7], "Posterior mean positive-area fraction")
    f1_close(mean(loaded$area_draws[, "negative"]), spatial$value[8], "Posterior mean negative-area fraction")
    fractions <- do.call(rbind, lapply(colnames(loaded$area_draws), function(nm) {
      x <- loaded$area_draws[, nm]; q <- f1_quantiles(x)
      data.frame(direction = nm, posterior_mean = mean(x), lower = q["lower"], median = q["median"], upper = q["upper"],
        unit = "fraction_of_eligible_land_area", estimand = "Sign of each location's supported-date mean delta_U; common posterior draws", row.names = NULL)
    }))
    f1_write_csv(fractions, file.path(out, "Figure1_area_fraction_posterior.csv"))
    saveRDS(list(draws = loaded$area_draws, paired_draw_ids = b$draws$paired_draw_ids,
      fingerprint = b$fingerprint, note = "Fractions use the sign of each cell's supported-date mean, not the fraction of individual dates with positive effects."),
      file.path(out, "Figure1_area_fraction_draws.rds"), compress = FALSE)
  }
  list(main = main, annual = annual, domain = domain, spatial = spatial, fractions = fractions)
}

f1_graphic_complete <- function(path, format) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 32) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  h <- as.integer(readBin(con, "raw", n = 8L))
  if (format == "pdf") return(identical(h[1:5], as.integer(charToRaw("%PDF-"))))
  if (format == "tiff") return(identical(h[1:2], c(73L, 73L)) || identical(h[1:2], c(77L, 77L)))
  if (format == "png") {
    if (!identical(h, c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))) return(FALSE)
    seek(con, where = file.info(path)$size - 12, origin = "start")
    tail <- as.integer(readBin(con, "raw", n = 12L))
    return(identical(tail[5:8], as.integer(charToRaw("IEND"))))
  }
  FALSE
}

f1_save_plot <- function(p, stem, width, height, config, tiff = FALSE) {
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else
    function(filename, ...) grDevices::pdf(file = filename, ..., useDingbats = FALSE)
  save_one <- function(format, device, ...) {
    target <- paste0(stem, ".", format)
    staging <- file.path(dirname(stem), "_figure_staging")
    if (!dir.exists(staging)) f1_assert(dir.create(staging), "Cannot create figure staging directory.")
    tmp <- tempfile(pattern = "figure_write_", tmpdir = staging, fileext = paste0(".", format))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    ggplot2::ggsave(tmp, plot = p, device = device, dpi = config$dpi,
      width = width, height = height, units = "mm", bg = "white", limitsize = FALSE, ...)
    if (!f1_graphic_complete(tmp, format) && format == "png") {
      warning("Incomplete PNG device write; retrying once with an explicit base bitmap device: ", basename(target), call. = FALSE)
      grDevices::png(filename = tmp, width = width / 25.4, height = height / 25.4,
        units = "in", res = config$dpi, bg = "white",
        type = if (isTRUE(capabilities("cairo"))) "cairo-png" else getOption("bitmapType", "windows"))
      dev <- grDevices::dev.cur()
      tryCatch(if (inherits(p, "grob")) grid::grid.draw(p) else print(p), finally = grDevices::dev.off(which = dev))
    }
    f1_assert(f1_graphic_complete(tmp, format), "Graphics device did not produce a complete file: ", target)
    f1_assert(!file.exists(target) && file.rename(tmp, target), "Cannot commit figure: ", target)
    invisible(target)
  }
  save_one("pdf", pdf_device)
  save_one("png", "png")
  if (tiff && isTRUE(config$save_tiff)) {
    if (isTRUE(capabilities("tiff"))) {
      save_one("tiff", "tiff", compression = "lzw")
    } else warning("TIFF is unavailable in this R build; PDF and PNG have been saved.", call. = FALSE)
  }
  invisible(stem)
}




f2_export_domain_results <- function(b, loaded, domains, spec, config, out) {
  original_dir <- file.path(out, "original_domain_tables"); dir.create(original_dir)
  original_tables <- f1_export_tables(b, loaded, original_dir)
  d <- loaded$data
  d$effective_eligible <- domains$membership[[config$main_domain]]
  d$main_domain <- config$main_domain
  def <- spec$definitions[spec$definitions$domain_id == config$main_domain, ]
  d$main_domain_reason <- ifelse(!d$eligible_for_area_day_summary, "Not eligible under original saved rules",
    ifelse(d$effective_eligible, "Retained in selected fixed domain", paste0("Excluded by fixed-domain rule: ", def$domain_id)))
  if (!is.null(spec$population)) {
    pi <- match(d$cell_id, spec$population$cell_id)
    for (nm in setdiff(names(spec$population), c("cell_id", "lon", "lat"))) d[[nm]] <- spec$population[[nm]][pi]
  }
  if (!is.null(spec$external)) d$external_reason <- spec$external$reason[match(d$cell_id, spec$external$cell_id)]
  if (!is.null(spec$elevation)) d$elevation_m <- spec$elevation$elevation_m[match(d$cell_id, spec$elevation$cell_id)]
  f1_write_csv(d, file.path(out, "Figure1_grid_data.csv"))
  f1_write_csv(domains$membership, file.path(out, "Figure1_grid_domain_membership.csv"))
  f1_write_csv(domains$coverage, file.path(out, "Figure1_domain_coverage.csv"))
  f1_write_csv(spec$definitions, file.path(out, "Figure1_domain_definitions.csv"))
  f1_write_csv(domains$summary, file.path(out, "Figure1_all_domain_summary.csv"))
  main_summary <- domains$summary[domains$summary$domain_id == config$main_domain, , drop = FALSE]
  f1_write_csv(main_summary, file.path(out, "Figure1_global_summary.csv"))
  f1_write_csv(domains$contrasts, file.path(out, "Figure1_domain_minus_original.csv"))
  f1_write_csv(domains$sign_summary, file.path(out, "Figure1_area_fraction_posterior.csv"))
  f1_write_csv(loaded$checks, file.path(out, "Figure1_tile_validation.csv"))
  saveRDS(list(schema = "Figure1-domain-posterior-v2", revision = FIGURE1_REVISION,
    fingerprint = b$fingerprint, period = b$period, main_domain = config$main_domain,
    domain_mean_draws = domains$draws, sign_fraction_draws = domains$sign_fraction_draws,
    paired_draw_ids = b$draws$paired_draw_ids, definitions = spec$definitions, coverage = domains$coverage,
    note = "All posterior means are recomputed on fixed selected locations using area times original supported days. No refitting, daily rescreening or joint SET-by-history validation."),
    file.path(out, "Figure1_domain_posterior_draws.rds"), compress = FALSE)
  templates <- file.path(out, "optional_input_templates"); dir.create(templates)
  td <- d[c("cell_id", "lon", "lat")]; td$run_fingerprint <- b$fingerprint
  td$include <- NA; td$reason <- ""
  f1_write_csv(td, file.path(templates, "external_domain_TEMPLATE.csv"))
  td$include <- NULL; td$reason <- NULL; td$elevation_m <- NA_real_
  f1_write_csv(td, file.path(templates, "elevation_TEMPLATE.csv"))
  f1_write_text(c("Optional input templates are intentionally incomplete and are never automatically applied.",
    "Populate every cell_id and provide the provenance/selection basis in the script configuration.",
    "external_domain: TRUE/FALSE in include, and a nonempty reason for every row. Do not use effect magnitude or significance alone to define inclusion.",
    "elevation: finite elevation_m for every row. Keep the original coordinate and run_fingerprint columns for checks.",
    "An external mask is user-defined evidence, not automatic certification of AOA. Elevation thresholds are domain sensitivity assumptions."),
    file.path(templates, "READ_FIRST.txt"))
  list(grid = d, main = main_summary, original_tables = original_tables)
}

f2_write_reports <- function(b, domains, spec, config, plots, out) {
  if (!is.null(spec$population)) return(p1_write_reports(b, domains, spec, config, plots, out))
  id <- config$main_domain; cov <- domains$coverage[domains$coverage$domain_id == id, ]
  s <- domains$summary[domains$summary$domain_id == id, ]
  text_interval <- function(m) {
    x <- s[s$metric == m, ]
    sprintf("%.3f（95%%可信区间 %.3f 至 %.3f）", x$median, x$lower, x$upper)
  }
  note <- c("Figure 1 更新版：数值结果与解释边界",
    paste0("研究期：", b$period, "；主图固定域：", id, "；", cov$domain_label, "。"),
    sprintf("主图保留 %d 个陆地点，占原有效分析域面积权重的 %.1f%%、面积×支持日分母的 %.1f%%。",
      cov$n_cells, 100 * cov$retained_area_fraction, 100 * cov$retained_area_day_fraction),
    paste0("ΔU = ", text_interval("delta_U"), " ℃ SET；ΔM = ", text_interval("delta_M"),
      " ℃ SET；ΔN = ", text_interval("delta_N"), " ℃ SET。"),
    paste0("样本P50固定SET下的辅助ΔTSV = ", text_interval("delta_TSV_P50"), "。P50不定义暖侧边界U。"),
    "U为标准化平均TSV=+1对应的SET。N=1/A为正向半宽，ΔN不另加负号；ΔU=ΔM+ΔN在每个共同抽样位置成立。",
    "各地图与森林图使用同一网格集合。每格仍使用原分析所有热史支持日期；未重新按每天暖侧根是否落入SET范围筛选。",
    "root_*为事后范围敏感性：先计算每个后验抽样中两模型暖侧根同时在训练SET边界内的支持日比例，再使用这个比例的后验中位数筛选固定网格。",
    sprintf("训练SET边界使用原分析保存的全样本P1–P99：%.3f 至 %.3f ℃ SET；不是样本最小值/最大值，也不是给定热史下的条件支持边界。", b$meta$set_bounds[1L], b$meta$set_bounds[2L]),
    "这不是单日后验概率、不是完整SET×热史联合适用域，也不是高海拔验证。相应可信区间条件于已选定的网格集合，未传播筛选规则本身的不确定性。",
    "history_*改变网格纳入所需的支持日覆盖比例，不改变入选网格内部原有支持日期。",
    if (is.null(spec$elevation)) "本次未提供海拔数据：没有按海拔排除任何地区，也没有验证高海拔人群迁移。" else
      paste0("已读取海拔数据；海拔筛选仅按配置阈值进行固定域敏感性。依据：", config$elevation_basis),
    if (is.null(spec$external)) "本次未提供外部适用域掩膜：没有新增训练资料证明的SET×热史适用域。" else
      paste0("外部固定域由用户提供；其科学有效性依赖提供的依据，而非本脚本自动认证。依据：", config$external_domain_basis),
    "当前保存数据中的热史支持是B21_raw与TV_lag1–3的边际范围及联合近邻支持；热史支持与SET根范围诊断是不同层面的检查。",
    "四指标受限域后验均值直接从检查点共同抽样重汇总。受限域的逐日绝对差异、根范围比例完整后验以及年度结果，不能由这四个全期间均值重建。",
    "original_domain_tables保存原域8指标及年度结果；不能将其混用为受限域结果。年度补图也明确只对应原域。",
    "没有重新拟合、没有加载Stan、没有读取ERA5、没有读取逐日能耗缓存。未支持或被排除区域不视为零效应，能耗端应使用一致范围。",
    "三张地图使用跨域可比较的共同完整线性色标；没有修改原值、裁剪极值或按显著性隐藏结果。",
    "EMPTY敏感性域保留空结果，不降低阈值、不回退原域。主图指定空域则报错。"
  )
  f1_write_text(note, file.path(out, "Figure1_results_and_limitations_CN.txt"))
  methods <- c("Fixed-domain posterior postprocessing, version 2",
    "The original numerical run is preserved. Main and scenario maps share the same full linear symmetric colour range.",
    "For domain D, draw s: mean_D(s) = sum_{i in D}(area_i * n_original_supported_i * cell_mean_i(s)) / sum_{i in D}(area_i * n_original_supported_i).",
    "Cells use identical original supported dates for M1 and M0. Static domain selection changes eligible locations only, never the date-specific support mask.",
    "Each grid mask is applied to both the maps and draw-wise aggregate. Intervals use quantile(type=8) at 0.025 and 0.975. The original fixed working product-posterior M1/M0 pairing is reused.",
    "Saved history support uses B21_raw and TV_lag1-3, with marginal training-range and joint k-nearest-neighbour checks. It does not establish joint SET-by-history support.",
    "The root-range fraction is first calculated per posterior draw across history-supported dates. The SAVED POSTERIOR MEDIAN of that date fraction is used for root_* location selection.",
    "root_* thresholds are post hoc diagnostic sensitivity assumptions, not data-validated AOA boundaries. Saved global training SET P1-P99 bounds (quantile type 7) are only a marginal check, not empirical minimum/maximum or conditional support bounds. Remaining locations can still lack conditional support or geographical transportability.",
    sprintf("Saved marginal training SET P1-P99: %.6f to %.6f degrees C SET.", b$meta$set_bounds[1L], b$meta$set_bounds[2L]),
    "Intervals for root-selected domains are conditional on the fixed selected set; uncertainty about the selection procedure is not propagated. Posterior draws outside the SET range are NOT discarded.",
    "External CSV masks are auditable user-defined domains; the software checks identifiers and data integrity, not the scientific validity of their inclusion criteria.",
    "No high-elevation exclusion is applied unless complete elevation data and explicit thresholds are supplied. Altitude alone is not proof of interpolation or extrapolation.",
    "Only delta_M, delta_N, delta_U and delta_TSV_P50 can be directly reaggregated from the saved full-period joint-grid draw means. The original-domain tables retain the additional unmodified metrics.",
    "Posterior sign-area fractions use the sign of location supported-date means. They are not daily frequencies, population shares or energy impacts.",
    "No restricted-domain annual posterior is fabricated. The annual figure, if written, remains explicitly on the ORIGINAL eligible domain.",
    "All coverage shares use the original eligible land-point domain as denominator; numerical area weights are not relabelled as km2.",
    "Source identification includes run fingerprints, paired draw IDs, partition-to-checkpoint equality, original aggregate reproduction, optional spatial-input provenance and file MD5 records.",
    paste("Period:", b$period, "; selected main domain:", id, "; source fingerprint:", b$fingerprint),
    "Domain definitions:", capture.output(print(spec$definitions, row.names = FALSE)),
    "Applicable methodological examples (they do not prescribe these diagnostic thresholds):",
    "French et al. Nature Communications 2023. https://doi.org/10.1038/s41467-023-40936-0 (non-analogue environmental masking).",
    "Aramburu-Merlos et al. Nature Food 2024. https://doi.org/10.1038/s43016-024-01029-3 (area of applicability).",
    "Meyer & Pebesma. Nature Communications 2022, Comment. https://doi.org/10.1038/s41467-022-29838-9 (limits of global-map assessment).")
  f1_write_text(methods, file.path(out, "Figure1_methods_and_estimands_EN.txt"))
  for (nm in names(plots$captions)) f1_write_text(plots$captions[[nm]], file.path(out, paste0("Figure1_caption_", nm, "_EN.txt")))
  f1_write_text(c("运行完成后，优先审阅：", "Figure1_thermal_response.pdf", "FigureS1_support_and_uncertainty.pdf",
    "FigureS2_domain_sensitivity.pdf", "Figure1_all_domain_summary.csv", "Figure1_domain_coverage.csv",
    "domain_figures/ 保存各非空敏感性域对应的完整地图与同步森林图；它们是探索性敏感性结果。",
    "仅更换main_domain即可从相同已完成结果重绘，不必重新运行暖侧模型。",
    "需要真正的SET×热史联合支持或按每日根范围重筛日期时，应另做逐日缓存后处理；本绘图程序不能替代这一步。"),
    file.path(out, "READ_ME_CN.txt"))
}

f2_main <- function(config = FIGURE1_CONFIG) {
  config <- utils::modifyList(FIGURE1_CONFIG, config)
  f1_assert(isTRUE(config$verify_checkpoints), "Version 2 requires checkpoint verification to recompute domain posteriors consistently.")
  f1_assert(f2_nonempty(config$main_domain), "main_domain must be one domain ID.")
  f1_assert(is.numeric(config$dpi) && length(config$dpi) == 1L && is.finite(config$dpi) && config$dpi >= 72, "Invalid dpi.")
  f1_assert(all(is.finite(c(config$figure_width_mm, config$figure_height_mm))) &&
    config$figure_width_mm > 0 && config$figure_height_mm > 0, "Invalid figure dimensions.")
  started <- Sys.time(); stage <- "resolve_completed_run"; warnings_seen <- character()
  message("Figure 1 updated postprocessing: ", FIGURE1_REVISION)
  message("Fixed-domain diagnostics only; no refitting, global prediction or daily-cache reading.")
  run <- f1_resolve_run(config)
  out <- file.path(run, paste0("Figure1_v2_", format(started, "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  if (!is.null(config$prepared_output_dir)) {
    out <- normalizePath(config$prepared_output_dir, winslash = "/", mustWork = TRUE)
    f1_assert(!file.exists(file.path(out, "Figure1_run_status.rds")), "Prepared output already has a figure run status; use a new directory.")
  } else {
    if (dir.exists(out)) out <- tempfile(pattern = paste0(basename(out), "_"), tmpdir = run)
    f1_assert(dir.create(out), "Cannot create output directory: ", out)
  }
  set_stage <- function(s) { stage <<- s; message("Stage: ", s) }
  result <- tryCatch(withCallingHandlers({
    set_stage("check_plotting_packages"); f1_packages(config)
    set_stage("validate_original_numerical_results"); b <- f1_load_identity(run)
    f1_verify_half_degree_run(b)
    set_stage("prepare_explicit_domain_definitions"); spec <- f2_prepare_domains(b, config)
    f1_write_csv(spec$definitions, file.path(out, "Figure1_domain_definitions.csv"))
    message("Main domain: ", config$main_domain, "; fixed domains: ", nrow(spec$definitions))
    if (is.null(spec$elevation)) message("No elevation CSV supplied: no automatic high-elevation exclusion.")
    message("Root-range sensitivities are diagnostics, not a validated SET-by-history applicability domain.")
    reader_config <- config; reader_config$domain_spec <- spec
    set_stage("stream_checkpoints_and_reaggregate_domains"); loaded <- f1_load_grid(b, reader_config)
    set_stage("validate_domain_posteriors"); domains <- f2_finalize_domains(b, loaded, spec)
    f1_write_csv(domains$coverage, file.path(out, "Figure1_domain_coverage.csv"))
    f1_write_csv(domains$summary, file.path(out, "Figure1_all_domain_summary.csv"))
    empty <- domains$coverage$domain_id[domains$coverage$status == "EMPTY"]
    if (length(empty)) message("Empty sensitivity domains retained as NA: ", paste(empty, collapse = ", "))
    f1_assert(!config$main_domain %in% empty, "Selected main domain has no eligible cells: ", config$main_domain,
      ". No fallback was applied. See Figure1_domain_coverage.csv and choose a scientifically justified domain.")
    set_stage("export_numerical_results"); tab <- f2_export_domain_results(b, loaded, domains, spec, config, out)
    # Freeze one common scale across original and every restricted-domain figure.
    v <- unlist(loaded$data[loaded$data$eligible_for_area_day_summary,
      paste0(c("delta_U", "delta_M", "delta_N"), "_median")], use.names = FALSE)
    scale_limit <- max(abs(v[is.finite(v)]))
    if (!is.null(config$map_scale_limit)) {
      f1_assert(is.numeric(config$map_scale_limit) && length(config$map_scale_limit) == 1 &&
        is.finite(config$map_scale_limit) && config$map_scale_limit > 0, "Invalid map_scale_limit.")
      if (config$map_scale_limit < scale_limit) warning("map_scale_limit expanded to the full original-domain range; no values were clipped.", call. = FALSE)
      scale_limit <- max(scale_limit, config$map_scale_limit)
    }
    if (scale_limit == 0) scale_limit <- .01
    cfg <- config; cfg$map_scale_limit <- scale_limit; cfg$period_label <- b$period
    cfg$main_domain_label <- spec$definitions$domain_label[match(config$main_domain, spec$definitions$domain_id)]
    cfg$domain_labels <- setNames(spec$definitions$domain_label, spec$definitions$domain_id)
    cfg$annual_summary <- tab$original_tables$annual
    set_stage("construct_updated_figures")
    plots <- f2_make_plots(tab$grid, tab$main, b$meta, cfg, domains$membership, domains$summary)
    f1_write_csv(plots$plot_data$maps, file.path(out, "Figure1_displayed_map_values.csv"))
    f1_write_csv(plots$plot_data$forest, file.path(out, "Figure1_displayed_forest_values.csv"))
    f2_write_reports(b, domains, spec, config, plots, out)
    set_stage("save_main_and_support_figures")
    f1_save_plot(plots$main, file.path(out, "Figure1_thermal_response"), config$figure_width_mm, config$figure_height_mm, config, tiff = TRUE)
    f1_save_plot(plots$support, file.path(out, "FigureS1_support_and_uncertainty"), 180, 215, config)
    if (!is.null(plots$sensitivity)) f1_save_plot(plots$sensitivity, file.path(out, "FigureS2_domain_sensitivity"), 185,
      max(140, 95 + nrow(spec$definitions) * 10), config)
    if (!is.null(plots$domain_map)) f1_save_plot(plots$domain_map, file.path(out, "FigureS3_main_domain_membership"), 180, 88, config)
    if (!is.null(plots$annual)) f1_save_plot(plots$annual, file.path(out, "FigureS4_annual_original_domain"), 130, 145, config)
    if (isTRUE(config$save_all_domain_figures)) {
      set_stage("save_nonempty_domain_sensitivity_maps")
      scenario_dir <- file.path(out, "domain_figures"); dir.create(scenario_dir)
      ids <- domains$coverage$domain_id[domains$coverage$status == "AVAILABLE"]
      for (id in ids) {
        message("Saving fixed-domain figure: ", id)
        dd <- loaded$data; dd$effective_eligible <- domains$membership[[id]]
        cc <- cfg; cc$main_domain <- id; cc$main_domain_label <- cfg$domain_labels[[id]]; cc$annual_summary <- NULL
        ps <- f2_make_plots(dd, domains$summary[domains$summary$domain_id == id, ], b$meta, cc, domains$membership, NULL)
        f1_save_plot(ps$main, file.path(scenario_dir, paste0("Figure1_", id)), config$figure_width_mm, config$figure_height_mm, config)
        f1_write_text(ps$captions$main, file.path(scenario_dir, paste0("Figure1_", id, "_caption_EN.txt")))
        f1_write_csv(ps$plot_data$maps, file.path(scenario_dir, paste0("Figure1_", id, "_map_values.csv")))
        f1_write_csv(ps$plot_data$forest, file.path(scenario_dir, paste0("Figure1_", id, "_forest_values.csv")))
      }
    }
    if (isTRUE(config$save_individual_panels)) {
      panel_dir <- file.path(out, "individual_panels"); dir.create(panel_dir)
      for (nm in names(plots$panels)) f1_save_plot(plots$panels[[nm]], file.path(panel_dir, paste0("panel_", nm)), 180,
        if (nm == "posterior_summary") 85 else 100, config)
    }
    saveRDS(list(main = plots$main, support = plots$support, sensitivity = plots$sensitivity,
      domain_map = plots$domain_map, annual_original = plots$annual, main_domain = config$main_domain,
      fingerprint = b$fingerprint, revision = FIGURE1_REVISION), file.path(out, "Figure1_plot_objects.rds"))
    set_stage("verify_committed_outputs_and_provenance")
    paths <- c(b$input_files, loaded$partition_files)
    if (!is.null(spec$external)) paths <- c(paths, attr(spec$external, "source_path"))
    if (!is.null(spec$elevation)) paths <- c(paths, attr(spec$elevation, "source_path"))
    inputs <- rbind(f1_file_inventory(paths), f1_file_inventory(loaded$checks$checkpoint, hash = FALSE))
    f1_write_csv(inputs, file.path(out, "Figure1_input_inventory.csv"))
    f1_write_text(capture.output(utils::sessionInfo()), file.path(out, "Figure1_sessionInfo.txt"))
    f1_write_text(if (length(warnings_seen)) unique(warnings_seen) else "No warnings.", file.path(out, "Figure1_warnings.txt"))
    invisible(gc())
    all_files <- list.files(out, full.names = TRUE, recursive = TRUE)
    is_staging <- grepl("(^|[/\\\\])_figure_staging([/\\\\]|$)", all_files)
    if (any(is_staging)) unlink(all_files[is_staging])
    files <- all_files[!is_staging]
    f1_assert(length(files) > 15L && all(file.info(files)$size > 0), "One or more committed output files are empty.")
    for (path in files[tolower(tools::file_ext(files)) %in% c("pdf", "png", "tiff")])
      f1_assert(f1_graphic_complete(path, tolower(tools::file_ext(path))), "Incomplete committed figure: ", path)
    mandatory <- c("Figure1_thermal_response.pdf", "Figure1_thermal_response.png", "Figure1_global_summary.csv",
      "Figure1_all_domain_summary.csv", "Figure1_domain_coverage.csv", "Figure1_domain_posterior_draws.rds")
    f1_assert(all(file.exists(file.path(out, mandatory))), "Required artifacts missing.")
    inventory <- f1_file_inventory(files); f1_write_csv(inventory, file.path(out, "Figure1_output_inventory.csv"))
    finished <- Sys.time()
    manifest <- list(revision = FIGURE1_REVISION, status = "COMPLETED", started = started, finished = finished,
      elapsed_minutes = as.numeric(difftime(finished, started, units = "mins")), source_run = run,
      source_fingerprint = b$fingerprint, source_metadata = b$meta, config = config, domain_definitions = spec$definitions,
      main_domain = config$main_domain, coverage = domains$coverage, shared_map_limit_degC_SET = scale_limit,
      input_inventory = inputs, output_inventory = inventory, warnings = unique(warnings_seen), sessionInfo = utils::sessionInfo())
    saveRDS(manifest, file.path(out, "Figure1_run_manifest.rds"))
    saveRDS(list(status = "COMPLETED", stage = "finished", revision = FIGURE1_REVISION,
      started = started, finished = finished, output_dir = out), file.path(out, "Figure1_run_status.rds"))
    message("Updated Figure 1 completed: ", out)
    message("Elapsed ", round(manifest$elapsed_minutes, 1), " min; warnings ", length(unique(warnings_seen)), "; empty domains ", length(empty))
    message("Share the main PDF, support/sensitivity PDFs, all-domain summary CSV and domain coverage CSV for interpretation.")
    list(output_dir = out, manifest = manifest, summary = tab$main, domains = domains$coverage)
  }, warning = function(w) {
    warnings_seen <<- c(warnings_seen, paste0(stage, ": ", conditionMessage(w)))
    message("Warning [", stage, "]: ", conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {
    msg <- paste0("Figure 1 v2 failed in stage: ", stage, "\n", conditionMessage(e), "\nOutput directory: ", out)
    f1_write_text(c(msg, "Original numerical outputs and checkpoints were not modified."), file.path(out, "Figure1_error.txt"))
    saveRDS(list(status = "FAILED", revision = FIGURE1_REVISION, stage = stage, started = started,
      finished = Sys.time(), error = conditionMessage(e), output_dir = out), file.path(out, "Figure1_run_status.rds"))
    f1_stop(msg)
  })
  invisible(result)
}


# Population aggregation module: source counts -> sampling-centred geographic cells.
# Conservative areal allocation on true pixel boundaries; no nearest-neighbour
# extraction, bilinear interpolation, population imputation or source-edge wrapping.
# Source: GHS-POP R2023A, E2015, EPSG:4326, 30 arc-seconds, population counts.
# terra cell vectors are row-major: https://rspatial.github.io/terra/reference/values.html

p1_abort <- function(...) stop(..., call. = FALSE)

p1_assert_integer <- function(x, label, tolerance = 1e-5) {
  if (any(!is.finite(x)) || any(abs(x - round(x)) > tolerance))
    p1_abort(label, " do not align with source pixel boundaries; no implicit resampling is allowed.")
  as.integer(round(x))
}

# Exact WGS84 ellipsoid surface area bounded by geographic meridians/parallels.
# Integral a^2*(1-e^2)*cos(phi)/(1-e^2*sin(phi)^2)^2 dphi dlambda.
# This is the full geographic cell area, NOT its land-only area.
p1_wgs84_rectangle_km2 <- function(south, north, width_degrees) {
  if (any(!is.finite(c(south, north, width_degrees))) || any(south < -90) ||
      any(north > 90) || any(north <= south) || any(width_degrees <= 0) ||
      any(width_degrees > 360)) p1_abort("Invalid geographic rectangle for cell area.")
  a <- 6378137
  f <- 1 / 298.257223563
  e2 <- f * (2 - f)
  e <- sqrt(e2)
  integral <- function(phi) {
    s <- sin(phi * pi / 180)
    s / (2 * (1 - e2 * s * s)) + atanh(e * s) / (2 * e)
  }
  (width_degrees * pi / 180) * a * a * (1 - e2) *
    (integral(north) - integral(south)) / 1e6
}

# Build column indices from true geographic edges, not aggregate() block origins.
# For a representative cell centred on -180 degrees this wraps both raster edges.
p1_source_columns <- function(centres, width, xmin, resolution, ncols) {
  nx <- p1_assert_integer(width / resolution, "Target longitude width")
  if (nx < 1L || nx > ncols) p1_abort("Invalid number of source columns per target cell.")
  starts <- p1_assert_integer((centres - width / 2 - xmin) / resolution,
                             "Target longitude edges")
  lapply(starts, function(s) as.integer(((s + seq_len(nx) - 1L) %% ncols) + 1L))
}

# Small pure-R core, independently testable with synthetic row-major values.
# Each output vector is indexed by source longitude column.
p1_band_statistics <- function(values, nrows, ncols, row_pixel_area_km2,
                               known_nodata_codes = -200) {
  if (!is.numeric(values) || length(values) != as.double(nrows) * ncols ||
      length(row_pixel_area_km2) != nrows ||
      any(!is.finite(row_pixel_area_km2)) || any(row_pixel_area_km2 <= 0))
    p1_abort("Invalid population band dimensions or pixel areas.")
  unexpected <- unique(values[is.finite(values) & values < 0 &
                               !(values %in% known_nodata_codes)])
  if (length(unexpected)) p1_abort("Unexpected negative population values: ",
                                  paste(head(unexpected, 20L), collapse = ", "),
                                  ". These values were not treated as zero.")
  if (any(is.infinite(values))) p1_abort("Infinite values encountered in population raster.")
  # terra returns row1 col1..ncol, then row2; one MATRIX COLUMN is one source row.
  z <- matrix(values, nrow = ncols, ncol = nrows)
  known <- is.finite(z) & z >= 0
  known_n <- rowSums(known)
  known_area <- as.numeric(known %*% row_pixel_area_km2)
  negative_n <- rowSums(is.finite(z) & z < 0)
  z[!known] <- 0  # lower-bound SUM only; missing flags are retained separately.
  list(population_lower_bound = rowSums(z), n_known = known_n,
       known_area_km2 = known_area, n_negative_nodata = negative_n,
       n_missing = nrows - known_n,
       n_na_nodata = nrows - known_n - negative_n,
       total_pixel_area_per_column_km2 = sum(row_pixel_area_km2))
}

p1_collect_targets <- function(stats, columns, area_km2, n_source_rows) {
  if (length(columns) != length(area_km2)) p1_abort("Target column/area mismatch.")
  n_pixels <- lengths(columns) * n_source_rows
  sums <- function(x) vapply(columns, function(cols) sum(x[cols]), numeric(1))
  n_known <- sums(stats$n_known)
  n_missing <- n_pixels - n_known
  complete <- n_missing == 0
  lower <- sums(stats$population_lower_bound)
  density_lower <- lower / area_km2
  known_area <- sums(stats$known_area_km2)
  total_source_area <- lengths(columns) * stats$total_pixel_area_per_column_km2
  if (any(abs(total_source_area - area_km2) > pmax(1e-8, area_km2 * 1e-10)))
    p1_abort("Source pixel areas do not conserve target cell area.")
  valid_area_fraction <- known_area / area_km2
  if (any(valid_area_fraction < -1e-10 | valid_area_fraction > 1 + 1e-10))
    p1_abort("Invalid known-area fraction.")
  data.frame(population_n = ifelse(complete, lower, NA_real_),
             population_lower_bound = lower,
             cell_area_km2 = area_km2,
             population_density = ifelse(complete, density_lower, NA_real_),
             population_density_lower_bound = density_lower,
             complete = complete,
             valid_pixel_fraction = n_known / n_pixels,
             valid_area_fraction = pmin(1, pmax(0, valid_area_fraction)),
             n_source_pixels = n_pixels, n_known = n_known, n_missing = n_missing,
             n_negative_nodata = sums(stats$n_negative_nodata),
             n_na_nodata = sums(stats$n_na_nodata),
             nodata_flag = ifelse(complete, "COMPLETE", ifelse(n_known == 0,
                       "ALL_UNKNOWN", "PARTIAL_UNKNOWN")), stringsAsFactors = FALSE)
}


# Canonical target segments; source padding outside [-180,180] is cropped, never
# wrapped back as an additional copy of population. A dateline target is split.
p1_longitude_overlap <- function(centres, width, xmin, resolution, ncols) {
  if (length(width) != 1L || !is.finite(width) || width <= 0 || width > 360 ||
      !is.finite(xmin) || !is.finite(resolution) || resolution <= 0 || ncols < 1)
    p1_abort("Invalid longitude geometry.")
  lapply(centres, function(centre) {
    west <- centre - width / 2
    east <- centre + width / 2
    segments <- if (west < -180) list(c(west + 360, 180), c(-180, east)) else
      if (east > 180) list(c(west, 180), c(-180, east - 360)) else list(c(west, east))
    entries <- lapply(segments, function(seg) {
      west <- max(-180, seg[1L]); east <- min(180, seg[2L])
      first <- max(1, floor((west - xmin) / resolution) + 1)
      last <- min(ncols, ceiling((east - xmin) / resolution))
      if (first > last || east <= west) return(NULL)
      cols <- seq.int(as.integer(first), as.integer(last))
      source_west <- xmin + (cols - 1) * resolution
      source_east <- xmin + cols * resolution
      overlap <- pmax(0, pmin(east, source_east) - pmax(west, source_west))
      take <- overlap > 0
      if (!any(take)) return(NULL)
      data.frame(col = cols[take], fraction = pmin(1, overlap[take] / resolution),
                 overlap_degrees = overlap[take])
    })
    entries <- Filter(Negate(is.null), entries)
    if (!length(entries)) return(data.frame(col = integer(), fraction = numeric(), overlap_degrees = numeric()))
    ans <- do.call(rbind, entries)
    if (anyDuplicated(ans$col)) ans <- stats::aggregate(ans[c("fraction", "overlap_degrees")],
                                                        list(col = ans$col), sum)
    if (any(ans$fraction > 1 + 1e-10)) p1_abort("A target cell allocates a source column more than once.")
    ans$fraction <- pmin(1, ans$fraction)
    ans[order(ans$col), , drop = FALSE]
  })
}

# Latitude fractions use exact WGS84 area ratios rather than latitude degrees.
# Source rows touched by this study (-55..80 centres) remain within physical poles.
p1_latitude_overlap <- function(centre, height, ymax, resolution, nrows) {
  north <- centre + height / 2
  south <- centre - height / 2
  first <- max(1, floor((ymax - north) / resolution) + 1)
  last <- min(nrows, ceiling((ymax - south) / resolution))
  if (first > last) return(data.frame(row = integer(), fraction = numeric(),
                                     overlap_area_per_degree_km2 = numeric()))
  rows <- seq.int(as.integer(first), as.integer(last))
  source_north <- ymax - (rows - 1) * resolution
  source_south <- ymax - rows * resolution
  intersect_north <- pmin(north, source_north)
  intersect_south <- pmax(south, source_south)
  take <- intersect_north > intersect_south
  if (!any(take)) return(data.frame(row = integer(), fraction = numeric(),
                                    overlap_area_per_degree_km2 = numeric()))
  rows <- rows[take]; source_north <- source_north[take]; source_south <- source_south[take]
  intersect_north <- intersect_north[take]; intersect_south <- intersect_south[take]
  if (any(source_north > 90 | source_south < -90))
    p1_abort("A source pixel intersecting a target crosses a geographic pole; its count footprint is undefined.")
  whole_area <- p1_wgs84_rectangle_km2(source_south, source_north, 1)
  overlap_area <- p1_wgs84_rectangle_km2(intersect_south, intersect_north, 1)
  fraction <- overlap_area / whole_area
  if (any(fraction > 1 + 1e-10 | fraction <= 0)) p1_abort("Invalid source latitude overlap fraction.")
  data.frame(row = rows, fraction = pmin(1, fraction),
             overlap_area_per_degree_km2 = overlap_area)
}

# Row-major terra vector -> columnwise weighted population and known-area sums.
# This assumes population density is uniform within a source pixel wherever a
# target boundary divides it; the result is a conservative areal allocation.
p1_weighted_band_statistics <- function(values, nrows, ncols, row_fraction,
                                         row_overlap_area_km2,
                                         known_nodata_codes = -200) {
  if (!is.numeric(values) || length(values) != as.double(nrows) * ncols ||
      length(row_fraction) != nrows || length(row_overlap_area_km2) != nrows ||
      any(!is.finite(c(row_fraction, row_overlap_area_km2))) ||
      any(row_fraction <= 0 | row_fraction > 1 + 1e-10) || any(row_overlap_area_km2 <= 0))
    p1_abort("Invalid weighted population band dimensions or overlap areas.")
  unexpected <- unique(values[is.finite(values) & values < 0 &
                               !(values %in% known_nodata_codes)])
  if (length(unexpected)) p1_abort("Unexpected negative population values: ",
      paste(head(unexpected, 20L), collapse = ", "), ". These values were not treated as zero.")
  if (any(is.infinite(values))) p1_abort("Infinite values encountered in population raster.")
  z <- matrix(values, nrow = ncols, ncol = nrows)
  known <- is.finite(z) & z >= 0
  n_known <- rowSums(known)
  negative_n <- rowSums(is.finite(z) & z < 0)
  known_area <- as.numeric(known %*% row_overlap_area_km2)
  z[!known] <- 0 # zero is used only to form a known-contribution sum, not imputation.
  population <- if (all(row_fraction == 1)) rowSums(z) else as.numeric(z %*% row_fraction)
  list(population_lower_bound = population, n_known = n_known,
       known_area_km2 = known_area, n_negative_nodata = negative_n,
       n_missing = nrows - n_known, n_na_nodata = nrows - n_known - negative_n,
       total_pixel_area_per_column_km2 = sum(row_overlap_area_km2))
}

p1_collect_weighted_targets <- function(stats, columns, area_km2, n_source_rows,
                                         coverage_tolerance = 1e-10) {
  n_targets <- length(columns)
  if (n_targets != length(area_km2)) p1_abort("Target overlap/area mismatch.")
  sums <- function(x, weighted = FALSE) vapply(columns, function(mm) {
    if (!nrow(mm)) return(0)
    sum(x[mm$col] * if (weighted) mm$fraction else 1)
  }, numeric(1))
  n_pixels <- vapply(columns, nrow, integer(1)) * n_source_rows
  n_known <- sums(stats$n_known)
  n_missing <- n_pixels - n_known
  lower <- sums(stats$population_lower_bound, TRUE)
  known_area <- sums(stats$known_area_km2, TRUE)
  covered_area <- vapply(columns, function(mm) sum(mm$fraction), numeric(1)) *
    stats$total_pixel_area_per_column_km2
  coverage <- covered_area / area_km2
  known_fraction <- known_area / area_km2
  if (any(coverage < -coverage_tolerance | coverage > 1 + coverage_tolerance) ||
      any(known_fraction < -coverage_tolerance | known_fraction > coverage + coverage_tolerance))
    p1_abort("Population overlap area is negative, duplicated or larger than the target cell.")
  coverage <- pmin(1, pmax(0, coverage))
  known_fraction <- pmin(coverage, pmax(0, known_fraction))
  missing_source <- pmax(0, 1 - coverage)
  missing_nodata <- pmax(0, coverage - known_fraction)
  complete <- n_missing == 0 & coverage >= 1 - coverage_tolerance
  density_lower <- lower / area_km2
  data.frame(population_n = ifelse(complete, lower, NA_real_),
    population_lower_bound = lower, cell_area_km2 = area_km2,
    population_density = ifelse(complete, density_lower, NA_real_),
    population_density_lower_bound = density_lower, complete = complete,
    valid_pixel_fraction = ifelse(n_pixels > 0, n_known / pmax(1, n_pixels), 0),
    valid_area_fraction = known_fraction, source_coverage_fraction = coverage,
    missing_source_area_fraction = missing_source,
    missing_nodata_area_fraction = missing_nodata,
    n_source_pixels = n_pixels, n_known = n_known, n_missing = n_missing,
    n_negative_nodata = sums(stats$n_negative_nodata), n_na_nodata = sums(stats$n_na_nodata),
    nodata_flag = ifelse(complete, "COMPLETE", ifelse(n_known == 0,
                 "ALL_UNKNOWN", "PARTIAL_UNKNOWN")), stringsAsFactors = FALSE)
}

# Sequential latitude-band reads. Adjacent non-aligned targets can share one
# boundary source row/column, with the corresponding fractions summing to <=1.
# No original data are changed. Only small new audit/matching files are written.
p1_aggregate_population <- function(pop_path, points, dx = 1, dy = 1, out_dir,
                                    config = list()) {
  cfg <- utils::modifyList(list(strict_product = TRUE, progress_every = 20L,
    known_nodata_codes = -200, max_band_mb = 512, coverage_tolerance = 1e-10), config)
  if (!requireNamespace("terra", quietly = TRUE)) p1_abort("Package terra is required.")
  if (!file.exists(pop_path)) p1_abort("Population TIFF not found: ", pop_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!dir.exists(out_dir)) p1_abort("Cannot create population audit directory: ", out_dir)
  if (!all(c("cell_id", "lat") %in% names(points))) p1_abort("points must contain cell_id and lat.")
  lon_name <- if ("lon_plot" %in% names(points)) "lon_plot" else "lon"
  if (!lon_name %in% names(points)) p1_abort("points must contain lon_plot or lon.")
  if (!nrow(points) || anyNA(points$cell_id) || anyDuplicated(points$cell_id))
    p1_abort("Population target cell IDs must be non-missing and unique.")
  if (length(dx) != 1L || length(dy) != 1L || !is.finite(dx) || !is.finite(dy) ||
      dx <= 0 || dx > 360 || dy <= 0 || dy > 180)
    p1_abort("Invalid target resolution.")
  lon <- as.numeric(points[[lon_name]]); lat <- as.numeric(points$lat)
  if (any(!is.finite(c(lon, lat))) || any(lon < -180 | lon >= 180) ||
      any(lat - dy / 2 < -90 | lat + dy / 2 > 90)) p1_abort("Invalid target coordinates.")
  if (anyDuplicated(paste(lon, lat, sep = ":"))) p1_abort("Duplicate target coordinates.")
  ras <- terra::rast(pop_path)
  rs <- terra::res(ras); ex <- as.vector(terra::ext(ras)); nc <- ncol(ras); nr <- nrow(ras)
  wkt <- terra::crs(ras)
  crs_details <- tryCatch(terra::crs(ras, describe = TRUE), error = function(e) NULL)
  na_flag <- tryCatch(terra::NAflag(ras), error = function(e) NA_real_)
  # Save geometry BEFORE product/coverage checks so failures remain diagnosable.
  audit <- list(schema = "population-aggregation-v2", created = format(Sys.time(), tz = "UTC"),
    source = normalizePath(pop_path, winslash = "/", mustWork = TRUE),
    source_file_size_bytes = file.info(pop_path)$size,
    source_file_mtime = as.character(file.info(pop_path)$mtime), source_crs = wkt,
    source_crs_details = crs_details, source_resolution = rs, source_extent = ex,
    source_nrow = nr, source_ncol = nc, source_nodata_metadata = na_flag,
    target_dx = dx, target_dy = dy, target_n = nrow(points),
    known_negative_nodata_codes = cfg$known_nodata_codes,
    coverage_tolerance = cfg$coverage_tolerance, status = "GEOMETRY_READ")
  saveRDS(audit, file.path(out_dir, "population_source_audit.rds"))
  utils::write.csv(data.frame(nrow = nr, ncol = nc, resolution_x = rs[1], resolution_y = rs[2],
    xmin = ex[1], xmax = ex[2], ymin = ex[3], ymax = ex[4],
    longitude_width = ex[2] - ex[1]), file.path(out_dir, "population_source_geometry.csv"), row.names = FALSE)
  message(sprintf("Population source geometry: %d rows x %d columns; resolution %.12g x %.12g degrees; extent [%.12g, %.12g] x [%.12g, %.12g].",
    nr, nc, rs[1], rs[2], ex[1], ex[2], ex[3], ex[4]))
  if (terra::nlyr(ras) != 1L || !terra::is.lonlat(ras))
    p1_abort("Population input must be one geographic longitude/latitude raster band.")
  epsg_ok <- !is.null(crs_details) && "code" %in% names(crs_details) &&
    any(as.character(crs_details$code) == "4326", na.rm = TRUE)
  epsg_ok <- epsg_ok || grepl('ID["EPSG",4326]', wkt, fixed = TRUE) ||
    grepl('AUTHORITY["EPSG","4326"]', wkt, fixed = TRUE)
  if (!epsg_ok) p1_abort("Population raster must explicitly identify EPSG:4326 / WGS84.")
  if (any(!is.finite(c(rs, ex))) || any(rs <= 0) || nc < 1L || nr < 1L ||
      ex[2L] <= ex[1L] || ex[4L] <= ex[3L]) p1_abort("Invalid population raster geometry.")
  if (ex[1L] >= 0 && ex[2L] > 180)
    p1_abort("This script expects the -180..180 longitude convention, not a 0..360 raster. Source geometry has been saved.")
  if (isTRUE(cfg$strict_product) && any(abs(rs - 1 / 120) > 1e-8))
    p1_abort("Expected nominal 30 arc-second GHS-POP resolution (tolerance 1e-8 degrees). Source geometry has been saved.")
  if (!is.numeric(cfg$coverage_tolerance) || length(cfg$coverage_tolerance) != 1L ||
      !is.finite(cfg$coverage_tolerance) || cfg$coverage_tolerance <= 0 || cfg$coverage_tolerance > 1e-7)
    p1_abort("Invalid numerical area-coverage tolerance.")
  areas_all <- p1_wgs84_rectangle_km2(lat - dy / 2, lat + dy / 2, dx)
  columns_all <- p1_longitude_overlap(lon, dx, ex[1L], rs[1L], nc)
  band_centres <- sort(unique(lat), decreasing = TRUE)
  if (length(band_centres) > 1L && any(-diff(band_centres) < dy - 1e-9))
    p1_abort("Representative target latitude bands overlap.")
  row_maps <- lapply(band_centres, p1_latitude_overlap, height = dy,
                     ymax = ex[4L], resolution = rs[2L], nrows = nr)
  row_allocation <- numeric(nr)
  for (mm in row_maps) if (nrow(mm)) row_allocation[mm$row] <- row_allocation[mm$row] + mm$fraction
  if (any(row_allocation > 1 + cfg$coverage_tolerance))
    p1_abort("Adjacent target latitude bands allocate a source row more than once.")
  maximum_rows <- max(vapply(row_maps, nrow, integer(1)))
  approximate_mb <- maximum_rows * nc * 8 / 1024^2
  if (approximate_mb > cfg$max_band_mb) p1_abort("One raw population latitude band requires ",
    round(approximate_mb), " MB; exceeds max_band_mb. Input resolution may be wrong.")
  aligned <- abs(ex[1L] + 180) < 1e-10 && abs(ex[2L] - 180) < 1e-10 &&
    all(vapply(columns_all, function(mm) all(abs(mm$fraction - 1) < 1e-12), logical(1))) &&
    all(vapply(row_maps, function(mm) all(abs(mm$fraction - 1) < 1e-12), logical(1)))
  audit$aggregation_path <- if (aligned) "aligned_count_sum" else "conservative_fractional_overlap"
  audit$source_padding_cropped <- ex[1L] < -180 || ex[2L] > 180
  audit$source_max_latitude_allocation <- max(row_allocation)
  audit$missing_policy <- paste("NoData and areas outside raster coverage are unknown, not zero population or water.",
    "Known counts provide a lower bound CONDITIONAL on within-pixel proportional areal allocation.",
    "Complete product-derived allocated density is reported only for complete source/known coverage within numerical tolerance.")
  audit$area_definition <- "Whole geographic target-cell area on WGS84 ellipsoid; no land-fraction adjustment."
  audit$spatial_assumption <- paste("Point climate response represents the sampling-centred geographic cell;",
    "population is uniform within each source pixel when target boundaries divide it.")
  audit$method <- paste("Conservative longitude and WGS84 latitude overlap fractions;",
    "canonical longitude clipping, explicit dateline target splitting, sequential source latitude bands.")
  audit$raw_band_memory_mb <- approximate_mb
  audit$status <- "READING"
  saveRDS(audit, file.path(out_dir, "population_source_audit.rds"))
  band_records <- vector("list", length(band_centres)); rows_out <- band_records
  begun <- proc.time()[["elapsed"]]
  terra::readStart(ras)
  on.exit(try(terra::readStop(ras), silent = TRUE), add = TRUE)
  for (b in seq_along(band_centres)) {
    idx <- which(lat == band_centres[b]); rmapping <- row_maps[[b]]
    current_columns <- columns_all[idx]
    union_fx <- numeric(nc)
    for (mm in current_columns) if (nrow(mm)) union_fx[mm$col] <- union_fx[mm$col] + mm$fraction
    if (any(union_fx > 1 + cfg$coverage_tolerance))
      p1_abort("Representative target longitude cells allocate a source pixel more than once.")
    if (nrow(rmapping)) {
      first_row <- min(rmapping$row); nread <- nrow(rmapping)
      if (!identical(as.integer(rmapping$row), seq.int(first_row, length.out = nread)))
        p1_abort("Noncontiguous source latitude overlap rows.")
      raw <- terra::readValues(ras, row = first_row, nrows = nread, col = 1L, ncols = nc, mat = FALSE)
      statistics <- p1_weighted_band_statistics(raw, nread, nc, rmapping$fraction,
        rmapping$overlap_area_per_degree_km2 * rs[1L], cfg$known_nodata_codes)
      rm(raw)
    } else {
      first_row <- NA_integer_; nread <- 0L
      zero <- numeric(nc)
      statistics <- list(population_lower_bound = zero, n_known = zero, known_area_km2 = zero,
        n_negative_nodata = zero, n_missing = zero, n_na_nodata = zero,
        total_pixel_area_per_column_km2 = 0)
    }
    cell_stats <- p1_collect_weighted_targets(statistics, current_columns, areas_all[idx], nread,
                                             coverage_tolerance = cfg$coverage_tolerance)
    rows_out[[b]] <- cbind(data.frame(.input_row = idx, cell_id = points$cell_id[idx],
                                     lon = lon[idx], lat = lat[idx]), cell_stats)
    direct_observed <- sum(statistics$population_lower_bound * union_fx)
    difference <- sum(cell_stats$population_lower_bound) - direct_observed
    if (abs(difference) > max(1e-7, abs(direct_observed) * 1e-12))
      p1_abort("Population allocation did not conserve known contributions within latitude band.")
    elapsed <- proc.time()[["elapsed"]] - begun
    band_records[[b]] <- data.frame(band = b, latitude = band_centres[b], source_row_start = first_row,
      source_row_count = nread, n_target_cells = length(idx), n_complete_cells = sum(cell_stats$complete),
      n_partial_cells = sum(cell_stats$nodata_flag == "PARTIAL_UNKNOWN"),
      n_all_unknown_cells = sum(cell_stats$nodata_flag == "ALL_UNKNOWN"),
      observed_population_sum = direct_observed, conservation_difference = difference,
      source_max_longitude_allocation = max(union_fx),
      n_raw_missing_pixels = sum(statistics$n_missing),
      n_raw_negative_nodata_pixels = sum(statistics$n_negative_nodata),
      target_missing_source_area_km2 = sum(cell_stats$missing_source_area_fraction * areas_all[idx]),
      target_missing_nodata_area_km2 = sum(cell_stats$missing_nodata_area_fraction * areas_all[idx]),
      elapsed_seconds = elapsed)
    rm(statistics)
    if (b == 1L || b %% max(1L, cfg$progress_every) == 0L || b == length(band_centres))
      message(sprintf("Population aggregation: %d/%d latitude bands; %d/%d target cells; elapsed %.1f min; estimated remaining %.1f min.",
        b, length(band_centres), sum(vapply(rows_out[seq_len(b)], nrow, integer(1))), nrow(points),
        elapsed / 60, elapsed / b * (length(band_centres) - b) / 60))
  }
  result <- do.call(rbind, rows_out)
  result <- result[order(result$.input_row), , drop = FALSE]; result$.input_row <- NULL
  rownames(result) <- NULL
  if (!identical(as.character(result$cell_id), as.character(points$cell_id)))
    p1_abort("Population output row IDs do not match target IDs.")
  band_audit <- do.call(rbind, band_records)
  audit$status <- "COMPLETED"; audit$elapsed_seconds <- proc.time()[["elapsed"]] - begun
  audit$n_complete <- sum(result$complete)
  audit$n_partial_unknown <- sum(result$nodata_flag == "PARTIAL_UNKNOWN")
  audit$n_all_unknown <- sum(result$nodata_flag == "ALL_UNKNOWN")
  audit$total_known_allocated_population_selected_cells <- sum(result$population_lower_bound)
  audit$max_abs_population_conservation_difference <- max(abs(band_audit$conservation_difference))
  audit$source_max_longitude_allocation <- max(band_audit$source_max_longitude_allocation)
  denom <- sum(result$n_source_pixels)
  audit$selected_known_touched_pixel_fraction <- if (denom > 0) sum(result$n_known) / denom else 0
  audit$selected_known_area_fraction <- sum(result$valid_area_fraction * result$cell_area_km2) / sum(result$cell_area_km2)
  audit$selected_source_coverage_fraction <- sum(result$source_coverage_fraction * result$cell_area_km2) / sum(result$cell_area_km2)
  saveRDS(audit, file.path(out_dir, "population_source_audit.rds"))
  utils::write.csv(band_audit, file.path(out_dir, "population_latitude_band_audit.csv"), row.names = FALSE, na = "")
  utils::write.csv(result, file.path(out_dir, "population_2015_aggregated_cells.csv"), row.names = FALSE, na = "")
  list(table = result, audit = audit, band_audit = band_audit)
}

# Population-specific domain definitions and orchestration.

p1_include <- function(d, threshold = 1, strict = FALSE) {
  z <- d$population_density_lower_bound
  f1_assert(is.numeric(z) && all(is.finite(z) & z >= 0), "Invalid population lower bounds.")
  if (strict) z > threshold else z >= threshold
}

p1_classify <- function(d, threshold = 1) {
  f1_assert(is.logical(d$complete) && !anyNA(d$complete), "Missing population completeness flags.")
  include <- p1_include(d, threshold)
  ifelse(include, ifelse(d$complete, "eligible_complete", "eligible_lower_bound"),
    ifelse(d$complete, "low_density", "unknown_population"))
}

p1_validate_matching <- function(d, b) {
  required <- c("cell_id", "lon", "lat", "population_n", "population_lower_bound",
    "cell_area_km2", "population_density", "population_density_lower_bound", "complete",
    "valid_pixel_fraction", "valid_area_fraction", "n_source_pixels", "n_known", "n_missing",
    "source_coverage_fraction", "missing_source_area_fraction", "missing_nodata_area_fraction")
  f1_assert(is.data.frame(d) && all(required %in% names(d)), "Incomplete population matching table.")
  f1_assert(!anyDuplicated(d$cell_id) && setequal(d$cell_id, b$land$cell_id), "Population matching cell IDs differ from run.")
  d <- d[match(b$land$cell_id, d$cell_id), , drop = FALSE]
  f1_close(d$lon, b$land$lon_plot, "Population/grid longitude")
  f1_close(d$lat, b$land$lat, "Population/grid latitude")
  f1_assert(all(is.finite(d$cell_area_km2) & d$cell_area_km2 > 0), "Invalid matched cell areas.")
  res <- b$meta$grid_resolution
  expected_area <- p1_wgs84_rectangle_km2(d$lat - res[["lat"]] / 2,
    d$lat + res[["lat"]] / 2, res[["lon"]])
  f1_close(d$cell_area_km2, expected_area, "Matched WGS84 cell areas")
  f1_assert(is.logical(d$complete) && !anyNA(d$complete), "Invalid matched completeness flags.")
  f1_assert(all(is.finite(d$population_lower_bound) & d$population_lower_bound >= 0), "Invalid known population counts.")
  for (nm in c("valid_pixel_fraction", "valid_area_fraction", "source_coverage_fraction",
    "missing_source_area_fraction", "missing_nodata_area_fraction"))
    f1_assert(all(is.finite(d[[nm]]) & d[[nm]] >= 0 & d[[nm]] <= 1), "Invalid population coverage: ", nm)
  for (nm in c("n_source_pixels", "n_known", "n_missing"))
    f1_assert(all(is.finite(d[[nm]]) & d[[nm]] >= 0 & d[[nm]] == floor(d[[nm]])), "Invalid source-pixel counts: ", nm)
  f1_assert(all(d$n_source_pixels >= 0) && all(d$n_known + d$n_missing == d$n_source_pixels),
    "Population source-pixel bookkeeping differs.")
  f1_close(d$missing_source_area_fraction, 1 - d$source_coverage_fraction, "Population uncovered area")
  f1_close(d$missing_nodata_area_fraction + d$valid_area_fraction, d$source_coverage_fraction,
    "Population missing and known area partition")
  f1_assert(identical(d$complete, d$n_missing == 0 & d$source_coverage_fraction >= 1 - 1e-10),
    "Population completeness requires both geographic coverage and known source pixels.")
  f1_assert(all(!d$complete | d$valid_area_fraction >= 1 - 1e-10), "Complete cells must have complete known area.")
  f1_assert(all(d$n_source_pixels > 0 | d$source_coverage_fraction == 0), "Zero intersecting pixels cannot cover target area.")
  f1_assert(all(d$n_known > 0 | (d$population_lower_bound == 0 & d$valid_area_fraction == 0)),
    "Zero known pixels must have zero allocated known population and area.")
  f1_close(d$population_density_lower_bound, d$population_lower_bound / d$cell_area_km2, "Population density units")
  f1_close(d$population_n, ifelse(d$complete, d$population_lower_bound, NA_real_), "Population exact count completeness")
  f1_close(d$population_density, ifelse(d$complete, d$population_density_lower_bound, NA_real_), "Population exact density completeness")
  f1_close(d$valid_pixel_fraction, ifelse(d$n_source_pixels > 0, d$n_known / d$n_source_pixels, 0), "Population pixel coverage")
  d$population_status <- p1_classify(d, 1)
  d$population_include_ge1 <- p1_include(d, 1)
  d$run_fingerprint <- b$fingerprint
  rownames(d) <- NULL
  d
}

p1_prepare_domains <- function(b, config) {
  d <- p1_validate_matching(config$population_table, b)
  defs <- data.frame(domain_id = "original", domain_label = "Original supported domain",
    kind = "original", threshold = NA_real_, basis = "Original eligible points and original supported dates", stringsAsFactors = FALSE)
  add <- function(id, label, kind, threshold, basis) {
    defs <<- rbind(defs, data.frame(domain_id = id, domain_label = label, kind = kind,
      threshold = threshold, basis = basis, stringsAsFactors = FALSE))
  }
  pop_basis <- paste("GHS-POP 2015 counts allocated by WGS84 source/target overlap area into sampling-centred representative cells;",
    "whole-cell WGS84 area; the observed-count lower bound is conditional on uniform within-source-pixel allocation;",
    "human-relevant geography, not an applicability-domain validation")
  add("pop_gt0", "Population > 0 / km2", "population_gt", 0, pop_basis)
  add("pop_ge1", "Population >= 1 / km2", "population_ge", 1, pop_basis)
  add("pop_ge5", "Population >= 5 / km2", "population_ge", 5, pop_basis)
  for (z in c(.50, .75)) add(sprintf("pop_ge1_history_%03d", round(z * 100)),
    sprintf("Pop >= 1; history >= %g%%", 100 * z), "population_history", z,
    "Population >=1 intersected with fraction of original expected dates supported; fixed locations, unchanged dates")
  for (z in c(.50, .80, .95)) add(sprintf("pop_ge1_root_%03d", round(z * 100)),
    sprintf("Pop >= 1; root range >= %g%%", 100 * z), "population_root", z,
    "Population >=1 intersected with saved posterior median root-range fraction; post hoc diagnostic, NOT validated AOA")
  f1_assert(config$main_domain %in% defs$domain_id, "Unknown population main_domain: ", config$main_domain)
  ext <- d[c("cell_id", "lon", "lat", "run_fingerprint")]
  ext$include <- d$population_include_ge1; ext$reason <- d$population_status
  attr(ext, "source_path") <- config$population_mask_csv
  list(definitions = defs, population = d, external = ext, elevation = NULL)
}

p1_domain_membership <- function(d, spec) {
  base <- as.logical(d$eligible_for_area_day_summary)
  f1_assert(!anyNA(base), "Missing original eligibility.")
  ix <- match(d$cell_id, spec$population$cell_id)
  f1_assert(!anyNA(ix), "Population table lacks checkpoint cell IDs.")
  p <- spec$population[ix, , drop = FALSE]
  ge1 <- p1_include(p, 1)
  ans <- data.frame(cell_id = d$cell_id)
  for (j in seq_len(nrow(spec$definitions))) {
    def <- spec$definitions[j, ]; cut <- def$threshold
    take <- switch(def$kind,
      original = base,
      population_gt = base & p1_include(p, cut, TRUE),
      population_ge = base & p1_include(p, cut),
      population_history = base & ge1 & is.finite(d$common_support_fraction) & d$common_support_fraction >= cut,
      population_root = base & ge1 & is.finite(d$U_pair_in_SET_range_median) & d$U_pair_in_SET_range_median >= cut,
      f1_stop("Unknown population domain kind: ", def$kind))
    f1_assert(length(take) == nrow(d) && !anyNA(take), "Invalid population domain: ", def$domain_id)
    ans[[def$domain_id]] <- take
  }
  ans
}

p1_write_reports <- function(b, domains, spec, config, plots, out) {
  id <- config$main_domain
  cov <- domains$coverage[domains$coverage$domain_id == id, ]
  s <- domains$summary[domains$summary$domain_id == id, ]
  p <- spec$population
  summary_lines <- vapply(c("delta_U", "delta_M", "delta_N", "delta_TSV_P50"), function(m) {
    x <- s[s$metric == m, ]; sprintf("%s = %.6f [95%% CrI %.6f, %.6f]; P(>0) = %.4f; unit = %s",
      m, x$median, x$lower, x$upper, x$prob_positive, x$unit)
  }, character(1))
  note <- c("人口筛选后的热响应结果：可审计的固定空间域后处理", "",
    paste("原运行指纹：", b$fingerprint), paste("原运行目录：", b$run),
    paste("时期：", b$period, "；主域：", id),
    sprintf("主域保留 %d 个原本有估计的格点；占原可估计域面积权重 %.2f%%、受支持日期面积权重 %.2f%%。",
      cov$n_cells, 100 * cov$retained_area_fraction, 100 * cov$retained_area_day_fraction),
    summary_lines, "",
    "1. 人口产品为GHS-POP R2023A 2015人口数量，不是人口密度；按源像元与目标区域WGS84面积交叠比例分配，再除以整格面积(km2)。",
    "1a. 源栅格两轴可能与研究网格错位；边界分摊假定源像元内人口按面积均匀分布。population_lower_bound条件于此分配假设和人口产品，不是真实人口的无条件数学下界。",
    "1b. 目标跨180度经线时分成两段；源像元超出[-180,180]的部分不回卷重复计数。源范围外的目标面积保持未知，新增coverage字段与NoData分开记录。",
    "2. 当前0.5度代表区域以ERA5采样点为中心(例如经度0度对应-0.25至0.25度)。气候响应是原生0.25度网格每2个坐标抽样所得的点估计，不宣称是0.5度气候空间平均。",
    "3. 密度分母包括代表区域内陆地及水体面积，属于整格密度；原land_fraction是中心点二元陆地判断，不能当真实陆地比例。",
    paste0("4. 预定义人口主规则pop_ge1为密度>=1人/km2；本次主图实际选域为", id, "。其他域定义见Figure1_domain_definitions.csv。缺测格点仅在已知人口下界已达该域门槛时保留。"),
    "5. 匹配表的population_status针对>=1规则；部分缺测且下界<1标为unknown_population，不能解释为低人口或零人口。FigureS3则按实际主域的门槛标识证据不足。-200仅作为官方NoData，不解释成水域。",
    sprintf("本次原选点总匹配点数=%d；完整=%d；部分缺测=%d；全缺测=%d。",
      nrow(p), sum(p$complete), sum(!p$complete & p$n_known > 0), sum(p$n_known == 0)),
    paste(">=1的人口判定状态(尚未交集原模型可估计域)：", paste(capture.output(print(table(p$population_status))), collapse = " ")),
    "6. 地图和总体区间同步使用相同空间域。总体仍采用原area_weight×原受支持天数逐后验抽样加权，不改成人口加权，不平均格点中位数。",
    "7. 保留原域、>0和>=5人口门槛，以及>=1交集热史支持比例和根范围诊断的敏感性结果；原域结果必须与原保存后验复现一致。",
    "8. 人口限制定义研究相关人群的地理范围，不能证明剩余地区不存在SET×热史或高海拔外推。它不能保证青藏高原全部被遮蔽，也不能按结果幅度另行删格。",
    "9. 0.5度范围内城市人口可能使整个代表区域过门槛；这与每个30角秒像元都有人口不同。人口固定为2015年的静态空间参照。",
    "10. 模型未传播人口产品误差或阈值选择不确定性；CrI条件于既定人口产品和空间筛选。人口门槛敏感性不是这些误差的概率区间。",
    "11. U=M+N为平均TSV=+1的暖端，N=1/A>0。所有差值均M1-M0，delta_N不额外加负号；P50仅是辅助delta_TSV的固定SET，不是U的定义。",
    "12. root范围是原支持日期中两模型根同时落在训练SET边际P1-P99内的日期比例，再取其后验中位数。按格点固定筛选，不丢弃某些抽样或逐日根。",
    "13. 根范围选择是事后模型诊断，区间条件于选择后的格点集合；不可称为经数据验证的适用域。原热史支持也不证明SET×热史联合支持。",
    "14. 原始年度图单独标注original domain；本脚本不编造人口限定域的年度后验、逐日绝对变化或制冷能耗。",
    "15. 选点覆盖不等于全球全部陆地和人口；不得把本表已知人口和称为全球总人口，也不据此报告未经验证的全球人口覆盖百分比。",
    "16. 原域和敏感性地图使用共同线性对称色标，无删极值或色阶截断。筛选后的图若颜色变浅，应解释为数值分布变化。")
  f1_write_text(note, file.path(out, "Figure1_results_and_limitations_CN.txt"))
  methods <- c("Population-restricted thermal-response postprocessing",
    paste("Revision:", POPULATION_REVISION, "; original run fingerprint:", b$fingerprint),
    paste("Period:", b$period, "; main domain:", id),
    "GHS-POP R2023A E2015 EPSG:4326 30 arc-second population counts are conservatively allocated into non-overlapping representative cells centred on the saved regularly sampled climate points. A source pixel contributes its count multiplied by WGS84 intersection area divided by its full pixel area; whole pixels retain their full counts.",
    "Actual source edges and resolution are used. The source is not required to have exactly 43200 columns or exactly 360 degrees of longitude. Target regions crossing the antimeridian are split into two canonical longitude intervals. Source parts outside [-180,180] are clipped, never wrapped into duplicate population. Geographic source gaps remain unknown.",
    "Conservative fractional-overlap allocation assumes a uniform within-source-pixel population density on the WGS84 ellipsoid. No nearest-neighbour or bilinear interpolation of population counts is used. Counts and assigned fractions are audited to prevent over-allocation; aggregation is conditional on this subpixel allocation assumption.",
    "Density = summed population counts / entire target-cell area in km2. Area is integrated on the WGS84 ellipsoid. It is not land-only density. The source centre-in-polygon binary land mask cannot supply land fractions.",
    "Complete targets have complete product-derived allocated counts and densities. For NoData or uncovered source regions, known nonnegative allocated counts yield a lower bound conditional on the stated subpixel allocation. A partially missing target is retained if this conditional lower bound meets the threshold; otherwise its status is unknown and excluded separately. Geographic coverage gaps and within-source NoData are reported separately; neither is treated as zero population or water.",
    "The main threshold is >=1 person/km2. Alternative thresholds (>0 and >=5) and intersections with saved history/root diagnostics are reported as fixed-domain sensitivities. Population is a fixed 2015 reference, not annual population or an additional posterior model.",
    "Climate projections were evaluated at sampled grid points, not spatially averaged over the representative population regions. For this run, native spacing is 0.25 degrees and sampling stride is 2, giving 0.5 degree centre spacing.",
    "For domain D and paired posterior draw s: mean_D(s) = sum_i(area_weight_i * n_supported_i * cell_mean_i(s)) / sum_i(area_weight_i * n_supported_i), i in D.",
    "The density denominator and posterior weights serve different purposes: density uses exact full-cell km2, whereas posterior aggregation preserves the original cos(latitude)-based area weights and original supported dates. No population-weighted estimand is substituted.",
    "Joint grid draw means are reused with the original fixed working product-posterior M1/M0 pairing; the same draw positions are used across all points. Intervals use quantile type=8 at 0.025/0.975. Maps use marginal posterior medians of supported-date means. Medians of components need not add.",
    "The saved per-cell means permit domain-specific full-period delta_M, delta_N, delta_U and delta_TSV_P50 reaggregation. They do not yield restricted-domain annual posteriors or daily absolute shifts; additional original-domain metrics remain separately labelled.",
    "U=M+N=(B+1)/A is the standardized mean-TSV=+1 SET root; M=B/A and N=1/A>0. Differences are M1 minus M0, with no additional sign reversal for N. P50 is only the reference for the auxiliary TSV response.",
    "History support checks B21_raw and TV_lag1-3 using saved marginal ranges and joint nearest-neighbour criteria. It does not establish joint SET-by-history support.",
    sprintf("The marginal training SET P1-P99 bounds are %.6f to %.6f degrees C SET.", b$meta$set_bounds[1], b$meta$set_bounds[2]),
    "Root-fraction sensitivities use the saved posterior median of the fraction of supported dates whose two roots fall within those marginal bounds. These post hoc output-related masks preserve all original supported dates at retained points; intervals are conditional on the fixed selection.",
    "Population relevance is not proof of geographical or physiological transportability. No automatic high-elevation exclusion is applied. The framework neither forces the Tibetan Plateau to disappear nor removes points based on effect size or significance.",
    "Population product uncertainty and domain-selection uncertainty are not propagated into the thermal posterior. Population totals and coverage diagnostics refer only to selected representative cells, not the entire global population.",
    "Official dataset: https://data.jrc.ec.europa.eu/dataset/2ff68a52-5b5b-4a22-8f40-c41da8332cfe",
    "Official product report (population units and NoData -200): https://publications.jrc.ec.europa.eu/repository/bitstream/JRC133256/JRC133256_01.pdf",
    "Population-domain example, not a validation of the present extrapolations: Tian et al. Communications Earth & Environment 2024. https://doi.org/10.1038/s43247-024-01784-y",
    "Domain definitions:", capture.output(print(spec$definitions, row.names = FALSE)))
  f1_write_text(methods, file.path(out, "Figure1_methods_and_estimands_EN.txt"))
  for (nm in names(plots$captions)) {
    extra <- if (nm == "main") paste("Population is a 2015 static geographic criterion, based on counts allocated by WGS84 overlap area divided by full representative-cell area. Boundary allocation assumes uniform population density within each source pixel; missing-data lower bounds are conditional on that assumption.",
      paste0("The primary population domain specified for this postprocessing, pop_ge1, requires density >=1 person/km2; the selected figure domain is ", id, ". Partially missing cells are retained only when their observed-population lower bound meets the applicable threshold."),
      "The main figure and posterior aggregate use the identical intersection with original thermal-history eligibility. This is not a validation of model extrapolation.") else
      "Population-domain diagnostics are conditional on a static 2015 population product; annual results, if shown, remain explicitly on the original domain."
    f1_write_text(c(plots$captions[[nm]], extra), file.path(out, paste0("Figure1_caption_", nm, "_EN.txt")))
  }
  f1_write_text(c("正式主图：Figure1_thermal_response.pdf/png/tiff (默认pop_ge1)。",
    "正式数值：Figure1_global_summary.csv；全部敏感性：Figure1_all_domain_summary.csv。",
    "检查筛选规模：Figure1_domain_coverage.csv；逐格原因：Figure1_grid_data.csv。",
    "人口匹配：population_grid_matching.csv/rds；population_mask_2015.csv可复用。",
    "FigureS1 保留原域支持和不确定性诊断；FigureS2 同时比较人口和根范围/热史敏感性。",
    "FigureS3 区分主域保留、排除、原模型无估计和人口证据不足；FigureS4 仅为原域年度结果。",
    "灰色并不都表示无人居住；完整分类见人口表。人口筛选不能证明已消除外推。",
    "所有人口密度按整格面积计算；总体热响应仍按原面积×支持天数加权，不是人口加权。",
    "如需重新画图，matching_rds仅填V2生成的population_grid_matching.rds，避免再次汇总人口TIF。",
    "人口匹配使用真实边界的保守面积分配；边缘分配假设源像元内均匀，不等于具有更高空间分辨率的真实人口资料。",
    "population目录在I盘；结果目录由脚本自动新建。无需重新跑ERA5-warm02或Stan。"), file.path(out, "READ_ME_CN.txt"))
}

p1_resolve_tif <- function(config) {
  if (!is.null(config$population_tif)) {
    f1_assert(is.character(config$population_tif) && length(config$population_tif) == 1L && file.exists(config$population_tif), "population_tif does not exist.")
    path <- config$population_tif
  } else {
    f1_assert(dir.exists(config$population_dir), "Population directory not found: ", config$population_dir)
    candidates <- list.files(config$population_dir, pattern = "\\.tiff?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    path <- candidates[grepl("GHS_POP_E2015_GLOBE_R2023A_4326_30ss", basename(candidates), ignore.case = TRUE)]
    f1_assert(length(path) == 1L, "Expected ONE unpacked GHS_POP_E2015_GLOBE_R2023A_4326_30ss TIF in ",
      config$population_dir, "; found ", length(path), ". Set population_tif explicitly if needed.")
  }
  if (isTRUE(config$strict_product)) f1_assert(grepl("GHS_POP_E2015_GLOBE_R2023A_4326_30ss", basename(path), ignore.case = TRUE),
    "Expected the GHS-POP 2015 R2023A WGS84 30ss count product, not density/built-up/another year.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

p1_main <- function(config = POPULATION_CONFIG) {
  cfg <- utils::modifyList(POPULATION_CONFIG, config)
  started <- Sys.time(); stage <- "preflight"; out <- NULL
  stage_to <- function(s) { stage <<- s; message("Stage: ", s) }
  result <- tryCatch({
    message("Population Figure 1: ", POPULATION_REVISION)
    message("Saved posterior postprocessing only. No model fitting, ERA5 projection or daily-cache reading.")
    f1_assert(isTRUE(cfg$strict_product), "Formal 0.5-degree population matching requires strict_product=TRUE.")
    f1_assert(is.numeric(cfg$gdal_cache_mb) && length(cfg$gdal_cache_mb) == 1L && is.finite(cfg$gdal_cache_mb) && cfg$gdal_cache_mb > 0, "Invalid GDAL cache size.")
    run <- f1_resolve_run(cfg)
    b <- f1_load_identity(run)
    resolution <- f1_verify_half_degree_run(b)
    points <- data.frame(cell_id = b$land$cell_id, lon = b$land$lon_plot, lat = b$land$lat)
    out <- file.path(run, paste0("Figure1_population_", format(started, "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
    if (dir.exists(out)) out <- tempfile(pattern = paste0(basename(out), "_"), tmpdir = run)
    f1_assert(dir.create(out), "Cannot create new output directory: ", out)
    f1_write_text(c(paste("Revision:", POPULATION_REVISION), paste("Original run:", run),
      paste("Fingerprint:", b$fingerprint), paste("Start:", started)), file.path(out, "population_workflow_provenance.txt"))
    if (!is.null(cfg$matching_rds)) {
      stage_to("validate_existing_population_matching")
      match_record <- f1_read(cfg$matching_rds)
      f1_assert(identical(match_record$schema, "ERA5-population-matching-v2") && identical(match_record$fingerprint, b$fingerprint), "Matching RDS schema/run differs. V2 requires fresh matching; set matching_rds=NULL if it was produced by V1.")
      f1_close(match_record$grid_resolution, resolution, "Reused population matching resolution")
      f1_assert(identical(match_record$population_year, 2015L) && identical(match_record$density_area, "whole_cell_WGS84_ellipsoid_km2") &&
        identical(match_record$nodata_rule, "allocated_nonnegative_population_lower_bound") &&
        identical(match_record$allocation, "uniform_source_pixel_WGS84_overlap"), "Matching RDS population definition differs.")
      if (isTRUE(cfg$strict_product)) f1_assert(isTRUE(match_record$strict_product), "Synthetic matching RDS cannot be reused for a production run.")
      match_record$table <- p1_validate_matching(match_record$table, b)
      f1_write_csv(f1_file_inventory(cfg$matching_rds), file.path(out, "population_reused_matching_source.csv"))
      message("Reusing validated population matching; no population raster read.")
    } else {
      stage_to("check_population_raster_package")
      if (!requireNamespace("terra", quietly = TRUE) && isTRUE(cfg$install_missing))
        utils::install.packages("terra", repos = "https://cloud.r-project.org", type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType"))
      f1_assert(requireNamespace("terra", quietly = TRUE), "Install package terra before running population matching.")
      old_cache <- Sys.getenv("GDAL_CACHEMAX", unset = NA_character_)
      Sys.setenv(GDAL_CACHEMAX = as.character(cfg$gdal_cache_mb))
      on.exit(if (is.na(old_cache)) Sys.unsetenv("GDAL_CACHEMAX") else Sys.setenv(GDAL_CACHEMAX = old_cache), add = TRUE)
      tif <- p1_resolve_tif(cfg)
      message("Population TIF: ", tif)
      stage_to("hash_population_source")
      source_inventory <- f1_file_inventory(tif)
      f1_write_csv(source_inventory, file.path(out, "population_raster_input_inventory.csv"))
      stage_to("aggregate_population_sequential_latitude_bands")
      agg <- p1_aggregate_population(tif, points, dx = unname(resolution["lon"]), dy = unname(resolution["lat"]), out_dir = out,
        config = list(strict_product = cfg$strict_product, progress_every = cfg$progress_every, known_nodata_codes = -200, max_band_mb = cfg$max_band_mb))
      table <- p1_validate_matching(agg$table, b)
      match_record <- list(schema = "ERA5-population-matching-v2", revision = POPULATION_REVISION,
        fingerprint = b$fingerprint, population_year = 2015L, strict_product = cfg$strict_product,
        density_area = "whole_cell_WGS84_ellipsoid_km2", nodata_rule = "allocated_nonnegative_population_lower_bound",
        allocation = "uniform_source_pixel_WGS84_overlap", grid_resolution = resolution, source_inventory = source_inventory, source_audit = agg$audit,
        band_audit = agg$band_audit, table = table, created = Sys.time())
    }
    stage_to("save_population_matching_and_masks")
    d <- match_record$table
    saveRDS(match_record, file.path(out, "population_grid_matching.rds"))
    f1_write_csv(d, file.path(out, "population_grid_matching.csv"))
    f1_write_csv(match_record$source_inventory, file.path(out, "population_product_provenance.csv"))
    f1_write_csv(match_record$band_audit, file.path(out, "population_conservation_audit.csv"))
    mask <- d[c("cell_id", "lon", "lat", "run_fingerprint")]
    mask$include <- d$population_include_ge1; mask$reason <- d$population_status
    mask_path <- file.path(out, "population_mask_2015.csv")
    f1_write_csv(mask, mask_path)
    message("Population threshold states (before original thermal eligibility):")
    print(table(d$population_status))
    plot_cfg <- utils::modifyList(FIGURE1_CONFIG, cfg)
    plot_cfg$RUN_DIR <- run; plot_cfg$population_table <- d
    plot_cfg$population_mask_csv <- mask_path; plot_cfg$prepared_output_dir <- out
    stage_to("population_domain_posterior_and_figures")
    answer <- f2_main(plot_cfg)
    finished <- Sys.time()
    saveRDS(list(status = "COMPLETED", revision = POPULATION_REVISION, started = started, finished = finished,
      elapsed_minutes = as.numeric(difftime(finished, started, units = "mins")), output_dir = out,
      fingerprint = b$fingerprint, matching_file = file.path(out, "population_grid_matching.rds")),
      file.path(out, "population_workflow_status.rds"))
    message("Population matching + Figure 1 COMPLETE in ", round(as.numeric(difftime(finished, started, units = "mins")), 1), " min.")
    message("Output: ", out)
    answer
  }, error = function(e) {
    msg <- paste0("Population Figure 1 failed in stage: ", stage, "\n", conditionMessage(e))
    if (!is.null(out) && dir.exists(out)) {
      f1_write_text(c(msg, paste("Output:", out), "Completed population matching, if present, can be reused via matching_rds."), file.path(out, "population_workflow_error.txt"))
      saveRDS(list(status = "FAILED", revision = POPULATION_REVISION, stage = stage, error = conditionMessage(e),
        started = started, finished = Sys.time(), output_dir = out), file.path(out, "population_workflow_status.rds"))
    }
    f1_stop(msg)
  })
  invisible(result)
}


# V5 plot-only entry point. Existing numerical outputs are inputs, never rewritten.

f3_read_bundle <- function(directory, config, identity) {
  f1_assert(is.character(directory) && length(directory) == 1L && !is.na(directory) &&
    dir.exists(directory), "Figure directory not found: ", directory)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  f1_assert(identical(dirname(directory), identity$run),
    "The population figure must be a direct child of the selected 0.5-degree run: ", identity$run)
  paths <- file.path(directory, c("Figure1_run_status.rds", "Figure1_run_manifest.rds", "Figure1_grid_data.csv"))
  f1_assert(all(file.exists(paths)), "Incomplete saved figure directory: ", directory)
  status <- f1_read(paths[1]); manifest <- f1_read(paths[2])
  f1_assert(identical(status$status, "COMPLETED") && identical(manifest$status, "COMPLETED"), "Saved figure run is not completed.")
  f1_assert(identical(manifest$source_fingerprint, identity$fingerprint), "Saved population figure fingerprint differs from the selected 0.5-degree warm02 run.")
  f1_assert(identical(manifest$main_domain, config$main_domain), "Saved figure main domain differs: ", manifest$main_domain)
  f1_assert(is.numeric(manifest$shared_map_limit_degC_SET) && length(manifest$shared_map_limit_degC_SET) == 1L &&
    is.finite(manifest$shared_map_limit_degC_SET) && manifest$shared_map_limit_degC_SET > 0, "Missing original shared colour-scale limit.")
  list(directory = directory, source_identity = identity, manifest = manifest, status = status, paths = paths)
}

f3_resolve_source <- function(config) {
  run <- f1_resolve_run(config)
  identity <- f1_load_identity(run)
  f1_verify_half_degree_run(identity)
  if (!is.null(config$FIGURE_DIR)) return(f3_read_bundle(config$FIGURE_DIR, config, identity))
  directories <- list.dirs(run, full.names = TRUE, recursive = FALSE)
  directories <- directories[grepl("^Figure1_population_", basename(directories))]
  bundles <- lapply(directories, function(path) tryCatch(f3_read_bundle(path, config, identity), error = function(e) NULL))
  bundles <- Filter(Negate(is.null), bundles)
  f1_assert(length(bundles) > 0, "No completed population Figure 1 output matches this run/domain. Set FIGURE_DIR to the folder containing Figure1_grid_data.csv and Figure1_run_manifest.rds.")
  stamps <- vapply(bundles, function(x) {
    stamp <- suppressWarnings(as.numeric(x$manifest$finished))
    if (length(stamp) == 1 && is.finite(stamp)) stamp else as.numeric(file.info(x$paths[2])$mtime)
  }, numeric(1))
  # Most recently completed matching run wins; incomplete/failed folders are ignored.
  bundles[[order(stamps, vapply(bundles, function(x) x$directory, character(1)), decreasing = TRUE)[1L]]]
}

f3_packages <- function(config) {
  need <- c("ggplot2", "maps", "scales", "gtable")
  missing <- need[!vapply(need, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) && isTRUE(config$install_missing)) utils::install.packages(missing,
    repos = "https://cloud.r-project.org", type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType"))
  missing <- need[!vapply(need, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  f1_assert(!length(missing), "Missing plotting packages: ", paste(missing, collapse = ", "))
  f1_assert(utils::packageVersion("ggplot2") >= "3.4.0", "ggplot2 >=3.4.0 is required.")
}

f3_read_grid <- function(bundle, config) {
  path <- bundle$paths[3]
  if (isTRUE(config$verify_saved_grid_hash)) {
    inventory <- bundle$manifest$output_inventory
    f1_assert(is.data.frame(inventory) && all(c("path", "md5") %in% names(inventory)), "Saved output inventory lacks grid hash.")
    canonical <- function(x) sub("/+$", "", gsub("\\\\", "/", x))
    original_output <- bundle$status$output_dir
    if (is.null(original_output)) original_output <- bundle$manifest$config$prepared_output_dir
    f1_assert(is.character(original_output) && length(original_output) == 1L, "Saved status lacks the original figure output directory.")
    expected_grid <- paste0(canonical(original_output), "/Figure1_grid_data.csv")
    item <- which(canonical(as.character(inventory$path)) == expected_grid)
    f1_assert(length(item) == 1L && !is.na(inventory$md5[item]) && nzchar(inventory$md5[item]), "Cannot identify one saved grid checksum.")
    f1_assert(identical(tolower(unname(tools::md5sum(path))), tolower(as.character(inventory$md5[item]))),
      "Saved grid CSV was changed after the completed run; its checksum differs.")
  }
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("cell_id", "lon", "lat", "effective_eligible", "eligible_for_area_day_summary", "main_domain",
    "delta_U_median", "delta_M_median", "delta_N_median")
  f1_assert(all(required %in% names(d)), "Saved grid lacks: ", paste(setdiff(required, names(d)), collapse = ", "))
  land <- bundle$source_identity$land
  f1_assert(nrow(d) == nrow(land) && !anyNA(d$cell_id) && !anyDuplicated(d$cell_id) &&
    setequal(d$cell_id, land$cell_id), "Saved population grid differs from the selected 0.5-degree run's land-point inventory.")
  f1_assert(all(is.finite(d$lon) & d$lon >= -180 & d$lon <= 180 & is.finite(d$lat) & d$lat >= -90 & d$lat <= 90), "Invalid saved map coordinates.")
  ix <- match(d$cell_id, land$cell_id)
  f1_close(d$lon, land$lon_plot[ix], "Saved grid longitude and 0.5-degree run inventory")
  f1_close(d$lat, land$lat[ix], "Saved grid latitude and 0.5-degree run inventory")
  f1_assert(!anyNA(d$main_domain) && all(d$main_domain == config$main_domain), "Saved grid mixes domains.")
  for (nm in c("effective_eligible", "eligible_for_area_day_summary"))
    f1_assert(is.logical(d[[nm]]) && !anyNA(d[[nm]]), "Invalid saved inclusion field: ", nm)
  f1_assert(any(d$effective_eligible) && !any(d$effective_eligible & !d$eligible_for_area_day_summary), "Invalid selected domain.")
  if ("run_fingerprint" %in% names(d)) f1_assert(!anyNA(d$run_fingerprint) && all(d$run_fingerprint == bundle$source_identity$fingerprint), "Grid fingerprint differs.")
  for (nm in c("delta_U_median", "delta_M_median", "delta_N_median"))
    f1_assert(is.numeric(d[[nm]]) && all(is.finite(d[[nm]][d$effective_eligible])), "Nonfinite displayed estimates: ", nm)
  limit <- bundle$manifest$shared_map_limit_degC_SET
  observed <- max(abs(unlist(d[d$effective_eligible, c("delta_U_median", "delta_M_median", "delta_N_median")], use.names = FALSE)))
  f1_assert(observed <= limit + 1e-10, "Original scale does not contain all selected map values. No clipping was applied.")
  meta <- bundle$manifest$source_metadata
  resolution <- meta$grid_resolution[c("lon", "lat")]
  f1_assert(length(resolution) == 2 && all(is.finite(resolution) & resolution > 0), "Invalid saved map resolution.")
  f1_close(resolution, bundle$source_identity$meta$grid_resolution[c("lon", "lat")], "Population figure/source grid resolution")
  list(grid = d, meta = meta, scale_limit = limit)
}

f3_make_maps <- function(d, meta, limit, config) {
  res <- meta$grid_resolution
  lat_range <- meta$global_config$lat_range
  if (is.null(lat_range)) lat_range <- range(d$lat) + c(-1, 1) * res[["lat"]] / 2
  f1_assert(length(lat_range) == 2 && all(is.finite(lat_range)) && lat_range[2] > lat_range[1], "Invalid latitude limits.")
  metrics <- c("delta_M", "delta_N", "delta_U")
  metric_titles <- c(delta_M = "Neutral position, \u0394M", delta_N = "Neutral-to-warm half-width, \u0394N", delta_U = "Warm boundary, \u0394U")
  world <- ggplot2::map_data("world")
  grey <- "#D9D9D9"
  # Five evenly spaced tick positions; limits and colours are unchanged.
  br <- seq(-limit, limit, length.out = 5)
  tick_labels <- function(z) {
    labels <- scales::label_number(accuracy = if (limit >= 1) .1 else .01)(z)
    labels[abs(z) < 1e-12] <- "0"
    labels
  }
  letter_space <- if (isTRUE(config$show_panel_letters)) config$panel_letter_space_mm else 0
  th <- ggplot2::theme_minimal(base_size = config$base_font_size, base_family = config$font_family) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = if (isTRUE(config$draw_map_frame)) ggplot2::element_rect(fill = NA,
        colour = config$frame_colour, linewidth = config$frame_linewidth_mm) else ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(config$map_margin_mm + letter_space, config$map_margin_mm, config$map_margin_mm, config$map_margin_mm, unit = "mm"),
      plot.title = ggplot2::element_text(size = config$base_font_size, hjust = 0),
      plot.tag = ggplot2::element_text(size = config$panel_letter_size_pt, face = "bold", hjust = 0, vjust = 1),
      plot.tag.position = c(0.01, 1),
      legend.position = "bottom", legend.box = "vertical", legend.box.just = "center", legend.justification = "center",
      legend.title = ggplot2::element_text(size = config$base_font_size),
      legend.text = ggplot2::element_text(size = config$base_font_size - .5),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.height = grid::unit(2.3, "mm"), legend.key.width = grid::unit(3, "mm"),
      legend.box.spacing = grid::unit(0, "mm"), legend.spacing.x = grid::unit(1, "mm"),
      legend.spacing.y = grid::unit(1, "mm"),
      legend.margin = ggplot2::margin(0, 1, 0, 1, unit = "mm"))
  plots <- list(); tables <- list()
  for (j in seq_along(metrics)) {
    metric <- metrics[j]
    x <- d[c("cell_id", "lon", "lat", "effective_eligible")]
    x$metric <- metric; x$saved_posterior_median <- d[[paste0(metric, "_median")]]
    x$value <- ifelse(x$effective_eligible, x$saved_posterior_median, NA_real_)
    tables[[metric]] <- x
    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data = world, ggplot2::aes(x = long, y = lat, group = group),
        inherit.aes = FALSE, fill = grey, colour = NA, show.legend = FALSE) +
      ggplot2::geom_tile(data = x, ggplot2::aes(x = lon, y = lat, fill = value),
        width = res[["lon"]], height = res[["lat"]], inherit.aes = FALSE) +
      ggplot2::geom_polygon(data = world, ggplot2::aes(x = long, y = lat, group = group),
        inherit.aes = FALSE, fill = NA, colour = "#777777", linewidth = .09, show.legend = FALSE) +
      ggplot2::geom_point(data = data.frame(lon = 0, lat = 0, category = "No estimate / outside analysis domain"),
        ggplot2::aes(x = lon, y = lat, shape = category), inherit.aes = FALSE,
        colour = grey, alpha = 0, size = 2, show.legend = TRUE) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D55E00", midpoint = 0,
        limits = c(-limit, limit), breaks = br, na.value = grey,
        labels = tick_labels,
        name = "Posterior median shift (\u00b0C SET)",
        guide = ggplot2::guide_colourbar(order = 1, direction = "horizontal", title.position = "top", title.hjust = .5,
          barwidth = grid::unit(config$colourbar_width_mm, "mm"), barheight = grid::unit(config$colourbar_height_mm, "mm"))) +
      ggplot2::scale_shape_manual(values = 15, name = NULL,
        guide = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1, size = 2.2))) +
      ggplot2::coord_quickmap(xlim = c(-180, 180), ylim = lat_range, expand = FALSE, clip = "on") + th +
      ggplot2::labs(title = if (isTRUE(config$show_metric_titles)) metric_titles[[metric]] else NULL,
        subtitle = NULL, caption = NULL, tag = if (isTRUE(config$show_panel_letters)) letters[j] else NULL)
    plots[[metric]] <- p
  }
  # Use identically sized viewports to guarantee that the lower map is centred
  # and has exactly the same physical size as each of the two upper maps.
  width <- config$figure_width_mm
  map_width <- (width - 2 * config$outer_margin_mm - config$column_gap_mm) / 2
  f1_assert(map_width > 30, "Figure is too narrow for two global maps.")
  aspect <- diff(lat_range) / 360 / cos(mean(lat_range) * pi / 180)
  title_space <- letter_space + if (isTRUE(config$show_metric_titles)) 5 else 0
  map_height <- (map_width - 2 * config$map_margin_mm) * aspect + 2 * config$map_margin_mm + title_space
  height <- 2 * config$outer_margin_mm + 2 * map_height + config$row_gap_mm + config$legend_height_mm
  layout <- data.frame(metric = metrics, panel_letter = if (isTRUE(config$show_panel_letters)) letters[seq_along(metrics)] else rep("", length(metrics)),
    x_mm = c(config$outer_margin_mm + map_width / 2, width - config$outer_margin_mm - map_width / 2, width / 2),
    y_mm = c(rep(height - config$outer_margin_mm - map_height / 2, 2),
      height - config$outer_margin_mm - 1.5 * map_height - config$row_gap_mm),
    width_mm = rep(map_width, 3), height_mm = rep(map_height, 3), stringsAsFactors = FALSE)
  all_grobs <- ggplot2::ggplotGrob(plots[[1]])
  legend_indices <- which(grepl("^guide-box", all_grobs$layout$name))
  legend_indices <- legend_indices[!vapply(all_grobs$grobs[legend_indices], inherits, logical(1), what = "zeroGrob")]
  f1_assert(length(legend_indices) == 1L, "Cannot identify exactly one shared legend.")
  legend <- all_grobs$grobs[[legend_indices]]
  map_grobs <- lapply(seq_along(metrics), function(i) {
    g <- ggplot2::ggplotGrob(plots[[i]] + ggplot2::theme(legend.position = "none"))
    grid::editGrob(g, name = paste0("map_", metrics[i]), vp = grid::viewport(
      x = grid::unit(layout$x_mm[i], "mm"), y = grid::unit(layout$y_mm[i], "mm"),
      width = grid::unit(map_width, "mm"), height = grid::unit(map_height, "mm")))
  })
  legend <- grid::editGrob(legend, name = "shared_legend", vp = grid::viewport(
    x = grid::unit(width / 2, "mm"), y = grid::unit(config$outer_margin_mm + config$legend_height_mm / 2, "mm"),
    width = grid::unit(width - 2 * config$outer_margin_mm, "mm"), height = grid::unit(config$legend_height_mm, "mm")))
  figure <- grid::grobTree(children = do.call(grid::gList, c(map_grobs, list(legend))), name = "three_global_maps")
  list(figure = figure, plots = plots, layout = layout, width_mm = width, height_mm = height,
    data = do.call(rbind, tables), scale_limit = limit, legend_breaks = br, legend_box = "vertical")
}

f3_main <- function(config = THREE_MAP_CONFIG) {
  config <- utils::modifyList(THREE_MAP_CONFIG, config)
  started <- Sys.time(); out <- NULL; stage <- "configuration"; warnings_seen <- character()
  set_stage <- function(s) { stage <<- s; message("Stage: ", s) }
  result <- tryCatch(withCallingHandlers({
    message("Three-map plot-only revision: ", THREE_MAP_REVISION)
    message("Reading saved figure tables only; no population raster, posterior checkpoints or model computation.")
    for (nm in c("figure_width_mm", "frame_linewidth_mm", "base_font_size", "dpi", "legend_height_mm", "panel_letter_size_pt", "colourbar_width_mm", "colourbar_height_mm"))
      f1_assert(is.numeric(config[[nm]]) && length(config[[nm]]) == 1L && is.finite(config[[nm]]) && config[[nm]] > 0, "Invalid config: ", nm)
    for (nm in c("outer_margin_mm", "column_gap_mm", "row_gap_mm", "map_margin_mm", "panel_letter_space_mm"))
      f1_assert(is.numeric(config[[nm]]) && length(config[[nm]]) == 1L && is.finite(config[[nm]]) && config[[nm]] >= 0, "Invalid spacing: ", nm)
    f1_assert(config$dpi >= 72 && config$base_font_size >= 6, "Use dpi>=72 and base_font_size>=6.")
    f1_assert(config$colourbar_width_mm <= config$figure_width_mm - 2 * config$outer_margin_mm - 10, "Colour bar is wider than the available figure space.")
    set_stage("find_completed_population_figure"); bundle <- f3_resolve_source(config)
    message("Source figure directory: ", bundle$directory)
    set_stage("check_plot_packages"); f3_packages(config)
    set_stage("validate_saved_grid_and_scale"); saved <- f3_read_grid(bundle, config)
    out <- file.path(bundle$directory, paste0("three_maps_", format(started, "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
    if (dir.exists(out)) out <- tempfile(pattern = paste0(basename(out), "_"), tmpdir = bundle$directory)
    f1_assert(dir.create(out), "Cannot create new figure output: ", out)
    set_stage("construct_three_maps"); maps <- f3_make_maps(saved$grid, saved$meta, saved$scale_limit, config)
    f1_write_csv(maps$data, file.path(out, "Figure1_three_maps_values.csv"))
    f1_write_csv(maps$layout, file.path(out, "Figure1_three_maps_layout_mm.csv"))
    set_stage("save_pdf_png_tiff")
    f1_save_plot(maps$figure, file.path(out, "Figure1_thermal_response"), maps$width_mm, maps$height_mm, config, tiff = TRUE)
    if (isTRUE(config$save_individual_maps)) for (metric in names(maps$plots))
      f1_save_plot(maps$plots[[metric]] + ggplot2::theme(legend.position = "none"), file.path(out, paste0("Map_", metric)),
        maps$layout$width_mm[1], maps$layout$height_mm[1], config)
    saveRDS(list(revision = THREE_MAP_REVISION, figure = maps$figure, plots = maps$plots,
      layout = maps$layout, width_mm = maps$width_mm, height_mm = maps$height_mm,
      fingerprint = bundle$source_identity$fingerprint, scale_limit = maps$scale_limit,
      legend_breaks = maps$legend_breaks, legend_box = maps$legend_box), file.path(out, "Figure1_three_maps_plot_objects.rds"))
    period <- paste(range(saved$meta$years), collapse = "-" )
    caption <- paste0("Figure 4. Global thermal-sensation shifts associated with short-term thermal history. ",
      "a, Neutral-position shift (\u0394M; predicted mean thermal sensation vote, TSV = 0). ",
      "b, Neutral-to-warm half-width change (\u0394N). ",
      "c, Warm-side boundary shift (\u0394U; predicted mean TSV = +1), where \u0394U = \u0394M + \u0394N. ",
      "Colours show posterior medians of changes averaged over supported dates during ", period, " using ERA5 thermal histories sampled every two native 0.25-degree coordinates (0.5-degree centre spacing, not spatial averaging). ",
      "Changes are expressed in \u00b0C standard effective temperature (SET), calculated as the model with thermal history minus the model without thermal history. ",
      "Each 0.5-degree display tile shows its sampled point, not a tile-area mean. ",
      "Maps retain representative cells whose mean population density is \u22651 person km\u207b\u00b2 (GHS-POP 2015). ",
      "Orange and blue indicate positive and negative changes, respectively; grey denotes excluded regions or unavailable estimates.")
    f1_write_text(caption, file.path(out, "Figure1_three_maps_caption_EN.txt"))
    f1_write_text(c("仅重绘版完成。", "上左a=ΔM中性位置；上右b=ΔN正半宽；下方居中c=ΔU暖侧边界。",
      "上排呈现位置与半宽两个组成，下方呈现合成的暖侧边界变化。ΔU=ΔM+ΔN在逐后验抽样层面成立，各面板的边际后验中位数不必精确相加。",
      "N沿用正值半宽定义，ΔN不添加负号。论文图号为Figure 4；图像文件名沿用Figure1_thermal_response，英文图注已同步更新。",
      "恢复a/b/c图号；取消顶部注释及英文子图标题。色带位于c图下方居中，灰色说明另起一行居中。",
      "图框为每张地图的细灰矩形边界，非整个复合图外边框。draw_map_frame可设为FALSE。",
      "主图文件：Figure1_thermal_response.pdf/png/tiff。各指标的位置对应关系必须保留在论文图注中。",
      "全部格点后验中位数、人口掩膜、原研究范围及原共同色标均沿用；未计算新的总体均值。",
      paste("来源目录：", bundle$directory), paste("本次目录：", out),
      sprintf("画布 %.2f x %.2f mm；三图各 %.2f x %.2f mm。", maps$width_mm, maps$height_mm, maps$layout$width_mm[1], maps$layout$height_mm[1]),
      "参考图框实例（实际图件观察，不是期刊强制格式）：",
      "Tian et al. 2024 Communications Earth & Environment, Fig.1：细灰矩形框。https://doi.org/10.1038/s43247-024-01784-y",
      "Fan & McColl 2024 Communications Earth & Environment, Fig.4：细矩形框。https://doi.org/10.1038/s43247-024-01930-6",
      "Li et al. 2024 Nature Communications, Fig.2a：沿投影外缘的曲线框。https://doi.org/10.1038/s41467-024-51355-0"), file.path(out, "READ_ME_three_maps_CN.txt"))
    f1_write_csv(f1_file_inventory(bundle$paths), file.path(out, "Figure1_three_maps_inputs.csv"))
    f1_write_text(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
    f1_write_text(if (length(warnings_seen)) unique(warnings_seen) else "No warnings.", file.path(out, "warnings.txt"))
    set_stage("verify_figures")
    files <- list.files(out, full.names = TRUE)
    files <- files[!file.info(files)$isdir]
    f1_assert(all(file.info(files)$size > 0), "Empty output file.")
    for (path in files[tolower(tools::file_ext(files)) %in% c("pdf", "png", "tiff")])
      f1_assert(f1_graphic_complete(path, tolower(tools::file_ext(path))), "Incomplete figure: ", path)
    f1_write_csv(f1_file_inventory(files), file.path(out, "Figure1_three_maps_output_inventory.csv"))
    finished <- Sys.time()
    status <- list(status = "COMPLETED", revision = THREE_MAP_REVISION, source_dir = bundle$directory,
      fingerprint = bundle$source_identity$fingerprint, main_domain = config$main_domain,
      started = started, finished = finished, elapsed_minutes = as.numeric(difftime(finished, started, units = "mins")),
      map_layout = maps$layout, shared_scale_limit = maps$scale_limit, legend_breaks = maps$legend_breaks, legend_box = maps$legend_box,
      config = config, warnings = unique(warnings_seen), output_dir = out)
    saveRDS(status, file.path(out, "three_maps_status.rds"))
    message("Three-map figure completed: ", out)
    message("Elapsed ", round(status$elapsed_minutes, 2), " min; warnings ", length(unique(warnings_seen)))
    list(output_dir = out, status = status, source_dir = bundle$directory)
  }, warning = function(w) {
    warnings_seen <<- c(warnings_seen, paste(stage, conditionMessage(w), sep = ": "))
    message("Warning [", stage, "]: ", conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {
    msg <- paste0("Three-map replot failed at ", stage, ": ", conditionMessage(e))
    if (!is.null(out) && dir.exists(out)) {
      f1_write_text(msg, file.path(out, "three_maps_error.txt"))
      saveRDS(list(status = "FAILED", revision = THREE_MAP_REVISION, stage = stage,
        error = conditionMessage(e), output_dir = out), file.path(out, "three_maps_status.rds"))
    }
    f1_stop(msg)
  })
  invisible(result)
}

if (!isTRUE(getOption("ERA5.HalfDegree.no_autorun", FALSE))) {
  f1_assert(HALF_DEGREE_MODE %in% c("population_then_maps", "population", "three_maps"),
    "HALF_DEGREE_MODE must be 'population_then_maps', 'population', or 'three_maps'.")
  f1_assert(is.character(HALF_DEGREE_RUN_DIR) && length(HALF_DEGREE_RUN_DIR) == 1L &&
    !is.na(HALF_DEGREE_RUN_DIR) && nzchar(HALF_DEGREE_RUN_DIR),
    "Set HALF_DEGREE_RUN_DIR to your completed NEW stride=2 warm02 run_... directory before sourcing this script.")
  if (HALF_DEGREE_MODE %in% c("population_then_maps", "population"))
    ERA5_HALF_DEGREE_POPULATION_RESULT <- p1_main(POPULATION_CONFIG)
  if (HALF_DEGREE_MODE %in% c("population_then_maps", "three_maps"))
    ERA5_HALF_DEGREE_THREE_MAPS_RESULT <- f3_main(THREE_MAP_CONFIG)
}
