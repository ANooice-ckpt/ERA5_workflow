# 0.25°_ERA5-warm02_J10TB_native.R -- standalone warm-boundary analysis + daily energy cache.
# Source this file in RStudio; it starts a clean R process automatically.
# Inputs are under E:/ERA5; generated analysis outputs go under J:/202610_Thermal history. Separate native-grid root preserves prior runs.
# KEEP old warm01 outputs: old daily responses were not saved and cannot be recovered from means.
# Full 2010-2019, stride1 on verified native 0.25-degree ERA5 axes: every native grid centre;
# weather is evaluated at its native grid centre, not interpolated or spatially averaged. 2000 paired draws; lag1-3; U=M+N at TSV +1.
# Exact double model calculations; supported daily delta_U cached as checked float32.
# The ALL-DATES cache upper bound is roughly 7.1 TB. Exact support is counted first,
# then the plan is checked against output-drive free space with margin and reserve.
# The 9.0-TB cap is an upper limit, not an estimate of actual free space on J:.
# Keeping the paired daily cache preserves later CDD scenarios.
# No compression of the large binary payload; readback precision QC + per-chunk digests.
# Raw daily temperature is embedded if the complete configured inventory is present.
# If raw T is absent, the warm cache still runs and explicitly reports not_embedded.
# Energy parameters can be decided later using energy_cache_reader.R, without rstan/model prediction.
# 8/12 PSOCK benchmark tiles are retained; measured ETA includes cache I/O. No 8-hour guarantee.
# Modes: global (default), plots (no model calculation), selftest (synthetic data only).
# Example: Rscript ERA5-warm02.R selftest
# Not executed against the user's Windows/RDS/ERA5 files in the delivery environment.
# ---------------------------------------------------------------------------
WARM02_PROGRAM <- quote({
# Full scientific core, spatial projection, diagnostics and figures are embedded below.
MODE <- "global"
ANALYSIS_REVISION <- "ERA5-warm02-2026-10-11-native-quarter-degree-J10TB"
PLOT_REVISION <- "ERA5-warm02-2026-10-11-native-quarter-degree-boundary-maps"
.warm_session_started <- Sys.time()
message("Analysis revision: ", ANALYSIS_REVISION)

# Windows startup: use the library created by Repair_RStan_Windows.R.
# Resolve the ACTUAL R major/minor version; do not hard-code another R version.
# On other systems, or without this optional library, use current .libPaths().
activate_thermal_history_library <- function() {
  if (.Platform$OS.type != "windows") return(invisible(NULL))
  local_app <- Sys.getenv("LOCALAPPDATA")
  if (!nzchar(local_app)) return(invisible(NULL))
  r_series <- paste(R.version$major,
                    strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L],
                    sep = ".")
  repaired_lib <- file.path(local_app, "R", "thermal-history-clean", r_series)
  if (!dir.exists(repaired_lib)) {
    message("No repaired library for R ", r_series,
            "; using the current .libPaths().")
    return(invisible(NULL))
  }
  repaired_lib <- normalizePath(repaired_lib, winslash = "/", mustWork = TRUE)
  core <- c("Rcpp", "RcppEigen", "RcppParallel", "QuickJSR", "StanHeaders", "rstan")
  installed_here <- core[file.exists(file.path(repaired_lib, core, "DESCRIPTION"))]
  already_loaded <- intersect(installed_here, loadedNamespaces())
  conflicting <- already_loaded[vapply(already_loaded, function(pkg) {
    actual <- normalizePath(getNamespaceInfo(getNamespace(pkg), "path"),
                            winslash = "/", mustWork = TRUE)
    expected <- normalizePath(file.path(repaired_lib, pkg),
                              winslash = "/", mustWork = TRUE)
    tolower(actual) != tolower(expected)
  }, logical(1))]
  if (length(conflicting)) {
    stop(paste0(
      "Packages are already loaded from another library: ",
      paste(conflicting, collapse = ", "), ".\n",
      "Restart R without restoring the workspace, then run this script FIRST.\n",
      "The repaired library will be activated automatically: ", repaired_lib
    ), call. = FALSE)
  }
  .libPaths(c(repaired_lib, .libPaths()))
  message("Using repaired package library: ", repaired_lib)
  invisible(repaired_lib)
}
activate_thermal_history_library()

# Edit only these input/output paths first. Source objects remain unchanged.
# Only P50 TSV is evaluated. Warm boundaries do not depend on a fixed SET.
# The global default is a full configured-domain run with 2000 posterior draws.
# Inf uses every available draw in the smaller fit; this can be expensive.
config <- list(
  m1_file = "E:/ERA5/20260815_180420_01-50-lag3/fit_01-50-lag3.rds",
  m0_file = "E:/ERA5/20260831_190333_01-M0-base/fit_01-M0-base.rds",
  data_file = "E:/ERA5/20260815_180420_01-50-lag3/data_01-50-lag3.rds",
  pred_file = "E:/ERA5/20260815_180420_01-50-lag3/pred_01-50-lag3.rds",
  output_dir = "J:/202610_Thermal history/ERA5-warm02_0p25deg_daily_cache",
  fast = list(
    backend = "auto", # exact C++ if a compatible compiler is available; otherwise exact R
    workers = 8L,
    worker_candidates = c(8L, 12L), # measured using complete production tiles
    autotune = TRUE,
    resume = "auto", # resume latest incomplete FAST run; "new" starts a separate run
    runtime_selftest = TRUE,
    target_total_hours = 8 # reporting target only; never truncates computation
  ),
  daily_cache = list(
    enabled = TRUE,
    max_abs_error = 5e-5, # absolute SET-degree tolerance of float32 cache roundtrip
    available_space_gb = 9000, # safety cap only; live free space is measured on J: at preflight
    reserve_space_gb = 500,
    temperature_mode = "auto", # auto embeds only a complete, valid raw-T inventory
    temperature_dir = "E:/ERA5/01_raw_ERA5",
    temperature_checksum_inputs = TRUE,
    verify_payload_on_resume = FALSE # TRUE rereads entire committed binary cache
  ),
  n_draws = Inf,
  seed = 20261002L,
  epred_qc = TRUE,
  strict_priors = TRUE,
  joint_support = TRUE,
  joint_k = 20L,
  joint_quantile = 0.95,
  joint_deduplicate = TRUE,
  scenario_weight_scheme = "observation", # alternative: "building"
  figure_dpi = 300L,
  global = list(
    input_dir = "E:/ERA5/02_Thermal history",
    years = 2010:2019,
    n_draws = 2000L, # complete coefficient rows; retain fixed pairing across all dates/cells
    weight_scheme = "observation", # alternative: "building"
    stride = 1L, # use every native 0.25-degree ERA5 centre; no interpolation
    tile_side = 6L, # at most 6 x 6 sampled grid centers per spatial tile
    day_chunk = 14L,
    lat_range = c(-60, 80),
    land_fraction_file = NULL, # preferred: file on the identical source grid
    land_fraction_variable = "lsm",
    land_fraction_threshold = 0.5,
    min_coverage = 0.95,
    min_supported_fraction = 0, # no whole-cell 80% date filter
    min_supported_days = 30L,
    save_annual_grid = FALSE, # annual area-day posterior summaries are still saved
    support_fraction_sensitivity = c(0.5, 0.8, 0.95),
    delta_threshold = 0.20, # descriptive TSV threshold, NOT validated comfort cutoff
    joint_support = TRUE,
    max_cells = NULL, # full configured spatial domain; no pilot subset
    preprocessing_acknowledged = FALSE # only TRUE after auditing legacy NC provenance
  )
)

# Small diagnostics only: never serialize model objects or weather arrays.
.analysis_progress <- new.env(parent = emptyenv())
.analysis_progress$stage <- "initialization"
.analysis_progress$context <- list()
.set_analysis_stage <- function(stage, ..., announce = FALSE) {
  .analysis_progress$stage <- stage
  .analysis_progress$context <- list(...)
  if (isTRUE(announce)) message("Stage: ", stage)
  invisible(NULL)
}

.analysis_error_details <- function(e, calls) {
  short_call <- function(x) {
    txt <- deparse(x, width.cutoff = 120L)
    paste(utils::head(txt, 8L), collapse = " ")
  }
  list(error_call = short_call(conditionCall(e)),
       traceback = vapply(calls, short_call, character(1)),
       stage = .analysis_progress$stage,
       context = .analysis_progress$context)
}

# ======================================================================
# MODULE: model_core.R
# ======================================================================
# Exact fixed-effect prediction for the supplied brms M1 and M0 fits.
# This file defines functions only. It does not install packages or start a run.
# A/B contract: standardized expected TSV(x) = A * x - B, excluding building
# effects. All nonlinear transforms occur within covariate stratum and draw.
# Between-model uncertainty uses a working product-posterior coupling. Separate
# fits to the same observations do NOT identify a joint posterior for their
# difference. Independent draw pairing is a declared sensitivity convention,
# not a claim that the fits or their estimation errors are independent.

history_columns <- function() c("B21_raw", "TV_lag1_raw", "TV_lag2_raw", "TV_lag3_raw")

require_core_packages <- function(joint = FALSE) {
  # Saved fits contain a stanfit object. Changing a brms backend option does
  # not convert that object; verify its namespace before deserializing it.
  tryCatch(loadNamespace("rstan"), error = function(e) {
    stop(paste0(
      "RStan cannot load, so the saved stanfit objects cannot be used.\n",
      "Original error: ", conditionMessage(e), "\n",
      "On Windows, run Repair_RStan_Windows.R in a clean R session first."
    ), call. = FALSE)
  })
  pkgs <- c("brms", "posterior", "dlnm", "splines")
  if (joint) pkgs <- c(pkgs, "FNN")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages explicitly before running: ",
                            paste(missing, collapse = ", "))
  invisible(pkgs)
}

.with_seed <- function(seed, fun) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  fun()
}

.config_value <- function(config, names, default = NULL) {
  for (n in names) if (!is.null(config[[n]])) return(config[[n]])
  default
}

.history_frame <- function(history) {
  history <- as.data.frame(history)
  wanted <- history_columns()
  aliases <- c(B21 = "B21_raw", TV1 = "TV_lag1_raw", TV2 = "TV_lag2_raw", TV3 = "TV_lag3_raw")
  for (a in names(aliases)) if (!(aliases[[a]] %in% names(history)) && a %in% names(history))
    names(history)[names(history) == a] <- aliases[[a]]
  missing <- setdiff(wanted, names(history))
  if (length(missing)) stop("Missing history columns: ", paste(missing, collapse = ", "))
  out <- history[wanted]
  if (!all(vapply(out, is.numeric, logical(1)))) stop("History columns must be numeric Celsius values.")
  out
}

.training_history <- function(data, pred) {
  if (!("B21_raw" %in% names(data))) stop("Training data has no B21_raw column.")
  if (!all(pred$tv_cols %in% names(data))) stop("Missing training TV columns specified by pred$tv_cols.")
  choose <- function(candidates) candidates[candidates %in% names(data)][1L]
  sources <- c(choose(c("B21_raw_unwinsorized", "B21_unwinsorized", "B21_raw")),
    vapply(seq_len(3L), function(k) choose(c(paste0(pred$tv_cols[k], "_unwinsorized"),
      paste0("TV_lag", k, "_unwinsorized"), paste0("TV_lag", k, "_raw_unwinsorized"),
      pred$tv_cols[k])), character(1)))
  out <- as.data.frame(data[sources], check.names = FALSE)
  names(out) <- history_columns()
  raw_available <- grepl("unwinsorized$", sources)
  if (!all(raw_available)) warning("Unwinsorized training histories are unavailable for: ",
    paste(names(out)[!raw_available], collapse = ", "),
    ". Density support will use already clipped model columns for those dimensions; inspect this limitation.")
  attr(out, "history_source_metadata") <- data.frame(dimension = names(out), source_column = sources,
                                                      unwinsorized = raw_available)
  out
}

.check_prediction_metadata <- function(pred) {
  required <- c("B21_ref", "B21_prediction_P1", "B21_prediction_P99", "b21_clip_to_boundaries",
    "tv_cols", "var_low", "var_high", "var_clip_to_boundaries",
    "b21_neutral_names", "b21_neutral_knots", "b21_neutral_boundary_knots",
    "b21_neutral_intercept", "b21_neutral_ref_basis", "b21_neutral_scale",
    "b21_slope_names", "b21_slope_knots", "b21_slope_boundary_knots",
    "b21_slope_intercept", "b21_slope_ref_basis", "b21_slope_scale",
    "argvar_neutral", "zero_basis_neutral", "cb_neutral_center", "cb_neutral_scale", "cb_neutral_names",
    "argvar_slope", "zero_basis_slope", "cb_slope_center", "cb_slope_scale", "cb_slope_names",
    "arglag", "lag_start", "lag_end", "factor_levels", "model")
  missing <- setdiff(required, names(pred))
  if (length(missing)) stop("Incomplete prediction metadata: ", paste(missing, collapse = ", "))
  if (!isTRUE(pred$b21_clip_to_boundaries) || !isTRUE(pred$var_clip_to_boundaries))
    stop("This implementation expects the supplied P1/P99 winsorized exposure operator.")
  bounds <- as.numeric(c(pred$B21_prediction_P1, pred$B21_prediction_P99, pred$var_low, pred$var_high))
  if (any(!is.finite(bounds)) || bounds[1] >= bounds[2] || bounds[3] >= bounds[4]) stop("Invalid exposure bounds.")
  if (abs(pred$B21_ref - mean(bounds[1:2])) > 1e-8) stop("Unexpected B21 reference.")
  if (pred$lag_start != 1 || pred$lag_end != 3 || length(pred$tv_cols) != 3L) stop("Expected lag1-lag3 history.")
  if (!identical(pred$model$slope_parameterization, "bounded") ||
      !isTRUE(all.equal(as.numeric(c(pred$model$slope_floor, pred$model$slope_range,
         pred$model$neutral_floor, pred$model$neutral_range)), c(.005, .545, 10, 30), tolerance = 1e-12)))
    stop("Unexpected S/N parameterization; check both fits before adapting this code.")
  invisible(TRUE)
}

.b21_basis <- function(x, pred, parameter) {
  x <- pmin(pred$B21_prediction_P99, pmax(pred$B21_prediction_P1, x)) - pred$B21_ref
  prefix <- paste0("b21_", parameter, "_")
  scale <- as.numeric(pred[[paste0(prefix, "scale")]])
  if (any(!is.finite(scale) | scale <= 0)) stop("Invalid saved B21 scaling: do not silently replace it.")
  basis <- as.matrix(splines::ns(x, knots = pred[[paste0(prefix, "knots")]],
    Boundary.knots = pred[[paste0(prefix, "boundary_knots")]],
    intercept = isTRUE(pred[[paste0(prefix, "intercept")]])))
  ref <- as.numeric(pred[[paste0(prefix, "ref_basis")]])
  names <- pred[[paste0(prefix, "names")]]
  if (ncol(basis) != length(ref) || ncol(basis) != length(scale) || ncol(basis) != length(names))
    stop("B21 basis dimension mismatch.")
  basis <- sweep(sweep(basis, 2, ref, "-"), 2, scale, "/")
  colnames(basis) <- names
  basis
}

.tv_basis <- function(v, pred, parameter) {
  v <- as.matrix(v)
  storage.mode(v) <- "double"
  v[v < pred$var_low] <- pred$var_low
  v[v > pred$var_high] <- pred$var_high
  argvar <- pred[[paste0("argvar_", parameter)]]
  # Explicit boundary/knots are saved; the template also supports a one-row input.
  template <- sort(unique(c(pred$var_low, as.numeric(argvar$knots), 0, pred$var_high)))
  template <- template[is.finite(template)]
  augmented <- rbind(v, matrix(template, nrow = length(template), ncol = ncol(v)))
  cb <- as.matrix(dlnm::crossbasis(augmented, lag = c(pred$lag_start, pred$lag_end),
    argvar = argvar, arglag = pred$arglag))[seq_len(nrow(v)), , drop = FALSE]
  zero <- as.numeric(pred[[paste0("zero_basis_", parameter)]])
  center <- as.numeric(pred[[paste0("cb_", parameter, "_center")]])
  scale <- as.numeric(pred[[paste0("cb_", parameter, "_scale")]])
  names <- pred[[paste0("cb_", parameter, "_names")]]
  if (any(!is.finite(scale) | scale <= 0)) stop("Invalid saved TV scaling.")
  if (any(c(length(zero), length(center), length(scale), length(names)) != ncol(cb)))
    stop("TV cross-basis dimension mismatch.")
  cb <- sweep(sweep(sweep(cb, 2, zero, "-"), 2, center, "-"), 2, scale, "/")
  colnames(cb) <- names
  cb
}

make_history_design <- function(history, pred) {
  h <- .history_frame(history)
  if (!nrow(h) || any(!is.finite(as.matrix(h)))) stop("Prediction history must be nonempty and finite; filter with evaluate_support first.")
  vn <- .tv_basis(as.matrix(h[2:4]), pred, "neutral")
  vs <- .tv_basis(as.matrix(h[2:4]), pred, "slope")
  as.data.frame(cbind(.b21_basis(h$B21_raw, pred, "neutral"),
    .b21_basis(h$B21_raw, pred, "slope"), vn, vs), check.names = FALSE)
}

.strip_formula <- function(f) gsub("[[:space:]]", "", paste(deparse(f, width.cutoff = 500L), collapse = ""))

.parameter_formula <- function(fit, name) {
  f <- fit$formula$pforms[[name]]
  if (is.null(f)) stop("Fit lacks nonlinear parameter formula: ", name)
  if (inherits(f, "brmsformula") || (is.list(f) && !is.null(f$formula))) f <- f$formula
  f
}

.check_models <- function(m1, m0, data, pred, strict_priors = TRUE) {
  if (!inherits(m1, "brmsfit") || !inherits(m0, "brmsfit")) stop("Both model files must contain brmsfit objects.")
  expected <- "thermal_sensation~(0.005+0.545*inv_logit(logslope))*(SET_raw-(10+30*inv_logit(etaneutral)))+etabuilding"
  main1 <- .strip_formula(m1$formula$formula)
  main0 <- .strip_formula(m0$formula$formula)
  if (!identical(main1, main0) || !identical(main1, expected)) stop("Main outcome formulas differ or do not match the supplied bounded model.")
  if (m1$family$family != "gaussian" || m0$family$family != "gaussian" ||
      m1$family$link != "identity" || m0$family$link != "identity") stop("Expected Gaussian identity-response models.")
  covars <- c("operation_mode_fe", "weekend_num", "aligned_season_fe")
  for (parameter in c("etaneutral", "logslope")) {
    history <- if (parameter == "etaneutral") c(pred$b21_neutral_names, pred$cb_neutral_names) else
      c(pred$b21_slope_names, pred$cb_slope_names)
    labels1 <- attr(stats::terms(.parameter_formula(m1, parameter)), "term.labels")
    labels0 <- attr(stats::terms(.parameter_formula(m0, parameter)), "term.labels")
    if (!setequal(labels1, c(history, covars)) || !setequal(labels0, covars))
      stop("Unexpected predictor terms in ", parameter, "; M0 must omit only the history basis blocks.")
    if (attr(stats::terms(.parameter_formula(m1, parameter)), "intercept") != 1L ||
        attr(stats::terms(.parameter_formula(m0, parameter)), "intercept") != 1L)
      stop("Both nonlinear predictors must include their supplied intercept.")
  }
  group1 <- .strip_formula(.parameter_formula(m1, "etabuilding"))
  group0 <- .strip_formula(.parameter_formula(m0, "etabuilding"))
  if (!identical(group1, group0) || !grepl("\\(1\\|building_id\\)", group1)) stop("Building effect specifications differ.")
  shared <- c("thermal_sensation", "SET_raw", "building_id", covars)
  for (col in shared) {
    if (!all(vapply(list(data, m1$data, m0$data), function(d) col %in% names(d), logical(1))))
      stop("Missing shared model data column: ", col)
    a <- data[[col]]; b <- m1$data[[col]]; c <- m0$data[[col]]
    same <- function(x, y) if (is.numeric(x) && is.numeric(y))
      isTRUE(all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE)) else identical(as.character(x), as.character(y))
    if (!same(a, b) || !same(a, c)) stop("Shared training data differs row by row: ", col)
  }
  for (col in c("operation_mode_fe", "aligned_season_fe")) {
    lev <- pred$factor_levels[[col]]
    for (d in list(data, m1$data, m0$data)) if (!identical(levels(d[[col]]), lev))
      stop("Factor-level order differs: ", col)
  }
  canonical_priors <- function(fit) {
    p <- as.data.frame(brms::prior_summary(fit), stringsAsFactors = FALSE)
    cols <- intersect(c("class", "coef", "group", "resp", "dpar", "nlpar", "lb", "ub", "prior"), names(p))
    if ("coef" %in% names(p)) p <- p[!grepl("^(B21[NS]|CB[NS])_", p$coef), , drop = FALSE]
    p <- p[cols]
    text <- apply(p, 1, function(x) paste(ifelse(is.na(x), "", x), collapse = "|"))
    sort(unique(text))
  }
  priors_match <- identical(canonical_priors(m1), canonical_priors(m0))
  if (!priors_match && strict_priors) stop("Shared prior tables differ; examine the fits and declare the comparison before continuing.")
  if (!priors_match) warning("Shared prior tables differ; comparison contains a prior change.")
  list(shared_columns = shared, priors_match = priors_match, formula_m1 = main1,
       formula_m0 = main0, n_rows = nrow(data), n_buildings = length(unique(data$building_id)))
}

.coefficient_columns <- function(draws, nlpar, terms) {
  terms[terms == "(Intercept)"] <- "Intercept"
  target <- paste0("b_", nlpar, "_", terms)
  available <- colnames(draws)
  canon <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  idx <- vapply(target, function(t) {
    exact <- which(available == t)
    if (length(exact) == 1L) return(exact)
    approximate <- which(canon(available) == canon(t))
    if (length(approximate) != 1L) stop("Cannot uniquely resolve coefficient: ", t)
    approximate
  }, integer(1))
  draws[, idx, drop = FALSE]
}

.strata_weights <- function(data, pred) {
  op <- pred$factor_levels$operation_mode_fe
  season <- pred$factor_levels$aligned_season_fe
  weekend <- sort(as.numeric(pred$factor_levels$weekend_num))
  strata <- expand.grid(operation_mode_fe = op, weekend_num = weekend,
    aligned_season_fe = season, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key <- function(d) paste(as.character(d$operation_mode_fe), as.numeric(d$weekend_num),
                          as.character(d$aligned_season_fe), sep = "||")
  training_key <- key(data); strata_key <- key(strata)
  counts <- table(training_key)
  strata$n_rows <- as.integer(counts[strata_key]); strata$n_rows[is.na(strata$n_rows)] <- 0L
  if (sum(strata$n_rows) != nrow(data)) stop("Strata do not cover the training sample.")
  strata$observation <- strata$n_rows / nrow(data)
  building_n <- table(as.character(data$building_id))
  row_building_w <- 1 / (length(building_n) * as.numeric(building_n[as.character(data$building_id)]))
  building_w <- tapply(row_building_w, training_key, sum)
  strata$building <- as.numeric(building_w[strata_key]); strata$building[is.na(strata$building)] <- 0
  strata$row_equal_weight <- strata$observation
  strata$building_equal_weight <- strata$building
  strata$operation_mode_fe <- factor(strata$operation_mode_fe, levels = op)
  strata$aligned_season_fe <- factor(strata$aligned_season_fe, levels = season)
  strata$weekend_num <- as.numeric(strata$weekend_num)
  if (abs(sum(strata$observation) - 1) > 1e-12 || abs(sum(strata$building) - 1) > 1e-12) stop("Invalid standardization weights.")
  strata
}

.covariate_matrix <- function(data) {
  # Supplied fits use treatment-coded factors. Avoid dependence on the user's
  # global options("contrasts"); official_epred QC verifies this assumption.
  stats::model.matrix(~ operation_mode_fe + weekend_num + aligned_season_fe, data,
    contrasts.arg = list(operation_mode_fe = "contr.treatment", aligned_season_fe = "contr.treatment"))
}

.basis_qc <- function(data, pred, seed, n = 500L) {
  h <- .training_history(data, pred)
  idx <- .with_seed(seed, function() sample.int(nrow(data), min(n, nrow(data))))
  tails <- unique(unlist(lapply(h, function(x) c(which.min(x), which.max(x)))))
  idx <- sort(unique(c(idx, tails)))
  rebuilt <- make_history_design(h[idx, , drop = FALSE], pred)
  blocks <- list(B21N = pred$b21_neutral_names, B21S = pred$b21_slope_names,
                 CBN = pred$cb_neutral_names, CBS = pred$cb_slope_names)
  result <- do.call(rbind, lapply(names(blocks), function(name) {
    columns <- blocks[[name]]
    if (!all(columns %in% names(data))) stop("Training design lacks ", name)
    expected <- as.matrix(data[idx, columns, drop = FALSE])
    actual <- as.matrix(rebuilt[columns])
    err <- max(abs(actual - expected))
    if (!is.finite(err) || err > 1e-8) stop(name, " reconstruction failed: ", err)
    data.frame(block = name, n_rows = length(idx), max_abs_error = err, tolerance = 1e-8)
  }))
  result
}

.make_joint_reference <- function(training, k = 20L, quantile = .95, deduplicate = TRUE) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Install FNN explicitly to enable joint-density diagnostics.")
  x <- as.matrix(.history_frame(training))
  x <- x[rowSums(!is.finite(x)) == 0L, , drop = FALSE]
  n_original <- nrow(x)
  if (deduplicate) x <- unique(x)
  if (nrow(x) < 2L || length(k) != 1L || !is.finite(k) || k < 1L)
    stop("At least two finite, distinct histories and a positive k are required for kNN support.")
  if (length(quantile) != 1L || !is.finite(quantile) || quantile <= 0 || quantile >= 1)
    stop("Joint-support threshold quantile must be strictly between 0 and 1.")
  requested_k <- as.integer(k)
  k <- min(requested_k, nrow(x) - 1L)
  center <- colMeans(x); scale <- apply(x, 2, stats::sd)
  if (any(!is.finite(scale) | scale <= 0)) stop("Degenerate history dimension in support reference.")
  z <- sweep(sweep(x, 2, center, "-"), 2, scale, "/")
  distance <- FNN::get.knn(z, k = as.integer(k), algorithm = "kd_tree")$nn.dist[, k]
  threshold <- as.numeric(stats::quantile(distance, quantile, names = FALSE, type = 8))
  if (!is.finite(threshold) || threshold <= 0) stop("Degenerate kNN support threshold; inspect duplicate histories.")
  list(z = z, center = center, scale = scale, k = as.integer(k), threshold = threshold,
       requested_k = requested_k, threshold_quantile = quantile, deduplicated = deduplicate, n_original = n_original,
       n_reference = nrow(z), training_source = attr(training, "history_source_metadata"),
       method = "standardized raw four-dimensional history; training LOO kth-neighbor distance")
}

history_support <- function(history, training, pred, joint = FALSE, k = 20L,
                            quantile = .95, deduplicate = TRUE, joint_reference = NULL) {
  h <- .history_frame(history); x <- as.matrix(h)
  finite <- rowSums(!is.finite(x)) == 0L
  lo <- c(pred$B21_prediction_P1, rep(pred$var_low, 3L))
  hi <- c(pred$B21_prediction_P99, rep(pred$var_high, 3L))
  # ERA5 was stored as float32. Tolerance is numerical only, not extra support.
  tol <- 1e-5
  outside <- sweep(x, 2, lo - tol, "<") | sweep(x, 2, hi + tol, ">")
  outside[!is.finite(x)] <- NA
  any_outside <- rowSums(outside, na.rm = TRUE) > 0
  result <- data.frame(finite = finite, marginal_inside = finite & !any_outside,
    B21_outside = outside[, 1], TV_lag1_outside = outside[, 2],
    TV_lag2_outside = outside[, 3], TV_lag3_outside = outside[, 4])
  result$joint_inside <- NA
  result$joint_distance <- NA_real_
  if (joint) {
    ref <- joint_reference
    if (is.null(ref)) ref <- .make_joint_reference(training, k, quantile, deduplicate)
    if (any(finite)) {
      query <- sweep(sweep(x[finite, , drop = FALSE], 2, ref$center, "-"), 2, ref$scale, "/")
      # Queries identical to a training record need one extra neighbor so their
      # zero-distance self match is excluded, matching the training LOO rule.
      nn <- FNN::get.knnx(ref$z, query, k = ref$k + 1L, algorithm = "kd_tree")$nn.dist
      same <- nn[, 1] <= 1e-12
      distance <- nn[cbind(seq_len(nrow(nn)), ref$k + as.integer(same))]
      result$joint_distance[finite] <- distance
      result$joint_inside[finite] <- distance <= ref$threshold + 1e-12
    }
    attr(result, "joint_reference_metadata") <- ref[setdiff(names(ref), "z")]
  }
  attr(result, "support_interpretation") <- paste("Marginal bounds and kNN are exposure-density diagnostics,",
    "not external validation or guarantees of transportability; they do not include conditional SET/covariate support.")
  result
}

evaluate_support <- function(analysis, history, joint = !is.null(analysis$joint_reference)) {
  history_support(history, analysis$training_history, analysis$pred, joint = joint,
    joint_reference = analysis$joint_reference)
}

.draw_indices <- function(analysis, draw_indices) {
  if (is.null(draw_indices)) return(seq_len(nrow(analysis$draws_m1)))
  if (!length(draw_indices) || any(!is.finite(draw_indices)) || any(draw_indices < 1L) ||
      any(draw_indices > nrow(analysis$draws_m1)) || any(draw_indices != as.integer(draw_indices)))
    stop("Invalid paired-draw indices.")
  as.integer(draw_indices)
}

predict_ab <- function(analysis, history, weight_scheme = "observation", draw_indices = NULL) {
  h <- .history_frame(history)
  design <- make_history_design(h, analysis$pred)
  idx <- .draw_indices(analysis, draw_indices)
  if (is.numeric(weight_scheme)) w <- weight_scheme else {
    aliases <- c(observation = "observation", building = "building", row_equal_weight = "observation",
                 building_equal_weight = "building")
    if (length(weight_scheme) != 1L || !(weight_scheme %in% names(aliases))) stop("Unknown weight scheme.")
    w <- analysis$weights[[aliases[[weight_scheme]]]]
  }
  if (length(w) != nrow(analysis$strata) || any(!is.finite(w) | w < 0) || sum(w) <= 0)
    stop("Weights must be finite, nonnegative, and match standardization strata.")
  w <- w / sum(w)
  p <- analysis$pred
  dn <- c(p$b21_neutral_names, p$cb_neutral_names)
  ds <- c(p$b21_slope_names, p$cb_slope_names)
  beta1 <- analysis$draws_m1[idx, , drop = FALSE]
  beta0 <- analysis$draws_m0[idx, , drop = FALSE]
  zN <- as.matrix(design[dn]) %*% t(.coefficient_columns(beta1, "etaneutral", dn))
  zS <- as.matrix(design[ds]) %*% t(.coefficient_columns(beta1, "logslope", ds))
  cov1N <- analysis$Xcov %*% t(.coefficient_columns(beta1, "etaneutral", colnames(analysis$Xcov)))
  cov1S <- analysis$Xcov %*% t(.coefficient_columns(beta1, "logslope", colnames(analysis$Xcov)))
  cov0N <- analysis$Xcov %*% t(.coefficient_columns(beta0, "etaneutral", colnames(analysis$Xcov)))
  cov0S <- analysis$Xcov %*% t(.coefficient_columns(beta0, "logslope", colnames(analysis$Xcov)))
  A1 <- B1 <- matrix(0, nrow(h), length(idx))
  for (g in which(w > 0)) {
    n1 <- p$model$neutral_floor + p$model$neutral_range * plogis(sweep(zN, 2, cov1N[g, ], "+"))
    s1 <- p$model$slope_floor + p$model$slope_range * plogis(sweep(zS, 2, cov1S[g, ], "+"))
    A1 <- A1 + w[g] * s1
    B1 <- B1 + w[g] * s1 * n1
  }
  n0 <- p$model$neutral_floor + p$model$neutral_range * plogis(cov0N)
  s0 <- p$model$slope_floor + p$model$slope_range * plogis(cov0S)
  A0 <- as.numeric(crossprod(w, s0))
  B0 <- as.numeric(crossprod(w, s0 * n0))
  if (any(!is.finite(c(A1, B1, A0, B0)))) stop("Nonfinite standardized prediction.")
  list(A1 = A1, B1 = B1, A0 = A0, B0 = B0,
       paired_draw_indices = idx, draw_ids_m1 = analysis$draw_ids_m1[idx],
       draw_ids_m0 = analysis$draw_ids_m0[idx], weight_scheme = weight_scheme,
       estimand = "study-sample-standardized expected TSV; building effects excluded")
}

predict_draws <- function(analysis, history, set_values = analysis$fixed_set,
                          weight_scheme = "observation", draw_indices = NULL) {
  set_values <- as.numeric(set_values)
  if (!length(set_values) || any(!is.finite(set_values))) stop("SET scenarios must be finite.")
  ab <- predict_ab(analysis, history, weight_scheme, draw_indices)
  m1 <- delta <- array(NA_real_, c(nrow(ab$A1), ncol(ab$A1), length(set_values)))
  m0 <- outer(ab$A0, set_values) - ab$B0
  for (s in seq_along(set_values)) {
    mu <- ab$A1 * set_values[s] - ab$B1
    m1[, , s] <- mu
    delta[, , s] <- sweep(mu, 2, m0[, s], "-")
  }
  list(m1 = m1, m0 = m0, delta = delta, set_values = set_values,
       dimensions = "m1/delta: history x paired draw x SET; m0: paired draw x SET", ab = ab)
}

.plain_draw_matrix <- function(x) {
  if (length(dim(x)) != 2L || !is.numeric(x)) stop("Expected a numeric draw matrix.")
  # Build a base matrix before any random row selection. as.matrix() alone
  # need not remove the posterior draws_matrix class and its S3 methods.
  matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
}

.official_epred_by_id <- function(fit, newdata, draw_ids, epred_fun = NULL) {
  if (!length(draw_ids) || anyNA(draw_ids) || any(!is.finite(draw_ids)) ||
      any(draw_ids < 1L | draw_ids != as.integer(draw_ids)) || anyDuplicated(draw_ids))
    stop("Official prediction requires distinct positive integer draw IDs.")
  draw_ids <- as.integer(draw_ids)
  if (is.null(epred_fun)) epred_fun <- brms::posterior_epred
  # brms/posterior may repair and reorder draw indices during subsequent
  # variable selection. Request ascending IDs, then explicitly restore the
  # caller's original pairing order. Never sort M1/M0 pairing globally.
  sorted_ids <- sort(draw_ids)
  out <- epred_fun(fit, newdata = newdata, re_formula = NA, draw_ids = sorted_ids)
  if (length(dim(out)) != 2L || !identical(as.integer(dim(out)),
      c(length(sorted_ids), nrow(newdata))))
    stop("Unexpected official prediction dimensions; expected draws x observations.")
  out <- .plain_draw_matrix(out)[match(draw_ids, sorted_ids), , drop = FALSE]
  rownames(out) <- as.character(draw_ids)
  out
}

.official_epred_qc <- function(analysis, n = 12L, n_draws = 20L, seed = 8103L,
                               output_dir = NULL) {
  result <- list()
  record_check <- function(model, stage, manual, official, draw_ids, history) {
    if (!identical(dim(manual), dim(official))) stop("QC prediction dimensions differ.")
    err <- max(abs(manual - official))
    passed <- is.finite(err) && err <= 1e-8
    result[[length(result) + 1L]] <<- data.frame(model = model, stage = stage,
      n_rows = ncol(manual), n_draws = nrow(manual), max_abs_error = err,
      tolerance = 1e-8, passed = passed, history = history,
      draw_alignment = "ascending API request, then explicit match to original draw IDs")
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(do.call(rbind, result),
        file.path(output_dir, "posterior_epred_QC.csv"), row.names = FALSE)
      if (!passed) {
        cells <- expand.grid(draw_index = seq_len(nrow(manual)),
                             observation_index = seq_len(ncol(manual)))
        cells$original_draw_id <- draw_ids[cells$draw_index]
        cells$manual <- as.numeric(manual)
        cells$official <- as.numeric(official)
        cells$difference <- cells$manual - cells$official
        utils::write.csv(cells, file.path(output_dir,
          paste0("posterior_epred_FAILED_", model, "_", stage, ".csv")), row.names = FALSE)
      }
    }
    if (!passed) stop(model, " ", stage, " official posterior_epred check failed: ",
      format(err, digits = 16), "; original draw IDs have been aligned.",
      if (!is.null(output_dir)) paste0(" Inspect posterior_epred_QC.csv and posterior_epred_FAILED_* in ", output_dir),
      call. = FALSE)
    message("PASS: ", model, " ", stage, " posterior_epred QC; max error = ",
            format(err, digits = 6))
    invisible(err)
  }
  p <- analysis$pred
  idx <- .with_seed(seed, function() sample.int(nrow(analysis$data), min(n, nrow(analysis$data))))
  nd <- analysis$data[idx, , drop = FALSE]
  h <- .with_seed(seed + 1L, function() data.frame(
    B21_raw = stats::runif(nrow(nd), p$B21_prediction_P1 - 2, p$B21_prediction_P99 + 2),
    TV_lag1_raw = stats::runif(nrow(nd), p$var_low - 2, p$var_high + 2),
    TV_lag2_raw = stats::runif(nrow(nd), p$var_low - 2, p$var_high + 2),
    TV_lag3_raw = stats::runif(nrow(nd), p$var_low - 2, p$var_high + 2)))
  basis <- make_history_design(h, p)
  for (name in names(basis)) nd[[name]] <- basis[[name]]
  nd$B21_raw <- h$B21_raw
  for (k in 1:3) nd[[p$tv_cols[k]]] <- h[[k + 1L]]
  nd$operation_mode_fe <- factor(nd$operation_mode_fe, levels = p$factor_levels$operation_mode_fe)
  nd$aligned_season_fe <- factor(nd$aligned_season_fe, levels = p$factor_levels$aligned_season_fe)
  X <- .covariate_matrix(nd)
  qi <- seq_len(min(n_draws, nrow(analysis$draws_m1)))
  for (model in c("m1", "m0")) {
    draws <- analysis[[paste0("draws_", model)]][qi, , drop = FALSE]
    etaN <- X %*% t(.coefficient_columns(draws, "etaneutral", colnames(X)))
    etaS <- X %*% t(.coefficient_columns(draws, "logslope", colnames(X)))
    if (model == "m1") {
      dn <- c(p$b21_neutral_names, p$cb_neutral_names)
      ds <- c(p$b21_slope_names, p$cb_slope_names)
      etaN <- etaN + as.matrix(basis[dn]) %*% t(.coefficient_columns(draws, "etaneutral", dn))
      etaS <- etaS + as.matrix(basis[ds]) %*% t(.coefficient_columns(draws, "logslope", ds))
    }
    neutral <- p$model$neutral_floor + p$model$neutral_range * plogis(etaN)
    slope <- p$model$slope_floor + p$model$slope_range * plogis(etaS)
    manual <- slope * (nd$SET_raw - neutral)
    draw_ids <- analysis[[paste0("draw_ids_", model)]][qi]
    official <- .official_epred_by_id(analysis[[model]], nd, draw_ids)
    # A single-ID API call independently anchors a row in the batched output.
    # Pick the largest selected original ID, often far from the first paired row.
    anchor <- which.max(draw_ids)
    single <- as.matrix(brms::posterior_epred(analysis[[model]], newdata = nd,
      re_formula = NA, draw_ids = draw_ids[anchor]))
    record_check(model, "draw_id_anchor", official[anchor, , drop = FALSE], single,
      draw_ids[anchor], "batched prediction versus a separate single-draw API call")
    record_check(model, "rowwise", t(manual), official, draw_ids,
      "independent random history, including P1/P99 tails")
  }
  # Also validate predict_ab itself against official predictions aggregated over
  # all strata. This independently checks dimensions, weights, and S*N ordering.
  h_small <- h[seq_len(min(2L, nrow(h))), , drop = FALSE]
  n_strata <- nrow(analysis$strata)
  gi <- rep(seq_len(n_strata), times = nrow(h_small))
  hi <- rep(seq_len(nrow(h_small)), each = n_strata)
  nd_standardized <- analysis$data[rep(1L, length(gi)), , drop = FALSE]
  for (name in c("operation_mode_fe", "weekend_num", "aligned_season_fe"))
    nd_standardized[[name]] <- analysis$strata[[name]][gi]
  nd_standardized$SET_raw <- rep(analysis$fixed_set[["P50"]], length(gi))
  basis_standardized <- make_history_design(h_small[hi, , drop = FALSE], p)
  for (name in names(basis_standardized)) nd_standardized[[name]] <- basis_standardized[[name]]
  nd_standardized$B21_raw <- h_small$B21_raw[hi]
  for (k in 1:3) nd_standardized[[p$tv_cols[k]]] <- h_small[[k + 1L]][hi]
  ab <- predict_ab(analysis, h_small, draw_indices = qi)
  for (model in c("m1", "m0")) {
    draw_ids <- analysis[[paste0("draw_ids_", model)]][qi]
    official <- .official_epred_by_id(analysis[[model]], nd_standardized, draw_ids)
    averaged <- vapply(seq_len(nrow(h_small)), function(i) {
      columns <- which(hi == i)
      as.numeric(official[, columns, drop = FALSE] %*% analysis$strata$observation)
    }, numeric(length(qi)))
    averaged <- matrix(averaged, nrow = length(qi), ncol = nrow(h_small))
    manual <- if (model == "m1") t(ab$A1 * analysis$fixed_set[["P50"]] - ab$B1) else
      matrix(ab$A0 * analysis$fixed_set[["P50"]] - ab$B0, nrow = length(qi), ncol = nrow(h_small))
    record_check(model, "standardized", manual, averaged, draw_ids,
      "random history x all standardization strata, aggregated to history rows")
  }
  do.call(rbind, result)
}

load_analysis <- function(config) {
  joint <- isTRUE(.config_value(config, c("joint_support"), TRUE))
  require_core_packages(joint)
  run_dir <- .config_value(config, c("run_dir", "model_dir"), ".")
  paths <- c(m1 = .config_value(config, c("m1_file", "fit_m1", "fit_m1_file"), file.path(run_dir, "fit_01-50-lag3.rds")),
    m0 = .config_value(config, c("m0_file", "fit_m0", "fit_m0_file"), file.path(run_dir, "fit_01-M0-base.rds")),
    data = .config_value(config, c("data_file", "analysis_data_file"), file.path(run_dir, "data_01-50-lag3.rds")),
    pred = .config_value(config, c("pred_file", "prediction_objects_file"), file.path(run_dir, "pred_01-50-lag3.rds")))
  if (any(!file.exists(paths))) stop("Missing model input: ", paste(paths[!file.exists(paths)], collapse = ", "))
  m1 <- readRDS(paths["m1"]); m0 <- readRDS(paths["m0"])
  data <- readRDS(paths["data"]); pred <- readRDS(paths["pred"])
  .check_prediction_metadata(pred)
  checks <- .check_models(m1, m0, data, pred,
                         strict_priors = isTRUE(.config_value(config, c("strict_priors"), TRUE)))
  seed <- as.integer(.config_value(config, c("seed"), 20261002L))
  # Keep complete joint rows for the variables used in fixed-effect prediction;
  # excluding unused building coefficients avoids a second large matrix copy.
  variables <- c("^b_(etaneutral|logslope)_", "^sigma$")
  draws1 <- .plain_draw_matrix(posterior::as_draws_matrix(m1, variable = variables, regex = TRUE))
  draws0 <- .plain_draw_matrix(posterior::as_draws_matrix(m0, variable = variables, regex = TRUE))
  requested <- .config_value(config, c("n_draws"), 1000L)
  max_draws <- min(nrow(draws1), nrow(draws0))
  if (length(requested) != 1L || is.na(requested) || (!is.infinite(requested) &&
      (!is.numeric(requested) || requested != as.integer(requested))))
    stop("n_draws must be a single positive integer or Inf.")
  n_draws <- if (length(requested) == 1L && is.infinite(requested) && requested > 0) max_draws else as.integer(requested)
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 1L || n_draws > max_draws)
    stop("n_draws must be between 1 and the smaller fit draw count, or Inf.")
  ids <- .with_seed(seed, function() list(m1 = sample.int(nrow(draws1), n_draws),
                                         m0 = sample.int(nrow(draws0), n_draws)))
  select <- function(x, id) {
    keep <- grepl("^b_(etaneutral|logslope)_", colnames(x)) | colnames(x) == "sigma"
    x[id, keep, drop = FALSE]
  }
  strata <- .strata_weights(data, pred)
  Xcov <- .covariate_matrix(strata)
  fixed_set <- stats::quantile(data$SET_raw, .5, na.rm = TRUE, names = FALSE, type = 7)
  names(fixed_set) <- "P50"
  training_history <- .training_history(data, pred)
  analysis <- list(m1 = m1, m0 = m0, data = data, pred = pred, strata = strata, weights = strata,
    Xcov = Xcov, draws_m1 = select(draws1, ids$m1), draws_m0 = select(draws0, ids$m0),
    draw_ids_m1 = ids$m1, draw_ids_m0 = ids$m0, draw_ids = ids$m1,
    fixed_set = fixed_set, training_history = training_history, model_checks = checks,
    history_source_metadata = attr(training_history, "history_source_metadata"),
    input_paths = paths, input_md5 = unname(tools::md5sum(paths)),
    draw_metadata = list(seed = seed, n_draws = n_draws, available_m1 = nrow(draws1), available_m0 = nrow(draws0),
      pairing = "independent random pairing without replacement within each fit",
      interpretation = paste("working product-posterior coupling; marginal within-fit dependence retained;",
        "same-data cross-fit dependence is not identified; difference intervals are coupling-conditional")))
  analysis$basis_qc <- .basis_qc(data, pred, seed + 11L)
  analysis$joint_reference <- if (joint) .make_joint_reference(training_history,
    k = as.integer(.config_value(config, c("joint_k"), 20L)),
    quantile = .config_value(config, c("joint_quantile"), .95),
    deduplicate = isTRUE(.config_value(config, c("joint_deduplicate"), TRUE))) else NULL
  if (isTRUE(.config_value(config, c("epred_qc"), TRUE))) analysis$epred_qc <- .official_epred_qc(
    analysis, seed = seed + 23L, output_dir = .config_value(config, c("output_dir"), NULL))
  pkgs <- c("R", "brms", "posterior", "dlnm", "splines", if (joint) "FNN")
  analysis$package_manifest <- data.frame(package = pkgs, version = vapply(pkgs, function(p)
    if (p == "R") as.character(getRversion()) else as.character(utils::packageVersion(p)), character(1)))
  analysis
}

# Compatibility alias for runners that prefer a bundle-oriented name.
prepare_model_bundle <- load_analysis

# ======================================================================

era5_required_variables <- c("B21_raw", "TV_lag1_raw", "TV_lag2_raw", "TV_lag3_raw")

global_defaults <- function(config) {
  g <- config$global
  if (is.null(g)) g <- list()
  defaults <- list(
    input_dir = NULL, years = 2010:2019, spatial_stride = 1L, n_draws = 2000L,
    tile_side = 6L, day_chunk = 14L, lat_range = c(-60, 80),
    land_fraction_file = NULL, land_fraction_variable = "lsm",
    land_fraction_threshold = 0.5, joint_support = TRUE,
    min_date_coverage = 0.95, min_supported_days = 30L,
    min_supported_fraction = 0,
    support_fraction_sensitivity = c(0.5, 0.8, 0.95), max_cells = NULL,
    delta_threshold = 0.2, weight_scheme = "observation",
    allow_kelvin_background = FALSE, checksum_inputs = TRUE,
    save_annual_grid = FALSE, file_template = "thermal_history_%04d.nc",
    preprocessing_acknowledged = FALSE
  )
  for (nm in names(defaults)) if (is.null(g[[nm]])) g[[nm]] <- defaults[[nm]]
  # Alias accepted for configuration files written before the interface settled.
  if (!is.null(g$stride)) g$spatial_stride <- g$stride
  if (!is.null(g$time_chunk_days)) g$day_chunk <- g$time_chunk_days
  if (!is.null(g$min_coverage)) g$min_date_coverage <- g$min_coverage
  if (!is.null(g$era5_dir)) g$input_dir <- g$era5_dir
  if (is.null(g$input_dir)) stop("config$global$input_dir is required.")
  g$years <- as.integer(g$years)
  if (!length(g$years) || anyDuplicated(g$years) || any(diff(g$years) != 1L))
    stop("Global years must be unique, consecutive, and increasing.")
  for (nm in c("spatial_stride", "tile_side", "day_chunk", "min_supported_days")) {
    if (length(g[[nm]]) != 1L || !is.finite(g[[nm]]) ||
        g[[nm]] < 1 || g[[nm]] != as.integer(g[[nm]])) stop("Invalid global setting: ", nm)
    g[[nm]] <- as.integer(g[[nm]])
  }
  if (length(g$lat_range) != 2L || any(!is.finite(g$lat_range)) ||
      g$lat_range[1] >= g$lat_range[2] || any(abs(g$lat_range) > 90))
    stop("lat_range must be an increasing pair within [-90,90].")
  if (!is.finite(g$delta_threshold) || g$delta_threshold <= 0)
    stop("delta_threshold is a positive illustrative TSV threshold, not a validated cutoff.")
  if (!is.finite(g$min_date_coverage) || g$min_date_coverage <= 0 || g$min_date_coverage > 1)
    stop("min_date_coverage must lie in (0,1].")
  if (length(g$min_supported_fraction) != 1L || !is.finite(g$min_supported_fraction) ||
      g$min_supported_fraction < 0 || g$min_supported_fraction > 1)
    stop("min_supported_fraction must lie in [0,1].")
  g$support_fraction_sensitivity <- sort(unique(c(g$support_fraction_sensitivity, g$min_supported_fraction)))
  if (any(!is.finite(g$support_fraction_sensitivity)) || any(g$support_fraction_sensitivity < 0 | g$support_fraction_sensitivity > 1))
    stop("support_fraction_sensitivity must contain fractions in [0,1].")
  if (!is.null(g$max_cells) && (length(g$max_cells) != 1L || !is.finite(g$max_cells) ||
      g$max_cells < 1 || g$max_cells != as.integer(g$max_cells))) stop("max_cells must be NULL or a positive integer.")
  if (length(g$n_draws) != 1L || is.na(g$n_draws) || g$n_draws < 1 ||
      (!is.infinite(g$n_draws) && g$n_draws != as.integer(g$n_draws)))
    stop("global$n_draws must be a positive integer or Inf.")
  if (length(g$land_fraction_threshold) != 1L || !is.finite(g$land_fraction_threshold) ||
      g$land_fraction_threshold < 0 || g$land_fraction_threshold > 1)
    stop("land_fraction_threshold must lie in [0,1].")
  g
}

normal_axis <- function(x) {
  x <- tolower(x)
  x[x %in% c("longitude", "lon")] <- "lon"
  x[x %in% c("latitude", "lat")] <- "lat"
  x[x %in% c("valid_time", "time", "date")] <- "time"
  x
}

nc_attribute <- function(nc, varid, name) {
  a <- ncdf4::ncatt_get(nc, varid, name)
  if (isTRUE(a$hasatt)) a$value else NA_character_
}

nc_axis <- function(nc, axis) {
  n <- names(nc$dim)
  hit <- which(normal_axis(n) == axis)
  if (length(hit) != 1L) stop("Expected exactly one ", axis, " coordinate in NetCDF.")
  d <- nc$dim[[hit]]
  if (!is.numeric(d$vals)) stop("Non-numeric NetCDF coordinate: ", n[hit])
  # ncdf4 can retain dim=c(n) for a coordinate. A single-index subset of
  # such an array also retains dim, so diff(x)-diff(x)[1] is not safe.
  values <- as.numeric(d$vals)
  if (length(values) != d$len || any(!is.finite(values)))
    stop("Invalid coordinate length or non-finite values: ", n[hit])
  list(name = n[hit], values = values, units = d$units,
       len = length(values), source_dim = dim(d$vals))
}

regular_coordinate_step <- function(x, label, tolerance = 1e-8) {
  if (!is.numeric(x)) stop("Non-numeric coordinate: ", label)
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x)))
    stop("Coordinate requires at least two finite values: ", label)
  spacing <- diff(x)
  step <- spacing[[1L]]
  if (step == 0 || max(abs(spacing - step)) > tolerance)
    stop("This pipeline requires regular, ordered coordinates: ", label)
  step
}

coordinate_input_preflight <- function(path, output_dir) {
  # Read coordinates only, before expensive model loading/QC or file hashes.
  .set_analysis_stage("coordinate_preflight", file = path, announce = TRUE)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  rows <- lapply(c("lon", "lat", "time"), function(axis) {
    coordinate <- nc_axis(nc, axis)
    step <- NA_real_
    legacy <- "not_applicable"
    if (axis != "time") {
      step <- regular_coordinate_step(coordinate$values, axis)
      raw <- nc$dim[[coordinate$name]]$vals
      # This caught probe reports whether the previous expression fails on
      # the actual local file; it never supplies analysis values.
      legacy <- tryCatch({
        max(abs(diff(raw) - diff(raw)[1L]))
        "passed"
      }, error = function(e) conditionMessage(e))
    }
    data.frame(axis = axis, name = coordinate$name, length = coordinate$len,
      source_dim = if (is.null(coordinate$source_dim)) "none" else
        paste(coordinate$source_dim, collapse = "x"),
      normalized_dim = "none", regular_step = step,
      legacy_spacing_check = legacy, normalized_check = "passed")
  })
  out <- do.call(rbind, rows)
  utils::write.csv(out, file.path(output_dir, "coordinate_preflight.csv"), row.names = FALSE)
  if (any(grepl("non-conformable arrays", out$legacy_spacing_check, fixed = TRUE)))
    message("Confirmed on local coordinates: the legacy array expression fails; normalized coordinates pass.")
  message("Coordinate preflight passed: ", paste(out$axis, out$length, sep = "=", collapse = ", "))
  invisible(out)
}

decode_cf_daily <- function(values, units, calendar, year) {
  if (length(calendar) != 1L || is.na(calendar) || !nzchar(calendar)) calendar <- "standard"
  if (!tolower(calendar) %in% c("standard", "gregorian", "proleptic_gregorian"))
    stop("Unsupported CF calendar for ERA5 civil dates: ", calendar)
  if (length(units) != 1L || is.na(units)) stop("CF time units are missing.")
  m <- regexec("^(seconds?|minutes?|hours?|days?)\\s+since\\s+(.+)$", units, ignore.case = TRUE)
  hit <- regmatches(units, m)[[1]]
  if (length(hit) != 3L) stop("Unsupported CF time units: ", units)
  origin <- trimws(hit[3])
  origin <- sub("^([0-9]{4}-[0-9]{2}-[0-9]{2})[Tt]", "\\1 ", origin)
  origin <- sub("(Z| UTC| GMT| ?[+]00:00)$", "", origin)
  if (grepl("[+-][0-9]{2}:?[0-9]{2}$", origin))
    stop("Non-UTC CF origin is unsupported; normalize the input explicitly.")
  if (nchar(origin) == 10L) origin <- paste(origin, "00:00:00")
  epoch <- as.POSIXct(origin, tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M"))
  if (is.na(epoch)) stop("Cannot decode CF time origin: ", units)
  factor <- switch(sub("s$", "", tolower(hit[2])), second = 1, minute = 60, hour = 3600, day = 86400)
  if (any(!is.finite(values))) stop("Non-finite CF time coordinate.")
  timestamps <- epoch + as.numeric(values) * factor
  dates <- as.Date(timestamps, tz = "UTC")
  expected <- seq(as.Date(sprintf("%04d-01-01", year)), as.Date(sprintf("%04d-12-31", year)), by = "day")
  if (length(dates) != length(expected) || anyDuplicated(dates) ||
      any(diff(as.numeric(timestamps)) <= 0) || !identical(as.character(dates), as.character(expected)))
    stop("Year ", year, " has missing, duplicate, unsorted, or incorrectly assigned dates; expected ", length(expected), ".")
  # Daily coordinates must also refer to a consistent clock time.
  tod <- as.numeric(timestamps) %% 86400
  if (max(abs(tod - tod[1])) > 1e-4) stop("Daily NetCDF coordinates use inconsistent times of day.")
  list(dates = dates, calendar = calendar, time_units = units,
       time_of_day_seconds = tod[1])
}

temperature_unit <- function(unit, variable, allow_kelvin_background = FALSE) {
  if (length(unit) != 1L || is.na(unit) || !nzchar(unit)) stop("Missing units for ", variable)
  key <- tolower(gsub("[[:space:]_]", "", unit))
  celsius <- key %in% c("c", "degc", "degreecelsius", "degreescelsius", "celsius", "°c")
  kelvin <- key %in% c("k", "kelvin", "degreekelvin", "degreeskelvin")
  if (celsius) return(0)
  if (kelvin && variable != "B21_raw") return(0) # temperature differences: 1 K = 1 degC
  if (kelvin && isTRUE(allow_kelvin_background)) return(-273.15)
  if (kelvin) stop("B21_raw is in Kelvin. Set allow_kelvin_background=TRUE only after checking preprocessing; conversion will be audited.")
  stop("Unrecognized temperature units for ", variable, ": ", unit)
}

same_coordinate <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b)) return(FALSE)
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) > 0L && length(a) == length(b) &&
    all(is.finite(a)) && all(is.finite(b)) && max(abs(a - b)) < 1e-8
}

inspect_history_files <- function(g) {
  paths <- file.path(g$input_dir, sprintf(g$file_template, g$years))
  if (!all(file.exists(paths))) stop("Missing thermal-history NetCDF files: ", paste(paths[!file.exists(paths)], collapse = ", "))
  info <- vector("list", length(paths)); audit <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    .set_analysis_stage("inspect_ERA5_input", file = paths[i], year = g$years[i])
    message("Checking ERA5 input ", i, "/", length(paths), ": ", basename(paths[i]))
    nc <- ncdf4::nc_open(paths[i])
    item <- tryCatch({
      lon <- nc_axis(nc, "lon"); lat <- nc_axis(nc, "lat"); tm <- nc_axis(nc, "time")
      if (anyDuplicated(lon$values) || anyDuplicated(lat$values)) stop("Duplicate spatial coordinates.")
      regular_coordinate_step(lon$values, "longitude")
      regular_coordinate_step(lat$values, "latitude")
      lon_normal <- ((lon$values + 180) %% 360) - 180
      if (anyDuplicated(round(lon_normal, 8))) stop("Duplicate cyclic meridian (e.g. both -180 and +180). Remove one explicitly.")
      if (any(abs(lat$values) > 90) || any(lon$values < -180 | lon$values > 360)) stop("Invalid spatial coordinate range.")
      time_info <- decode_cf_daily(tm$values, tm$units, nc_attribute(nc, tm$name, "calendar"), g$years[i])
      missing <- setdiff(era5_required_variables, names(nc$var))
      if (length(missing)) stop("Missing history variables: ", paste(missing, collapse = ", "))
      offsets <- setNames(numeric(4), era5_required_variables)
      units <- character(4)
      for (j in seq_along(era5_required_variables)) {
        v <- era5_required_variables[j]; dims <- nc$var[[v]]$dim
        axes <- normal_axis(vapply(dims, function(x) x$name, character(1)))
        if (any(vapply(c("lon", "lat", "time"), function(z) sum(axes == z) != 1L, logical(1))))
          stop("Variable ", v, " does not have unique longitude, latitude, and time axes.")
        extra <- which(!axes %in% c("lon", "lat", "time"))
        if (length(extra) && any(vapply(dims[extra], function(x) as.integer(x$len), integer(1)) != 1L))
          stop("Non-singleton additional NetCDF dimensions in ", v)
        expected_len <- c(lon = lon$len, lat = lat$len, time = tm$len)
        for (z in names(expected_len)) if (dims[[which(axes == z)]]$len != expected_len[z]) stop("Variable coordinate length mismatch: ", v)
        units[j] <- nc_attribute(nc, v, "units")
        offsets[v] <- temperature_unit(units[j], v, g$allow_kelvin_background)
      }
      first_attribute <- function(names) {
        for (a in names) {
          x <- nc_attribute(nc, 0, a)
          if (length(x) && !is.na(x[1]) && nzchar(as.character(x[1])))
            return(paste(x, collapse = ";"))
        }
        "missing"
      }
      # The user's supplied Thermal history(1).py writes daily_time_zone,
      # B21_definition and TV_definition. Recognize these documented aliases.
      history_definition <- first_attribute("history_definition")
      if (history_definition == "missing") {
        bdef <- first_attribute("B21_definition")
        tdef <- first_attribute("TV_definition")
        if (bdef != "missing" && tdef != "missing")
          history_definition <- paste("B21:", bdef, "| TV:", tdef)
      }
      provenance <- c(source = first_attribute("source"),
        time_zone = first_attribute(c("time_zone", "daily_time_zone")),
        history_definition = history_definition,
        history_definitions = first_attribute("history_definitions"))
      if (provenance["time_zone"] != "missing" &&
          !tolower(gsub("[[:space:]]", "", provenance["time_zone"])) %in% c("utc", "utc+00:00", "utc+0", "gmt", "utc+00", "0"))
        stop("Input aggregation time_zone is not UTC; match the fitted exposure definition explicitly before running.")
      list(path = paths[i], year = g$years[i], lon = lon$values, lat = lat$values,
           dates = time_info$dates, time_info = time_info, offsets = offsets,
           units = units, provenance = provenance)
    }, finally = ncdf4::nc_close(nc))
    if (i > 1L && (!same_coordinate(item$lon, info[[1]]$lon) || !same_coordinate(item$lat, info[[1]]$lat)))
      stop("Spatial coordinates/order differ between annual files: ", paths[i])
    info[[i]] <- item
    file_meta <- file.info(paths[i])
    .set_analysis_stage("audit_and_checksum_ERA5_input", file = paths[i], year = g$years[i])
    audit[[i]] <- data.frame(file = normalizePath(paths[i], winslash = "/"), year = item$year,
      n_dates = length(item$dates), start = as.character(min(item$dates)), end = as.character(max(item$dates)),
      time_units = item$time_info$time_units, calendar = item$time_info$calendar,
      nlon = length(item$lon), nlat = length(item$lat), input_units = paste(item$units, collapse = "|"),
      B21_conversion_offset = item$offsets["B21_raw"], source = item$provenance["source"],
      aggregation_timezone = item$provenance["time_zone"],
      history_definition = paste(item$provenance[c("history_definition", "history_definitions")], collapse = "|"),
      bytes = file_meta$size, modification_time = as.character(file_meta$mtime),
      md5 = if (isTRUE(g$checksum_inputs)) unname(tools::md5sum(paths[i])) else "not_computed")
  }
  if (any(vapply(info, function(x) any(x$provenance[c("source", "time_zone")] == "missing") ||
                            all(x$provenance[c("history_definition", "history_definitions")] == "missing"), logical(1)))) {
    msg <- paste("Preprocessing provenance is incomplete. Daily dates alone cannot verify aggregation timezone,",
      "B21 exclusion of current day, TV lag formulas, Celsius conversion, or cross-year warm-up.",
      "Audit the preprocessing and fitted weather-source settings, then explicitly set",
      "global$preprocessing_acknowledged=TRUE. This does not supply missing evidence.")
    if (!isTRUE(g$preprocessing_acknowledged)) stop(msg)
    warning(msg)
  }
  list(files = info, audit = do.call(rbind, audit))
}

read_rectangular_history <- function(nc, variable, lon_indices, lat_indices, start_day, nday, offset = 0) {
  dims <- nc$var[[variable]]$dim
  axes <- normal_axis(vapply(dims, function(x) x$name, character(1)))
  start <- rep(1L, length(dims)); count <- vapply(dims, function(x) as.integer(x$len), integer(1))
  selections <- list(lon = lon_indices, lat = lat_indices, time = seq.int(start_day, length.out = nday))
  for (axis in names(selections)) {
    k <- which(axes == axis); ix <- selections[[axis]]
    start[k] <- min(ix); count[k] <- max(ix) - min(ix) + 1L
  }
  arr <- ncdf4::ncvar_get(nc, variable, start = start, count = count, collapse_degen = FALSE)
  order <- c(match(c("lon", "lat", "time"), axes), which(!axes %in% c("lon", "lat", "time")))
  arr <- aperm(arr, order)
  dim(arr) <- count[match(c("lon", "lat", "time"), axes)]
  arr <- arr[lon_indices - min(lon_indices) + 1L, lat_indices - min(lat_indices) + 1L, , drop = FALSE]
  matrix(as.vector(arr) + offset, nrow = length(lon_indices) * length(lat_indices), ncol = nday)
}

load_land_fraction <- function(g, lon, lat, selected_lon, selected_lat) {
  points <- expand.grid(lon = lon[selected_lon], lat = lat[selected_lat], KEEP.OUT.ATTRS = FALSE)
  if (is.null(g$land_fraction_file)) {
    if (!requireNamespace("maps", quietly = TRUE)) stop("Install maps for the explicit approximate land-mask fallback.")
    fraction <- as.numeric(!is.na(maps::map.where("world", x = ((points$lon + 180) %% 360) - 180, y = points$lat)))
    method <- "approximate maps::map.where world polygon; center-in-polygon; binary"
    warning("No land_fraction_file: using an approximate polygon land mask; coastlines and small islands require sensitivity checks.")
  } else {
    if (!file.exists(g$land_fraction_file)) stop("Land-fraction NetCDF file not found.")
    nc <- ncdf4::nc_open(g$land_fraction_file)
    fraction <- tryCatch({
      if (!same_coordinate(nc_axis(nc, "lon")$values, lon) || !same_coordinate(nc_axis(nc, "lat")$values, lat))
        stop("Land-fraction coordinates must exactly match annual history coordinates/order.")
      v <- g$land_fraction_variable
      if (!v %in% names(nc$var)) stop("Land-fraction variable is missing: ", v)
      dims <- nc$var[[v]]$dim; axes <- normal_axis(vapply(dims, function(x) x$name, character(1)))
      if (sum(axes == "lon") != 1L || sum(axes == "lat") != 1L) stop("Land fraction lacks lon/lat axes.")
      extra <- which(!axes %in% c("lon", "lat"))
      if (length(extra) && any(vapply(dims[extra], function(x) as.integer(x$len), integer(1)) != 1L)) stop("Land fraction must be time-independent or have singleton time.")
      a <- ncdf4::ncvar_get(nc, v, collapse_degen = FALSE)
      a <- aperm(a, c(match(c("lon", "lat"), axes), extra))
      dim(a) <- c(length(lon), length(lat))
      as.vector(a[selected_lon, selected_lat, drop = FALSE])
    }, finally = ncdf4::nc_close(nc))
    if (any(!is.finite(fraction)) || any(fraction < 0 | fraction > 1)) stop("Land fractions must be finite and in [0,1].")
    method <- paste("matched NetCDF land fraction:", normalizePath(g$land_fraction_file, winslash = "/"), g$land_fraction_variable)
  }
  points$cell_id <- seq_len(nrow(points))
  points$lon_plot <- ((points$lon + 180) %% 360) - 180
  points$land_fraction <- fraction
  points$land <- fraction >= g$land_fraction_threshold & fraction > 0
  # Equal-size angular point-sampling cells: cos(latitude) gives relative area.
  # This is not population weighting and not native-grid spatial aggregation.
  points$area_weight <- cos(points$lat * pi / 180) * fraction
  list(points = points, method = method)
}

sum_by_cell <- function(x, selected, cell_id, ncell) {
  out <- matrix(0, nrow = ncell, ncol = ncol(x))
  if (!any(selected)) return(out)
  grouped <- rowsum(x[selected, , drop = FALSE], cell_id[selected], reorder = FALSE)
  out[as.integer(rownames(grouped)), ] <- grouped
  out
}

posterior_rows <- function(x, denominator) {
  out <- matrix(NA_real_, nrow(x), 5L)
  good <- denominator > 0 & is.finite(denominator)
  if (any(good)) {
    y <- x[good, , drop = FALSE] / denominator[good]
    q <- t(apply(y, 1L, stats::quantile, probs = c(0.025, 0.5, 0.975), names = FALSE, type = 8))
    out[good, 1:3] <- q
    out[good, 4] <- rowMeans(y > 0)
    out[good, 5] <- rowMeans(y < 0)
  }
  colnames(out) <- c("lower", "median", "upper", "prob_positive", "prob_negative")
  as.data.frame(out)
}

# EXACT FAST GLOBAL MODULES (embedded; no external scripts needed)
# Exact arithmetic acceleration for the unchanged M1/M0 estimand.
# This file defines functions and source text only: no packages are installed,
# no models are read, and no run starts when it is sourced. The original
# make_history_design(), .coefficient_columns(), .history_frame(), and
# predict_ab() must also be available. Compile once in the master and call
# fast_init_backend() separately in each worker; never export an externalptr.

FAST_CPP_CODE <- paste(c(
  '#include <Rcpp.h>',
  '#include <cmath>',
  '#include <limits>',
  'using namespace Rcpp;',
  '',
  '// Positive/negative forms avoid overflow, including infinite eta.',
  'inline double fast_stable_eta(const double eta) {',
  '  if (eta >= 0.0) { const double e = std::exp(-eta); return 1.0 / (1.0 + e); }',
  '  const double e = std::exp(eta); return e / (1.0 + e);',
  '}',
  '',
  '// e = exp(-abs(z)); cached ep/en = exp(+/- c). Extreme z/c use',
  '// the direct stable form to avoid underflow followed by cancellation.',
  'inline double fast_shifted_probability(const double z, const double e,',
  '    const double c, const double ep, const double en) {',
  '  if (std::abs(z) > 700.0 || std::abs(c) > 700.0)',
  '    return fast_stable_eta(z + c);',
  '  if (z >= 0.0) return 1.0 / (1.0 + e * en);',
  '  const double q = e * ep; return q / (1.0 + q);',
  '}',
  '',
  '// [[Rcpp::export]]',
  'List fast_shift_ab_cpp(NumericMatrix zN, NumericMatrix zS,',
  '    NumericMatrix covN, NumericMatrix covS,',
  '    NumericMatrix exp_posN, NumericMatrix exp_negN,',
  '    NumericMatrix exp_posS, NumericMatrix exp_negS,',
  '    NumericVector w, double sfloor, double srange,',
  '    double nfloor, double nrange) {',
  '  const int nh = zN.nrow(), nd = zN.ncol(), ng = w.size();',
  '  if (nh < 1 || nd < 1 || ng < 1 || zS.nrow() != nh || zS.ncol() != nd)',
  '    stop("Invalid history/draw dimensions in fused kernel.");',
  '  if (covN.nrow() != ng || covN.ncol() != nd || covS.nrow() != ng || covS.ncol() != nd ||',
  '      exp_posN.nrow() != ng || exp_posN.ncol() != nd || exp_negN.nrow() != ng || exp_negN.ncol() != nd ||',
  '      exp_posS.nrow() != ng || exp_posS.ncol() != nd || exp_negS.nrow() != ng || exp_negS.ncol() != nd)',
  '    stop("Invalid cached stratum/draw dimensions in fused kernel.");',
  '  NumericMatrix A(nh, nd), B(nh, nd);',
  '  for (int d = 0; d < nd; ++d) {',
  '    if ((d & 15) == 0) checkUserInterrupt();',
  '    for (int h = 0; h < nh; ++h) {',
  '      const double zn = zN(h,d), zs = zS(h,d);',
  '      const double en = std::exp(-std::abs(zn)), es = std::exp(-std::abs(zs));',
  '      double a = 0.0, b = 0.0;',
  '      for (int g = 0; g < ng; ++g) {',
  '        const double pn = fast_shifted_probability(zn, en, covN(g,d), exp_posN(g,d), exp_negN(g,d));',
  '        const double ps = fast_shifted_probability(zs, es, covS(g,d), exp_posS(g,d), exp_negS(g,d));',
  '        const double n = nfloor + nrange * pn, s = sfloor + srange * ps;',
  '        a += w[g] * s;',
  '        b += (w[g] * s) * n;',
  '      }',
  '      A(h,d) = a; B(h,d) = b;',
  '    }',
  '  }',
  '  return List::create(Named("A1") = A, Named("B1") = B);',
  '}'
), collapse = "\n")

.fast_get_state <- function() {
  # Export functions/constants only. Each PSOCK worker creates its own state
  # here, in the functions' enclosing environment, after deserialization.
  target <- environment(fast_init_backend)
  if (!exists(".fast_kernel_state", envir = target, inherits = FALSE)) {
    state <- new.env(parent = emptyenv())
    state$backend <- "r"
    state$cpp <- NULL
    state$compile_error <- NULL
    assign(".fast_kernel_state", state, envir = target)
  }
  get(".fast_kernel_state", envir = target, inherits = FALSE)
}

fast_init_backend <- function(cache_dir, backend = "auto") {
  if (length(backend) != 1L || !(backend %in% c("auto", "cpp", "r")))
    stop("fast backend must be 'auto', 'cpp', or 'r'.")
  state <- .fast_get_state()
  state$backend <- "r"
  state$cpp <- NULL
  state$compile_error <- NULL
  if (backend == "r") return("r")
  compile <- tryCatch({
    if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is unavailable.")
    if (length(cache_dir) != 1L || is.na(cache_dir) || !nzchar(cache_dir))
      stop("A nonempty persistent fast cache directory is required.")
    if (!dir.exists(cache_dir) && !dir.create(cache_dir, recursive = TRUE))
      stop("Cannot create the fast kernel cache directory: ", cache_dir)
    source_file <- file.path(cache_dir, "thermal_history_shift_kernel.cpp")
    source_lines <- strsplit(FAST_CPP_CODE, "\n", fixed = TRUE)[[1L]]
    previous <- if (file.exists(source_file)) readLines(source_file, warn = FALSE) else character()
    if (!identical(previous, source_lines)) writeLines(source_lines, source_file, useBytes = TRUE)
    cpp_env <- new.env(parent = baseenv())
    Rcpp::sourceCpp(file = source_file, env = cpp_env, rebuild = FALSE,
      cacheDir = file.path(cache_dir, "rcpp"), showOutput = FALSE, verbose = FALSE)
    if (!exists("fast_shift_ab_cpp", envir = cpp_env, inherits = FALSE))
      stop("Compiled kernel did not export fast_shift_ab_cpp.")
    get("fast_shift_ab_cpp", envir = cpp_env, inherits = FALSE)
  }, error = function(e) e)
  if (inherits(compile, "error")) {
    state$compile_error <- conditionMessage(compile)
    if (backend == "cpp") stop("Requested C++ fast backend cannot load: ", state$compile_error)
    warning("C++ fast backend unavailable; using the stable R fallback. ", state$compile_error,
            call. = FALSE)
    return("r")
  }
  state$cpp <- compile
  state$backend <- "cpp"
  "cpp"
}

.fast_plain_matrix <- function(x) {
  if (!is.numeric(x) || length(dim(x)) != 2L) stop("Expected a numeric matrix in fast cache.")
  matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
}

.fast_cached_exp <- function(cov, sign) {
  # The dummy value 1 in unsafe columns is overwritten by the direct fallback
  # before any returned probability is used; no overflow warnings are emitted.
  safe <- is.finite(cov) & abs(cov) <= 700
  out <- matrix(1, nrow(cov), ncol(cov))
  out[safe] <- exp(sign * cov[safe])
  out
}

fast_make_cache <- function(analysis, weight_scheme = "observation") {
  required <- c("make_history_design", ".coefficient_columns", ".history_frame")
  if (!all(vapply(required, exists, logical(1), mode = "function", inherits = TRUE)))
    stop("Source the unchanged model core before constructing the fast cache.")
  if (is.numeric(weight_scheme)) w <- as.numeric(weight_scheme) else {
    aliases <- c(observation = "observation", building = "building", row_equal_weight = "observation",
                 building_equal_weight = "building")
    if (length(weight_scheme) != 1L || !(weight_scheme %in% names(aliases))) stop("Unknown weight scheme.")
    w <- as.numeric(analysis$weights[[aliases[[weight_scheme]]]])
  }
  if (length(w) != nrow(analysis$strata) || any(!is.finite(w) | w < 0) || sum(w) <= 0)
    stop("Invalid stratum weights in fast cache.")
  w <- w / sum(w)
  beta1 <- .fast_plain_matrix(analysis$draws_m1)
  beta0 <- .fast_plain_matrix(analysis$draws_m0)
  if (!nrow(beta1) || nrow(beta1) != nrow(beta0)) stop("M1/M0 fast cache must use the same paired draw count.")
  if (any(!is.finite(beta1)) || any(!is.finite(beta0))) stop("Nonfinite posterior coefficients in fast cache.")
  p <- analysis$pred
  dn <- c(p$b21_neutral_names, p$cb_neutral_names)
  ds <- c(p$b21_slope_names, p$cb_slope_names)
  Xcov <- .fast_plain_matrix(analysis$Xcov)
  cov1N <- Xcov %*% t(.coefficient_columns(beta1, "etaneutral", colnames(Xcov)))
  cov1S <- Xcov %*% t(.coefficient_columns(beta1, "logslope", colnames(Xcov)))
  cov0N <- Xcov %*% t(.coefficient_columns(beta0, "etaneutral", colnames(Xcov)))
  cov0S <- Xcov %*% t(.coefficient_columns(beta0, "logslope", colnames(Xcov)))
  if (any(!is.finite(cov1N)) || any(!is.finite(cov1S)) ||
      any(!is.finite(cov0N)) || any(!is.finite(cov0S))) stop("Nonfinite stratum predictors in fast cache.")
  model <- p$model
  pars <- c(slope_floor = model$slope_floor, slope_range = model$slope_range,
            neutral_floor = model$neutral_floor, neutral_range = model$neutral_range)
  if (length(pars) != 4L || any(!is.finite(pars))) stop("Invalid bounded model constants in fast cache.")
  # Preserve the original crossprod over all strata for the history-free model.
  n0 <- model$neutral_floor + model$neutral_range * plogis(cov0N)
  s0 <- model$slope_floor + model$slope_range * plogis(cov0S)
  active <- which(w > 0)
  cov1N <- cov1N[active, , drop = FALSE]
  cov1S <- cov1S[active, , drop = FALSE]
  nd <- nrow(beta1)
  ids1 <- analysis$draw_ids_m1
  ids0 <- analysis$draw_ids_m0
  if (length(ids1) != nd || length(ids0) != nd || anyNA(ids1) || anyNA(ids0))
    stop("Original model draw IDs are missing or do not match the fast cache.")
  structure(list(pred = p, dn = dn, ds = ds,
    history_coef_N = t(.coefficient_columns(beta1, "etaneutral", dn)),
    history_coef_S = t(.coefficient_columns(beta1, "logslope", ds)),
    covN = cov1N, covS = cov1S,
    exp_posN = .fast_cached_exp(cov1N, 1), exp_negN = .fast_cached_exp(cov1N, -1),
    exp_posS = .fast_cached_exp(cov1S, 1), exp_negS = .fast_cached_exp(cov1S, -1),
    w = as.numeric(w[active]), active_strata = active, model = model,
    A0 = as.numeric(crossprod(w, s0)), B0 = as.numeric(crossprod(w, s0 * n0)),
    paired_draw_indices = seq_len(nd), draw_ids_m1 = ids1, draw_ids_m0 = ids0,
    weight_scheme = weight_scheme, n_draws = nd, n_strata = length(active),
    estimand = "study-sample-standardized expected TSV; building effects excluded",
    arithmetic = "exact logistic-shift algebra with extreme-logit direct fallback; fixed paired draws"),
    class = "thermal_history_fast_cache")
}

.fast_shift_probability_r <- function(z, e, negative, bad_z, bad_z_col, cov, ep, en) {
  out <- 1 / (1 + sweep(e, 2, en, "*"))
  if (length(negative)) {
    q <- sweep(e, 2, ep, "*")
    out[negative] <- q[negative] / (1 + q[negative])
  }
  bad_cols <- which(abs(cov) > 700)
  if (length(bad_cols)) out[, bad_cols] <- plogis(sweep(z[, bad_cols, drop = FALSE], 2,
                                                      cov[bad_cols], "+"))
  if (length(bad_z)) out[bad_z] <- plogis(z[bad_z] + cov[bad_z_col])
  out
}

.fast_r_kernel <- function(cache, zN, zS) {
  en <- exp(-abs(zN)); es <- exp(-abs(zS))
  negN <- which(zN < 0); negS <- which(zS < 0)
  badN <- which(abs(zN) > 700); badS <- which(abs(zS) > 700)
  # Column-major linear indices retain singleton dimensions safely.
  badNcol <- if (length(badN)) ((badN - 1L) %/% nrow(zN)) + 1L else integer()
  badScol <- if (length(badS)) ((badS - 1L) %/% nrow(zS)) + 1L else integer()
  A <- B <- matrix(0, nrow(zN), ncol(zN))
  m <- cache$model
  for (g in seq_along(cache$w)) {
    pn <- .fast_shift_probability_r(zN, en, negN, badN, badNcol,
      cache$covN[g, ], cache$exp_posN[g, ], cache$exp_negN[g, ])
    ps <- .fast_shift_probability_r(zS, es, negS, badS, badScol,
      cache$covS[g, ], cache$exp_posS[g, ], cache$exp_negS[g, ])
    n <- m$neutral_floor + m$neutral_range * pn
    s <- m$slope_floor + m$slope_range * ps
    A <- A + cache$w[g] * s
    B <- B + cache$w[g] * s * n
  }
  list(A1 = A, B1 = B)
}

fast_predict_ab <- function(cache, history) {
  if (!inherits(cache, "thermal_history_fast_cache")) stop("Use fast_make_cache() before fast_predict_ab().")
  h <- .history_frame(history)
  design <- make_history_design(h, cache$pred)
  zN <- as.matrix(design[cache$dn]) %*% cache$history_coef_N
  zS <- as.matrix(design[cache$ds]) %*% cache$history_coef_S
  if (anyNA(zN) || anyNA(zS)) stop("NA/NaN history predictor in fast kernel.")
  state <- .fast_get_state()
  if (identical(state$backend, "cpp")) {
    if (!is.function(state$cpp)) stop("C++ backend was selected without a process-local compiled kernel.")
    m <- cache$model
    out <- state$cpp(zN, zS, cache$covN, cache$covS,
      cache$exp_posN, cache$exp_negN, cache$exp_posS, cache$exp_negS, cache$w,
      m$slope_floor, m$slope_range, m$neutral_floor, m$neutral_range)
  } else out <- .fast_r_kernel(cache, zN, zS)
  if (!identical(dim(out$A1), c(nrow(h), cache$n_draws)) ||
      !identical(dim(out$B1), c(nrow(h), cache$n_draws))) stop("Fast kernel returned incorrect history/draw dimensions.")
  if (any(!is.finite(out$A1)) || any(!is.finite(out$B1)) ||
      any(!is.finite(cache$A0)) || any(!is.finite(cache$B0))) stop("Nonfinite fast standardized prediction.")
  c(out, list(A0 = cache$A0, B0 = cache$B0,
    paired_draw_indices = cache$paired_draw_indices, draw_ids_m1 = cache$draw_ids_m1,
    draw_ids_m0 = cache$draw_ids_m0, weight_scheme = cache$weight_scheme,
    estimand = cache$estimand, backend = state$backend))
}

.fast_qc_histories <- function(analysis) {
  p <- analysis$pred
  train <- if (!is.null(analysis$training_history)) .history_frame(analysis$training_history) else
    .training_history(analysis$data, p)
  finite <- which(rowSums(!is.finite(as.matrix(train))) == 0L)
  if (!length(finite)) stop("No finite training history is available for fast-kernel QC.")
  selected <- finite[unique(as.integer(round(seq(1, length(finite), length.out = min(24L, length(finite))))))]
  for (column in names(train)) {
    x <- train[[column]][finite]
    selected <- c(selected, finite[which.min(x)], finite[which.max(x)])
  }
  representative <- unique(train[unique(selected), , drop = FALSE])
  lower <- c(p$B21_prediction_P1, rep(p$var_low, 3L))
  upper <- c(p$B21_prediction_P99, rep(p$var_high, 3L))
  mid <- (lower + upper) / 2
  corners <- expand.grid(rep(list(c(0, 1)), 4L), KEEP.OUT.ATTRS = FALSE)
  boundary <- as.data.frame(sweep(sweep(as.matrix(corners), 2, upper - lower, "*"), 2, lower, "+"))
  names(boundary) <- history_columns()
  ref <- as.data.frame(t(c(p$B21_ref, 0, 0, 0)))
  names(ref) <- history_columns()
  tails <- matrix(rep(mid, each = 8L), nrow = 8L)
  for (j in seq_len(4L)) {
    tails[2L * j - 1L, j] <- lower[j] - 10 * (upper[j] - lower[j])
    tails[2L * j, j] <- upper[j] + 10 * (upper[j] - lower[j])
  }
  tails <- as.data.frame(tails); names(tails) <- history_columns()
  list(training = representative, boundaries_and_reference = rbind(boundary, ref),
       clipped_tails = tails, singleton_training = representative[1L, , drop = FALSE],
       singleton_tail = tails[nrow(tails), , drop = FALSE])
}

fast_kernel_QC <- function(analysis, cache, output_dir, tolerance = 1e-8) {
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance > 1e-8)
    stop("Fast-kernel tolerance must be positive and no greater than 1e-8.")
  if (!identical(as.integer(cache$draw_ids_m1), as.integer(analysis$draw_ids_m1)) ||
      !identical(as.integer(cache$draw_ids_m0), as.integer(analysis$draw_ids_m0)))
    stop("Fast-kernel QC draw IDs differ from the unchanged original prediction IDs.")
  if (length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) stop("Fast QC needs an output directory.")
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create fast QC directory.")
  histories <- .fast_qc_histories(analysis)
  records <- list()
  append_check <- function(case, quantity, reference, accelerated, nh) {
    same_dimension <- identical(dim(reference), dim(accelerated)) && length(reference) == length(accelerated)
    err <- if (same_dimension && all(is.finite(reference)) && all(is.finite(accelerated)))
      max(abs(as.numeric(reference) - as.numeric(accelerated))) else Inf
    records[[length(records) + 1L]] <<- data.frame(case = case, quantity = quantity,
      backend = .fast_get_state()$backend, n_history = nh, n_draws = cache$n_draws,
      max_absolute_error = err, tolerance = tolerance, pass = is.finite(err) && err <= tolerance,
      stringsAsFactors = FALSE)
  }
  sets <- as.numeric(analysis$fixed_set)
  if (!length(sets) || any(!is.finite(sets))) sets <- c(15, 25, 35)
  for (case in names(histories)) {
    h <- histories[[case]]
    original <- predict_ab(analysis, h, weight_scheme = cache$weight_scheme)
    fast <- fast_predict_ab(cache, h)
    for (quantity in c("A1", "B1", "A0", "B0"))
      append_check(case, quantity, original[[quantity]], fast[[quantity]], nrow(h))
    for (s in sets) {
      ref_mu <- original$A1 * s - original$B1
      fast_mu <- fast$A1 * s - fast$B1
      append_check(case, paste0("M1_TSV_SET_", format(s, digits = 12)), ref_mu, fast_mu, nrow(h))
      append_check(case, paste0("delta_TSV_SET_", format(s, digits = 12)),
        sweep(ref_mu, 2, original$A0 * s - original$B0, "-"),
        sweep(fast_mu, 2, fast$A0 * s - fast$B0, "-"), nrow(h))
    }
  }
  result <- do.call(rbind, records)
  utils::write.csv(result, file.path(output_dir, "fast_kernel_qc.csv"), row.names = FALSE)
  saveRDS(list(checks = result, histories = histories, draw_ids_m1 = cache$draw_ids_m1,
    draw_ids_m0 = cache$draw_ids_m0, arithmetic = cache$arithmetic,
    note = "This additional check supplements rather than replaces the unchanged official posterior_epred QC."),
    file.path(output_dir, "fast_kernel_qc.rds"))
  if (!all(result$pass)) stop("Fast kernel differs from the unchanged predictor beyond tolerance; see fast_kernel_qc.csv.")
  message("Fast kernel QC PASS; backend=", .fast_get_state()$backend,
          "; max absolute error=", format(max(result$max_absolute_error), digits = 16))
  invisible(result)
}


# Independent exact ERA5 spatial-tile computation and immutable checkpoints.
# This module is sourced by the main process and by each PSOCK worker.

fast_atomic_publish <- function(path, writer, replace = FALSE) {
  directory <- dirname(path)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE))
    stop("Cannot create output directory: ", directory)
  temporary <- tempfile(paste0(".", basename(path), ".tmp_"), tmpdir = directory)
  backup <- NULL
  on.exit({
    if (file.exists(temporary)) unlink(temporary)
    if (!is.null(backup) && file.exists(backup) && !file.exists(path))
      file.rename(backup, path)
  }, add = TRUE)
  writer(temporary)
  if (!file.exists(temporary)) stop("Output writer did not create its temporary file: ", path)
  if (file.exists(path)) {
    if (!isTRUE(replace)) stop("Refusing to replace committed output: ", path)
    # Windows cannot rename onto an existing file. Only uncommitted partition
    # outputs, never checkpoints, use this recoverable replacement sequence.
    backup <- tempfile(paste0(".", basename(path), ".old_"), tmpdir = directory)
    if (!file.rename(path, backup)) stop("Cannot move uncommitted previous output: ", path)
  }
  if (!file.rename(temporary, path)) stop("Cannot commit output with same-directory rename: ", path)
  if (!is.null(backup) && file.exists(backup)) unlink(backup)
  invisible(path)
}

fast_atomic_save_rds <- function(object, path, replace = FALSE) {
  fast_atomic_publish(path, function(temporary) saveRDS(object, temporary, compress = "gzip"), replace)
}

fast_tile_signature <- function(task) {
  list(lon_indices = as.integer(task$lon_indices), lat_indices = as.integer(task$lat_indices),
       keep = as.integer(task$keep), cell_ids = as.integer(task$tile_points$cell_id))
}


fast_object_hash <- function(x) {
  path <- tempfile("era5_object_hash_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, compress = FALSE, version = 3L)
  value <- unname(tools::md5sum(path))
  if (length(value) != 1L || is.na(value) || !nzchar(value)) stop("Cannot hash serialized run inputs.")
  value
}

fast_local_host <- function() {
  host <- unname(Sys.info()["nodename"])
  if (length(host) != 1L || is.na(host) || !nzchar(host)) stop("Cannot determine the local host for the run lock.")
  tolower(host)
}

fast_process_state <- function(pid) {
  if (length(pid) != 1L || !is.numeric(pid) || !is.finite(pid) || pid < 1 ||
      pid > .Machine$integer.max || pid != floor(pid))
    return("unknown")
  pid <- as.integer(pid)
  if (identical(pid, as.integer(Sys.getpid()))) return("alive")
  if (.Platform$OS.type == "windows") {
    executable <- Sys.which("tasklist")
    if (!nzchar(executable)) return("unknown")
    output <- tryCatch(suppressWarnings(system2(executable, c("/FO", "CSV", "/NH"),
      stdout = TRUE, stderr = TRUE)), error = function(e) NULL)
    if (is.null(output) || (!is.null(attr(output, "status")) && attr(output, "status") != 0L))
      return("unknown")
    # Listing all PIDs avoids language-dependent 'no matching tasks' messages.
    # A complete CSV with our own PID is required before absence is accepted.
    output <- output[nzchar(trimws(output))]
    if (!length(output) || any(!startsWith(trimws(output), '"'))) return("unknown")
    tab <- tryCatch(utils::read.csv(text = paste(output, collapse = "\n"), header = FALSE,
      stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE), error = function(e) NULL)
    if (is.null(tab) || ncol(tab) != 5L || !nrow(tab)) return("unknown")
    if (any(!grepl("^[0-9]+$", trimws(tab[[2L]])))) return("unknown")
    pids <- suppressWarnings(as.integer(trimws(tab[[2L]])))
    if (anyNA(pids) || !Sys.getpid() %in% pids) return("unknown")
    return(if (pid %in% pids) "alive" else "dead")
  }
  kill <- Sys.which("kill")
  if (nzchar(kill)) {
    status <- tryCatch(suppressWarnings(system2(kill, c("-0", as.character(pid)),
      stdout = FALSE, stderr = FALSE)), error = function(e) NA_integer_)
    if (isTRUE(status == 0L)) return("alive")
    # A failed kill -0 can also mean insufficient permission. Do not infer death.
  }
  ps <- Sys.which("ps")
  if (!nzchar(ps)) return("unknown")
  output <- tryCatch(suppressWarnings(system2(ps, c("-p", as.character(pid), "-o", "pid="),
    stdout = TRUE, stderr = TRUE)), error = function(e) NULL)
  if (is.null(output)) return("unknown")
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  lines <- trimws(output)
  lines <- lines[nzchar(lines)]
  if (status == 0L && length(lines) && all(grepl("^[0-9]+$", lines)) &&
      pid %in% suppressWarnings(as.integer(lines))) return("alive")
  if (status == 1L && !length(lines)) return("dead")
  "unknown"
}

fast_lock_owner <- function(lock_dir) {
  owner_path <- file.path(lock_dir, "owner.rds")
  if (!file.exists(owner_path)) stop("Run lock has no readable owner record: ", lock_dir,
    ". Its owner cannot be verified, so the lock will not be removed automatically.", call. = FALSE)
  owner <- tryCatch(readRDS(owner_path), error = function(e) NULL)
  if (!is.list(owner) || !identical(owner$schema, "ERA5-global01-fast-lock-v1") ||
      length(owner$pid) != 1L || !is.numeric(owner$pid) || !is.finite(owner$pid) ||
      owner$pid < 1 || owner$pid > .Machine$integer.max || owner$pid != floor(owner$pid) ||
      !is.character(owner$host) || length(owner$host) != 1L || is.na(owner$host) || !nzchar(owner$host) ||
      !is.character(owner$token) || length(owner$token) != 1L || is.na(owner$token) || !nzchar(owner$token))
    stop("Run lock owner record is invalid: ", owner_path,
         ". The lock will not be removed automatically.", call. = FALSE)
  if (!is.null(owner$worker_pids) && (!is.numeric(owner$worker_pids) ||
      any(!is.finite(owner$worker_pids)) || any(owner$worker_pids < 1 | owner$worker_pids > .Machine$integer.max) ||
      any(owner$worker_pids != floor(owner$worker_pids)) || anyDuplicated(owner$worker_pids) ||
      any(owner$worker_pids == owner$pid)))
    stop("Run lock contains invalid worker PIDs: ", owner_path,
         ". The lock will not be removed automatically.", call. = FALSE)
  owner
}

fast_worker_process_states <- function(owner) {
  pids <- owner$worker_pids
  if (is.null(pids) || !length(pids)) return(setNames(character(), character()))
  setNames(vapply(as.integer(pids), fast_process_state, character(1)), as.character(pids))
}

fast_register_workers <- function(path, worker_audit, locktoken = NULL) {
  lock_dir <- file.path(path, ".run_lock")
  owner <- fast_lock_owner(lock_dir)
  if (!identical(as.integer(owner$pid), as.integer(Sys.getpid())) ||
      !identical(tolower(owner$host), fast_local_host()) ||
      (!is.null(locktoken) && !identical(owner$token, locktoken)))
    stop("Cannot register workers: this process does not own the run lock.", call. = FALSE)
  if (!is.list(worker_audit)) stop("Worker audit must be a list of worker records.")
  pids <- vapply(worker_audit, function(worker) {
    if (!is.list(worker)) stop("Worker audit contains a non-list record.")
    pid <- worker$pid
    if (length(pid) != 1L || !is.numeric(pid) || !is.finite(pid) ||
        pid < 1 || pid > .Machine$integer.max || pid != floor(pid))
      stop("Worker audit contains an invalid PID.")
    as.integer(pid)
  }, integer(1))
  # A single-process backend may report the master as its worker. Its lifetime
  # is already represented by owner$pid, and it must not block its own release.
  pids <- sort(unique(pids[pids != as.integer(owner$pid)]))
  owner$worker_pids <- sort(unique(c(as.integer(owner$worker_pids), pids)))
  owner$workers_registered <- Sys.time()
  current <- fast_lock_owner(lock_dir)
  if (!identical(current$token, owner$token) || !identical(current$pid, owner$pid) ||
      !identical(current$host, owner$host)) stop("Run lock changed before worker registration; no update was written.")
  fast_atomic_save_rds(owner, file.path(lock_dir, "owner.rds"), replace = TRUE)
  invisible(owner$worker_pids)
}

fast_acquire_run_lock <- function(path) {
  if (!dir.exists(path)) stop("Run directory does not exist: ", path)
  lock_dir <- file.path(path, ".run_lock")
  host <- fast_local_host()
  for (attempt in seq_len(5L)) {
    if (dir.create(lock_dir, showWarnings = FALSE)) {
      token <- paste0(Sys.getpid(), "_", basename(tempfile("era5_lock_")))
      owner <- list(schema = "ERA5-global01-fast-lock-v1", pid = Sys.getpid(),
                    host = host, token = token, started = Sys.time())
      owner_ok <- FALSE
      tryCatch({
        fast_atomic_save_rds(owner, file.path(lock_dir, "owner.rds"))
        owner_ok <- TRUE
      }, error = function(e) {
        if (!owner_ok) unlink(lock_dir, recursive = TRUE)
        stop(e)
      })
      return(token)
    }
    if (!dir.exists(lock_dir)) stop("Cannot create the run lock: ", lock_dir)
    owner <- fast_lock_owner(lock_dir)
    if (!identical(tolower(owner$host), host)) stop(
      "Run is locked by host ", owner$host, " (PID ", owner$pid,
      "). A process on another host cannot be verified locally; no lock was removed.", call. = FALSE)
    state <- fast_process_state(owner$pid)
    if (!identical(state, "dead")) stop("Run is locked by PID ", owner$pid,
      "; process state is ", state, ". No process was stopped and no lock was removed.", call. = FALSE)
    worker_states <- fast_worker_process_states(owner)
    blocked <- worker_states[worker_states != "dead"]
    if (length(blocked)) stop("The previous master is absent, but registered workers are still active or unverifiable: ",
      paste(paste0(names(blocked), "=", blocked), collapse = ", "),
      ". Leaving the run lock in place to prevent concurrent writes.", call. = FALSE)
    # The claim directory serializes stale-lock recovery. After claiming, check
    # the owner again so a delayed contender never removes a newly acquired lock.
    claim <- file.path(lock_dir, ".recovery_claim")
    if (!dir.create(claim, showWarnings = FALSE)) stop(
      "Stale-lock recovery is already claimed or was interrupted: ", claim,
      ". Inspect this lock before trying again.", call. = FALSE)
    moved <- FALSE
    stale <- tempfile(".recovered_run_lock_", tmpdir = path)
    tryCatch({
      verified <- fast_lock_owner(lock_dir)
      if (!identical(verified$token, owner$token) || !identical(tolower(verified$host), host) ||
          !identical(fast_process_state(verified$pid), "dead") ||
          any(fast_worker_process_states(verified) != "dead"))
        stop("Run lock ownership changed during recovery; no lock was removed.")
      if (!file.rename(lock_dir, stale)) stop("Cannot quarantine the verified stale run lock.")
      moved <- TRUE
      unlink(stale, recursive = TRUE)
      message("Recovered stale run lock from verified absent local PID ", owner$pid, ".")
    }, finally = {
      if (!moved && dir.exists(claim)) unlink(claim, recursive = TRUE)
    })
  }
  stop("Could not acquire the run lock after concurrent recovery attempts: ", path)
}

fast_release_run_lock <- function(path, locktoken) {
  lock_dir <- file.path(path, ".run_lock")
  if (!dir.exists(lock_dir)) return(invisible(FALSE))
  owner <- tryCatch(fast_lock_owner(lock_dir), error = function(e) NULL)
  if (is.null(owner) || !identical(owner$token, locktoken) ||
      !identical(as.integer(owner$pid), as.integer(Sys.getpid())) ||
      !identical(tolower(owner$host), fast_local_host())) {
    warning("Run lock ownership no longer matches this process; leaving lock untouched: ", lock_dir)
    return(invisible(FALSE))
  }
  worker_states <- fast_worker_process_states(owner)
  blocked <- worker_states[worker_states != "dead"]
  if (length(blocked)) {
    warning("Leaving run lock in place because registered workers are active or unverifiable: ",
      paste(paste0(names(blocked), "=", blocked), collapse = ", "),
      ". Close the worker cluster first; no process was stopped.")
    return(invisible(FALSE))
  }
  unlink(lock_dir, recursive = TRUE)
  invisible(!dir.exists(lock_dir))
}

fast_is_same_path <- function(a, b) {
  a <- normalizePath(a, winslash = "/", mustWork = TRUE)
  b <- normalizePath(b, winslash = "/", mustWork = TRUE)
  if (.Platform$OS.type == "windows") identical(tolower(a), tolower(b)) else identical(a, b)
}

fast_existing_run <- function(path, root) {
  if (!dir.exists(path)) stop("Requested fast run directory does not exist: ", path)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!fast_is_same_path(dirname(path), root)) stop(
    "A resumed fast run must be directly inside the configured fast output root: ", root,
    ". Refusing a latest-run pointer to another location.", call. = FALSE)
  identity_path <- file.path(path, "fast_run_identity.rds")
  identity <- if (file.exists(identity_path)) tryCatch(readRDS(identity_path), error = function(e) NULL) else NULL
  if (!is.list(identity) || !identical(identity$schema, "ERA5-warm02-run-v1"))
    stop("This directory lacks a valid fast-run identity: ", path,
         ". Existing serial-run summaries are not resumable fast checkpoints.", call. = FALSE)
  status_path <- file.path(path, "standalone_run_status.rds")
  status <- ""
  if (file.exists(status_path)) {
    audit <- tryCatch(readRDS(status_path), error = function(e) NULL)
    if (!is.list(audit) || !is.character(audit$status) || length(audit$status) != 1L || is.na(audit$status))
      stop("Cannot read the previous run status: ", status_path, ". Inspect it before resuming.")
    status <- audit$status
  }
  if (!status %in% c("", "RUNNING", "FAILED", "COMPLETED"))
    stop("Unrecognized previous fast-run status: ", status)
  list(path = path, status = status)
}

fast_select_run_dir <- function(root, fast_config = list()) {
  if (!is.character(root) || length(root) != 1L || is.na(root) || !nzchar(root))
    stop("Fast output root must be a nonempty path.")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(root)) stop("Cannot create fast output root: ", root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  resume <- fast_config$resume
  if (is.null(resume)) resume <- "auto"
  if (!is.character(resume) || length(resume) != 1L || is.na(resume) || !nzchar(resume))
    stop("config$fast$resume must be 'auto', 'new', or a fast run directory.")
  if (!resume %in% c("auto", "new")) {
    path <- path.expand(resume)
    absolute <- grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)
    if (!absolute) path <- file.path(root, path)
    prior <- fast_existing_run(path, root)
    return(list(path = prior$path, resumed = TRUE, previous_status = prior$status, root = root))
  }
  latest <- file.path(root, "latest_run.txt")
  if (resume == "auto" && file.exists(latest)) {
    lines <- trimws(readLines(latest, warn = FALSE, encoding = "UTF-8"))
    lines <- lines[nzchar(lines)]
    if (length(lines) != 1L) stop("latest_run.txt must identify exactly one fast run directory: ", latest)
    prior <- fast_existing_run(lines, root)
    if (prior$status != "COMPLETED")
      return(list(path = prior$path, resumed = TRUE, previous_status = prior$status, root = root))
    stop("Latest warm02 run is already complete. Use MODE='plots' or the saved energy_cache_reader.R. To deliberately start another run set fast$resume='new'.")
  }
  stem <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  for (attempt in seq_len(1000L)) {
    suffix <- if (attempt == 1L) "" else sprintf("_%03d", attempt)
    path <- file.path(root, paste0(stem, suffix))
    if (dir.create(path, showWarnings = FALSE))
      return(list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
                  resumed = FALSE, previous_status = "", root = root))
    if (!file.exists(path) && !dir.exists(path)) stop("Cannot create fast run directory: ", path)
  }
  stop("Cannot allocate a distinct fast run directory under: ", root)
}

fast_prepare_run <- function(root, fast_config = list()) {
  selected <- fast_select_run_dir(root, fast_config)
  token <- fast_acquire_run_lock(selected$path)
  completed <- FALSE
  on.exit(if (!completed) fast_release_run_lock(selected$path, token), add = TRUE)
  # A previous process can finish between directory selection and lock
  # acquisition. Recheck after acquiring the lock before automatic resumption.
  automatic <- is.null(fast_config$resume) || identical(fast_config$resume, "auto")
  if (isTRUE(selected$resumed) && automatic) {
    current <- fast_existing_run(selected$path, selected$root)
    if (identical(current$status, "COMPLETED")) {
      fast_release_run_lock(selected$path, token)
      completed <- TRUE
      stop("The selected run finished before lock acquisition; no new run was started. Use plots or the saved cache reader.")
    }
  }
  if (!selected$resumed) {
    identity <- list(schema = "ERA5-warm02-run-v1", created = Sys.time(),
      path = selected$path, description = "Fast run identity; the full manifest is validated before checkpoint use.")
    fast_atomic_save_rds(identity, file.path(selected$path, "fast_run_identity.rds"))
  }
  latest <- file.path(selected$root, "latest_run.txt")
  fast_atomic_publish(latest, function(temporary)
    writeLines(enc2utf8(selected$path), temporary, useBytes = TRUE), replace = TRUE)
  completed <- TRUE
  list(path = selected$path, locktoken = token, resumed = selected$resumed,
       previous_status = selected$previous_status, root = selected$root)
}



# Warm-side boundary tile engine. M = B/A; N = 1/A is a positive
# half-width on the SET scale, NOT the original model's neutral parameter N.
# U = M + N is the TSV = +1 boundary of the standardized response curve.

# Daily warm-boundary cache: raw little-endian IEEE float32 by date chunk.
# Delta_U rows are supported cell-days, columns are the fixed paired draws.
# Full date-by-cell history flags are retained; unsupported is never zero demand.
# This module reads/writes only one tile/chunk at a time, without RStan.

fast_daily_fail <- function(...) stop(..., call. = FALSE)

fast_daily_scalar <- function(x, lower = -Inf, integer = FALSE) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= lower && (!integer || x == floor(x))
}

fast_daily_digest <- function(bytes) {
  if (!requireNamespace("digest", quietly = TRUE))
    fast_daily_fail("Daily caching requires the digest package.")
  digest::digest(bytes, algo = "xxhash64", serialize = FALSE)
}

fast_daily_signature <- function(task) {
  list(lon_indices = as.integer(task$lon_indices),
       lat_indices = as.integer(task$lat_indices), keep = as.integer(task$keep),
       cell_ids = as.integer(task$tile_points$cell_id))
}

fast_daily_flag_pack <- function(flags) {
  flags <- as.matrix(flags)
  if (is.null(colnames(flags)) || anyDuplicated(colnames(flags)) ||
      ncol(flags) < 4L || ncol(flags) > 8L ||
      !all(c("finite", "marginal", "joint", "primary") %in% colnames(flags)) ||
      anyNA(flags) || !all(flags %in% c(FALSE, TRUE)))
    fast_daily_fail("Daily cache flags must have 4-8 named Boolean columns and no missing values.")
  storage.mode(flags) <- "logical"
  if (any(flags[, "primary"] & (!flags[, "finite"] | !flags[, "marginal"])) ||
      any((flags[, "marginal"] | flags[, "joint"]) & !flags[, "finite"]))
    fast_daily_fail("Inconsistent daily finite/marginal/joint/primary flags.")
  code <- integer(nrow(flags))
  for (j in seq_len(ncol(flags))) code <- code + as.integer(flags[, j]) * 2L^(j - 1L)
  as.raw(code)
}

fast_daily_flag_unpack <- function(bytes, flag_names) {
  if (!is.raw(bytes) || !is.character(flag_names) || length(flag_names) < 4L ||
      length(flag_names) > 8L || anyDuplicated(flag_names) || anyNA(flag_names) ||
      !all(c("finite", "marginal", "joint", "primary") %in% flag_names))
    fast_daily_fail("Invalid daily cache quality-mask encoding.")
  values <- as.integer(bytes)
  if (any(values >= 2^length(flag_names))) fast_daily_fail("Unknown bits in the daily quality mask.")
  result <- matrix(FALSE, length(bytes), length(flag_names), dimnames = list(NULL, flag_names))
  for (j in seq_along(flag_names)) result[, j] <- bitwAnd(values, bitwShiftL(1L, j - 1L)) != 0L
  result
}

fast_daily_dates <- function(x, label = "dates") {
  result <- tryCatch(as.Date(x), error = function(e) NULL)
  if (is.null(result) || !length(result) || anyNA(result) ||
      any(diff(as.numeric(result)) <= 0))
    fast_daily_fail("Daily cache ", label, " must be nonempty, unique, ascending calendar dates.")
  result
}

fast_daily_open <- function(task, bundle, fingerprint, cache_dir = bundle$daily_cache$dir) {
  if (!is.character(fingerprint) || length(fingerprint) != 1L || is.na(fingerprint) || !nzchar(fingerprint))
    fast_daily_fail("Daily cache needs one nonempty run fingerprint.")
  if (!requireNamespace("digest", quietly = TRUE)) fast_daily_fail("Daily caching requires the digest package.")
  if (!is.character(cache_dir) || length(cache_dir) != 1L || is.na(cache_dir) || !nzchar(cache_dir))
    fast_daily_fail("Daily cache directory is missing.")
  if (!fast_daily_scalar(task$tile_id, 1, TRUE) || !is.data.frame(task$tile_points) ||
      nrow(task$tile_points) < 1L ||
      !all(c("cell_id", "lon_plot", "lat", "area_weight") %in% names(task$tile_points)))
    fast_daily_fail("Daily cache spatial tile is invalid.")
  points <- task$tile_points
  if (anyNA(points$cell_id) || anyDuplicated(points$cell_id) ||
      any(!is.finite(points$lon_plot)) || any(!is.finite(points$lat)) ||
      any(!is.finite(points$area_weight) | points$area_weight <= 0))
    fast_daily_fail("Daily cache cell identifiers, coordinates or weights are invalid.")
  ndraw <- as.integer(bundle$ndraw)
  cache <- bundle$cache
  if (!fast_daily_scalar(bundle$ndraw, 1, TRUE) ||
      length(cache$A0) != ndraw || length(cache$B0) != ndraw ||
      any(!is.finite(cache$A0) | cache$A0 <= 0) || any(!is.finite(cache$B0)))
    fast_daily_fail("Daily cache needs the positive-slope M0 baseline for each retained draw.")
  ids <- data.frame(position = seq_len(ndraw), m1 = cache$draw_ids_m1, m0 = cache$draw_ids_m0)
  if (nrow(ids) != ndraw || anyNA(ids) || anyDuplicated(ids$m1) || anyDuplicated(ids$m0) ||
      length(cache$draw_ids_m1) != ndraw || length(cache$draw_ids_m0) != ndraw)
    fast_daily_fail("Daily cache needs unique, complete M1/M0 draw identifiers.")
  tolerance <- bundle$daily_cache$max_abs_error
  if (is.null(tolerance)) tolerance <- 5e-5
  if (!fast_daily_scalar(tolerance, 0)) fast_daily_fail("Invalid daily float32 error tolerance.")
  if (!dir.exists(cache_dir) && !dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE))
    fast_daily_fail("Cannot create daily cache directory: ", cache_dir)
  stem <- sprintf("tile_%05d", as.integer(task$tile_id))
  data_file <- paste0(stem, ".f32"); index_file <- paste0(stem, ".rds")
  final_data <- file.path(cache_dir, data_file); final_index <- file.path(cache_dir, index_file)
  if (file.exists(final_index) || file.exists(final_data))
    fast_daily_fail("Daily cache already exists for tile ", task$tile_id,
      ". Validate/reuse its committed output; do not overwrite cache files. Directory: ", cache_dir)
  expected_dates <- NULL
  if (length(bundle$files)) {
    date_parts <- lapply(bundle$files, function(x) as.character(x$dates))
    if (any(vapply(date_parts, length, integer(1)) == 0L))
      fast_daily_fail("Missing planned dates in the daily cache file manifest.")
    expected_dates <- fast_daily_dates(unlist(date_parts, use.names = FALSE), "planned dates")
  }
  state <- new.env(parent = emptyenv())
  state$schema <- "ERA5-warm02-daily-f32-v1"
  state$fingerprint <- fingerprint; state$tile_id <- as.integer(task$tile_id)
  state$task_signature <- fast_daily_signature(task); state$cells <- points
  state$ncells <- nrow(points); state$ndraw <- ndraw; state$paired_draw_ids <- ids
  state$U_M0 <- as.numeric((cache$B0 + 1) / cache$A0)
  if (any(!is.finite(state$U_M0))) fast_daily_fail("Nonfinite daily M0 boundary.")
  state$cache_dir <- cache_dir; state$data_file <- data_file; state$index_file <- index_file
  state$temp_data <- tempfile(paste0(stem, "_partial_"), tmpdir = cache_dir, fileext = ".f32")
  state$temp_index <- tempfile(paste0(stem, "_partial_"), tmpdir = cache_dir, fileext = ".rds")
  state$con <- file(state$temp_data, open = "wb")
  state$closed <- FALSE; state$finished <- FALSE
  state$chunks <- list(); state$flag_names <- NULL; state$dates <- as.Date(character())
  state$expected_dates <- expected_dates; state$nrows <- 0; state$nbytes <- 0
  state$max_abs_error <- 0; state$error_tolerance <- tolerance; state$elapsed_write_seconds <- 0
  state$created <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  state
}

fast_daily_cleanup <- function(state) {
  if (!is.environment(state)) return(invisible(NULL))
  if (!isTRUE(state$closed) && inherits(state$con, "connection")) {
    try(close(state$con), silent = TRUE); state$closed <- TRUE
  }
  # Only this writer's temporary files are eligible for cleanup. Published data
  # are preserved even if publication of the index or tile checkpoint failed.
  for (path in c(state$temp_data, state$temp_index))
    if (length(path) == 1L && file.exists(path)) try(unlink(path), silent = TRUE)
  invisible(NULL)
}

fast_daily_append <- function(state, dates, flags, delta_U = NULL, temperature = NULL) {
  if (!is.environment(state) || isTRUE(state$closed) || isTRUE(state$finished))
    fast_daily_fail("Daily cache writer is not open.")
  started <- proc.time()[["elapsed"]]
  dates <- fast_daily_dates(dates)
  if (length(state$dates) && dates[1L] <= tail(state$dates, 1L))
    fast_daily_fail("Daily cache dates overlap or are out of order.")
  if (!is.null(state$expected_dates)) {
    expected <- state$expected_dates[length(state$dates) + seq_along(dates)]
    if (anyNA(expected) || !identical(as.numeric(dates), as.numeric(expected)))
      fast_daily_fail("Daily cache chunk dates do not match the input date manifest.")
  }
  flags <- as.matrix(flags)
  if (nrow(flags) != state$ncells * length(dates))
    fast_daily_fail("Daily cache flags must use cell-fastest, date-second row order.")
  packed <- fast_daily_flag_pack(flags)
  if (is.null(state$flag_names)) state$flag_names <- colnames(flags)
  if (!identical(colnames(flags), state$flag_names)) fast_daily_fail("Daily flag columns changed between chunks.")
  if (!is.null(temperature) && (!is.numeric(temperature) || !is.matrix(temperature) ||
      !identical(as.integer(dim(temperature)), c(as.integer(state$ncells), as.integer(length(dates)))) ||
      any(is.infinite(temperature) | is.nan(temperature))))
    fast_daily_fail("Daily temperature must be a cell-by-date numeric Celsius matrix; NA is allowed.")
  has_temperature <- !is.null(temperature)
  if (is.null(state$temperature_available)) state$temperature_available <- has_temperature
  if (!identical(state$temperature_available, has_temperature))
    fast_daily_fail("Daily temperature availability changed between chunks.")
  selected <- which(as.logical(flags[, "primary"]))
  nr <- length(selected); nd <- state$ndraw
  if (is.null(delta_U)) {
    if (nr) fast_daily_fail("Supported dates require their daily posterior delta_U matrix.")
    delta_U <- matrix(numeric(), nrow = 0L, ncol = nd)
  }
  if (!is.numeric(delta_U) || !is.matrix(delta_U) ||
      !identical(as.integer(dim(delta_U)), c(as.integer(nr), nd)) || any(!is.finite(delta_U)))
    fast_daily_fail("Daily delta_U must contain finite supported-row by paired-draw values.")
  bytes_expected <- as.double(nr) * nd * 4
  if (bytes_expected > .Machine$integer.max)
    fast_daily_fail("Daily cache chunk exceeds the raw-vector limit; reduce global$day_chunk.")
  bytes <- writeBin(as.numeric(delta_U), raw(), size = 4L, endian = "little")
  if (length(bytes) != bytes_expected) fast_daily_fail("Unexpected daily float32 encoding length.")
  decoded <- readBin(bytes, what = "numeric", n = length(delta_U), size = 4L, endian = "little")
  error <- if (length(decoded)) max(abs(decoded - as.numeric(delta_U))) else 0
  if (any(!is.finite(decoded)) || !is.finite(error) || error > state$error_tolerance)
    fast_daily_fail("Daily float32 round-trip error exceeds tolerance: ", error,
      " > ", state$error_tolerance, ". Preserve precision; do not silently accept this cache.")
  hash <- fast_daily_digest(bytes)
  if (length(bytes)) writeBin(bytes, state$con)
  chunk <- list(dates = dates, offset_bytes = state$nbytes, nbytes = bytes_expected,
    nrows = nr, supported_row_indices = as.integer(selected), flags = packed,
    hash_xxhash64 = hash, max_abs_error = error, temperature = temperature)
  state$chunks[[length(state$chunks) + 1L]] <- chunk
  state$dates <- c(state$dates, dates); state$nrows <- state$nrows + nr
  state$nbytes <- state$nbytes + bytes_expected
  state$max_abs_error <- max(state$max_abs_error, error)
  state$elapsed_write_seconds <- state$elapsed_write_seconds + proc.time()[["elapsed"]] - started
  invisible(chunk)
}

fast_daily_finish <- function(state) {
  if (!is.environment(state) || isTRUE(state$closed) || isTRUE(state$finished))
    fast_daily_fail("Daily cache writer is not open for finalization.")
  if (!length(state$chunks)) fast_daily_fail("A daily cache must record dates, including unsupported dates.")
  if (!is.null(state$expected_dates) &&
      !identical(as.numeric(state$dates), as.numeric(state$expected_dates)))
    fast_daily_fail("Daily cache is incomplete: expected dates are missing.")
  started <- proc.time()[["elapsed"]]
  flush(state$con); close(state$con); state$closed <- TRUE
  actual_bytes <- unname(file.info(state$temp_data)$size)
  if (length(actual_bytes) != 1L || is.na(actual_bytes) || actual_bytes != state$nbytes)
    fast_daily_fail("Daily cache file size differs from its chunk manifest.")
  index <- list(schema = state$schema, fingerprint = state$fingerprint,
    tile_id = state$tile_id, task_signature = state$task_signature,
    data_file = state$data_file, index_file = state$index_file,
    storage = "IEEE-754 float32 little-endian", byte_order = "little", value_size = 4L,
    layout = "within each chunk: supported cell-days fastest, paired draws second (R column-major)",
    row_order = "cell fastest, calendar date second", quantity = "delta_U = U_M1 - U_M0",
    units = "degree_C_SET", boundary = "study-standardized mean TSV = +1; U = M + positive_N",
    ncells = state$ncells, ndraw = state$ndraw, cells = state$cells,
    paired_draw_ids = state$paired_draw_ids, U_M0 = state$U_M0,
    flag_names = state$flag_names, dates = state$dates, chunks = state$chunks,
    temperature_available = isTRUE(state$temperature_available), temperature_units = "degree_C",
    nrows = state$nrows, nbytes = state$nbytes, max_abs_error = state$max_abs_error,
    max_allowed_abs_error = state$error_tolerance, elapsed_write_seconds = state$elapsed_write_seconds,
    created = state$created,
    interpretation = paste("Only raw-history-supported rows have posterior values.",
      "Unsupported or unavailable dates are missing, not zero demand.",
      "Use the same paired_draw_ids at all cells and dates; reconstruct U_M1 as delta_U + U_M0.",
      "The cross-model pairing is the declared working product-posterior coupling."))
  if (!is.null(state$recovery_checkpoint)) index$recovery_checkpoint <- state$recovery_checkpoint
  index$metadata_sha256 <- digest::digest(index, algo = "sha256", serialize = TRUE)
  saveRDS(index, state$temp_index, compress = FALSE, version = 3L)
  index_md5 <- unname(tools::md5sum(state$temp_index))
  if (length(index_md5) != 1L || is.na(index_md5)) fast_daily_fail("Cannot checksum daily cache index.")
  final_data <- file.path(state$cache_dir, state$data_file)
  final_index <- file.path(state$cache_dir, state$index_file)
  if (file.exists(final_data) || file.exists(final_index))
    fast_daily_fail("Refusing to overwrite a daily cache during publication.")
  # Completion order is deliberate: binary first, index second. The enclosing
  # tile's numerical checkpoint must be published after this function succeeds.
  if (!file.rename(state$temp_data, final_data)) fast_daily_fail("Cannot publish daily binary cache: ", final_data)
  if (!file.rename(state$temp_index, final_index))
    fast_daily_fail("Daily binary published but index publication failed. Preserve the orphan binary: ", final_data)
  state$finished <- TRUE
  state$elapsed_write_seconds <- state$elapsed_write_seconds + proc.time()[["elapsed"]] - started
  list(schema = state$schema, fingerprint = state$fingerprint, tile_id = state$tile_id,
    index_file = state$index_file, data_file = state$data_file, index_md5 = index_md5,
    nrows = state$nrows, ndraw = state$ndraw, nbytes = state$nbytes,
    n_chunks = length(state$chunks), n_days = length(state$dates),
    max_abs_error = state$max_abs_error, elapsed_write_seconds = state$elapsed_write_seconds)
}

fast_daily_descriptor <- function(index, cache_dir) {
  if (!is.list(index) || !identical(index$schema, "ERA5-warm02-daily-f32-v1") ||
      !fast_daily_relative_file(index$index_file, "rds"))
    fast_daily_fail("Cannot reconstruct a descriptor from an invalid daily index.")
  path <- file.path(cache_dir, index$index_file)
  if (!file.exists(path)) fast_daily_fail("Daily index is absent while recovering its descriptor: ", path)
  hash <- unname(tools::md5sum(path))
  if (length(hash) != 1L || is.na(hash)) fast_daily_fail("Cannot checksum the daily cache recovery index.")
  list(schema = index$schema, fingerprint = index$fingerprint, tile_id = index$tile_id,
    index_file = index$index_file, data_file = index$data_file, index_md5 = hash,
    nrows = index$nrows, ndraw = index$ndraw, nbytes = index$nbytes,
    n_chunks = length(index$chunks), n_days = length(index$dates),
    max_abs_error = index$max_abs_error, elapsed_write_seconds = index$elapsed_write_seconds)
}

fast_daily_relative_file <- function(x, suffix) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x) &&
    identical(basename(x), x) && !grepl("[/\\\\:]", x) &&
    grepl(paste0("[.]", suffix, "$"), x)
}

fast_daily_validate <- function(descriptor, cache_dir, fingerprint, task = NULL,
                                ndraw = NULL, verify_payload = FALSE) {
  needed <- c("schema", "fingerprint", "tile_id", "index_file", "data_file", "index_md5",
              "nrows", "ndraw", "nbytes", "max_abs_error")
  if (!is.list(descriptor) || !all(needed %in% names(descriptor)) ||
      !identical(descriptor$schema, "ERA5-warm02-daily-f32-v1") ||
      !identical(descriptor$fingerprint, fingerprint) ||
      !fast_daily_relative_file(descriptor$index_file, "rds") ||
      !fast_daily_relative_file(descriptor$data_file, "f32"))
    fast_daily_fail("Invalid daily cache descriptor or run fingerprint.")
  index_path <- file.path(cache_dir, descriptor$index_file)
  data_path <- file.path(cache_dir, descriptor$data_file)
  if (!file.exists(index_path) || !file.exists(data_path)) fast_daily_fail("A committed daily cache file is missing.")
  md5 <- unname(tools::md5sum(index_path))
  if (!identical(md5, descriptor$index_md5)) fast_daily_fail("Daily cache index checksum mismatch: ", index_path)
  index <- readRDS(index_path)
  fields <- c("schema", "fingerprint", "tile_id", "task_signature", "data_file", "index_file", "byte_order",
    "value_size", "ncells", "ndraw", "cells", "paired_draw_ids", "U_M0", "flag_names", "dates", "chunks",
    "nrows", "nbytes", "max_abs_error", "max_allowed_abs_error", "metadata_sha256", "temperature_available")
  if (!is.list(index) || !all(fields %in% names(index))) fast_daily_fail("Daily cache index fields are missing.")
  metadata_hash <- index$metadata_sha256; metadata_copy <- index; metadata_copy$metadata_sha256 <- NULL
  if (!is.character(metadata_hash) || length(metadata_hash) != 1L || is.na(metadata_hash) ||
      !identical(metadata_hash, digest::digest(metadata_copy, algo = "sha256", serialize = TRUE)))
    fast_daily_fail("Daily cache metadata checksum mismatch, including its recovery checkpoint.")
  for (nm in c("schema", "fingerprint", "tile_id", "index_file", "data_file", "nrows", "ndraw", "nbytes", "max_abs_error"))
    if (!identical(index[[nm]], descriptor[[nm]])) fast_daily_fail("Daily descriptor/index mismatch: ", nm)
  if (!identical(index$byte_order, "little") || !identical(index$value_size, 4L) ||
      !is.logical(index$temperature_available) || length(index$temperature_available) != 1L ||
      is.na(index$temperature_available) ||
      !fast_daily_scalar(index$ncells, 1, TRUE) || !fast_daily_scalar(index$ndraw, 1, TRUE) ||
      !fast_daily_scalar(index$nrows, 0, TRUE) || !fast_daily_scalar(index$nbytes, 0, TRUE) ||
      !fast_daily_scalar(index$max_abs_error, 0) || !fast_daily_scalar(index$max_allowed_abs_error, 0) ||
      index$max_abs_error > index$max_allowed_abs_error ||
      index$nbytes != 4 * as.double(index$nrows) * index$ndraw ||
      unname(file.info(data_path)$size) != index$nbytes)
    fast_daily_fail("Daily cache precision, dimensions or file size are invalid.")
  if (!is.null(ndraw) && !identical(as.integer(index$ndraw), as.integer(ndraw)))
    fast_daily_fail("Daily cache draw count differs from the analysis.")
  if (!is.data.frame(index$cells) || nrow(index$cells) != index$ncells ||
      !all(c("cell_id", "lon_plot", "lat", "area_weight") %in% names(index$cells)) ||
      anyNA(index$cells$cell_id) || anyDuplicated(index$cells$cell_id) ||
      any(!is.finite(index$cells$lon_plot)) || any(!is.finite(index$cells$lat)) ||
      any(!is.finite(index$cells$area_weight) | index$cells$area_weight <= 0))
    fast_daily_fail("Invalid daily cache cell metadata.")
  if (!is.null(task) && (!identical(index$task_signature, fast_daily_signature(task)) ||
      !identical(index$cells, task$tile_points) || as.integer(index$tile_id) != as.integer(task$tile_id)))
    fast_daily_fail("Daily cache spatial tile differs from the analysis.")
  ids <- index$paired_draw_ids
  if (!is.data.frame(ids) || !all(c("position", "m1", "m0") %in% names(ids)) ||
      nrow(ids) != index$ndraw || anyNA(ids) ||
      !identical(as.integer(ids$position), seq_len(index$ndraw)) ||
      anyDuplicated(ids$m1) || anyDuplicated(ids$m0) || length(index$U_M0) != index$ndraw ||
      any(!is.finite(index$U_M0))) fast_daily_fail("Invalid daily cache draw pairing or M0 baseline.")
  dates <- fast_daily_dates(index$dates, "indexed dates")
  if (!length(index$chunks)) fast_daily_fail("Daily cache has no date chunks.")
  offset <- 0; rows <- 0; chunk_dates <- as.Date(character()); observed_error <- 0
  for (j in seq_along(index$chunks)) {
    ch <- index$chunks[[j]]
    if (!is.list(ch) || !all(c("dates", "offset_bytes", "nbytes", "nrows", "supported_row_indices",
         "flags", "hash_xxhash64", "max_abs_error") %in% names(ch))) fast_daily_fail("Malformed daily chunk index.")
    cd <- fast_daily_dates(ch$dates, "chunk dates")
    if (length(chunk_dates) && cd[1L] <= tail(chunk_dates, 1L)) fast_daily_fail("Overlapping daily chunk dates.")
    masks <- fast_daily_flag_unpack(ch$flags, index$flag_names)
    if (nrow(masks) != index$ncells * length(cd) ||
        !identical(as.integer(which(masks[, "primary"])), ch$supported_row_indices) ||
        !identical(fast_daily_flag_pack(masks), ch$flags) ||
        !fast_daily_scalar(ch$nrows, 0, TRUE) || ch$nrows != length(ch$supported_row_indices) ||
        !fast_daily_scalar(ch$offset_bytes, 0, TRUE) || ch$offset_bytes != offset ||
        !fast_daily_scalar(ch$nbytes, 0, TRUE) || ch$nbytes != 4 * as.double(ch$nrows) * index$ndraw ||
        ch$nbytes > .Machine$integer.max ||
        !is.character(ch$hash_xxhash64) || length(ch$hash_xxhash64) != 1L ||
        is.na(ch$hash_xxhash64) || !grepl("^[0-9a-f]{16}$", ch$hash_xxhash64) ||
        !fast_daily_scalar(ch$max_abs_error, 0) || ch$max_abs_error > index$max_allowed_abs_error)
      fast_daily_fail("Inconsistent daily chunk rows, flags, offsets, digest or precision.")
    if (!identical(!is.null(ch$temperature), isTRUE(index$temperature_available)) ||
        (!is.null(ch$temperature) && (!is.numeric(ch$temperature) || !is.matrix(ch$temperature) ||
          !identical(as.integer(dim(ch$temperature)), c(as.integer(index$ncells), as.integer(length(cd)))) ||
          any(is.infinite(ch$temperature) | is.nan(ch$temperature)))))
      fast_daily_fail("Daily cache temperature availability, dimensions or values differ from its manifest.")
    offset <- offset + ch$nbytes; rows <- rows + ch$nrows
    chunk_dates <- c(chunk_dates, cd); observed_error <- max(observed_error, ch$max_abs_error)
    if (isTRUE(verify_payload)) fast_daily_read_chunk(index, j, cache_dir, verify = TRUE)
  }
  if (offset != index$nbytes || rows != index$nrows || observed_error != index$max_abs_error ||
      !identical(as.numeric(chunk_dates), as.numeric(dates)))
    fast_daily_fail("Daily chunk totals or dates differ from the index.")
  index
}

fast_daily_read_chunk <- function(index, chunk_i, cache_dir, verify = TRUE) {
  if (!is.list(index) || !identical(index$schema, "ERA5-warm02-daily-f32-v1") ||
      !fast_daily_scalar(chunk_i, 1, TRUE) || chunk_i > length(index$chunks) ||
      !fast_daily_relative_file(index$data_file, "f32")) fast_daily_fail("Invalid daily cache read request.")
  ch <- index$chunks[[chunk_i]]
  if (!fast_daily_scalar(ch$nbytes, 0, TRUE) || ch$nbytes > .Machine$integer.max ||
      !fast_daily_scalar(ch$offset_bytes, 0, TRUE)) fast_daily_fail("Invalid daily cache byte range.")
  path <- file.path(cache_dir, index$data_file)
  con <- file(path, open = "rb"); on.exit(close(con), add = TRUE)
  seek(con, where = ch$offset_bytes, origin = "start", rw = "read")
  bytes <- if (ch$nbytes > 0) readBin(con, what = "raw", n = as.integer(ch$nbytes)) else raw()
  if (length(bytes) != ch$nbytes) fast_daily_fail("Truncated daily cache payload: ", path)
  if (isTRUE(verify) && !identical(fast_daily_digest(bytes), ch$hash_xxhash64))
    fast_daily_fail("Daily cache chunk checksum mismatch: ", path, ", chunk ", chunk_i)
  values <- readBin(bytes, what = "numeric", n = as.double(ch$nrows) * index$ndraw, size = 4L, endian = "little")
  if (length(values) != as.double(ch$nrows) * index$ndraw || any(!is.finite(values)))
    fast_daily_fail("Daily cache payload has invalid posterior values.")
  delta_U <- matrix(values, nrow = ch$nrows, ncol = index$ndraw)
  rows <- ch$supported_row_indices
  row_cell <- ((rows - 1L) %% index$ncells) + 1L
  row_day <- ((rows - 1L) %/% index$ncells) + 1L
  flags <- fast_daily_flag_unpack(ch$flags, index$flag_names)
  temperature_supported <- if (is.null(ch$temperature)) NULL else as.vector(ch$temperature)[rows]
  list(dates = as.Date(ch$dates), cells = index$cells, delta_U = delta_U,
    U_M0 = index$U_M0, U0 = index$U_M0, paired_draw_ids = index$paired_draw_ids,
    flags = flags, quality = flags, quality_packed = ch$flags,
    supported_row_indices = rows, supported_rows = rows, row_cell_index = row_cell, row_day_index = row_day,
    row_cell_id = index$cells$cell_id[row_cell], row_date = as.Date(ch$dates)[row_day],
    cell_index = row_cell, row_dates = as.Date(ch$dates)[row_day],
    temperature_full = ch$temperature, temperature_supported = temperature_supported,
    units = index$units, max_abs_error = ch$max_abs_error)
}

fast_daily_selftest <- function() {
  folder <- tempfile("warm02_daily_selftest_"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  task <- list(tile_id = 1L, lon_indices = 1:2, lat_indices = 1L, keep = 1:2,
    tile_points = data.frame(cell_id = c(11L, 12L), lon_plot = c(0, 1), lat = c(30, 30),
      land_fraction = c(1, 1), area_weight = c(1, 1)))
  dates <- as.Date("2010-01-01") + 0:2
  bundle <- list(ndraw = 3L, cache = list(A0 = c(.1, .2, .3), B0 = c(2, 3, 4),
    draw_ids_m1 = c(7L, 8L, 9L), draw_ids_m0 = c(17L, 18L, 19L)),
    files = list(list(dates = dates)), daily_cache = list(dir = folder, max_abs_error = 5e-5))
  state <- fast_daily_open(task, bundle, "synthetic-daily-v1", folder)
  on.exit(fast_daily_cleanup(state), add = TRUE)
  flags <- cbind(finite = rep(TRUE, 4), marginal = c(TRUE, FALSE, TRUE, TRUE),
    joint = c(TRUE, FALSE, TRUE, TRUE), primary = c(TRUE, FALSE, TRUE, TRUE))
  truth <- matrix(c(.123456789, -.456789123, 1.98765432, 2.12345678,
    -.987654321, 1.00000001, -.1, .7, .9), 3L, 3L)
  temp <- matrix(c(20, 21, NA, 23), nrow = 2L)
  fast_daily_append(state, dates[1:2], flags, truth, temperature = temp)
  none <- cbind(finite = c(TRUE, FALSE), marginal = FALSE, joint = FALSE, primary = FALSE)
  fast_daily_append(state, dates[3], none, temperature = matrix(c(24, NA), nrow = 2L))
  state$recovery_checkpoint <- list(schema = "synthetic-recovery", value = 42L)
  descriptor <- fast_daily_finish(state)
  index <- fast_daily_validate(descriptor, folder, "synthetic-daily-v1", task, 3L, TRUE)
  got <- fast_daily_read_chunk(index, 1L, folder)
  stopifnot(identical(dim(got$delta_U), c(3L, 3L)), max(abs(got$delta_U - truth)) <= 5e-5,
    descriptor$nbytes == 36, descriptor$nrows == 3, descriptor$n_days == 3,
    identical(as.integer(got$row_cell_id), c(11L, 11L, 12L)),
    identical(as.numeric(got$row_date), as.numeric(dates[c(1, 2, 2)])),
    identical(got$flags, flags), identical(index$paired_draw_ids$m1, c(7L, 8L, 9L)),
    identical(got$temperature_full, temp), identical(got$temperature_supported, c(20, NA, 23)),
    identical(index$recovery_checkpoint$value, 42L))
  recovered <- fast_daily_descriptor(index, folder)
  stopifnot(identical(recovered$index_md5, descriptor$index_md5))
  empty <- fast_daily_read_chunk(index, 2L, folder)
  stopifnot(identical(dim(empty$delta_U), c(0L, 3L)), !length(empty$row_date), identical(empty$flags, none))
  failure <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = identity), "error")
  stopifnot(failure(fast_daily_validate(descriptor, folder, "wrong-fingerprint")),
    failure(fast_daily_open(task, bundle, "synthetic-daily-v1", folder)))
  # Deliberately corrupt one byte only within the disposable test cache.
  path <- file.path(folder, descriptor$data_file)
  con <- file(path, "r+b"); original <- readBin(con, "raw", 1L)
  seek(con, where = 0, origin = "start", rw = "write")
  writeBin(as.raw(bitwXor(as.integer(original), 1L)), con); close(con)
  stopifnot(failure(fast_daily_read_chunk(index, 1L, folder, verify = TRUE)))
  # A tile with zero supported dates still has a valid zero-byte payload and masks.
  task$tile_id <- 2L
  state2 <- fast_daily_open(task, bundle, "synthetic-daily-v1", folder)
  on.exit(fast_daily_cleanup(state2), add = TRUE)
  zero <- matrix(FALSE, 6L, 4L, dimnames = list(NULL, c("finite", "marginal", "joint", "primary")))
  fast_daily_append(state2, dates, zero)
  d2 <- fast_daily_finish(state2)
  i2 <- fast_daily_validate(d2, folder, "synthetic-daily-v1", task, 3L, TRUE)
  stopifnot(d2$nbytes == 0, d2$nrows == 0, length(i2$dates) == 3L)
  list(passed = TRUE, checks = c("float32 precision", "posterior pairing", "cell/date mapping",
    "byte counts", "zero-supported chunk", "zero-supported tile", "fingerprint rejection",
    "immutable publication", "payload corruption detection", "daily temperature alignment",
    "recoverable tile checkpoint"), max_abs_error = descriptor$max_abs_error)
}


# Daily outdoor-temperature input for ERA5-warm02 energy-ready caches.
# The raw monthly t2m files contain absolute temperatures; Kelvin conversion is
# therefore -273.15 (unlike Kelvin-valued history differences).

fast_daily_temperature_dates <- function(values, units, calendar, month) {
  if (length(calendar) != 1L || is.na(calendar) || !nzchar(calendar)) calendar <- "standard"
  if (!tolower(calendar) %in% c("standard", "gregorian", "proleptic_gregorian"))
    stop("Unsupported temperature-file calendar: ", calendar)
  if (length(units) != 1L || is.na(units)) stop("Temperature CF time units are missing.")
  hit <- regmatches(units, regexec("^(seconds?|minutes?|hours?|days?)\\s+since\\s+(.+)$",
                                  units, ignore.case = TRUE))[[1L]]
  if (length(hit) != 3L) stop("Unsupported temperature CF time units: ", units)
  origin <- sub("^([0-9]{4}-[0-9]{2}-[0-9]{2})[Tt]", "\\1 ", trimws(hit[3L]))
  origin <- sub("(Z| UTC| GMT| ?[+]00:00)$", "", origin)
  if (grepl("[+-][0-9]{2}:?[0-9]{2}$", origin))
    stop("Temperature time origin must be UTC.")
  if (nchar(origin) == 10L) origin <- paste(origin, "00:00:00")
  epoch <- as.POSIXct(origin, tz = "UTC",
                     tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M"))
  if (is.na(epoch)) stop("Cannot decode temperature CF time origin: ", units)
  multiplier <- switch(sub("s$", "", tolower(hit[2L])),
                       second = 1, minute = 60, hour = 3600, day = 86400)
  if (any(!is.finite(values))) stop("Non-finite temperature time coordinate.")
  timestamps <- epoch + as.numeric(values) * multiplier
  dates <- as.Date(timestamps, tz = "UTC")
  first_day <- as.Date(paste0(month, "-01"))
  next_month <- seq(first_day, by = "month", length.out = 2L)[2L]
  expected <- seq(first_day, next_month - 1, by = "day")
  if (anyDuplicated(dates) || any(diff(as.numeric(timestamps)) <= 0) ||
      !identical(as.character(dates), as.character(expected)))
    stop("Raw temperature file does not contain exactly one ordered daily value for every day of ", month, ".")
  time_of_day <- as.numeric(timestamps) %% 86400
  if (max(abs(time_of_day - time_of_day[1L])) > 1e-4)
    stop("Temperature time coordinates do not use a consistent UTC time of day.")
  list(dates = dates, calendar = calendar, time_units = units,
       time_of_day_seconds = time_of_day[1L])
}

fast_daily_temperature_offset <- function(units) {
  if (length(units) != 1L || is.na(units) || !nzchar(units))
    stop("Raw absolute t2m units are missing.")
  key <- tolower(gsub("[[:space:]_]", "", units))
  if (key %in% c("k", "kelvin", "degreekelvin", "degreeskelvin")) return(-273.15)
  if (key %in% c("c", "degc", "degreecelsius", "degreescelsius", "celsius", "\u00b0c")) return(0)
  stop("Unrecognized absolute t2m units: ", units)
}

fast_daily_coordinate_index <- function(wanted, available, axis) {
  wanted <- as.numeric(wanted); available <- as.numeric(available)
  if (any(!is.finite(wanted)) || any(!is.finite(available)))
    stop("Non-finite ", axis, " coordinate in temperature alignment.")
  if (axis == "lon") {
    wanted <- ((wanted + 180) %% 360) - 180
    available <- ((available + 180) %% 360) - 180
  }
  if (anyDuplicated(round(available, 8L)))
    stop("Duplicate ", axis, " coordinate after normalization in temperature input.")
  # Exact grid matching within coordinate precision, with no interpolation or
  # nearest-neighbour substitution. Order and longitude convention may differ.
  ix <- match(round(wanted, 8L), round(available, 8L))
  if (anyNA(ix) || any(abs(available[ix] - wanted) > 1e-7))
    stop("Raw temperature grid does not contain every requested ", axis,
         " coordinate. Regrid explicitly before running; index substitution is prohibited.")
  as.integer(ix)
}

fast_daily_inspect_temperature <- function(config, inspected_history) {
  settings <- config$daily_cache
  mode <- if (is.null(settings$temperature_mode)) "auto" else settings$temperature_mode
  if (length(mode) != 1L || is.na(mode) || !mode %in% c("auto", "required", "off"))
    stop("daily_cache$temperature_mode must be 'auto', 'required', or 'off'.")
  absent <- function(reason, missing_files = character()) {
    if (mode == "required") stop(reason)
    if (mode != "off") warning(reason, " Daily warm-boundary caches will still be produced; ",
      "outdoor temperature is NOT embedded. Later energy calculations need the raw daily t2m files.",
      call. = FALSE)
    list(status = "not_embedded", mode = mode, files = list(), audit = data.frame(),
         reason = reason, missing_files = missing_files,
         identity_method = "no_temperature_embedded")
  }
  if (mode == "off") return(absent("Embedding of daily outdoor temperature was explicitly disabled."))
  if (!length(inspected_history$files)) stop("History input inspection must precede temperature inspection.")
  history_dates <- as.Date(unlist(lapply(inspected_history$files, function(x) as.character(x$dates))))
  if (anyNA(history_dates) || anyDuplicated(history_dates)) stop("Invalid or duplicate history dates.")
  months <- sort(unique(format(history_dates, "%Y-%m")))
  directory <- settings$temperature_dir
  if (is.null(directory) || length(directory) != 1L || is.na(directory) || !nzchar(directory))
    return(absent("No daily_cache$temperature_dir was configured."))
  paths <- file.path(directory, paste0("ERA5_t2m_daily_mean_", gsub("-", "_", months), ".nc"))
  names(paths) <- months
  present <- file.exists(paths)
  if (!any(present)) return(absent(paste0("No expected monthly t2m files found in ", directory, "."), paths))
  checksum <- if (is.null(settings$temperature_checksum_inputs)) TRUE else settings$temperature_checksum_inputs
  if (length(checksum) != 1L || is.na(checksum) || !is.logical(checksum))
    stop("temperature_checksum_inputs must be TRUE or FALSE.")
  if (!checksum) warning("Raw temperature identity uses path, size, mtime and coordinate/date metadata only; ",
    "content checksums are disabled. Same-size content replacement preserving mtime cannot be detected.", call. = FALSE)
  first_history <- inspected_history$files[[1L]]
  files <- list(); audit <- list(); previous_clock <- NULL
  for (month in months[present]) {
    path <- paths[[month]]
    if (exists(".set_analysis_stage", mode = "function"))
      .set_analysis_stage("inspect_daily_outdoor_temperature", file = path, month = month)
    message("Checking daily outdoor temperature ", match(month, months), "/", length(months), ": ", basename(path))
    nc <- ncdf4::nc_open(path)
    item <- tryCatch({
      lon <- nc_axis(nc, "lon"); lat <- nc_axis(nc, "lat"); tm <- nc_axis(nc, "time")
      if (any(abs(lat$values) > 90) || any(lon$values < -180 | lon$values > 360))
        stop("Invalid raw-temperature spatial coordinate range: ", path)
      fast_daily_coordinate_index(first_history$lon, lon$values, "lon")
      fast_daily_coordinate_index(first_history$lat, lat$values, "lat")
      time <- fast_daily_temperature_dates(tm$values, tm$units,
                                           nc_attribute(nc, tm$name, "calendar"), month)
      variable <- settings$temperature_variable
      if (is.null(variable)) {
        candidates <- intersect(c("t2m", "2m_temperature", "temperature_2m"), names(nc$var))
        if (length(candidates) != 1L)
          stop("Expected one t2m variable in ", path, "; set daily_cache$temperature_variable explicitly.")
        variable <- candidates[[1L]]
      }
      if (length(variable) != 1L || is.na(variable) || !variable %in% names(nc$var))
        stop("Configured raw temperature variable is missing in ", path)
      dims <- nc$var[[variable]]$dim
      axes <- normal_axis(vapply(dims, function(d) d$name, character(1)))
      for (axis in c("lon", "lat", "time")) if (sum(axes == axis) != 1L)
        stop("Raw temperature variable lacks a unique ", axis, " dimension: ", path)
      extra <- which(!axes %in% c("lon", "lat", "time"))
      if (length(extra) && any(vapply(dims[extra], function(d) as.integer(d$len), integer(1)) != 1L))
        stop("Raw temperature has a non-singleton extra dimension; select it explicitly: ", path)
      expected_len <- c(lon = lon$len, lat = lat$len, time = tm$len)
      for (axis in names(expected_len)) if (dims[[which(axes == axis)]]$len != expected_len[[axis]])
        stop("Raw temperature dimension/coordinate length mismatch: ", path)
      units <- nc_attribute(nc, variable, "units")
      offset <- fast_daily_temperature_offset(units)
      tz_att <- nc_attribute(nc, 0, "daily_time_zone")
      if (is.na(tz_att[1L])) tz_att <- nc_attribute(nc, 0, "time_zone")
      if (length(tz_att) && !is.na(tz_att[1L]) &&
          !tolower(gsub("[[:space:]]", "", tz_att[1L])) %in% c("utc", "utc+00:00", "utc+0", "gmt", "utc+00", "0"))
        stop("Raw daily temperature declares a non-UTC aggregation timezone: ", path)
      list(path = normalizePath(path, winslash = "/", mustWork = TRUE), month = month,
           variable = variable, lon = lon$values, lat = lat$values, dates = time$dates,
           time_info = time, units = units, offset = offset,
           aggregation_timezone_attribute = if (length(tz_att)) as.character(tz_att[1L]) else NA_character_)
    }, finally = ncdf4::nc_close(nc))
    if (!is.null(previous_clock) && abs(item$time_info$time_of_day_seconds - previous_clock) > 1e-4)
      stop("Daily raw temperature files use different UTC clock times across months.")
    previous_clock <- item$time_info$time_of_day_seconds
    meta <- file.info(path)
    item$bytes <- as.numeric(meta$size)
    item$mtime <- as.numeric(meta$mtime)
    item$md5 <- if (checksum) unname(tools::md5sum(path)) else "not_computed"
    if (checksum && (length(item$md5) != 1L || is.na(item$md5)))
      stop("Failed to checksum raw temperature input: ", path)
    files[[month]] <- item
    audit[[month]] <- data.frame(file = item$path, month = month, variable = item$variable,
      n_dates = length(item$dates), first_date = as.character(min(item$dates)),
      last_date = as.character(max(item$dates)), nlon = length(item$lon), nlat = length(item$lat),
      units = item$units, offset_to_celsius = item$offset, calendar = item$time_info$calendar,
      time_units = item$time_info$time_units, bytes = item$bytes, modification_time = item$mtime,
      md5 = item$md5, aggregation_timezone_attribute = item$aggregation_timezone_attribute,
      stringsAsFactors = FALSE)
  }
  # An incomplete inventory is never silently treated as complete. Existing
  # malformed files have already failed above even in auto mode.
  if (!all(present)) return(absent(paste0("Incomplete raw daily-temperature inventory: ",
    sum(!present), " of ", length(paths), " monthly files are missing. Example: ",
    paste(head(paths[!present], 3L), collapse = "; ")), paths[!present]))
  all_dates <- as.Date(unlist(lapply(files, function(x) as.character(x$dates))))
  if (anyDuplicated(all_dates) || anyNA(match(history_dates, all_dates)))
    stop("Raw daily-temperature dates do not uniquely cover every history date.")
  message("Daily outdoor temperature is available for embedding: ", length(files), " monthly files; ",
          length(history_dates), " matched UTC dates; absolute temperatures converted to degrees C.")
  list(status = "embedded", mode = mode, files = files, audit = do.call(rbind, audit),
       reason = "Complete validated raw t2m inventory; UTC date and coordinate matching required for every read.",
       missing_files = character(), identity_method = if (checksum) "md5+path+size+mtime+coordinates+dates"
       else "path+size+mtime+coordinates+dates_without_content_checksum",
       timezone_provenance = "CF timestamps use UTC; daily aggregation timezone follows the supplied UTC downloader unless explicitly declared in NetCDF.")
}

fast_daily_temperature_chunk <- function(temp_spec, task, dates) {
  dates <- as.Date(dates)
  if (anyNA(dates) || anyDuplicated(dates)) stop("Temperature chunk dates must be valid and unique.")
  points <- task$tile_points
  if (!all(c("lon_plot", "lat") %in% names(points))) stop("Temperature chunk task lacks point coordinates.")
  out <- matrix(NA_real_, nrow = nrow(points), ncol = length(dates))
  if (!identical(temp_spec$status, "embedded") || !length(dates) || !nrow(points)) return(out)
  keys <- format(dates, "%Y-%m")
  for (month in unique(keys)) {
    item <- temp_spec$files[[month]]
    if (is.null(item)) stop("Temperature inventory lacks requested month: ", month)
    columns <- which(keys == month)
    day_indices <- match(as.character(dates[columns]), as.character(item$dates))
    if (anyNA(day_indices)) stop("Requested true date absent from raw temperature input: ", item$path)
    lon_indices <- fast_daily_coordinate_index(points$lon_plot, item$lon, "lon")
    lat_indices <- fast_daily_coordinate_index(points$lat, item$lat, "lat")
    meta <- file.info(item$path)
    if (!is.finite(meta$size) || meta$size != item$bytes || as.numeric(meta$mtime) != item$mtime)
      stop("Raw temperature file changed or disappeared after inspection: ", item$path)
    nc <- ncdf4::nc_open(item$path)
    block <- tryCatch({
      unique_lon <- sort(unique(lon_indices)); unique_lat <- sort(unique(lat_indices))
      # Split widely separated longitude indices, including a cyclic seam, so a
      # small tile cannot accidentally read almost the whole globe.
      lon_groups <- split(unique_lon, cumsum(c(1L, as.integer(diff(unique_lon) > 16L))))
      values <- matrix(NA_real_, nrow = nrow(points), ncol = length(columns))
      first_day <- min(day_indices); nday <- max(day_indices) - first_day + 1L
      for (selection in lon_groups) {
        rows <- which(lon_indices %in% selection)
        rect <- read_rectangular_history(nc, item$variable, selection, unique_lat,
                                         first_day, nday, offset = item$offset)
        rectangle_rows <- match(lon_indices[rows], selection) +
          (match(lat_indices[rows], unique_lat) - 1L) * length(selection)
        values[rows, ] <- rect[rectangle_rows, day_indices - first_day + 1L, drop = FALSE]
      }
      values[!is.finite(values)] <- NA_real_
      values
    }, finally = ncdf4::nc_close(nc))
    out[, columns] <- block
  }
  out
}

fast_daily_temperature_selftest <- function(tempdir) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("ncdf4 is required for the temperature selftest.")
  directory <- tempfile("daily_temperature_test_", tmpdir = tempdir)
  dir.create(directory, recursive = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  origin <- as.Date("2010-01-01")
  expected_value <- function(lon, lat, dates) (((lon + 180) %% 360) - 180) / 100 +
    lat / 10 + as.numeric(dates - origin)
  make_month <- function(month, lon, lat, units, extra = FALSE) {
    first_day <- as.Date(paste0(month, "-01"))
    dates <- seq(first_day, seq(first_day, by = "month", length.out = 2L)[2L] - 1, by = "day")
    dims <- list(ncdf4::ncdim_def("longitude", "degrees_east", lon),
                 ncdf4::ncdim_def("latitude", "degrees_north", lat),
                 ncdf4::ncdim_def("valid_time", "days since 2010-01-01 00:00:00", as.numeric(dates - origin)))
    if (extra) dims <- c(dims, list(ncdf4::ncdim_def("member", "1", 1:2)))
    path <- file.path(directory, paste0("ERA5_t2m_daily_mean_", gsub("-", "_", month), ".nc"))
    variable <- ncdf4::ncvar_def("t2m", units, dims, missval = -9999, prec = "double")
    nc <- ncdf4::nc_create(path, variable)
    tryCatch({
      lengths <- vapply(dims, function(d) as.integer(d$len), integer(1))
      x <- array(NA_real_, dim = lengths)
      for (k in seq_along(dates)) for (j in seq_along(lat)) for (i in seq_along(lon)) {
        v <- expected_value(lon[i], lat[j], dates[k]) - fast_daily_temperature_offset(units)
        if (extra) x[i, j, k, ] <- v else x[i, j, k] <- v
      }
      ncdf4::ncvar_put(nc, "t2m", x)
      ncdf4::ncatt_put(nc, "valid_time", "calendar", "standard")
      ncdf4::ncatt_put(nc, 0, "daily_time_zone", "UTC+00:00")
    }, finally = ncdf4::nc_close(nc))
    invisible(path)
  }
  make_month("2010-01", c(0, 90, 180, 270), c(20, 10), "K")
  make_month("2010-02", c(-180, -90, 0, 90), c(10, 20), "degree_Celsius")
  config <- list(daily_cache = list(temperature_mode = "required", temperature_dir = directory,
                                    temperature_checksum_inputs = TRUE))
  history <- list(files = list(list(lon = c(-180, -90, 0, 90), lat = c(10, 20),
                                    dates = seq(origin, as.Date("2010-02-28"), by = "day"))))
  spec <- fast_daily_inspect_temperature(config, history)
  task <- list(tile_points = data.frame(lon_plot = c(-90, 0, -180, 90), lat = c(10, 20, 10, 20)))
  dates <- as.Date(c("2010-01-31", "2010-02-01", "2010-02-02"))
  actual <- fast_daily_temperature_chunk(spec, task, dates)
  expected <- vapply(dates, function(d) expected_value(task$tile_points$lon_plot,
                                                     task$tile_points$lat, as.Date(d, origin = "1970-01-01")),
                     numeric(nrow(task$tile_points)))
  if (!identical(dim(actual), dim(expected)) || max(abs(actual - expected)) > 1e-10)
    stop("Daily temperature selftest failed: unit, date, coordinate or cross-month alignment.")
  single <- fast_daily_temperature_chunk(spec, list(tile_points = task$tile_points[1L, , drop = FALSE]), dates[2L])
  if (!identical(dim(single), c(1L, 1L)) || abs(single[1L] - expected[1L, 2L]) > 1e-10)
    stop("Daily temperature selftest failed: singleton matrix dimensions.")
  missing_nc <- ncdf4::nc_open(spec$files[["2010-01"]]$path, write = TRUE)
  tryCatch(ncdf4::ncvar_put(missing_nc, "t2m", NA_real_, start = c(1L, 1L, 31L),
                           count = c(1L, 1L, 1L)), finally = ncdf4::nc_close(missing_nc))
  missing_spec <- fast_daily_inspect_temperature(config, history)
  missing_chunk <- fast_daily_temperature_chunk(missing_spec, task, dates)
  if (!is.na(missing_chunk[2L, 1L]) || sum(is.na(missing_chunk)) != 1L)
    stop("Daily temperature selftest failed: missing input temperature was not preserved as NA.")
  if (!inherits(try(fast_daily_temperature_offset("m/s"), silent = TRUE), "try-error"))
    stop("Daily temperature selftest failed: invalid units accepted.")
  if (!inherits(try(fast_daily_coordinate_index(12.5, c(10, 20), "lat"), silent = TRUE), "try-error"))
    stop("Daily temperature selftest failed: nonmatching coordinate accepted.")
  if (!inherits(try(fast_daily_temperature_dates(c(0, 0), "days since 2010-01-01", "standard", "2010-01"),
                   silent = TRUE), "try-error"))
    stop("Daily temperature selftest failed: duplicate/incomplete dates accepted.")
  make_month("2010-02", c(-180, -90, 0, 90), c(10, 20), "degree_Celsius", extra = TRUE)
  if (!inherits(try(fast_daily_inspect_temperature(config, history), silent = TRUE), "try-error"))
    stop("Daily temperature selftest failed: extra non-singleton dimension accepted.")
  unlink(file.path(directory, "ERA5_t2m_daily_mean_2010_02.nc"))
  config$daily_cache$temperature_mode <- "auto"
  incomplete <- suppressWarnings(fast_daily_inspect_temperature(config, history))
  if (!identical(incomplete$status, "not_embedded") || length(incomplete$missing_files) != 1L)
    stop("Daily temperature selftest failed: incomplete inventory was called embedded.")
  message("Daily outdoor-temperature selftest passed: K/C, true dates across months, longitude wrap, reversed latitude, singleton reads and input rejection.")
  invisible(TRUE)
}


# The daily cache is a reusable result, not a new energy model.
fast_daily_defaults <- function(config) {
  d <- list(enabled = TRUE, max_abs_error = 5e-5,
    temperature_mode = "auto", temperature_dir = "E:/ERA5/01_raw_ERA5",
    temperature_checksum_inputs = TRUE, available_space_gb = 9000,
    reserve_space_gb = 500, verify_payload_on_resume = FALSE)
  if (!is.null(config$daily_cache)) d <- utils::modifyList(d, config$daily_cache)
  if (!isTRUE(d$enabled)) stop("ERA5-warm02 requires daily caching; use warm01 for summaries only.")
  if (!d$temperature_mode %in% c("auto", "required", "off")) stop("Invalid daily_cache$temperature_mode.")
  for (nm in c("max_abs_error", "available_space_gb", "reserve_space_gb"))
    if (length(d[[nm]]) != 1L || !is.finite(d[[nm]]) || d[[nm]] <= 0) stop("Invalid daily_cache$", nm)
  d
}

fast_daily_identity_settings <- function(d) {
  d[setdiff(names(d), c("available_space_gb", "reserve_space_gb", "verify_payload_on_resume"))]
}

fast_daily_live_free_gb <- function(out_dir) {
  # On Windows, query the OUTPUT volume rather than confusing disk capacity
  # with the space still available after ERA5 inputs and prior results.
  if (.Platform$OS.type != "windows") return(NA_real_)
  normalized <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  if (!grepl("^[A-Za-z]:/", normalized)) return(NA_real_)
  drive <- substr(normalized, 1L, 1L)
  command <- sprintf("(Get-PSDrive -Name '%s' -ErrorAction Stop).Free", drive)
  result <- tryCatch(suppressWarnings(system2("powershell.exe",
    c("-NoProfile", "-NonInteractive", "-Command", shQuote(command, type = "cmd")),
    stdout = TRUE, stderr = FALSE)), error = function(e) character())
  numeric <- suppressWarnings(as.double(trimws(result)))
  numeric <- numeric[is.finite(numeric) & numeric >= 0]
  if (length(numeric) != 1L) return(NA_real_)
  numeric / 1e9
}

fast_daily_support_preflight <- function(d, tasks, files, support_analysis, g,
                                         input_audit, ndraw, out_dir) {
  # The float32 writer stores exactly one row for each primary-supported
  # cell/date, regardless of the number of dates or cells in the tile.
  # Count those rows with the SAME reader and support rule as fast_compute_tile.
  checkpoint_dir <- file.path(out_dir, "daily_capacity_support_checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  starting_free_gb <- fast_daily_live_free_gb(out_dir)
  if (.Platform$OS.type == "windows" && !is.finite(starting_free_gb))
    stop("Could not measure free space on the output drive. Check drive availability ",
      "and PowerShell Get-PSDrive; do not use disk capacity as available_space_gb.", call. = FALSE)
  initial_budget_gb <- if (is.finite(starting_free_gb)) {
    min(d$available_space_gb, starting_free_gb)
  } else d$available_space_gb
  fingerprint <- fast_object_hash(list(schema = "ERA5-warm02-capacity-support-v1",
    files = files, input_audit = input_audit, support = support_analysis,
    joint_support = isTRUE(g$joint_support), day_chunk = g$day_chunk,
    tasks = lapply(tasks, fast_tile_signature),
    functions = lapply(c("read_rectangular_history", "evaluate_support", "history_support"),
      function(nm) deparse(body(get(nm, mode = "function"))))))
  counts <- numeric(length(tasks))
  for (ti in seq_along(tasks)) {
    task <- tasks[[ti]]
    checkpoint <- file.path(checkpoint_dir, sprintf("tile_%05d.rds", task$tile_id))
    old <- if (file.exists(checkpoint)) tryCatch(readRDS(checkpoint), error = function(e) NULL) else NULL
    valid <- is.list(old) && identical(old$schema, "ERA5-warm02-capacity-support-v1") &&
      identical(old$fingerprint, fingerprint) && identical(old$task_signature, fast_tile_signature(task)) &&
      length(old$supported_rows) == 1L && is.finite(old$supported_rows) &&
      old$supported_rows >= 0 && old$supported_rows <= nrow(task$tile_points) * sum(vapply(files,
        function(fi) length(fi$dates), integer(1)))
    if (valid) {
      counts[ti] <- old$supported_rows
    } else {
      ncs <- list(); supported <- 0
      tryCatch({
        for (yy in seq_along(files)) {
          fi <- files[[yy]]
          nc <- ncdf4::nc_open(fi$path)
          ncs[[yy]] <- nc
          nt <- length(fi$dates)
          for (ds in seq.int(1L, nt, by = g$day_chunk)) {
            nday <- min(g$day_chunk, nt - ds + 1L)
            mats <- lapply(era5_required_variables, function(v)
              read_rectangular_history(nc, v, task$lon_indices, task$lat_indices,
                ds, nday, fi$offsets[v])[task$keep, , drop = FALSE])
            history <- as.data.frame(setNames(lapply(mats, as.vector), era5_required_variables))
            support <- evaluate_support(support_analysis, history, joint = isTRUE(g$joint_support))
            if (!is.data.frame(support) || nrow(support) != nrow(history) ||
                !all(c("finite", "marginal_inside", "joint_inside") %in% names(support)))
              stop("Support evaluator returned invalid rows or columns.")
            finite <- as.logical(support$finite); finite[is.na(finite)] <- FALSE
            marginal <- finite & !is.na(support$marginal_inside) & support$marginal_inside
            joint <- finite & !is.na(support$joint_inside) & support$joint_inside
            primary <- marginal & (if (isTRUE(g$joint_support)) joint else TRUE)
            supported <- supported + sum(primary)
          }
          ncdf4::nc_close(nc); ncs[[yy]] <- NULL
        }
      }, finally = invisible(lapply(ncs, function(nc) try(ncdf4::nc_close(nc), silent = TRUE))))
      counts[ti] <- supported
      fast_atomic_save_rds(list(schema = "ERA5-warm02-capacity-support-v1",
        fingerprint = fingerprint, task_signature = fast_tile_signature(task),
        supported_rows = supported), checkpoint, replace = file.exists(checkpoint))
    }
    if (ti %% 100L == 0L || ti == length(tasks)) {
      message(sprintf("Capacity support scan: %d/%d tiles; %.1f GB paired draws counted so far.",
        ti, length(tasks), sum(counts[seq_len(ti)]) * ndraw * 4 / 1e9))
      partial_gb <- sum(counts[seq_len(ti)]) * ndraw * 4 / 1e9 * 1.10 + d$reserve_space_gb
      if (partial_gb > initial_budget_gb)
        stop(sprintf(paste("Exact support scan stopped after %d/%d tiles:",
          "the already-counted cache requires at least %.1f GB including margin and reserve,",
          "exceeding the measured/declared %.1f GB budget. Completed scan checkpoints remain reusable."),
          ti, length(tasks), partial_gb, initial_budget_gb), call. = FALSE)
    }
  }
  ncells <- sum(vapply(tasks, function(task) nrow(task$tile_points), integer(1)))
  ndates <- sum(vapply(files, function(fi) length(fi$dates), integer(1)))
  maximum_gb <- as.double(ncells) * ndates * ndraw * 4 / 1e9
  raw_gb <- sum(counts) * ndraw * 4 / 1e9
  need <- raw_gb * 1.10 + d$reserve_space_gb
  live_free_gb <- fast_daily_live_free_gb(out_dir)
  if (.Platform$OS.type == "windows" && !is.finite(live_free_gb))
    stop("Could not measure free space on the output drive. Check drive availability ",
      "and PowerShell Get-PSDrive; do not use disk capacity as available_space_gb.", call. = FALSE)
  budget_gb <- if (is.finite(live_free_gb)) min(d$available_space_gb, live_free_gb) else d$available_space_gb
  tab <- data.frame(land_points = ncells, dates = ndates, draws = ndraw,
    supported_cell_dates = sum(counts), bytes_per_value = 4,
    full_support_upper_bound_GB = maximum_gb, exact_delta_U_GB = raw_gb,
    planning_total_GB_including_margin = need,
    user_reported_available_GB = d$available_space_gb,
    measured_available_GB = live_free_gb, effective_available_GB = budget_gb,
    capacity_source = if (is.finite(live_free_gb)) "Minimum of reported cap and live output-volume free space"
      else "Declared free space; live measurement unavailable on this OS")
  utils::write.csv(tab, file.path(out_dir, "daily_cache_capacity_plan.csv"), row.names = FALSE)
  message(sprintf("Daily delta_U cache after exact support scan: %.1f GB binary; %.1f GB incl. 10%% margin and %.1f GB reserve; measured output-drive free space %.1f GB; effective budget %.1f GB.",
    raw_gb, need, d$reserve_space_gb, live_free_gb, budget_gb))
  if (need > budget_gb)
    stop(sprintf(paste("Exact support-based capacity check stopped safely:",
      "the plan needs %.1f GB including margin and reserve but effective free space is %.1f GB.",
      "Do not increase available_space_gb without increasing actual free space;",
      "a streaming CDD implementation would be needed for this drive."),
      need, budget_gb), call. = FALSE)
  list(table = tab, supported_rows_by_tile = counts)
}

fast_daily_recover_warm <- function(task, bundle, fingerprint, checkpoint_path, part_dir) {
  if (!isTRUE(bundle$daily_cache$enabled) || file.exists(checkpoint_path)) return(invisible(FALSE))
  cache_dir <- bundle$daily_cache$dir
  index_path <- file.path(cache_dir, sprintf("tile_%05d.rds", task$tile_id))
  binary_path <- file.path(cache_dir, sprintf("tile_%05d.f32", task$tile_id))
  if (file.exists(index_path)) {
    index <- readRDS(index_path)
    descriptor <- fast_daily_descriptor(index, cache_dir)
    fast_daily_validate(descriptor, cache_dir, fingerprint, task, bundle$ndraw,
      verify_payload = isTRUE(bundle$daily_cache$verify_payload_on_resume))
    cp <- index$recovery_checkpoint
    if (!is.list(cp) || !identical(cp$fingerprint, fingerprint) ||
        !identical(cp$task_signature, fast_tile_signature(task)))
      stop("Committed daily cache has no matching warm-checkpoint recovery record: ", index_path)
    cp$daily_cache <- descriptor
    fast_atomic_publish(checkpoint_path, function(temporary) {
      saveRDS(cp, temporary, compress = "gzip")
      fast_read_checkpoint(temporary, fingerprint, task, bundle)
    })
    fast_restore_partition(cp, part_dir)
    message("Recovered warm checkpoint from committed daily cache for tile ", task$tile_id, "; no model prediction.")
    return(invisible(TRUE))
  }
  if (file.exists(binary_path)) {
    # Only an uncommitted binary with NO index and NO checkpoint reaches here.
    orphan_dir <- file.path(cache_dir, "uncommitted_orphans")
    dir.create(orphan_dir, recursive = TRUE, showWarnings = FALSE)
    target <- tempfile(paste0(basename(binary_path), "_"), tmpdir = orphan_dir)
    if (!file.rename(binary_path, target)) stop("Cannot preserve uncommitted binary: ", binary_path)
    message("Preserved uncommitted cache binary: ", target, ". This unfinished tile will be recalculated.")
  }
  invisible(FALSE)
}

fast_daily_snapshot <- function(config, analysis, bundle, tasks, fingerprint, out_dir) {
  path <- bundle$daily_cache$dir
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  metadata <- list(schema = "ERA5-warm02-cache-manifest-v1", fingerprint = fingerprint,
    ndraw = bundle$ndraw, dates = do.call(c, lapply(bundle$files, function(x) as.Date(x$dates))),
    tile_ids = vapply(tasks, function(x) as.integer(x$tile_id), integer(1)),
    tile_cell_ids = setNames(lapply(tasks, function(x) as.integer(x$tile_points$cell_id)),
      vapply(tasks, function(x) as.character(x$tile_id), character(1))),
    grid = do.call(rbind, lapply(tasks, function(x) x$tile_points)),
    paired_draw_ids = data.frame(position = seq_len(bundle$ndraw),
      m1 = bundle$cache$draw_ids_m1, m0 = bundle$cache$draw_ids_m0),
    U_M0 = (bundle$cache$B0 + 1) / bundle$cache$A0,
    M_M0 = bundle$cache$B0 / bundle$cache$A0, N_M0 = 1 / bundle$cache$A0,
    A_M0 = bundle$cache$A0, B_M0 = bundle$cache$B0,
    temperature = bundle$temperature, cache_settings = fast_daily_identity_settings(bundle$daily_cache),
    SET_bounds = bundle$set_bounds, history_settings = bundle$g,
    storage = "Little-endian IEEE float32; each chunk is an R column-major row-by-draw matrix. Separate digest per chunk.",
    row_order = "Each date: cells in saved tile order. Only primary-supported rows are stored; all-date masks retained.",
    boundary = "Expected TSV=+1; U=M+N, M=B/A, N=1/A positive. delta_U=M1-M0.",
    inference = "Fixed working product-posterior M1/M0 coupling across every date and grid; not an identified paired-data joint posterior.",
    scope = "Supports future SET-to-Ta conversions and cooling scenarios at the cached TSV +1 boundary and the same support domain. Changing the model, response threshold or expanding support can require new predictions.",
    energy_status = "No SET-to-Ta relationship, energy threshold or observed electricity demand is inferred by this cache.",
    completion = "Manifest is an input contract, not evidence of completed tiles. Check cache_completion.rds and tile inventory.")
  target <- file.path(path, "cache_manifest.rds")
  if (file.exists(target)) {
    old <- readRDS(target)
    if (!identical(old, metadata)) stop("Daily cache manifest differs. Preserve the previous run and use a separate output.")
  } else fast_atomic_save_rds(metadata, target)
  training <- list(schema = "ERA5-warm02-bridge-training-v1", fingerprint = fingerprint,
    data = analysis$data, model_input_md5 = analysis$input_md5,
    note = "Original analysis data only, without fitted model/S4 objects. SET-to-Ta calibration and warm-range validation are not performed here.")
  training_path <- file.path(path, "energy_bridge_training_data.rds")
  if (!file.exists(training_path)) fast_atomic_save_rds(training, training_path)
  fast_daily_write_reader(path)
  invisible(metadata)
}

fast_daily_write_reader <- function(cache_dir) {
  # Export ordinary functions only. Source this file without brms/rstan/models.
  env <- environment(fast_daily_write_reader)
  symbols <- ls(envir = env, all.names = TRUE)
  chosen <- symbols[grepl("^fast_daily_", symbols) &
    vapply(symbols, function(x) is.function(get(x, envir = env)), logical(1))]
  chosen <- unique(c(chosen, "nc_axis", "nc_attribute", "normal_axis", "read_rectangular_history"))
  lines <- c("# Model-free reader for ERA5-warm02 daily posterior cache.",
    "# Required for integrity checking: digest. No rstan/brms model loading.",
    "# source('energy_cache_reader.R'); fast_daily_visit(cache_dir, function(x, index) { ... })",
    "# Each x contains every cached posterior draw; never average across dates before cooling truncation.")
  for (nm in sort(chosen)) lines <- c(lines, "", paste0(nm, " <- "),
    deparse(get(nm, envir = env), width.cutoff = 120L))
  fast_atomic_publish(file.path(cache_dir, "energy_cache_reader.R"), function(tmp)
    writeLines(enc2utf8(lines), tmp, useBytes = TRUE), replace = TRUE)
}

fast_daily_visit <- function(cache_dir, FUN, tile_ids = NULL, verify = TRUE) {
  if (!is.function(FUN)) stop("FUN must be a function(x, index).")
  manifest <- readRDS(file.path(cache_dir, "cache_manifest.rds"))
  if (!identical(manifest$schema, "ERA5-warm02-cache-manifest-v1")) stop("Unknown cache manifest.")
  if (is.null(tile_ids)) tile_ids <- manifest$tile_ids
  if (anyDuplicated(tile_ids) || any(!tile_ids %in% manifest$tile_ids)) stop("Invalid or repeated tile IDs.")
  for (id in tile_ids) {
    path <- file.path(cache_dir, sprintf("tile_%05d.rds", id))
    if (!file.exists(path)) stop("Tile is not committed: ", path)
    index <- readRDS(path)
    if (!identical(as.integer(index$tile_id), as.integer(id)) ||
        !identical(index$index_file, basename(path)) ||
        !identical(as.integer(index$cells$cell_id), manifest$tile_cell_ids[[as.character(id)]]))
      stop("Daily index identity or tile membership differs from its filename/manifest.")
    reference <- manifest$grid[match(index$cells$cell_id, manifest$grid$cell_id), , drop = FALSE]
    for (nm in c("lon_plot", "lat", "area_weight", "land_fraction"))
      if (!identical(as.numeric(index$cells[[nm]]), as.numeric(reference[[nm]])))
        stop("Daily cache grid metadata differs from its manifest: ", nm)
    fast_daily_validate(fast_daily_descriptor(index, cache_dir), cache_dir, manifest$fingerprint,
      ndraw = manifest$ndraw, verify_payload = FALSE)
    if (!identical(index$paired_draw_ids, manifest$paired_draw_ids) ||
        !identical(as.numeric(index$dates), as.numeric(manifest$dates)) ||
        max(abs(index$U_M0 - manifest$U_M0)) > 1e-10)
      stop("Daily tile dates, posterior pairing or baseline differ from the cache manifest.")
    for (j in seq_along(index$chunks)) FUN(fast_daily_read_chunk(index, j, cache_dir, verify), index)
  }
  invisible(manifest)
}

fast_daily_cdd_chunk <- function(x, slope, theta0, kappa = 1) {
  # A reusable daily calculation, not a fitted/validated energy model.
  # theta0 is the explicitly chosen cooling balance threshold, scalar or per draw.
  # Returned baseline_all_weather permits cooling-demand support-coverage audits.
  nd <- ncol(x$delta_U)
  expand <- function(v, name, positive = FALSE) {
    if (!is.numeric(v) || !length(v) || !length(v) %in% c(1L, nd) ||
        any(!is.finite(v)) || (positive && any(v <= 0))) stop("Invalid ", name)
    rep(v, length.out = nd)
  }
  slope <- expand(slope, "SET-on-Ta slope", TRUE)
  theta0 <- expand(theta0, "baseline threshold")
  kappa <- expand(kappa, "threshold transfer")
  if (any(kappa < 0)) stop("kappa must be nonnegative.")
  if (is.null(x$temperature_full)) stop("Same-day temperature was not embedded; join it by date/cell before calculating CDD.")
  all_t <- as.vector(x$temperature_full)
  if (length(all_t) != length(x$quality_packed)) stop("Temperature and mask alignment differs.")
  t <- all_t[x$supported_rows]
  good <- is.finite(t)
  q <- matrix(0, nrow(x$delta_U), nd)
  q0 <- q
  if (any(good)) {
    delta_theta <- sweep(x$delta_U[good, , drop = FALSE], 2L, kappa / slope, "*")
    baseline <- outer(t[good], theta0, "-")
    q0[good, ] <- pmax(baseline, 0)
    q[good, ] <- pmax(baseline - delta_theta, 0)
  }
  all_finite <- is.finite(all_t)
  total0 <- if (any(all_finite)) colSums(pmax(outer(all_t[all_finite], theta0, "-"), 0)) else rep(0, nd)
  # Rows excluded for temperature are returned explicitly, never counted as observed zero demand.
  q[!good, ] <- NA_real_; q0[!good, ] <- NA_real_
  list(CDD_M0 = q0, CDD_M1 = q, delta_CDD = q - q0,
    valid_supported_rows = good, cell_index = x$cell_index, row_dates = x$row_dates,
    baseline_all_weather = total0, n_valid_weather_rows = sum(all_finite),
    support_baseline_sum = colSums(q0, na.rm = TRUE),
    statement = "Unit: degree-C days for daily input. These are conditional climatic cooling indicators, not kWh/TWh.")
}

fast_daily_cdd_selftest <- function() {
  # Constant mean delta but different timing must give different cooling sums.
  x <- list(delta_U = matrix(c(2, -2), 2, 1), temperature_full = matrix(c(30, 20), 1, 2),
    quality_packed = as.raw(c(15, 15)), supported_rows = 1:2, cell_index = c(1L, 1L),
    row_dates = as.Date(c("2010-01-01", "2010-01-02")))
  a <- fast_daily_cdd_chunk(x, 1, 25)
  x$delta_U <- -x$delta_U
  b <- fast_daily_cdd_chunk(x, 1, 25)
  stopifnot(sum(a$CDD_M0) == 5, sum(a$CDD_M1) == 3, sum(b$CDD_M1) == 7)
  x$delta_U[,] <- 1
  c <- fast_daily_cdd_chunk(x, 1, 25)
  stopifnot(all(c$CDD_M1 <= c$CDD_M0))
  x$temperature_full[1] <- NA_real_
  z <- fast_daily_cdd_chunk(x, 1, 25)
  stopifnot(is.na(z$CDD_M1[1]), !z$valid_supported_rows[1])
  invisible(TRUE)
}

fast_daily_finalize_inventory <- function(bundle, tasks, fingerprint, out_dir) {
  path <- bundle$daily_cache$dir
  rows <- lapply(tasks, function(task) {
    index <- readRDS(file.path(path, sprintf("tile_%05d.rds", task$tile_id)))
    d <- fast_daily_descriptor(index, path)
    fast_daily_validate(d, path, fingerprint, task, bundle$ndraw, verify_payload = FALSE)
    data.frame(tile_id = task$tile_id, cached_rows = d$nrows, posterior_draws = d$ndraw,
      bytes = d$nbytes, max_float32_absolute_error = d$max_abs_error,
      cache_processing_write_seconds = d$elapsed_write_seconds,
      index_file = d$index_file, data_file = d$data_file)
  })
  tab <- do.call(rbind, rows)
  utils::write.csv(tab, file.path(path, "cache_inventory.csv"), row.names = FALSE)
  completion <- list(schema = "ERA5-warm02-cache-completion-v1", status = "COMPLETED",
    fingerprint = fingerprint, tiles = nrow(tab), total_bytes = sum(tab$bytes),
    supported_cell_days = sum(tab$cached_rows), max_float32_absolute_error = max(tab$max_float32_absolute_error),
    temperature_status = bundle$temperature$status,
    payload_verification = "Per-chunk digest recorded at write; byte counts and metadata checked at finalization. Payload digests checked when chunks are read.",
    finished = Sys.time())
  fast_atomic_save_rds(completion, file.path(path, "cache_completion.rds"), replace = TRUE)
  message(sprintf("Daily cache complete: %.1f GB; %s supported cell-days; same-day temperature: %s.",
    completion$total_bytes / 1e9, format(completion$supported_cell_days, scientific = FALSE), bundle$temperature$status))
  invisible(completion)
}


fast_daily_all_selftests <- function(output_dir = NULL) {
  for (pkg in c("digest", "ncdf4")) if (!requireNamespace(pkg, quietly = TRUE))
    stop("Daily-cache selftests need package: ", pkg)
  tmp <- tempfile("warm02_temperature_qc_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  a <- fast_daily_selftest()
  b <- fast_daily_temperature_selftest(tmp)
  fast_daily_cdd_selftest()
  report <- list(status = "PASSED", time = Sys.time(), cache = a, temperature = b,
    daily_CDD_ordering = "Two-day timing example, positive-shift monotonicity and missing-T handling passed.")
  if (!is.null(output_dir)) fast_atomic_save_rds(report,
    file.path(output_dir, "daily_cache_runtime_QC.rds"), replace = TRUE)
  invisible(report)
}


fast_warm_metrics <- function(ab, fixed_set, set_bounds) {
  # These are coefficients of the ALREADY standardized response curve.
  # Transform each day/draw before any temporal or spatial aggregation.
  M1 <- ab$B1 / ab$A1; N1 <- 1 / ab$A1; U1 <- M1 + N1
  M0 <- ab$B0 / ab$A0; N0 <- 1 / ab$A0; U0 <- M0 + N0
  delta_M <- sweep(M1, 2L, M0, "-")
  delta_N <- sweep(N1, 2L, N0, "-")
  delta_U <- delta_M + delta_N
  delta_TSV <- sweep(ab$A1, 2L, ab$A0, "-") * unname(fixed_set) -
    sweep(ab$B1, 2L, ab$B0, "-")
  U0_matrix <- matrix(rep(U0, each = nrow(U1)), nrow = nrow(U1), ncol = ncol(U1))
  # This diagnostic is not used to discard a draw or alter its denominator.
  in_range <- (U1 >= set_bounds[1L] & U1 <= set_bounds[2L] &
               U0_matrix >= set_bounds[1L] & U0_matrix <= set_bounds[2L]) * 1
  list(delta_M = delta_M, delta_N = delta_N, delta_U = delta_U,
       delta_TSV_P50 = delta_TSV, U_M1 = U1, U_M0 = U0_matrix,
       U_pair_in_SET_range = in_range, absolute_delta_U = abs(delta_U))
}

fast_read_checkpoint <- function(path, fingerprint, task, bundle) {
  fail <- function(reason) stop("Invalid committed warm-boundary checkpoint: ", path,
    "\n", reason, "\nPreserve this checkpoint; do not silently mix or overwrite it.", call. = FALSE)
  if (!file.exists(path)) fail("File does not exist.")
  cp <- tryCatch(readRDS(path), error = function(e) fail(conditionMessage(e)))
  fields <- c("schema", "fingerprint", "tile_id", "task_signature", "ncells", "ndraw", "periods",
    "planned", "metric_names", "aggregate_sums", "aggregate_denominator", "aggregate_cell_counts",
    "support_sums", "support_denominator", "support_cell_counts", "summary_table",
    "joint_draw_means", "joint_draw_denominator", "joint_draw_cell_ids", "elapsed_seconds")
  if (!is.list(cp) || !all(fields %in% names(cp))) fail("Required fields are missing.")
  if (!identical(cp$schema, "ERA5-warm02-tile-v1")) fail("Checkpoint schema is incompatible.")
  if (!identical(cp$fingerprint, fingerprint)) fail("Run fingerprint differs.")
  if (length(cp$tile_id) != 1L || !identical(as.integer(cp$tile_id), as.integer(task$tile_id)))
    fail("Tile ID differs.")
  if (!identical(cp$task_signature, fast_tile_signature(task))) fail("Tile membership differs.")
  ncells <- nrow(task$tile_points); ndraw <- as.integer(bundle$ndraw)
  np <- length(bundle$periods); nsupport <- length(bundle$g$support_fraction_sensitivity)
  metrics <- c("delta_M", "delta_N", "delta_U", "delta_TSV_P50", "U_M1", "U_M0",
               "U_pair_in_SET_range", "absolute_delta_U")
  if (!identical(cp$metric_names, metrics)) fail("Metric definition/order differs.")
  if (!identical(as.integer(cp$ncells), as.integer(ncells)) ||
      !identical(as.integer(cp$ndraw), ndraw)) fail("Cell/draw count differs.")
  if (!identical(as.character(cp$periods), as.character(bundle$periods)) ||
      !identical(as.numeric(cp$planned), as.numeric(bundle$planned))) fail("Period definitions differ.")
  expected <- list(aggregate_sums = c(np, 8L, ndraw),
    support_sums = c(np, 8L, nsupport, ndraw), support_denominator = c(np, nsupport),
    support_cell_counts = c(np, nsupport))
  for (nm in names(expected)) {
    x <- cp[[nm]]
    if (!is.numeric(x) || !identical(as.integer(dim(x)), as.integer(expected[[nm]])) ||
        any(!is.finite(x))) fail(paste("Invalid dimensions/type/values in", nm))
  }
  for (nm in c("aggregate_denominator", "aggregate_cell_counts")) {
    x <- cp[[nm]]
    if (!is.numeric(x) || !is.null(dim(x)) || length(x) != np || any(!is.finite(x)))
      fail(paste("Invalid vector", nm))
  }
  for (nm in c("aggregate_denominator", "aggregate_cell_counts", "support_denominator", "support_cell_counts")) {
    if (any(cp[[nm]] < 0)) fail(paste("Negative count/denominator in", nm))
    if (grepl("cell_counts$", nm) && any(cp[[nm]] > ncells | cp[[nm]] != floor(cp[[nm]])))
      fail(paste("Invalid eligible cell count in", nm))
  }
  if (length(cp$elapsed_seconds) != 1L || !is.finite(cp$elapsed_seconds) || cp$elapsed_seconds < 0)
    fail("Invalid elapsed time.")
  saved_periods <- if (isTRUE(bundle$g$save_annual_grid)) bundle$periods else tail(bundle$periods, 1L)
  tab <- cp$summary_table; nblocks <- length(saved_periods)
  needed <- c("cell_id", "lon", "lat", "land_fraction", "area_weight", "period", "fixed_set_name",
    "fixed_set", "n_expected_days", "n_valid_history_days", "n_common_support_days", "n_estimation_days",
    "date_coverage", "common_support_fraction", "eligible_for_area_day_summary", "eligible_for_supported_day_map")
  quantiles <- as.vector(outer(metrics, c("lower", "median", "upper", "prob_positive", "prob_negative"),
                               paste, sep = "_"))
  if (!is.data.frame(tab) || !all(c(needed, quantiles) %in% names(tab)) || nrow(tab) != ncells * nblocks)
    fail("Grid summary columns/row count differ.")
  if (!identical(as.integer(tab$cell_id), rep(as.integer(task$tile_points$cell_id), nblocks)) ||
      !identical(as.character(tab$period), rep(as.character(saved_periods), each = ncells)) ||
      !all(tab$fixed_set_name == "P50") || !all(tab$fixed_set == unname(bundle$fixed_set)))
    fail("Grid cell/period/P50 metadata differs.")
  for (nm in names(tab)[vapply(tab, is.numeric, logical(1))])
    if (any(is.infinite(tab[[nm]]) | is.nan(tab[[nm]]))) fail(paste("Invalid numeric summary:", nm))
  for (nm in c("cell_id", "lon", "lat", "land_fraction", "area_weight", "fixed_set", "n_expected_days",
               "n_valid_history_days", "n_common_support_days", "n_estimation_days", "date_coverage", "common_support_fraction"))
    if (!is.numeric(tab[[nm]]) || any(!is.finite(tab[[nm]]))) fail(paste("Nonfinite grid metadata:", nm))
  for (nm in c("lon", "lat", "land_fraction", "area_weight")) {
    values <- if (nm == "lon") task$tile_points$lon_plot else task$tile_points[[nm]]
    if (!identical(as.numeric(tab[[nm]]), rep(as.numeric(values), nblocks)))
      fail(paste("Grid coordinate/weight differs:", nm))
  }
  if (any(tab$n_common_support_days > tab$n_valid_history_days | tab$n_valid_history_days > tab$n_expected_days) ||
      any(tab$n_estimation_days != tab$n_common_support_days)) fail("Invalid support/date counts.")
  joint <- cp$joint_draw_means
  if (!is.numeric(joint) || !identical(as.integer(dim(joint)), c(as.integer(ncells), ndraw, 4L)) ||
      !identical(dimnames(joint)[[3L]], metrics[1:4]) || any(is.infinite(joint) | is.nan(joint)))
    fail("Invalid retained joint draws.")
  denominator <- cp$joint_draw_denominator
  if (!is.numeric(denominator) || length(denominator) != ncells ||
      any(!is.finite(denominator) | denominator < 0 | denominator != floor(denominator)) ||
      !identical(as.integer(cp$joint_draw_cell_ids), as.integer(task$tile_points$cell_id)))
    fail("Invalid joint-draw cell IDs/denominators.")
  full <- tab[tab$period == tail(bundle$periods, 1L), , drop = FALSE]
  if (!identical(as.numeric(denominator), as.numeric(full$n_common_support_days)))
    fail("Joint-draw denominator differs from the full-period summary.")
  valid <- denominator > 0
  if (any(valid) && any(!is.finite(joint[valid, , , drop = FALSE]))) fail("Missing draws for a supported cell.")
  if (any(!valid) && !all(is.na(joint[!valid, , , drop = FALSE]))) fail("No-data joint draws must be NA.")
  if (any(valid)) {
    error <- joint[valid, , 3L, drop = FALSE] - joint[valid, , 1L, drop = FALSE] - joint[valid, , 2L, drop = FALSE]
    if (max(abs(error)) > 1e-9) fail("Joint posterior identity delta_U = delta_M + delta_N failed.")
  }
  if (isTRUE(bundle$daily_cache$enabled)) {
    if (!is.list(cp$daily_cache)) fail("Missing required daily energy cache descriptor.")
    daily_index <- fast_daily_validate(cp$daily_cache, bundle$daily_cache$dir, fingerprint, task, ndraw,
      verify_payload = isTRUE(bundle$daily_cache$verify_payload_on_resume))
    if (!identical(as.integer(daily_index$paired_draw_ids$m1), as.integer(bundle$cache$draw_ids_m1)) ||
        !identical(as.integer(daily_index$paired_draw_ids$m0), as.integer(bundle$cache$draw_ids_m0)) ||
        !identical(as.character(daily_index$dates),
          unlist(lapply(bundle$files, function(x) as.character(x$dates)), use.names = FALSE)))
      fail("Daily cache draw identities or dates differ from the model/weather manifest.")
  }
  cp
}

fast_restore_partition <- function(cp, part_dir) {
  filename <- file.path(part_dir, sprintf("tile_%05d.rds", cp$tile_id))
  csv <- sub("[.]rds$", ".csv", filename)
  if (file.exists(filename)) {
    existing <- tryCatch(readRDS(filename), error = function(e) NULL)
    if (!identical(existing, cp$summary_table)) stop("Partition differs from committed checkpoint: ", filename,
      ". Preserve both files and inspect before resuming.", call. = FALSE)
  } else fast_atomic_save_rds(cp$summary_table, filename)
  if (!file.exists(csv)) fast_atomic_publish(csv, function(temporary)
    utils::write.csv(cp$summary_table, temporary, row.names = FALSE))
  invisible(filename)
}

fast_compute_tile <- function(task, bundle, part_dir, checkpoint_dir, fingerprint) {
  started <- Sys.time(); context <- list(stage = "initialize_tile", tile = task$tile_id)
  set_context <- function(stage, ...) context <<- c(list(stage = stage, tile = task$tile_id), list(...))
  ncs <- list(); daily_state <- NULL
  on.exit(if (!is.null(daily_state)) fast_daily_cleanup(daily_state), add = TRUE)
  on.exit(invisible(lapply(ncs, function(nc) try(ncdf4::nc_close(nc), silent = TRUE))), add = TRUE)
  tryCatch({
    dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    tile_id <- as.integer(task$tile_id)
    checkpoint_path <- file.path(checkpoint_dir, sprintf("tile_%05d.rds", tile_id))
    fast_daily_recover_warm(task, bundle, fingerprint, checkpoint_path, part_dir)
    if (file.exists(checkpoint_path)) {
      set_context("validate_existing_checkpoint")
      cp <- fast_read_checkpoint(checkpoint_path, fingerprint, task, bundle)
      filename <- fast_restore_partition(cp, part_dir)
      return(list(tile_id = tile_id, ncells = cp$ncells, elapsed_seconds = cp$elapsed_seconds,
        checkpoint = checkpoint_path, partition = filename, resumed = TRUE))
    }
    g <- bundle$g; fixed_set <- bundle$fixed_set; ndraw <- as.integer(bundle$ndraw)
    periods <- bundle$periods; planned <- bundle$planned; np <- length(periods)
    nsupport <- length(g$support_fraction_sensitivity); tile_points <- task$tile_points
    ncells <- nrow(tile_points); keep <- as.integer(task$keep); set_bounds <- bundle$set_bounds
    if (length(tile_id) != 1L || is.na(tile_id) || tile_id < 1L || ncells < 1L ||
        length(keep) != ncells || anyNA(keep) || any(keep < 1L) || anyDuplicated(keep) ||
        any(keep > length(task$lon_indices) * length(task$lat_indices))) stop("Invalid spatial tile.")
    if (np != length(bundle$files) + 1L || length(planned) != np ||
        any(!is.finite(planned) | planned <= 0) || ndraw < 1L ||
        length(fixed_set) != 1L || !identical(names(fixed_set), "P50") || !is.finite(fixed_set) ||
        length(set_bounds) != 2L || any(!is.finite(set_bounds)) || set_bounds[1L] >= set_bounds[2L] ||
        nsupport < 1L || any(!is.finite(g$support_fraction_sensitivity)) ||
        any(g$support_fraction_sensitivity < 0 | g$support_fraction_sensitivity > 1))
      stop("Invalid periods, draws, P50, SET bounds, or support thresholds.")
    if (any(!is.finite(tile_points$area_weight) | tile_points$area_weight <= 0)) stop("Invalid area weights.")
    metrics <- c("delta_M", "delta_N", "delta_U", "delta_TSV_P50", "U_M1", "U_M0",
                 "U_pair_in_SET_range", "absolute_delta_U")
    quality_names <- c("finite", "marginal", "joint", "primary", "B21_outside", "TV_lag1_outside",
                       "TV_lag2_outside", "TV_lag3_outside")
    aggregate_sums <- array(0, c(np, 8L, ndraw))
    aggregate_denominator <- numeric(np); aggregate_cell_counts <- integer(np)
    support_sums <- array(0, c(np, 8L, nsupport, ndraw))
    support_denominator <- matrix(0, np, nsupport); support_cell_counts <- matrix(0L, np, nsupport)
    sums <- array(0, c(ncells, ndraw, np, 8L))
    quality <- array(0L, c(ncells, np, 8L), dimnames = list(NULL, periods, quality_names))
    if (isTRUE(bundle$daily_cache$enabled))
      daily_state <- fast_daily_open(task, bundle, fingerprint, bundle$daily_cache$dir)
    set_context("open_annual_history_files")
    for (i in seq_along(bundle$files)) ncs[[i]] <- ncdf4::nc_open(bundle$files[[i]]$path)
    for (yy in seq_along(ncs)) {
      fi <- bundle$files[[yy]]; nt <- length(fi$dates)
      for (ds in seq.int(1L, nt, by = g$day_chunk)) {
        nday <- min(g$day_chunk, nt - ds + 1L)
        info <- list(year = fi$year, start_day = ds, n_days = nday, n_cells = ncells, n_draws = ndraw)
        set_context("read_history_chunk", chunk = info)
        mats <- lapply(era5_required_variables, function(v)
          read_rectangular_history(ncs[[yy]], v, task$lon_indices, task$lat_indices,
            ds, nday, fi$offsets[v])[keep, , drop = FALSE])
        history <- as.data.frame(setNames(lapply(mats, as.vector), era5_required_variables))
        cells <- rep(seq_len(ncells), times = nday)
        set_context("evaluate_support", chunk = info)
        support <- evaluate_support(bundle$support_analysis, history, joint = isTRUE(g$joint_support))
        needed <- c("finite", "marginal_inside", "joint_inside", quality_names[5:8])
        if (!is.data.frame(support) || nrow(support) != nrow(history) || !all(needed %in% names(support)))
          stop("Support evaluator returned invalid rows or columns.")
        finite <- as.logical(support$finite); finite[is.na(finite)] <- FALSE
        marginal <- finite & !is.na(support$marginal_inside) & support$marginal_inside
        joint <- finite & !is.na(support$joint_inside) & support$joint_inside
        primary <- marginal & (if (isTRUE(g$joint_support)) joint else TRUE)
        flags <- cbind(finite, marginal, joint, primary, as.matrix(support[, quality_names[5:8], drop = FALSE]))
        flags[is.na(flags)] <- FALSE; flags[!finite, 5:8] <- FALSE
        for (q in seq_along(quality_names)) {
          counts <- tabulate(cells[as.logical(flags[, q])], nbins = ncells)
          quality[, yy, q] <- quality[, yy, q] + counts
          quality[, np, q] <- quality[, np, q] + counts
        }
        chunk_dates <- as.Date(fi$dates[seq.int(ds, length.out = nday)])
        temperature <- NULL
        if (!is.null(daily_state) && identical(bundle$temperature$status, "embedded")) {
          set_context("read_same_day_temperature", chunk = info)
          temperature <- fast_daily_temperature_chunk(bundle$temperature, task, chunk_dates)
        }
        # Record ALL dates/masks, including entirely unsupported chunks.
        if (!any(primary)) {
          if (!is.null(daily_state)) fast_daily_append(daily_state, chunk_dates,
            flags[, 1:4, drop = FALSE], NULL, temperature = temperature)
          next
        }
        set_context("predict_supported_ab", chunk = info, supported_rows = sum(primary))
        ab <- fast_predict_ab(bundle$cache, history[primary, , drop = FALSE])
        expected <- c(sum(primary), ndraw)
        if (!identical(as.integer(dim(ab$A1)), as.integer(expected)) ||
            !identical(as.integer(dim(ab$B1)), as.integer(expected)) ||
            length(ab$A0) != ndraw || length(ab$B0) != ndraw ||
            any(!is.finite(ab$A1)) || any(!is.finite(ab$B1)) ||
            any(!is.finite(ab$A0)) || any(!is.finite(ab$B0)) ||
            any(ab$A1 <= 0) || any(ab$A0 <= 0)) stop("Invalid positive-slope paired prediction output.")
        set_context("accumulate_warm_posterior", chunk = info)
        values <- fast_warm_metrics(ab, fixed_set, set_bounds)
        if (!is.null(daily_state)) {
          set_context("write_daily_energy_cache", chunk = info)
          fast_daily_append(daily_state, chunk_dates, flags[, 1:4, drop = FALSE],
            values$delta_U, temperature = temperature)
        }
        selected <- rep(TRUE, sum(primary)); selected_cells <- cells[primary]
        for (metric in seq_along(metrics)) {
          if (any(!is.finite(values[[metric]]))) stop("Nonfinite transformed metric; coverage must not vary by draw.")
          increment <- sum_by_cell(values[[metric]], selected, selected_cells, ncells)
          sums[, , yy, metric] <- matrix(sums[, , yy, metric], ncells, ndraw) + increment
          sums[, , np, metric] <- matrix(sums[, , np, metric], ncells, ndraw) + increment
        }
      }
    }
    if (!is.null(bundle$capacity_supported_rows_by_tile)) {
      expected_supported <- bundle$capacity_supported_rows_by_tile[tile_id]
      actual_supported <- sum(quality[, np, "primary"])
      if (length(expected_supported) != 1L || !is.finite(expected_supported) ||
          actual_supported != expected_supported)
        stop("Exact capacity scan and production support counts disagree in tile ", tile_id,
          ": expected ", expected_supported, "; observed ", actual_supported,
          ". Halt before committing this tile; inspect input stability.")
    }
    set_context("summarize_warm_tile")
    rows <- list(); area <- tile_points$area_weight
    for (p in seq_along(periods)) {
      finite_n <- quality[, p, "finite"]; denominator <- quality[, p, "primary"]
      date_coverage <- finite_n / planned[p]; fraction <- denominator / planned[p]
      eligible <- date_coverage >= g$min_date_coverage & denominator >= g$min_supported_days &
                  fraction >= g$min_supported_fraction
      aggregate_denominator[p] <- sum(area[eligible] * denominator[eligible])
      aggregate_cell_counts[p] <- sum(eligible)
      save_grid <- p == np || isTRUE(g$save_annual_grid)
      qtables <- if (save_grid) vector("list", length(metrics)) else NULL
      for (metric in seq_along(metrics)) {
        mat <- matrix(sums[, , p, metric], ncells, ndraw)
        if (save_grid) qtables[[metric]] <- posterior_rows(mat, denominator)
        aggregate_sums[p, metric, ] <- colSums(mat[eligible, , drop = FALSE] * area[eligible])
        for (cut in seq_len(nsupport)) {
          chosen <- date_coverage >= g$min_date_coverage & denominator >= g$min_supported_days &
                    fraction >= g$support_fraction_sensitivity[cut]
          support_sums[p, metric, cut, ] <- colSums(mat[chosen, , drop = FALSE] * area[chosen])
          if (metric == 1L) {
            support_denominator[p, cut] <- sum(area[chosen] * denominator[chosen])
            support_cell_counts[p, cut] <- sum(chosen)
          }
        }
      }
      if (!save_grid) next
      row <- data.frame(cell_id = tile_points$cell_id, lon = tile_points$lon_plot, lat = tile_points$lat,
        land_fraction = tile_points$land_fraction, area_weight = area, period = periods[p],
        fixed_set_name = "P50", fixed_set = unname(fixed_set),
        fixed_SET_applies_only_to = "delta_TSV_P50; warm boundaries do not depend on P50",
        n_expected_days = planned[p], n_valid_history_days = finite_n,
        n_marginal_support_days = quality[, p, "marginal"], n_joint_support_days = quality[, p, "joint"],
        n_common_support_days = denominator, n_estimation_days = denominator,
        date_coverage = date_coverage, common_support_fraction = fraction,
        common_support_fraction_of_valid = ifelse(finite_n > 0, denominator / finite_n, NA_real_),
        eligible_for_area_day_summary = eligible, eligible_for_supported_day_map = eligible,
        required_common_support_fraction = g$min_supported_fraction,
        SET_range_lower = set_bounds[1L], SET_range_upper = set_bounds[2L],
        pilot_not_global = !is.null(g$max_cells))
      for (q in quality_names[5:8]) row[[paste0(q, "_fraction_of_valid")]] <-
        ifelse(finite_n > 0, quality[, p, q] / finite_n, NA_real_)
      for (metric in seq_along(metrics)) for (field in names(qtables[[metric]]))
        row[[paste0(metrics[metric], "_", field)]] <- qtables[[metric]][[field]]
      for (metric in metrics[1:4]) row[[paste0(metric, "_CrI_excludes_zero")]] <-
        row[[paste0(metric, "_lower")]] > 0 | row[[paste0(metric, "_upper")]] < 0
      rows[[length(rows) + 1L]] <- row
    }
    table <- do.call(rbind, rows); rownames(table) <- NULL
    joint_denominator <- as.numeric(quality[, np, "primary"])
    joint_draw_means <- array(NA_real_, c(ncells, ndraw, 4L), dimnames = list(NULL, NULL, metrics[1:4]))
    valid <- joint_denominator > 0
    for (metric in 1:4) if (any(valid)) {
      mat <- matrix(sums[, , np, metric], ncells, ndraw)
      joint_draw_means[valid, , metric] <- mat[valid, , drop = FALSE] / joint_denominator[valid]
    }
    set_context("commit_warm_tile")
    filename <- file.path(part_dir, sprintf("tile_%05d.rds", tile_id)); csv <- sub("[.]rds$", ".csv", filename)
    cp <- list(schema = "ERA5-warm02-tile-v1", fingerprint = fingerprint, tile_id = tile_id,
      task_signature = fast_tile_signature(task), ncells = ncells, ndraw = ndraw,
      periods = periods, planned = planned, metric_names = metrics,
      aggregate_sums = aggregate_sums, aggregate_denominator = aggregate_denominator,
      aggregate_cell_counts = aggregate_cell_counts, support_sums = support_sums,
      support_denominator = support_denominator, support_cell_counts = support_cell_counts,
      summary_table = table, joint_draw_means = joint_draw_means,
      joint_draw_denominator = joint_denominator, joint_draw_cell_ids = as.integer(tile_points$cell_id),
      joint_draw_note = "Full-period supported-date mean per cell and fixed paired draw; delta_U = delta_M + delta_N.",
      history_row_counts = c(all = ncells * planned[np], finite = sum(quality[, np, "finite"]),
        predicted = sum(quality[, np, "primary"]),
        finite_skipped_before_prediction = sum(quality[, np, "finite"] - quality[, np, "primary"])),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), created = Sys.time())
    if (!is.null(daily_state)) {
      set_context("commit_daily_energy_cache")
      daily_state$recovery_checkpoint <- cp
      cp$daily_cache <- fast_daily_finish(daily_state)
      cp$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    }
    fast_atomic_save_rds(table, filename, replace = TRUE)
    fast_atomic_publish(csv, function(temporary) utils::write.csv(table, temporary, row.names = FALSE), replace = TRUE)
    fast_atomic_publish(checkpoint_path, function(temporary) {
      saveRDS(cp, temporary, compress = "gzip")
      invisible(fast_read_checkpoint(temporary, fingerprint, task, bundle))
    })
    list(tile_id = tile_id, ncells = ncells, elapsed_seconds = cp$elapsed_seconds,
      checkpoint = checkpoint_path, partition = filename, resumed = FALSE)
  }, error = function(e) {
    report <- list(schema = "ERA5-warm02-tile-error-v1", fingerprint = fingerprint,
      tile_id = task$tile_id, time = Sys.time(), context = context, error = conditionMessage(e),
      call = paste(deparse(conditionCall(e), width.cutoff = 120L), collapse = " "))
    path <- file.path(checkpoint_dir, sprintf("tile_%05d_error.rds", as.integer(task$tile_id)))
    try(fast_atomic_save_rds(report, path, replace = TRUE), silent = TRUE)
    stop("Tile ", task$tile_id, " failed at ", context$stage, ": ", conditionMessage(e), call. = FALSE)
  })
}

# Deterministic regression checks run in a temporary directory, using the real
# NetCDF reader, tile accumulator, serializer, resume checks and summaries.
# Synthetic inputs never enter the production run or any scientific figure.
fast_tile_selftest <- function(report_path = NULL) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("ncdf4 is required for warm tile regression QC.")
  test_dir <- tempfile("ERA5_warm01_regression_")
  if (!dir.create(test_dir)) stop("Cannot create tile-test directory.")
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  g <- list(day_chunk = 1L, joint_support = FALSE, max_cells = NULL,
    min_date_coverage = .5, min_supported_days = 1L, min_supported_fraction = 0,
    support_fraction_sensitivity = c(0, .5, .8, .95), save_annual_grid = TRUE)
  files <- vector("list", 2L)
  for (year in 1:2) {
    path <- file.path(test_dir, paste0("test_", year, ".nc"))
    dims <- list(ncdf4::ncdim_def("lon", "degrees_east", c(0, 1)),
      ncdf4::ncdim_def("lat", "degrees_north", 10),
      ncdf4::ncdim_def("time", "days since 2000-01-01", 0:1))
    definitions <- lapply(era5_required_variables, function(v)
      ncdf4::ncvar_def(v, "degC", dims, missval = -9999, prec = "double"))
    nc <- ncdf4::nc_create(path, definitions)
    tryCatch({
      for (v in era5_required_variables) {
        x <- if (v == "B21_raw") c(year, NA_real_, year + 2, NA_real_) else c(0, NA_real_, 0, NA_real_)
        ncdf4::ncvar_put(nc, v, array(x, c(2L, 1L, 2L)))
      }
    }, finally = ncdf4::nc_close(nc))
    files[[year]] <- list(path = path, year = 2000L + year,
      dates = as.Date(sprintf("%04d-01-01", 2000L + year)) + 0:1, offsets = setNames(rep(0, 4), era5_required_variables))
  }
  nd <- 5L
  task <- list(tile_id = 1L, lon_indices = 1:2, lat_indices = 1L, keep = 1:2,
    tile_points = data.frame(cell_id = 1:2, lon_plot = c(0, 1), lat = c(10, 10),
      land_fraction = c(1, 1), area_weight = c(.8, .9)))
  bundle <- list(g = g, files = files, fixed_set = c(P50 = 24.3), set_bounds = c(15.5, 19),
    cache = list(n_draws = nd), support_analysis = list(), ndraw = nd,
    periods = c("2001", "2002", "2001-2002"), planned = c(2L, 2L, 4L))
  env <- new.env(parent = environment(fast_compute_tile)); env$predicted_rows <- 0L
  env$evaluate_support <- function(analysis, history, joint) {
    finite <- rowSums(!is.finite(as.matrix(history))) == 0L
    supported <- finite & history$B21_raw <= 2
    data.frame(finite = finite, marginal_inside = supported, joint_inside = supported,
      B21_outside = !supported, TV_lag1_outside = FALSE, TV_lag2_outside = FALSE, TV_lag3_outside = FALSE)
  }
  env$fast_predict_ab <- function(cache, history) {
    if (!nrow(history) || any(!is.finite(history$B21_raw) | history$B21_raw > 2))
      stop("QC: predictor received an unsupported or empty row.")
    env$predicted_rows <- env$predicted_rows + nrow(history)
    draws <- seq_len(cache$n_draws)
    list(A1 = outer(history$B21_raw, draws, function(h, d) .1 + .01 * h * d),
      B1 = outer(history$B21_raw, draws, function(h, d) .4 + .3 * h * d),
      A0 = .08 + draws * .001, B0 = .3 + draws * .002)
  }
  scan <- fast_daily_support_preflight
  environment(scan) <- env
  capacity_test_dir <- file.path(test_dir, "capacity_preflight")
  scanned <- scan(list(available_space_gb = 1, reserve_space_gb = .01), list(task),
    files, bundle$support_analysis, g, list(test_input = TRUE), nd, capacity_test_dir)
  scanned_again <- scan(list(available_space_gb = 1, reserve_space_gb = .01), list(task),
    files, bundle$support_analysis, g, list(test_input = TRUE), nd, capacity_test_dir)
  if (!identical(scanned$supported_rows_by_tile, c(2)) ||
      !identical(scanned$supported_rows_by_tile, scanned_again$supported_rows_by_tile))
    stop("Support-based disk preflight and restart checkpoint QC failed.")
  compute <- fast_compute_tile; environment(compute) <- env
  run_case <- function(name, changed_task = task, chunk = 1L, annual = TRUE, n_draws = nd, daily = FALSE) {
    changed <- bundle; changed$g$day_chunk <- chunk; changed$g$save_annual_grid <- annual
    changed$ndraw <- as.integer(n_draws); changed$cache$n_draws <- as.integer(n_draws)
    directory <- file.path(test_dir, name); before <- env$predicted_rows
    if (daily) {
      changed$daily_cache <- list(enabled = TRUE, dir = file.path(directory, "daily"),
        max_abs_error = 5e-5, verify_payload_on_resume = FALSE)
      d <- seq_len(changed$ndraw)
      changed$cache$A0 <- .08 + d * .001; changed$cache$B0 <- .3 + d * .002
      changed$cache$draw_ids_m1 <- d; changed$cache$draw_ids_m0 <- d + 1000L
      changed$temperature <- list(status = "embedded")
      env$fast_daily_temperature_chunk <- function(temp_spec, task, dates)
        matrix(30, nrow(task$tile_points), length(dates))
    }
    result <- compute(changed_task, changed, file.path(directory, "parts"),
      file.path(directory, "checkpoints"), "test-only-warm-fingerprint")
    cp <- fast_read_checkpoint(result$checkpoint, "test-only-warm-fingerprint", changed_task, changed)
    list(cp = cp, result = result, bundle = changed, task = changed_task,
      rows_predicted = env$predicted_rows - before)
  }
  full <- run_case("full"); chunked <- run_case("chunk2", chunk = 2L)
  guarded <- bundle; guarded$capacity_supported_rows_by_tile <- scanned$supported_rows_by_tile
  guarded_dir <- file.path(test_dir, "capacity_guarded")
  guarded_result <- compute(task, guarded, file.path(guarded_dir, "parts"),
    file.path(guarded_dir, "checkpoints"), "test-only-warm-fingerprint")
  if (readRDS(guarded_result$checkpoint)$history_row_counts["predicted"] != 2)
    stop("Production tile did not match exact capacity support count.")
  guarded$capacity_supported_rows_by_tile <- 3
  mismatch <- tryCatch({
    compute(task, guarded, file.path(test_dir, "capacity_mismatch", "parts"),
      file.path(test_dir, "capacity_mismatch", "checkpoints"), "test-only-warm-fingerprint")
    NULL
  }, error = function(e) e)
  if (!inherits(mismatch, "error") || !grepl("Exact capacity scan", conditionMessage(mismatch), fixed = TRUE))
    stop("Production tile did not reject a changed support count.")
  if (full$rows_predicted != 2L || chunked$rows_predicted != 2L)
    stop("Support-skip QC failed: exactly two supported histories must be predicted.")
  fields <- c("aggregate_sums", "aggregate_denominator", "aggregate_cell_counts", "support_sums",
    "support_denominator", "support_cell_counts", "summary_table", "joint_draw_means", "joint_draw_denominator")
  for (field in fields) if (!isTRUE(all.equal(full$cp[[field]], chunked$cp[[field]], tolerance = 1e-12)))
    stop("Chunk-invariance QC failed: ", field)
  metrics <- c("delta_M", "delta_N", "delta_U", "delta_TSV_P50", "U_M1", "U_M0",
               "U_pair_in_SET_range", "absolute_delta_U")
  # Independent scalar reference, including U computed directly as (B+1)/A.
  reference <- array(0, c(3L, 8L, nd))
  histories <- list(1, 2, c(1, 2))
  reference_days <- c(1L, 1L, 2L)
  for (p in 1:3) for (d in seq_len(nd)) for (h in histories[[p]]) {
    a1 <- .1 + .01 * h * d; b1 <- .4 + .3 * h * d
    a0 <- .08 + d * .001; b0 <- .3 + d * .002
    m1 <- b1 / a1; m0 <- b0 / a0; n1 <- 1 / a1; n0 <- 1 / a0
    u1 <- (b1 + 1) / a1; u0 <- (b0 + 1) / a0
    expected_values <- c(m1 - m0, n1 - n0, u1 - u0,
      (a1 * 24.3 - b1) - (a0 * 24.3 - b0), u1, u0,
      as.numeric(u1 >= 15.5 && u1 <= 19 && u0 >= 15.5 && u0 <= 19), abs(u1 - u0))
    reference[p, , d] <- reference[p, , d] + expected_values
  }
  expected_support <- array(0, c(3L, 8L, 4L, nd))
  expected_support[, , 1L, ] <- reference * .8
  expected_support[, , 2L, ] <- reference * .8
  expected_support_den <- cbind(.8 * reference_days, .8 * reference_days, rep(0, 3L), rep(0, 3L))
  checks <- list(aggregate_sums = reference * .8, aggregate_denominator = .8 * reference_days,
    aggregate_cell_counts = rep(1L, 3L), support_sums = expected_support,
    support_denominator = expected_support_den,
    support_cell_counts = cbind(rep(1L, 3L), rep(1L, 3L), rep(0L, 3L), rep(0L, 3L)))
  errors <- vapply(names(checks), function(field) {
    if (!identical(dim(checks[[field]]), dim(full$cp[[field]]))) stop("Reference dimension QC failed: ", field)
    max(abs(full$cp[[field]] - checks[[field]]))
  }, numeric(1))
  if (any(!is.finite(errors) | errors > 1e-10)) stop("Independent scalar reference QC failed.")
  tab <- full$cp$summary_table
  for (r in which(tab$cell_id == 1L)) {
    p <- match(tab$period[r], bundle$periods)
    for (m in seq_along(metrics)) {
      estimates <- reference[p, m, ] / reference_days[p]
      expected_quantiles <- c(stats::quantile(estimates, c(.025, .5, .975), type = 8, names = FALSE),
                             mean(estimates > 0), mean(estimates < 0))
      columns <- paste0(metrics[m], "_", c("lower", "median", "upper", "prob_positive", "prob_negative"))
      actual <- as.numeric(unlist(tab[r, columns], use.names = FALSE))
      if (max(abs(actual - expected_quantiles)) > 1e-10) stop("Posterior-quantile QC failed.")
    }
  }
  joint <- full$cp$joint_draw_means
  if (max(abs(joint[1L, , 3L] - joint[1L, , 1L] - joint[1L, , 2L])) > 1e-10 ||
      max(abs(joint[1L, , ] - t(reference[3L, 1:4, ]) / 2)) > 1e-10)
    stop("Retained joint posterior or positive N sign QC failed.")
  # Verify daily nonlinear transformation is not replaced by transforming mean A/B.
  wrong_u <- (mean(.4 + .3 * c(1, 2)) + 1) / mean(.1 + .01 * c(1, 2))
  correct_u <- reference[3L, 5L, 1L] / 2
  if (abs(wrong_u - correct_u) < 1e-4) stop("Test fixture does not distinguish transform-before-average.")
  metric_columns <- as.vector(outer(metrics, c("lower", "median", "upper", "prob_positive", "prob_negative"), paste, sep = "_"))
  if (!all(is.na(as.matrix(tab[tab$cell_id == 2L, metric_columns]))) || !all(is.na(joint[2L, , ])))
    stop("No-data summaries/draws must be NA.")
  one <- task; one$keep <- 1L; one$tile_points <- task$tile_points[1L, , drop = FALSE]
  single <- run_case("single", one)
  for (field in names(checks)) if (!isTRUE(all.equal(single$cp[[field]], checks[[field]], tolerance = 1e-12)))
    stop("Single-cell regression QC failed: ", field)
  missing <- task; missing$keep <- 2L; missing$tile_points <- task$tile_points[2L, , drop = FALSE]
  nodata <- run_case("nodata", missing)
  for (field in names(checks)) if (any(nodata$cp[[field]] != 0)) stop("No-data totals must be zero: ", field)
  if (nodata$rows_predicted != 0L) stop("No-data tile must not call predictor.")
  singleton <- run_case("single_draw", one, n_draws = 1L)
  if (!identical(dim(singleton$cp$joint_draw_means), c(1L, 1L, 4L))) stop("Singleton dimensions were lost.")
  compact <- run_case("no_annual_grid", annual = FALSE)
  if (nrow(compact$cp$summary_table) != 2L ||
      !isTRUE(all.equal(compact$cp$aggregate_sums, full$cp$aggregate_sums, tolerance = 1e-12)))
    stop("Disabling annual grid saving changed numerical totals.")
  before_hash <- unname(tools::md5sum(full$result$checkpoint)); before_rows <- env$predicted_rows
  # Also exercise restoring a missing non-authoritative partition from its checkpoint.
  unlink(c(full$result$partition, sub("[.]rds$", ".csv", full$result$partition)))
  resumed <- compute(task, bundle, dirname(full$result$partition), dirname(full$result$checkpoint),
                     "test-only-warm-fingerprint")
  if (!isTRUE(resumed$resumed) || !identical(before_hash, unname(tools::md5sum(full$result$checkpoint))) ||
      env$predicted_rows != before_rows || !identical(readRDS(resumed$partition), full$cp$summary_table))
    stop("Immutable resume/partition restore QC failed.")
  rejection <- tryCatch({
    fast_read_checkpoint(full$result$checkpoint, "wrong-fingerprint", task, bundle); NULL
  }, error = function(e) e)
  if (!inherits(rejection, "error") || !grepl("Run fingerprint differs", conditionMessage(rejection), fixed = TRUE))
    stop("Checkpoint fingerprint mismatch was not rejected.")
  cached <- run_case("daily_cache", chunk = 2L, daily = TRUE)
  for (field in fields) if (!isTRUE(all.equal(cached$cp[[field]], full$cp[[field]], tolerance = 1e-12)))
    stop("Enabling daily caching changed warm scientific output: ", field)
  idx <- fast_daily_validate(cached$cp$daily_cache, cached$bundle$daily_cache$dir,
    "test-only-warm-fingerprint", task, nd, TRUE)
  z <- lapply(seq_along(idx$chunks), function(j)
    fast_daily_read_chunk(idx, j, cached$bundle$daily_cache$dir))
  daily_draws <- do.call(rbind, lapply(z, function(x) x$delta_U))
  if (nrow(daily_draws) != 2L || ncol(daily_draws) != nd ||
      max(abs(colMeans(daily_draws) - full$cp$joint_draw_means[1L, , 3L])) > 5e-5 ||
      length(idx$dates) != 4L || !all(vapply(z, function(x) all(x$temperature_full == 30), logical(1))))
    stop("Integrated daily cache does not match the independent tile reference.")
  before_rows <- env$predicted_rows
  unlink(c(cached$result$checkpoint, cached$result$partition,
    sub("[.]rds$", ".csv", cached$result$partition)))
  restored <- compute(task, cached$bundle, dirname(cached$result$partition),
    dirname(cached$result$checkpoint), "test-only-warm-fingerprint")
  if (!isTRUE(restored$resumed) || env$predicted_rows != before_rows ||
      !identical(readRDS(restored$partition), cached$cp$summary_table))
    stop("Committed daily cache recovery reran the model or changed a warm summary.")
  tests <- c("exact_support_capacity_scan_and_restart", "capacity_count_matches_tile",
    "capacity_count_mismatch_rejected", "integrated_daily_cache", "cache_commit_gap_recovery_without_prediction",
    "independent_scalar_reference", "daily_transform_before_average", "positive_N_joint_identity",
    "supported_rows_only", "day_chunk_invariance", "single_cell", "no_data", "single_draw",
    "annual_grid_optional", "immutable_resume_and_partition_restore", "fingerprint_rejection")
  report <- data.frame(test = tests, status = "PASS", max_abs_error = c(rep(0, 3L),
    idx$max_abs_error, 0, max(errors), rep(0, length(tests) - 6L)))
  if (!is.null(report_path)) utils::write.csv(report, report_path, row.names = FALSE)
  message("PASS: sixteen warm-boundary/cache/checkpoint regression checks.")
  invisible(report)
}

fast_runtime_selftest <- function(output_dir) {
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE))
    stop("Cannot create self-test directory: ", output_dir)
  path <- file.path(output_dir, "warm_tile_regression_QC.csv")
  tryCatch(fast_tile_selftest(path), error = function(e) {
    utils::write.csv(data.frame(test = "warm_tile_regression", status = "FAILED", max_abs_error = NA_real_,
      detail = conditionMessage(e)), path, row.names = FALSE)
    stop("Warm-boundary tile/checkpoint regression failed: ", conditionMessage(e), call. = FALSE)
  })
}

fast_defaults <- function(config) {
  f <- config$fast
  if (is.null(f)) f <- list()
  defaults <- list(backend = "auto", workers = 8L, worker_candidates = c(8L, 12L),
    autotune = TRUE, resume = "auto", runtime_selftest = TRUE, target_total_hours = 8)
  for (nm in names(defaults)) if (is.null(f[[nm]])) f[[nm]] <- defaults[[nm]]
  if (!is.character(f$backend) || length(f$backend) != 1L || is.na(f$backend) ||
      !f$backend %in% c("auto", "cpp", "r")) stop("fast$backend must be auto, cpp, or r.")
  for (nm in c("autotune", "runtime_selftest"))
    if (!is.logical(f[[nm]]) || length(f[[nm]]) != 1L || is.na(f[[nm]]))
      stop("fast$", nm, " must be TRUE or FALSE.")
  for (nm in c("workers", "worker_candidates")) {
    x <- f[[nm]]
    if (!is.numeric(x) || !length(x) || any(!is.finite(x) | x < 1 | x != as.integer(x)))
      stop("Invalid fast setting: ", nm)
    f[[nm]] <- as.integer(x)
  }
  if (length(f$workers) != 1L) stop("fast$workers must be one integer.")
  if (length(f$target_total_hours)!=1L || !is.finite(f$target_total_hours) || f$target_total_hours<=0) stop("fast$target_total_hours must be positive.")
  f$worker_candidates <- sort(unique(f$worker_candidates))
  f
}

fast_worker_task <- function(task) {
  fast_compute_tile(task, .fast_worker_bundle, .fast_worker_part_dir,
                    .fast_worker_checkpoint_dir, .fast_worker_fingerprint)
}

fast_worker_initialize <- function(paths, cache_dir, backend, bundle, part_dir,
                                   checkpoint_dir, fingerprint, qc_history, qc_expected) {
  .libPaths(paths)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
  for (p in c("ncdf4", "dlnm", "splines", if (isTRUE(bundle$g$joint_support)) "FNN"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Worker is missing package: ", p)
  .analysis_progress <- new.env(parent = emptyenv())
  .analysis_progress$stage <- "worker_initialization"
  .analysis_progress$context <- list()
  assign(".analysis_progress", .analysis_progress, envir = .GlobalEnv)
  selected <- fast_init_backend(cache_dir, backend)
  if (!identical(selected, backend)) stop("Worker backend differs from master.")
  actual <- fast_predict_ab(bundle$cache, qc_history)
  error <- max(vapply(c("A1", "B1", "A0", "B0"), function(nm) {
    if (!identical(dim(actual[[nm]]), dim(qc_expected[[nm]])) ||
        length(actual[[nm]]) != length(qc_expected[[nm]])) return(Inf)
    max(abs(actual[[nm]] - qc_expected[[nm]]))
  }, numeric(1)))
  if (!is.finite(error) || error > 1e-8) stop("Worker prediction QC failed: ", error)
  assign(".fast_worker_bundle", bundle, envir = .GlobalEnv)
  assign(".fast_worker_part_dir", part_dir, envir = .GlobalEnv)
  assign(".fast_worker_checkpoint_dir", checkpoint_dir, envir = .GlobalEnv)
  assign(".fast_worker_fingerprint", fingerprint, envir = .GlobalEnv)
  list(pid = Sys.getpid(), backend = selected, max_absolute_error = error)
}

fast_make_tasks <- function(grid, ilon, ilat, tile_side) {
  xb <- split(seq_along(ilon), ceiling(seq_along(ilon) / tile_side))
  yb <- split(seq_along(ilat), ceiling(seq_along(ilat) / tile_side))
  tasks <- list()
  for (y in yb) for (x in xb) {
    local <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE)
    ids <- (local$y - 1L) * length(ilon) + local$x
    keep <- which(grid$land[ids])
    if (!length(keep)) next
    id <- length(tasks) + 1L
    tasks[[id]] <- list(tile_id = id, lon_indices = ilon[x], lat_indices = ilat[y],
      keep = keep, tile_points = grid[ids[keep], , drop = FALSE])
  }
  tasks
}

run_global_era5 <- function(config, analysis = NULL) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Install ncdf4 before global analysis.")
  g <- global_defaults(config); f <- fast_defaults(config)
  if (is.null(analysis)) analysis <- load_analysis(config)
  if (isTRUE(g$joint_support) && is.null(analysis$joint_reference))
    stop("Global joint support requires the original cached training reference.")
  fixed_set <- analysis$fixed_set["P50"]
  if (length(fixed_set) != 1L || !identical(names(fixed_set), "P50")) stop("This analysis requires P50 only.")
  set_bounds <- as.numeric(stats::quantile(analysis$data$SET_raw, c(.01,.99), na.rm=TRUE, type=7))
  if (length(set_bounds)!=2L || any(!is.finite(set_bounds)) || diff(set_bounds)<=0) stop("Invalid training SET bounds.")
  if (!length(fixed_set) || any(!is.finite(fixed_set))) stop("Invalid fixed SET values.")
  if (is.null(names(fixed_set))) names(fixed_set) <- paste0("SET", seq_along(fixed_set))
  available <- nrow(analysis$draws_m1)
  if (nrow(analysis$draws_m0) != available) stop("M1/M0 posterior row counts differ.")
  requested <- if (is.infinite(g$n_draws)) available else as.integer(g$n_draws)
  ndraw <- min(requested, available)
  if (requested > available) warning("Using all available paired draws: ", ndraw)
  if (ndraw < 50L) stop("At least 50 paired posterior draws are required.")
  positions <- seq_len(ndraw)
  analysis$draws_m1 <- analysis$draws_m1[positions, , drop = FALSE]
  analysis$draws_m0 <- analysis$draws_m0[positions, , drop = FALSE]
  analysis$draw_ids_m1 <- analysis$draw_ids_m1[positions]
  analysis$draw_ids_m0 <- analysis$draw_ids_m0[positions]
  analysis$draw_ids <- analysis$draw_ids_m1
  analysis$draw_metadata$global_draws_requested <- g$n_draws
  analysis$draw_metadata$global_draws_used <- ndraw
  analysis$draw_metadata$paired_draws_before_global_subset <- available
  analysis$draw_metadata$n_draws <- ndraw
  out_dir <- file.path(config$output_dir, "global")
  part_dir <- file.path(out_dir, "grid_partitions")
  checkpoint_dir <- file.path(out_dir, "checkpoints")
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  .set_analysis_stage("inspect_global_inputs", announce = TRUE)
  inspected <- inspect_history_files(g)
  daily_settings <- fast_daily_defaults(config)
  config$daily_cache <- daily_settings
  .set_analysis_stage("inspect_same_day_temperature_inputs", announce = TRUE)
  temperature <- fast_daily_inspect_temperature(config, inspected)
  daily_settings$dir <- file.path(out_dir, "energy_cache")
  dir.create(daily_settings$dir, recursive = TRUE, showWarnings = FALSE)
  first <- inspected$files[[1L]]
  # Enforce the native 0.25-degree spacing and stride=1.
  # ERA5 sample centres are retained; no interpolation or cell-area averaging is done.
  if (g$spatial_stride != 1L ||
      length(first$lon) < 2L || length(first$lat) < 2L ||
      any(!is.finite(diff(first$lon))) || any(!is.finite(diff(first$lat))) ||
      any(abs(abs(diff(first$lon)) - 0.25) > 1e-6) ||
      any(abs(abs(diff(first$lat)) - 0.25) > 1e-6))
    stop("This 0.25-degree script requires stride=1 on regular native 0.25-degree ",
         "longitude and latitude axes; check the ERA5 input coordinates.", call. = FALSE)
  ilon <- seq.int(1L, length(first$lon), by = g$spatial_stride)
  latitude_candidates <- which(first$lat >= g$lat_range[1L] & first$lat <= g$lat_range[2L])
  if (!length(latitude_candidates)) stop("No latitude coordinates in requested domain.")
  ilat <- latitude_candidates[seq.int(1L, length(latitude_candidates), by = g$spatial_stride)]
  .set_analysis_stage("land_mask", announce = TRUE)
  land <- load_land_fraction(g, first$lon, first$lat, ilon, ilat)
  grid <- land$points
  if (!any(grid$land)) stop("No land grid points selected.")
  if (!is.null(g$max_cells)) {
    chosen <- which(grid$land)
    chosen <- chosen[unique(round(seq(1, length(chosen), length.out = min(length(chosen), g$max_cells))))]
    grid$land <- seq_len(nrow(grid)) %in% chosen
    warning("PILOT subset: these are not full-domain global estimates.")
  }
  grid$selected_lon_index <- rep(seq_along(ilon), times = length(ilat))
  grid$selected_lat_index <- rep(seq_along(ilat), each = length(ilon))
  tasks <- fast_make_tasks(grid, ilon, ilat, g$tile_side)
  cell_counts <- vapply(tasks, function(x) nrow(x$tile_points), integer(1))
  planned <- vapply(inspected$files, function(x) length(x$dates), integer(1))
  periods <- c(as.character(g$years), paste0(min(g$years), "-", max(g$years)))
  planned <- c(planned, sum(planned))
  .set_analysis_stage("scan_exact_daily_cache_support", announce = TRUE)
  capacity <- fast_daily_support_preflight(daily_settings, tasks, inspected$files,
    analysis[c("pred", "training_history", "joint_reference")], g, inspected$audit, ndraw, out_dir)
  message("Formal global run: ", length(tasks), " spatial tiles; ", sum(cell_counts),
    " selected land points; ", sum(head(planned, -1L)), " daily dates; ", ndraw,
    " fixed paired posterior draws; lag1-3 retained. Only P50 TSV; M, positive N=1/S and U=M+N computed analytically.")
  .set_analysis_stage("initialize_fast_kernel", announce = TRUE)
  cache <- fast_make_cache(analysis, g$weight_scheme)
  cache_dir <- file.path(config$output_dir, "kernel_cache", paste0("R-", getRversion(), "-", R.version$arch))
  # A resume retains its original backend even if a compiler is installed later.
  manifest_path <- file.path(out_dir, "fast_manifest.rds")
  old_manifest <- if (file.exists(manifest_path)) readRDS(manifest_path) else NULL
  requested_backend <- f$backend
  if (!is.null(old_manifest) && identical(requested_backend, "auto")) requested_backend <- old_manifest$backend
  backend <- fast_init_backend(cache_dir, requested_backend)
  fast_kernel_QC(analysis, cache, config$output_dir)
  fast_boundary_preflight(cache, analysis, set_bounds, out_dir)
  bundle <- list(g = g, files = inspected$files, fixed_set = fixed_set, cache = cache,
    support_analysis = analysis[c("pred", "training_history", "joint_reference")],
    ndraw = ndraw, periods = periods, planned = planned, set_bounds = set_bounds,
    daily_cache = daily_settings, temperature = temperature,
    capacity_supported_rows_by_tile = capacity$supported_rows_by_tile)
  package_names <- c("dlnm", "FNN", "ncdf4", "maps", "Rcpp", "digest")
  packages <- vapply(package_names, function(p) if (requireNamespace(p, quietly = TRUE))
    as.character(utils::packageVersion(p)) else "unavailable", character(1))
  scientific_functions <- c("fast_compute_tile", "fast_make_cache", "fast_predict_ab", ".fast_r_kernel",
    ".fast_shift_probability_r", ".fast_cached_exp", "make_history_design", ".b21_basis", ".tv_basis",
    "history_support", "evaluate_support", "read_rectangular_history", "sum_by_cell", "posterior_rows", "fast_warm_metrics")
  all_symbols <- ls(envir = environment(), all.names = TRUE)
  daily_env <- environment(fast_daily_open)
  daily_symbols <- ls(envir = daily_env, all.names = TRUE)
  scientific_functions <- unique(c(scientific_functions,
    daily_symbols[grepl("^fast_daily_", daily_symbols) &
      vapply(daily_symbols, function(nm) is.function(get(nm, envir = daily_env)), logical(1))]))
  implementation <- lapply(scientific_functions, function(nm) deparse(body(get(nm, mode = "function"))))
  names(implementation) <- scientific_functions
  signature <- list(schema = "ERA5-warm02-run-v1", revision = ANALYSIS_REVISION,
    implementation = implementation, cpp_source = FAST_CPP_CODE,
    global_config = g, daily_cache = fast_daily_identity_settings(daily_settings), temperature = temperature,
    model_input_md5 = analysis$input_md5,
    model_cache = cache, support_reference = bundle$support_analysis, fixed_set = fixed_set,
    set_bounds = set_bounds, input_audit = inspected$audit, files = inspected$files, grid = grid,
    land_method = land$method, backend = backend, packages = packages, R = as.character(getRversion()))
  fingerprint <- fast_object_hash(signature)
  rm(signature)
  if (!is.null(old_manifest) && (!identical(old_manifest$schema, "ERA5-warm02-run-v1") ||
      !identical(old_manifest$fingerprint, fingerprint)))
    stop("Resume fingerprint differs (inputs, scientific settings, software, backend, or draw IDs). ",
         "Set config$fast$resume = 'new' for a separate run; committed results are unchanged.")
  if (is.null(old_manifest)) {
    if (length(list.files(checkpoint_dir, pattern = "^tile_[0-9]+[.]rds$")))
      stop("Checkpoints exist without an immutable fast manifest; cannot establish draw pairing.")
    manifest <- list(schema = "ERA5-warm02-run-v1", revision = ANALYSIS_REVISION,
      fingerprint = fingerprint, backend = backend, created = Sys.time(),
      total_tiles = length(tasks), total_land_points = sum(cell_counts), g = g,
      paired_draw_ids = data.frame(m1 = cache$draw_ids_m1, m0 = cache$draw_ids_m0),
      fixed_set = fixed_set, set_bounds = set_bounds, packages = packages, input_audit = inspected$audit,
      model_input_md5 = analysis$input_md5)
    fast_atomic_save_rds(manifest, manifest_path)
  }
  fast_daily_snapshot(config, analysis, bundle, tasks, fingerprint, out_dir)
  utils::write.csv(inspected$audit, file.path(out_dir, "input_audit.csv"), row.names = FALSE)
  utils::write.csv(grid[grid$land, ], file.path(out_dir, "selected_land_points.csv"), row.names = FALSE)
  checkpoint_paths <- file.path(checkpoint_dir, sprintf("tile_%05d.rds", seq_along(tasks)))
  done <- logical(length(tasks))
  .set_analysis_stage("validate_resume_checkpoints", announce = TRUE)
  for (i in seq_along(tasks))
    fast_daily_recover_warm(tasks[[i]], bundle, fingerprint, checkpoint_paths[i], part_dir)
  for (i in which(file.exists(checkpoint_paths))) {
    cp <- fast_read_checkpoint(checkpoint_paths[i], fingerprint, tasks[[i]], bundle)
    fast_restore_partition(cp, part_dir); done[i] <- TRUE
  }
  initial_done <- sum(done)
  message("Validated resume: ", initial_done, "/", length(tasks), " committed tiles; ",
          sum(cell_counts[done]), "/", sum(cell_counts), " land points.")
  start <- Sys.time(); benchmark <- data.frame(); progress <- list(); chosen_workers <- f$workers
  rate <- NA_real_; elapsed_new <- 0; cells_new <- 0L
  record_progress <- function(label, workers, seconds = NA_real_, new_cells = 0L) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    remaining <- sum(cell_counts[!done])
    eta <- if (remaining == 0) 0 else if (is.finite(rate) && rate > 0) remaining / rate else NA_real_
    row <- data.frame(time = as.character(Sys.time()), phase = label, backend = backend, workers = workers,
      completed_tiles = sum(done), total_tiles = length(tasks), completed_land_points = sum(cell_counts[done]),
      total_land_points = sum(cell_counts), seconds_this_wave = seconds, land_points_this_wave = new_cells,
      elapsed_this_session_seconds = elapsed, land_points_per_second = rate,
      estimated_remaining_compute_hours = eta / 3600,
      estimated_remaining_recent_hours = if (remaining==0) 0 else if (is.finite(seconds) && seconds>0 && new_cells>0) remaining/(new_cells/seconds)/3600 else NA_real_,
      elapsed_from_launch_hours = as.numeric(difftime(Sys.time(), .warm_session_started, units="hours")),
      projected_total_session_hours_excluding_finalize = as.numeric(difftime(Sys.time(), .warm_session_started, units="hours")) + eta/3600)
    progress[[length(progress) + 1L]] <<- row
    utils::write.csv(do.call(rbind, progress), file.path(out_dir, "progress.csv"), row.names = FALSE)
    message("Progress ", sum(done), "/", length(tasks), " tiles, ", sum(cell_counts[done]), "/", sum(cell_counts),
      " land points; ", workers, " workers; remaining compute ~",
      if (is.finite(eta)) sprintf("%.2f h", eta / 3600) else "pending benchmark",
      ". ETA INCLUDES cache processing/write time in measured waves; excludes final aggregation/plots and may change with geography/I/O.")
    if (is.finite(eta) && remaining>0) {
      forecast <- row$projected_total_session_hours_excluding_finalize
      message(sprintf("Measured forecast for this session: %.2f h including elapsed setup, excluding final aggregation/plots.", forecast))
      if (is.finite(f$target_total_hours) && forecast > f$target_total_hours)
        message("Current measured speed exceeds the ", f$target_total_hours, " h planning target; calculation continues without reducing draws or dates.")
    }
    invisible(row)
  }
  cl <- NULL
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (any(!done)) {
    candidates <- if (isTRUE(f$autotune)) f$worker_candidates else f$workers
    physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (is.finite(physical) && physical > 0) candidates <- pmin(candidates, as.integer(physical))
    candidates <- sort(unique(pmin(candidates, sum(!done))))
    chosen_workers <- candidates[1L]
    .set_analysis_stage("start_PSOCK_workers", workers = max(candidates), announce = TRUE)
    cl <- parallel::makePSOCKcluster(max(candidates), rscript_args = "--vanilla",
      outfile = file.path(out_dir, "worker_log.txt"))
    source_env <- environment(run_global_era5)
    symbols <- ls(envir = source_env, all.names = TRUE)
    fast_functions <- symbols[grepl("^(fast_|[.]fast_)", symbols) &
      vapply(symbols, function(x) is.function(get(x, envir = source_env)), logical(1))]
    exports <- unique(c(fast_functions, "FAST_CPP_CODE", "era5_required_variables", "history_columns",
      ".history_frame", ".b21_basis", ".tv_basis", "make_history_design", ".make_joint_reference",
      "history_support", "evaluate_support", "normal_axis", "read_rectangular_history", "sum_by_cell",
      "posterior_rows", ".set_analysis_stage", ".analysis_error_details"))
    parallel::clusterExport(cl, exports, envir = source_env)
    qc_history <- .fast_qc_histories(analysis)$training
    qc_history <- qc_history[seq_len(min(2L, nrow(qc_history))), , drop = FALSE]
    qc_expected <- fast_predict_ab(cache, qc_history)[c("A1", "B1", "A0", "B0")]
    worker_audit <- vector("list", length(cl))
    # Serial initialization prevents competing sourceCpp cache writes.
    for (i in seq_along(cl)) worker_audit[[i]] <- parallel::clusterCall(cl[i], fast_worker_initialize,
      .libPaths(), cache_dir, backend, bundle, part_dir, checkpoint_dir, fingerprint, qc_history, qc_expected)[[1L]]
    fast_register_workers(config$output_dir, worker_audit)
    saveRDS(worker_audit, file.path(out_dir, "worker_QC.rds"))
    message("All ", length(cl), " workers passed prediction QC; backend=", backend,
            ". Complete production tiles will be used for speed measurement.")
    run_wave <- function(ids, workers, phase) {
      .set_analysis_stage(phase, tile_ids = ids, workers = workers, announce = TRUE)
      message("Computing ", length(ids), " tiles (", sum(cell_counts[ids]), " land points) with ", workers,
              " workers. Individual tile commits are recorded in global/checkpoints.")
      t0 <- Sys.time()
      result <- parallel::clusterApplyLB(cl[seq_len(workers)], tasks[ids], fast_worker_task)
      seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      returned <- vapply(result, function(x) as.integer(x$tile_id), integer(1))
      if (!identical(returned, as.integer(ids))) stop("Worker result tile order/IDs differ.")
      for (i in ids) if (!file.exists(checkpoint_paths[i])) stop("Worker returned without a committed checkpoint.")
      done[ids] <<- TRUE
      list(seconds = seconds, cells = sum(cell_counts[ids]), rate = sum(cell_counts[ids]) / max(seconds, .001))
    }
    benchmark_file <- file.path(out_dir, "worker_benchmark.rds")
    saved_benchmark <- if (file.exists(benchmark_file)) readRDS(benchmark_file) else NULL
    reuse_benchmark <- isTRUE(f$autotune) && !is.null(saved_benchmark) &&
      identical(saved_benchmark$fingerprint, fingerprint) && identical(saved_benchmark$candidates, candidates)
    if (reuse_benchmark) {
      benchmark <- saved_benchmark$table
      chosen_workers <- saved_benchmark$chosen_workers
      rate <- saved_benchmark$chosen_rate
      message("Reusing completed 8/12-worker selection from this unchanged run: ", chosen_workers, " workers.")
    } else if (isTRUE(f$autotune) && length(candidates) > 1L && sum(!done) >= sum(candidates)) {
      for (workers in candidates) {
        remaining <- which(!done)
        full <- remaining[cell_counts[remaining] == max(cell_counts[remaining])]
        pool <- if (length(full) >= workers) full else remaining
        ids <- pool[unique(as.integer(round(seq(1, length(pool), length.out = workers))))]
        wave <- run_wave(ids, workers, "benchmark_production_tiles")
        row <- data.frame(workers = workers, tiles = length(ids), land_points = wave$cells,
          wall_seconds = wave$seconds, land_points_per_second = wave$rate)
        benchmark <- rbind(benchmark, row); rate <- wave$rate
        utils::write.csv(benchmark, file.path(out_dir, "worker_benchmark.csv"), row.names = FALSE)
        record_progress("benchmark", workers, wave$seconds, wave$cells)
      }
      best <- which.max(benchmark$land_points_per_second)
      chosen_workers <- benchmark$workers[best]; rate <- benchmark$land_points_per_second[best]
      fast_atomic_save_rds(list(fingerprint = fingerprint, candidates = candidates, table = benchmark,
        chosen_workers = chosen_workers, chosen_rate = rate), benchmark_file, replace = TRUE)
      message("Measured worker selection: ", chosen_workers,
        ". Benchmark tiles are retained as completed scientific output; speeds can vary by region.")
    }
    record_progress("production_start", chosen_workers)
    while (any(!done)) {
      remaining <- which(!done)
      ids <- head(remaining, 2L * chosen_workers)
      wave <- run_wave(ids, chosen_workers, "parallel_global_projection")
      elapsed_new <- elapsed_new + wave$seconds; cells_new <- cells_new + wave$cells
      rate <- cells_new / max(elapsed_new, .001)
      record_progress("production", chosen_workers, wave$seconds, wave$cells)
    }
    parallel::stopCluster(cl); cl <- NULL
  }
  .set_analysis_stage("combine_committed_draw_sums", announce = TRUE)
  totals <- NULL; partition_files <- character(length(tasks))
  fields <- c("aggregate_sums", "aggregate_denominator", "aggregate_cell_counts",
              "support_sums", "support_denominator", "support_cell_counts")
  # Fixed tile-ID order makes final sums independent of worker scheduling.
  for (i in seq_along(tasks)) {
    cp <- fast_read_checkpoint(checkpoint_paths[i], fingerprint, tasks[[i]], bundle)
    partition_files[i] <- fast_restore_partition(cp, part_dir)
    if (is.null(totals)) totals <- cp[fields] else
      for (nm in fields) totals[[nm]] <- totals[[nm]] + cp[[nm]]
  }
  cache_completion <- fast_daily_finalize_inventory(bundle, tasks, fingerprint, out_dir)
  execution <- list(backend = backend, selected_workers = chosen_workers, benchmark = benchmark,
    daily_cache_completion = cache_completion,
    fingerprint = fingerprint, checkpoint_schema = "ERA5-warm02-tile-v1",
    resumed_tiles = initial_done, total_tiles = length(tasks), total_land_points = sum(cell_counts),
    elapsed_projection_session_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
    numerical_QC = "Original official posterior_epred checks plus fast kernel and every worker checked at tolerance <=1e-8.",
    resume = "Only validated per-draw checkpoints of this exact input/configuration/draw fingerprint are combined.",
    reproducibility = "No worker RNG; fixed draw IDs; final tile-ID summation order.")
  fast_finalize_global(config, analysis, bundle, inspected, land, ilon, ilat,
                       totals, partition_files, out_dir, execution)
}


# ======================================================================

# Warm-boundary outputs. N is the POSITIVE SET half-width: N = 1/A.
warm_output_metrics <- function() c("delta_M", "delta_N", "delta_U", "delta_TSV_P50",
  "U_M1", "U_M0", "U_pair_in_SET_range", "absolute_delta_U")

warm_draw_summary <- function(x) {
  if (!length(x) || any(!is.finite(x))) return(c(lower = NA_real_, median = NA_real_,
    upper = NA_real_, prob_positive = NA_real_, prob_negative = NA_real_))
  c(stats::setNames(stats::quantile(x, c(.025, .5, .975), names = FALSE, type = 8),
    c("lower", "median", "upper")), prob_positive = mean(x > 0), prob_negative = mean(x < 0))
}

fast_finalize_global <- function(config, analysis, bundle, inspected, land, ilon, ilat,
                                 totals, partition_files, out_dir, execution) {
  g <- bundle$g; metrics <- warm_output_metrics(); periods <- bundle$periods
  np <- length(periods); ndraw <- as.integer(bundle$ndraw)
  cuts <- g$support_fraction_sensitivity; nc <- length(cuts)
  if (!np || !nc || length(ndraw) != 1L || is.na(ndraw) || ndraw < 1L ||
      length(bundle$planned) != np || any(!is.finite(bundle$planned) | bundle$planned <= 0))
    stop("Invalid warm-boundary finalizer dimensions.")
  if (length(analysis$draw_ids) != ndraw || length(analysis$draw_ids_m0) != ndraw ||
      anyNA(analysis$draw_ids) || anyNA(analysis$draw_ids_m0)) stop("Posterior draw-ID mismatch.")
  dims <- list(aggregate_sums = c(np, 8L, ndraw), support_sums = c(np, 8L, nc, ndraw),
    support_denominator = c(np, nc), support_cell_counts = c(np, nc))
  for (nm in names(dims)) if (!is.numeric(totals[[nm]]) ||
      !identical(as.integer(dim(totals[[nm]])), as.integer(dims[[nm]])) ||
      any(!is.finite(totals[[nm]]))) stop("Invalid warm-boundary total: ", nm)
  for (nm in c("aggregate_denominator", "aggregate_cell_counts")) if (!is.numeric(totals[[nm]]) ||
      length(totals[[nm]]) != np || any(!is.finite(totals[[nm]]))) stop("Invalid total: ", nm)
  for (nm in c("aggregate_denominator", "aggregate_cell_counts", "support_denominator", "support_cell_counts"))
    if (any(totals[[nm]] < 0)) stop("Negative denominator/count: ", nm)
  if (!length(partition_files) || any(!file.exists(partition_files))) stop("Missing committed grid partitions.")
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory.")
  .set_analysis_stage("summarize_global_warm_boundary_posterior", announce = TRUE)
  draws <- array(NA_real_, c(np, 8L, ndraw), dimnames = list(period = periods, metric = metrics,
    paired_draw_position = as.character(seq_len(ndraw))))
  sensitivity <- array(NA_real_, c(np, 8L, nc, ndraw), dimnames = list(period = periods,
    metric = metrics, min_supported_fraction = as.character(cuts),
    paired_draw_position = as.character(seq_len(ndraw))))
  rows <- list(); sens_rows <- list(); ir <- 0L; is <- 0L
  pilot <- !is.null(g$max_cells)
  for (p in seq_len(np)) for (m in seq_along(metrics)) {
    den <- totals$aggregate_denominator[p]
    if (den > 0) draws[p, m, ] <- totals$aggregate_sums[p, m, ] / den
    ir <- ir + 1L; q <- warm_draw_summary(draws[p, m, ])
    rows[[ir]] <- data.frame(period = periods[p], scenario = "common_support", metric = metrics[m],
      fixed_set_name = "P50", fixed_set = unname(bundle$fixed_set[1L]),
      lower = q["lower"], median = q["median"], upper = q["upper"],
      prob_positive = q["prob_positive"], prob_negative = q["prob_negative"],
      weighted_supported_cell_day_denominator = den, eligible_grid_points = totals$aggregate_cell_counts[p],
      pilot_not_global = pilot, row.names = NULL)
    for (j in seq_len(nc)) {
      den_j <- totals$support_denominator[p, j]
      if (den_j > 0) sensitivity[p, m, j, ] <- totals$support_sums[p, m, j, ] / den_j
      is <- is + 1L; qj <- warm_draw_summary(sensitivity[p, m, j, ])
      sens_rows[[is]] <- data.frame(period = periods[p], scenario = "common_support", metric = metrics[m],
        min_supported_fraction = cuts[j], lower = qj["lower"], median = qj["median"], upper = qj["upper"],
        prob_positive = qj["prob_positive"], prob_negative = qj["prob_negative"],
        weighted_supported_cell_day_denominator = den_j,
        eligible_grid_points = totals$support_cell_counts[p, j], pilot_not_global = pilot, row.names = NULL)
    }
  }
  summary <- do.call(rbind, rows); sensitivity_summary <- do.call(rbind, sens_rows)
  utils::write.csv(summary, file.path(out_dir, "warm_boundary_area_day_posterior_summary.csv"), row.names = FALSE)
  utils::write.csv(sensitivity_summary, file.path(out_dir, "warm_boundary_support_sensitivity_summary.csv"), row.names = FALSE)
  paired_ids <- data.frame(position = seq_len(ndraw), m1 = analysis$draw_ids, m0 = analysis$draw_ids_m0)
  draw_output <- list(estimates = draws, sensitivity_estimates = sensitivity, paired_draw_ids = paired_ids,
    aggregate_denominator = totals$aggregate_denominator, aggregate_cell_counts = totals$aggregate_cell_counts,
    support_denominator = totals$support_denominator, support_cell_counts = totals$support_cell_counts,
    note = "Metrics retain common posterior positions; joint uncertainty must use these paired draws, not combine marginal quantiles.")
  saveRDS(draw_output, file.path(out_dir, "warm_boundary_area_day_posterior_draws.rds"), compress = "gzip")
  saveRDS(summary, file.path(out_dir, "warm_boundary_area_day_posterior_summary.rds"))
  first <- inspected$files[[1L]]
  methods <- list(global_config = g, input_audit = inspected$audit, fixed_set = bundle$fixed_set,
    set_bounds = bundle$set_bounds, periods = periods, years = g$years,
    n_expected_total_days = bundle$planned[np], global_draws_actual = ndraw,
    paired_draw_ids = paired_ids, draw_metadata = analysis$draw_metadata,
    metrics = metrics, definitions = c(M = "B/A: neutral SET", N = "1/A: positive SET half-width",
      U = "(B+1)/A = M+N: model-derived SET at TSV +1", delta = "M1 minus M0"),
    fixed_SET_statement = "Only delta_TSV_P50 uses fixed SET. M, N and U are analytic curve parameters/roots and do not depend on P50.",
    root_statement = "Roots are calculated within every day and posterior position before temporal/area averaging; roots are never clipped to the SET range.",
    standardization = paste("Study-sample covariate standardization:", g$weight_scheme),
    posterior_statement = paste("Independent M1 and M0 posterior draws paired once and reused in all space/time:",
      "a working product-posterior coupling, not an identified joint paired-data posterior."),
    primary_estimand = paste("Area-weighted supported grid-cell days, retaining configured finite-weather coverage",
      "and minimum supported-day requirements, with no primary whole-grid support-fraction cutoff.",
      "Each grid mean is conditional on its supported dates, not necessarily an annual climate mean."),
    support_interpretation = paste("History marginal/kNN support is an exposure-density diagnostic.",
      "Both roots falling inside the marginal training SET P1-P99 range is an extrapolation audit,",
      "not validation of conditional SET/history/covariate support; this audit does not filter the primary estimates."),
    pair_SET_fraction_definition = "For each posterior draw, fraction of supported dates on which BOTH U_M1 and U_M0 lie within the saved SET bounds; map shows its posterior median.",
    boundary_interpretation = "TSV +1 is a declared warm-sensation response boundary, not a validated thermal acceptability cutoff or an observed thermostat setting.",
    energy_statement = "Daily paired warm-boundary posterior cache retained in energy_cache for model-free energy scenarios. No SET-to-Ta conversion or actual energy estimate is fitted here.",
    daily_cache_directory = bundle$daily_cache$dir,
    same_day_temperature_status = bundle$temperature$status,
    spatial_method = paste("Point sampling every", g$spatial_stride, "native coordinate indices; not spatial block averaging."),
    grid_resolution = c(lon = abs(diff(first$lon)[1L]) * g$spatial_stride,
      lat = abs(diff(first$lat)[1L]) * g$spatial_stride), land_mask_method = land$method,
    pilot_not_global = pilot, partition_files = partition_files, output_dir = out_dir,
    fast_execution = execution, checkpoint_schema = "ERA5-warm02-tile-v1")
  saveRDS(methods, file.path(out_dir, "global_run_metadata.rds"))
  dput(methods[setdiff(names(methods), c("input_audit", "partition_files", "draw_metadata"))],
    file = file.path(out_dir, "global_run_methods.txt"))
  list(output_dir = out_dir, partition_files = partition_files, aggregate_summary = summary,
    aggregate_draws = draw_output, metadata = methods)
}

warm_read_period <- function(result, period = NULL) {
  directory <- if (is.list(result)) result$output_dir else as.character(result)
  meta <- if (is.list(result)) result$metadata else readRDS(file.path(directory, "global_run_metadata.rds"))
  if (is.null(period)) period <- tail(meta$periods, 1L)
  if (length(period) != 1L || !period %in% meta$periods) stop("Unknown warm-boundary output period.")
  paths <- meta$partition_files
  if (!length(paths) || any(!file.exists(paths))) stop("Grid partitions are missing.")
  pieces <- lapply(paths, function(path) {
    x <- readRDS(path)
    if (!is.data.frame(x) || !"period" %in% names(x)) stop("Invalid grid partition: ", path)
    x[x$period == period, , drop = FALSE]
  })
  d <- do.call(rbind, pieces)
  if (!nrow(d)) stop("No saved grid summaries for ", period)
  if (anyDuplicated(d$cell_id)) stop("Duplicate grid IDs within warm-boundary period.")
  list(data = d, metadata = meta, period = period, directory = directory)
}

warm_plot_preflight <- function() {
  needed <- c("ggplot2", "maps", "scales")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing plotting packages: ", paste(absent, collapse = ", "))
  world <- ggplot2::map_data("world")
  p <- ggplot2::ggplot(data.frame(lon = c(0, .5), lat = c(0, .5), z = c(-1, 1),
      panel = c("M", "N")), ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(ggplot2::aes(fill = z), width = .5, height = .5) +
    ggplot2::geom_polygon(data = world, mapping = ggplot2::aes(x = long, y = lat, group = group),
      inherit.aes = FALSE, fill = NA, colour = "grey40", linewidth = .1) + ggplot2::facet_wrap(~panel)
  invisible(ggplot2::ggplot_build(p)); invisible(p)
}

fast_plot_preflight <- function(output_dir) {
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create preflight directory.")
  p <- warm_plot_preflight()
  stem <- tempfile("warm_plot_preflight_", tmpdir = output_dir)
  files <- paste0(stem, c(".png", ".pdf"))
  on.exit(unlink(files), add = TRUE)
  ggplot2::ggsave(files[1L], p, width = 4, height = 2.5, dpi = 72, bg = "white")
  ggplot2::ggsave(files[2L], p, width = 4, height = 2.5, bg = "white")
  if (any(!file.exists(files)) || any(file.info(files)$size <= 0)) stop("Plot device preflight failed.")
  message("PASS: faceted world-outline plot and PNG/PDF devices.")
  invisible(TRUE)
}

warm_map_base <- function(d, meta, caption, world) {
  res <- meta$grid_resolution
  if (length(res) != 2L || any(!is.finite(res) | res <= 0)) stop("Invalid map grid resolution.")
  ggplot2::ggplot(d, ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(ggplot2::aes(fill = map_value), width = res[["lon"]], height = res[["lat"]]) +
    ggplot2::geom_polygon(data = world, mapping = ggplot2::aes(x = long, y = lat, group = group),
      inherit.aes = FALSE, fill = NA, colour = "grey35", linewidth = .14, show.legend = FALSE) +
    ggplot2::coord_quickmap(xlim = c(-180, 180), ylim = meta$global_config$lat_range, expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, caption = paste(strwrap(caption, width = 145), collapse = "\n")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), legend.position = "bottom",
      plot.caption = ggplot2::element_text(size = 7.5, hjust = 0), strip.text = ggplot2::element_text(face = "bold"))
}

warm_signed_scale <- function(x, title) {
  lim <- max(abs(x[is.finite(x)]), 0)
  if (lim <= 0) lim <- .01
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-lim, lim), na.value = "#E8E8E8", name = title)
}

plot_global_maps <- function(result, period = NULL, output_dir = NULL, width = 11, height = 5.4) {
  warm_plot_preflight(); b <- warm_read_period(result, period); d <- b$data; meta <- b$metadata
  if (is.null(output_dir)) output_dir <- file.path(b$directory, "figures")
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create figures directory.")
  world <- ggplot2::map_data("world"); plots <- list()
  required <- c("eligible_for_area_day_summary", "delta_U_median", "delta_M_median", "delta_N_median",
    "delta_U_lower", "delta_U_upper", "delta_TSV_P50_median", "U_pair_in_SET_range_median",
    "date_coverage", "common_support_fraction", "lon", "lat")
  if (!all(required %in% names(d))) stop("Saved grid summary lacks warm-boundary columns.")
  eligible <- !is.na(d$eligible_for_area_day_summary) & d$eligible_for_area_day_summary
  pilot <- if (isTRUE(meta$pilot_not_global)) "PILOT SUBSET; no global claim. " else ""
  caption <- paste0(pilot, "M1 minus M0; study-sample standardization. Each grid averages its supported dates,",
    " not necessarily a representative full year. Posterior medians from fixed working product-posterior pairing.",
    " Grey: unavailable or insufficient dates. Warm boundaries are model-derived SET roots; extrapolation is audited separately.",
    " Every native 0.25-degree ERA5 centre is used; display tiles represent point estimates, not area-averaged weather.")
  label <- if (isTRUE(meta$pilot_not_global)) "PILOT - " else ""
  for (field in c("delta_U", "delta_TSV_P50", "absolute_delta_U")) {
    x <- d; x$map_value <- x[[paste0(field, "_median")]]; x$map_value[!eligible] <- NA_real_
    if (!any(is.finite(x$map_value))) { warning("No eligible estimates for ", field); next }
    titles <- c(delta_U = "Warm-sensation SET boundary shift", delta_TSV_P50 = "Thermal-sensation contrast at P50",
      absolute_delta_U = "Magnitude of daily warm-boundary shifts")
    units <- if (field == "delta_TSV_P50") "Mean delta TSV\n(posterior median)" else
      if (field == "absolute_delta_U") "Mean |delta U| (degrees C SET)\n(posterior median)" else
        "Mean delta U (degrees C SET)\n(posterior median)"
    note <- if (field == "delta_TSV_P50") paste0("Fixed SET = ", meta$fixed_set[[1L]], " degrees C; P50 reference only") else
      if (field == "absolute_delta_U") "Absolute daily differences averaged within each draw" else
        "TSV = +1; U = M + N; positive delta U means a higher SET boundary"
    scale <- if (field == "absolute_delta_U") ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B",
      na.value = "#E8E8E8", name = units) else warm_signed_scale(x$map_value, units)
    plots[[field]] <- warm_map_base(x, meta, caption, world) + scale +
      ggplot2::labs(title = paste0(label, titles[field], ", ", b$period), subtitle = note)
  }
  components <- do.call(rbind, lapply(c("delta_M", "delta_N"), function(nm) {
    x <- d; x$map_value <- x[[paste0(nm, "_median")]]; x$map_value[!eligible] <- NA_real_
    x$component <- if (nm == "delta_M") "Delta M: neutral SET position" else "Delta N: positive SET half-width"
    x
  }))
  plots$components_M_N <- warm_map_base(components, meta, caption, world) +
    ggplot2::facet_wrap(~component, nrow = 1) + warm_signed_scale(components$map_value, "Mean change (degrees C SET)\n(posterior median)") +
    ggplot2::labs(title = paste0(label, "Warm-boundary components, ", b$period),
      subtitle = "M = B/A; N = 1/A > 0; increased M and increased N both shift U to higher SET")
  for (field in c("date_coverage", "common_support_fraction", "U_pair_in_SET_range_median")) {
    x <- d; x$map_value <- x[[field]]
    if (field == "U_pair_in_SET_range_median") x$map_value[!eligible] <- NA_real_
    title <- switch(field, date_coverage = "Finite-history date coverage",
      common_support_fraction = "History-supported dates / expected dates",
      U_pair_in_SET_range_median = "Both warm boundaries inside the training SET range")
    note <- if (field == "U_pair_in_SET_range_median") paste0("Posterior median fraction of supported dates; marginal SET P1-P99: ",
      paste(format(as.numeric(meta$set_bounds), digits = 4), collapse = " to "), " degrees C.",
      " A diagnostic, not conditional-support validation or a primary mask.") else
      "Denominator: all expected dates. Coverage does not establish seasonally representative support."
    plots[[field]] <- warm_map_base(x, meta, paste0(pilot, note), world) +
      ggplot2::scale_fill_gradient(low = "#FFF5F0", high = "#08519C", limits = c(0, 1),
        labels = scales::label_percent(), na.value = "#E8E8E8", name = "Fraction") +
      ggplot2::labs(title = paste0(label, title, ", ", b$period))
  }
  x <- d; valid <- eligible & is.finite(x$delta_U_lower) & is.finite(x$delta_U_upper)
  classes <- c("95% CrI below zero", "95% CrI includes zero", "95% CrI above zero")
  x$map_value <- factor(ifelse(x$delta_U_upper < 0, classes[1L],
    ifelse(x$delta_U_lower > 0, classes[3L], classes[2L])), levels = classes)
  x$map_value[!valid] <- NA
  plots$delta_U_interval_direction <- warm_map_base(x, meta, paste0(caption,
    " Pointwise intervals only; no multiple-comparison or field-significance claim."), world) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#2166AC", "#BDBDBD", "#B2182B"), classes),
      drop = FALSE, na.value = "#F2F2F2", name = NULL) +
    ggplot2::labs(title = paste0(label, "Uncertainty in warm-boundary direction, ", b$period))
  x <- d[eligible & is.finite(d$delta_M_median) & is.finite(d$delta_N_median) & is.finite(d$delta_U_median), ]
  if (nrow(x)) {
    lim <- max(abs(x$delta_U_median)); if (lim <= 0) lim <- .01
    plots$joint_M_N <- ggplot2::ggplot(x, ggplot2::aes(delta_M_median, delta_N_median, colour = delta_U_median)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = .3) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = .3) +
      ggplot2::geom_abline(intercept = 0, slope = -1, linetype = "dashed", colour = "grey35", linewidth = .4) +
      ggplot2::geom_point(size = .65, alpha = .5) + ggplot2::coord_equal() +
      ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        limits = c(-lim, lim), name = "Median delta U\n(degrees C SET)") +
      ggplot2::labs(x = "Median delta M (degrees C SET)", y = "Median delta N (degrees C SET)",
        title = paste0(label, "Joint position and half-width shifts, ", b$period),
        caption = paste(strwrap(paste0(pilot, "Each point is a grid's pair of marginal posterior medians, not an independent observation or a posterior draw.",
          " Dashed line: delta M + delta N = 0 as a coordinate reference. Median(delta U) need not equal median(delta M) + median(delta N)."), 110), collapse = "\n")) +
      ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "right", plot.caption = ggplot2::element_text(hjust = 0, size = 8))
  }
  tag <- gsub("[^A-Za-z0-9_-]", "_", b$period); paths <- list()
  for (nm in names(plots)) {
    stem <- file.path(output_dir, paste0("global_", tag, "_", nm))
    w <- if (nm == "components_M_N") 13 else if (nm == "joint_M_N") 8 else width
    h <- if (nm == "joint_M_N") 6.5 else height
    ggplot2::ggsave(paste0(stem, ".png"), plots[[nm]], width = w, height = h, dpi = 400, bg = "white")
    ggplot2::ggsave(paste0(stem, ".pdf"), plots[[nm]], width = w, height = h, bg = "white")
    paths[[nm]] <- c(png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf"))
  }
  saveRDS(plots, file.path(output_dir, paste0("global_", tag, "_plot_objects.rds")))
  utils::write.csv(d, file.path(b$directory, paste0("warm_boundary_grid_summary_", tag, ".csv")), row.names = FALSE)
  saveRDS(d, file.path(b$directory, paste0("warm_boundary_grid_summary_", tag, ".rds")))
  invisible(list(plots = plots, paths = paths, output_dir = output_dir))
}

plot_global_annual_summary <- function(result, output_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing ggplot2.")
  directory <- if (is.list(result)) result$output_dir else as.character(result)
  d <- if (is.list(result)) result$aggregate_summary else utils::read.csv(file.path(directory, "warm_boundary_area_day_posterior_summary.csv"))
  d <- d[grepl("^[0-9]{4}$", d$period) & d$metric == "delta_U", , drop = FALSE]
  if (!nrow(d)) { message("No annual warm-boundary summaries; annual figure skipped."); return(invisible(NULL)) }
  d$year <- as.integer(d$period)
  p <- ggplot2::ggplot(d, ggplot2::aes(year, median)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = .3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#416788", alpha = .2) +
    ggplot2::geom_line(colour = "#416788", linewidth = .7) + ggplot2::geom_point(colour = "#416788", size = 1.7) +
    ggplot2::labs(x = "Year", y = "Area-day mean delta U (degrees C SET)", title = "Annual warm-sensation boundary shifts",
      caption = paste0(if (any(d$pilot_not_global)) "PILOT SUBSET. " else "", "M1 minus M0; posterior median and pointwise 95% CrI.\n",
        "Area-weighted supported dates; yearly support domains may differ. No trend or population-weighted inference.")) +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8))
  if (is.null(output_dir)) output_dir <- file.path(directory, "figures")
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create figure directory.")
  for (ext in c("png", "pdf")) ggplot2::ggsave(file.path(output_dir, paste0("global_annual_warm_boundary_summary.", ext)),
    p, width = 9, height = 5.5, dpi = 400, bg = "white")
  invisible(p)
}

# A small diagnostic BEFORE the long run. This is not a global estimate.
fast_boundary_preflight <- function(cache, analysis, set_bounds, out_dir) {
  ids <- .with_seed(20261004L, function()
    sample.int(nrow(analysis$training_history), min(96L, nrow(analysis$training_history))))
  h <- analysis$training_history[ids, , drop=FALSE]
  ab <- fast_predict_ab(cache, h)
  if (any(ab$A1<=0) || any(ab$A0<=0)) stop("Warm-boundary inversion requires positive slopes.")
  u1 <- (ab$B1+1)/ab$A1
  u0 <- (ab$B0+1)/ab$A0
  good1 <- u1>=set_bounds[1L] & u1<=set_bounds[2L]
  good0 <- u0>=set_bounds[1L] & u0<=set_bounds[2L]
  pair <- good1 & matrix(rep(good0, each=nrow(u1)), nrow=nrow(u1))
  q <- t(apply(u1,1L,stats::quantile,probs=c(.025,.5,.975),type=8,names=FALSE))
  diagnostic <- data.frame(training_row=ids,U_M1_lower=q[,1L],U_M1_median=q[,2L],U_M1_upper=q[,3L],
    probability_both_roots_within_marginal_SET_P1_P99=rowMeans(pair),
    SET_P1=set_bounds[1L],SET_P99=set_bounds[2L])
  utils::write.csv(diagnostic,file.path(out_dir,"boundary_preflight_training_histories.csv"),row.names=FALSE)
  message(sprintf("Boundary diagnostic: global training SET P1/P99 = %.2f / %.2f C; sampled-history mean paired root-in-range probability = %.3f.",
    set_bounds[1L],set_bounds[2L],mean(pair)))
  message("This marginal SET check does not establish conditional support. Roots outside it remain labelled model-derived; no clipping is applied.")
  invisible(diagnostic)
}

run_standalone <- function() {
  cli <- commandArgs(trailingOnly=TRUE)
  mode <- MODE
  if (length(cli) && cli[1L] %in% c("global","plots","selftest")) {
    mode <- cli[1L]
    if (length(cli)>=2L) sys.source(cli[2L],envir=.GlobalEnv)
  }
  if (!mode %in% c("global","plots","selftest")) stop("MODE must be global, plots or selftest.")
  config <- config
  if (mode == "selftest") {
    fast_tile_selftest()
    fast_daily_all_selftests()
    message("All daily cache / weather / energy-ordering selftests passed. No model or global projection was run.")
    return(invisible(TRUE))
  }
  if (mode=="plots") {
    root <- normalizePath(config$output_dir,winslash="/",mustWork=TRUE)
    explicit <- fast_defaults(config)$resume
    if (!explicit %in% c("auto","new")) {
      path <- if (grepl("^([A-Za-z]:[/\\\\]|/)",explicit)) explicit else file.path(root,explicit)
    } else {
      pointer <- file.path(root,"latest_run.txt")
      if (!file.exists(pointer)) stop("No warm-analysis latest_run.txt exists. Run MODE='global' first.")
      path <- trimws(readLines(pointer,warn=FALSE,encoding="UTF-8"))
      path <- path[nzchar(path)]
      if (length(path)!=1L) stop("Latest-run pointer must identify one warm run.")
    }
    prior <- fast_existing_run(path,root)
    token <- fast_acquire_run_lock(prior$path)
    on.exit(fast_release_run_lock(prior$path,token),add=TRUE)
    directory <- file.path(prior$path,"global")
    metadata_path <- file.path(directory,"global_run_metadata.rds")
    complete_path <- file.path(directory,"numerical_completion.rds")
    if (!file.exists(metadata_path) || !file.exists(complete_path))
      stop("Warm numerical finalization is incomplete. Resume global mode; committed tiles will be reused.")
    meta <- readRDS(metadata_path)
    completion <- readRDS(complete_path)
    manifest <- readRDS(file.path(directory,"fast_manifest.rds"))
    if (!identical(meta$checkpoint_schema,"ERA5-warm02-tile-v1") ||
        !identical(completion$status,"COMPLETED") ||
        !identical(completion$fingerprint,manifest$fingerprint) ||
        !identical(meta$fast_execution$fingerprint,manifest$fingerprint))
      stop("Saved warm-run identity/completion fingerprint mismatch.")
    result <- list(output_dir=directory,metadata=meta,
      aggregate_summary=utils::read.csv(file.path(directory,"warm_boundary_area_day_posterior_summary.csv")))
    fast_plot_preflight(prior$path)
    plot_global_maps(result)
    plot_global_annual_summary(result)
    saveRDS(list(status="COMPLETED",revision=ANALYSIS_REVISION,finished=Sys.time(),
      numerical_outputs_changed=FALSE),file.path(directory,"plot_recovery_status.rds"))
    message("Plot-only recovery complete. No model or ERA5 prediction was run. Figures: ",file.path(directory,"figures"))
    return(invisible(result))
  }
  destination <- fast_prepare_run(config$output_dir,fast_defaults(config))
  on.exit(fast_release_run_lock(destination$path,destination$locktoken),add=TRUE)
  config$output_root <- destination$root
  config$output_dir <- destination$path
  message(if (destination$resumed) "Resuming warm-boundary run: " else "New warm-boundary run: ",config$output_dir)
  message("M increases toward warmer neutrality. N=1/S is POSITIVE half-width; increasing N means a wider response interval. U=M+N.")
  message("Only P50 is evaluated for the supplementary fixed-SET TSV contrast. Boundary U is independent of P50.")
  audit <- list(status="RUNNING",revision=ANALYSIS_REVISION,started=Sys.time(),mode=mode,config=config)
  audit_file <- file.path(config$output_dir,"standalone_run_status.rds")
  fast_atomic_save_rds(audit,audit_file,replace=TRUE)
  details <- NULL
  warning_messages <- character()
  on.exit({
    writeLines(unique(warning_messages),file.path(config$output_dir,"analysis_warnings.txt"))
    writeLines(capture.output(sessionInfo()),file.path(config$output_dir,"sessionInfo.txt"))
  },add=TRUE)
  tryCatch(withCallingHandlers({
    .set_analysis_stage("check_inputs_and_plotting",announce=TRUE)
    needed <- c("ncdf4","ggplot2","maps","scales","digest")
    absent <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
    if (length(absent)) stop("Missing packages: ",paste(absent,collapse=", "),". Install them in the repaired library before restarting this script.")
    g <- global_defaults(config)
    required <- file.path(g$input_dir,sprintf(g$file_template,g$years))
    if (any(!file.exists(required))) stop("Missing ERA5 thermal-history files: ",paste(required[!file.exists(required)],collapse=", "))
    coordinate_input_preflight(required[1L],config$output_dir)
    fast_plot_preflight(config$output_dir)
    if (isTRUE(fast_defaults(config)$runtime_selftest)) {
      .set_analysis_stage("warm_boundary_accumulator_QC",announce=TRUE)
      fast_runtime_selftest(config$output_dir)
      fast_daily_all_selftests(config$output_dir)
    }
    .set_analysis_stage("load_models_and_official_prediction_QC",announce=TRUE)
    analysis <- load_analysis(config)
    if (!identical(names(analysis$fixed_set),"P50")) stop("P50-only contract violated.")
    utils::write.csv(analysis$basis_qc,file.path(config$output_dir,"basis_QC.csv"),row.names=FALSE)
    saveRDS(list(model_checks=analysis$model_checks,epred_qc=analysis$epred_qc,
      draw_metadata=analysis$draw_metadata,input_md5=analysis$input_md5,
      package_manifest=analysis$package_manifest,fixed_set=analysis$fixed_set,
      history_source_metadata=analysis$history_source_metadata),file.path(config$output_dir,"model_audit.rds"))
    .set_analysis_stage("global_warm_boundary_projection",announce=TRUE)
    result <- run_global_era5(config,analysis)
    fast_atomic_save_rds(list(status="COMPLETED",revision=ANALYSIS_REVISION,
      finished=Sys.time(),fingerprint=result$metadata$fast_execution$fingerprint),
      file.path(result$output_dir,"numerical_completion.rds"),replace=TRUE)
    .set_analysis_stage("plot_warm_boundary_maps",announce=TRUE)
    plot_global_maps(result)
    .set_analysis_stage("plot_warm_boundary_annual_summary",announce=TRUE)
    plot_global_annual_summary(result)
    audit$status <- "COMPLETED"
    audit$finished <- Sys.time()
    audit$total_wall_hours <- as.numeric(difftime(Sys.time(),.warm_session_started,units="hours"))
    fast_atomic_save_rds(audit,audit_file,replace=TRUE)
    message(sprintf("Completed warm-boundary analysis in %.2f h this session. Output: %s",audit$total_wall_hours,config$output_dir))
    invisible(result)
  },warning=function(w) {
    warning_messages <<- c(warning_messages,conditionMessage(w))
  },error=function(e) {
    details <<- .analysis_error_details(e,sys.calls())
  }),error=function(e) {
    audit$status <- "FAILED"; audit$finished <- Sys.time(); audit$error <- conditionMessage(e)
    if (!is.null(details)) audit <- c(audit,details)
    report <- c(paste("Revision:",ANALYSIS_REVISION),paste("Stage:",audit$stage),
      paste("Error:",audit$error),paste("Call:",audit$error_call),
      "Context:",capture.output(str(audit$context)),"Traceback:",audit$traceback)
    path <- file.path(config$output_dir,"analysis_error.txt")
    tryCatch(writeLines(enc2utf8(report),path,useBytes=TRUE),error=function(w)message(conditionMessage(w)))
    tryCatch(fast_atomic_save_rds(audit,audit_file,replace=TRUE),error=function(w)message(conditionMessage(w)))
    message("Analysis failed in stage: ",audit$stage,"\nOriginal error: ",audit$error,"\nError report: ",path)
    stop(e)
  })
}
run_standalone()

})

run_warm02_embedded <- function(program) {
  thread_names <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  old_threads <- Sys.getenv(thread_names, unset = NA_character_)
  on.exit({
    missing <- is.na(old_threads)
    if (any(missing)) Sys.unsetenv(thread_names[missing])
    if (any(!missing)) do.call(Sys.setenv, as.list(old_threads[!missing]))
  }, add = TRUE)
  do.call(Sys.setenv, setNames(as.list(rep("1", length(thread_names))), thread_names))
  exe_name <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  candidates <- unique(c(file.path(R.home("bin"), exe_name),
                         file.path(R.home(), "bin", exe_name),
                         file.path(R.home(), "bin", "x64", exe_name)))
  available <- candidates[file.exists(candidates)]
  if (!length(available)) stop("Rscript was not found in the current R installation: ",
                                R.home(), call. = FALSE)
  rscript_exe <- normalizePath(available[1L], winslash = "/", mustWork = TRUE)
  child_script <- tempfile(pattern = "warm02_clean_", fileext = ".R")
  on.exit(unlink(child_script), add = TRUE)
  script_lines <- enc2utf8(deparse(program, width.cutoff = 500L))
  invisible(parse(text = script_lines, keep.source = FALSE))
  writeLines(script_lines, con = child_script, useBytes = TRUE)
  child_script <- normalizePath(child_script, winslash = "/", mustWork = TRUE)
  quote_type <- if (.Platform$OS.type == "windows") "cmd" else "sh"
  cli <- commandArgs(trailingOnly = TRUE)
  # Forward only the explicit modes supported by the embedded runner.
  # With Source/paste there are no such arguments; MODE="global" is used.
  forwarded <- character()
  if (length(cli) && cli[1L] %in% c("global", "plots", "selftest")) {
    forwarded <- cli[seq_len(min(length(cli), 2L))]
  }
  args <- c("--vanilla",
            "--default-packages=methods,datasets,utils,grDevices,graphics,stats",
            shQuote(child_script, type = quote_type),
            shQuote(forwarded, type = quote_type))
  message("ERA5-warm02: automatically starting a clean R process (no file dialog).")
  message("R executable: ", rscript_exe)
  message("Progress and errors from the analysis will appear below.")
  status <- system2(command = rscript_exe, args = args,
                    stdout = "", stderr = "", wait = TRUE, timeout = 0)
  if (!isTRUE(status == 0L)) stop(
    "ERA5-warm02 child process failed (exit status ", status,
    "). Inspect the FIRST analysis error above.", call. = FALSE)
  message("ERA5-warm02 completed successfully (child exit status 0).")
  invisible(status)
}

run_warm02_embedded(WARM02_PROGRAM)
