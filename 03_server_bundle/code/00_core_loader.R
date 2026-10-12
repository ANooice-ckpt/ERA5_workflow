# Load scientific definitions from the supplied warm02 script WITHOUT executing its runner.
# Evaluates only top-level assignments within WARM02_PROGRAM, excluding old config.
# Source this file; call stream_load_warm02(script_path).
stream_this_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE)[1L])
  if (length(script) != 1L || is.na(script)) stop("Launch via Rscript --vanilla <script> <config.R>.")
  dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
}

stream_load_warm02 <- function(path, target = .GlobalEnv) {
  stopifnot(file.exists(path))
  exprs <- as.list(parse(path, encoding = "UTF-8", keep.source = FALSE))
  assignment_name <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.symbol(e[[2L]]))
      as.character(e[[2L]]) else ""
  }
  keys <- vapply(exprs, assignment_name, character(1))
  hit <- which(keys == "WARM02_PROGRAM")
  if (length(hit) != 1L) stop("Expected exactly one WARM02_PROGRAM in reference source.")
  wrapped <- exprs[[hit]][[3L]]
  if (!is.call(wrapped) || !identical(wrapped[[1L]], as.name("quote")))
    stop("Reference WARM02_PROGRAM must be quoted; never source it as a runner.")
  block <- wrapped[[2L]]
  if (!is.call(block) || !identical(block[[1L]], as.name("{"))) stop("Malformed quoted core.")
  defs <- as.list(block)[-1L]
  allowed_values <- c("MODE", "ANALYSIS_REVISION", "PLOT_REVISION", ".warm_session_started",
                      ".analysis_progress", "prepare_model_bundle", "era5_required_variables", "FAST_CPP_CODE")
  loaded <- character()
  for (e in defs) {
    nm <- assignment_name(e)
    if (!nzchar(nm) || identical(nm, "config")) next
    if ((is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) || nm %in% allowed_values) {
      eval(e, envir = target)
      loaded <- c(loaded, nm)
    }
  }
  if (!all(c("load_analysis", "fast_predict_ab", "fast_compute_tile", "run_global_era5",
             "fast_finalize_global") %in% loaded)) stop("Reference definitions incomplete.")
  stage <- get(".analysis_progress", envir = target)
  stage$stage <- "stream_startup"
  stage$context <- list()
  invisible(loaded)
}

stream_read_config <- function(path) {
  e <- new.env(parent = .GlobalEnv)
  sys.source(path, envir = e, keep.source = FALSE)
  if (!exists("ERA5_CONFIG", envir = e, inherits = FALSE))
    stop("The configuration must define ERA5_CONFIG.")
  get("ERA5_CONFIG", envir = e, inherits = FALSE)
}
