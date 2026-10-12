# Shared local/server configuration. Launchers set ERA5_BUNDLE_ROOT automatically.
# Scientific settings must be explicit; no default cooling scenario is assumed.
bundle_root <- Sys.getenv("ERA5_BUNDLE_ROOT",unset="/data/era5")
local_project <- "E:/ERA5/ERA5_workflow"
ERA5_CONFIG <- list(
  local = list(
    root = local_project,
    mask_rds = file.path(local_project,"02_local/mask/mask_0p25deg_pop_gt1.rds"),
    history_dir = file.path(local_project,"01_source_data/02_Thermal history"),
    temperature_dir = file.path(local_project,"01_source_data/01_raw_ERA5"),
    pack_dir = file.path(bundle_root,"weather")
  ),
  remote = list(
    root = bundle_root,
    pack_dir = file.path(bundle_root,"weather"),
    output_dir = file.path(bundle_root,"results"),
    m1_file = file.path(bundle_root,"models/fit_01-50-lag3.rds"),
    m0_file = file.path(bundle_root,"models/fit_01-M0-base.rds"),
    data_file = file.path(bundle_root,"models/data_01-50-lag3.rds"),
    pred_file = file.path(bundle_root,"models/pred_01-50-lag3.rds")
  ),
  pack = list(
    years = 2010:2019, tile_side = 12L, latitude_range = c(-60, 80),
    expected_selected_cells = 104344L, latitude_band = 90L,
    netcdf_compression_level = 1L,
    # Hash original input bytes once locally; can take minutes for a large inventory.
    checksum_original_inputs = TRUE
  ),
  run = list(
    n_draws = 2000L,
    workers = 32L,
    scheduling = "lpt", # Continuous queue; longest estimated tile first.
    progress_seconds = 30, # Update on completed tiles, at most this often.
    worker_candidates = c(16L, 32L, 48L),
    autotune = FALSE, # First benchmark with a fixed worker count.
    backend = "cpp", resume = "auto", day_chunk = 31L, max_cells = NULL,
    min_date_coverage = 0.95, min_supported_days = 30L,
    min_supported_fraction = 0, support_fraction_sensitivity = c(0.5, 0.8, 0.95),
    joint_support = TRUE, weight_scheme = "observation",
    save_annual_warm_grid = FALSE, make_plots = FALSE,
    # Fill after auditing the preprocessing metadata in source NetCDF files.
    preprocessing_acknowledged = FALSE
  ),
  cdd = list(
    # Required before ECS production. Example FIELD NAMES only (not recommended values):
    # scenarios = list(list(name = "declared_scenario", theta0 = ..., slope = ..., kappa = ...))
    # theta0: degree C Ta balance point; slope: positive SET-on-Ta slope;
    # kappa: >=0 degree SET shift transfer multiplier.
    # Each can be a scalar or a 2000-value posterior vector; do not resample per grid.
    scenarios = list()
  )
)
