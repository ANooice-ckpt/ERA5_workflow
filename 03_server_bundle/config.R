# Shared local/server configuration. Launchers set ERA5_BUNDLE_ROOT automatically.
# Scientific settings are explicit; the CDD scenario below was declared on 2026-10-12.
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
    # Declared 2026-10-12: baseline 18 C; supplied image gives delta_Ta = 0.9683 * delta_U.
    # The engine uses delta_theta = (kappa / slope) * delta_U, so slope is 1/lambda.
    # kappa=1 explicitly assumes equal indoor-shift to outdoor-CDD-threshold transfer.
    # The supplied intercept 0.583 cancels in differences; it is not added to theta0.
    # theta0: outdoor CDD balance point (C); slope: positive SET-on-indoor-Ta slope.
    # Each can be a scalar or a 2000-value posterior vector; do not resample per grid.
    scenarios = list(list(
      name = "CDD18_lambda0.9683_kappa1",
      theta0 = 18,
      slope = 1 / 0.9683,
      kappa = 1
    ))
  )
)
