#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export ERA5_BUNDLE_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mode="${1:-check}"
config="${2:-config_memory.R}"
mkdir -p logs results
case "$mode" in
  check)
    sha256sum -c SHA256SUMS.txt
    Rscript --vanilla code/06_test_regressions.R
    Rscript --vanilla code/10_test_parallel_runtime.R
    Rscript --vanilla code/12_test_memory_scheduler.R
    Rscript --vanilla code/04_verify_outputs.R "$config" pack
    Rscript --vanilla code/07_check_server.R "$config"
    ;;
  pilot)
    Rscript --vanilla code/03_run_ecs.R "${2:-config_pilot.R}"
    Rscript --vanilla code/04_verify_outputs.R "${2:-config_pilot.R}" run
    ;;
  full)
    Rscript --vanilla code/03_run_ecs.R "$config"
    Rscript --vanilla code/04_verify_outputs.R "$config" run
    ;;
  verify)
    Rscript --vanilla code/04_verify_outputs.R "$config" run
    ;;
  *) printf 'Usage: bash run_server.sh [check|pilot|full|verify] [config.R]\n' >&2; exit 2 ;;
esac
