#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export ERA5_BUNDLE_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mode="${1:-check}"
mkdir -p logs results
case "$mode" in
  check)
    sha256sum -c SHA256SUMS.txt
    Rscript --vanilla code/06_test_regressions.R
    Rscript --vanilla code/10_test_parallel_runtime.R
    Rscript --vanilla code/04_verify_outputs.R config.R pack
    Rscript --vanilla code/07_check_server.R config.R
    ;;
  pilot)
    Rscript --vanilla code/03_run_ecs.R config_pilot.R
    Rscript --vanilla code/04_verify_outputs.R config_pilot.R run
    ;;
  full)
    Rscript --vanilla code/03_run_ecs.R config.R
    Rscript --vanilla code/04_verify_outputs.R config.R run
    ;;
  verify)
    Rscript --vanilla code/04_verify_outputs.R config.R run
    ;;
  *) printf 'Usage: bash run_server.sh [check|pilot|full|verify]\n' >&2; exit 2 ;;
esac
