#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export ERA5_BUNDLE_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mode="${1:-check}"
config="${2:-config_memory.R}"
case "$mode" in
  check|pilot|full|verify|env) ;;
  *) printf 'Usage: bash run_server.sh [env|check|pilot|full|verify] [config.R]\n' >&2; exit 2 ;;
esac
# Noninteractive activation also works under nohup and in a fresh SSH session.
# Prefer the existing ECS environment even when a different system R is on PATH.
env_mode="${ERA5_R_ENV:-auto}"
case "$env_mode" in
  auto|current|micromamba) ;;
  *) printf 'ERA5_R_ENV must be auto, current or micromamba.\n' >&2; exit 2 ;;
esac
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/data/era5/micromamba}"
env_prefix="$MAMBA_ROOT_PREFIX/envs/${ERA5_MAMBA_ENV:-era5}"
if [[ "${ERA5_ENV_READY:-0}" != 1 && "$env_mode" != current ]]; then
  if [[ "$env_mode" == micromamba || -d "$env_prefix" ]]; then
    mamba_bin="${MAMBA_EXE:-$HOME/.local/bin/micromamba}"
    if [[ ! -x "$mamba_bin" ]]; then mamba_bin="$(command -v micromamba || true)"; fi
    if [[ ! -x "$mamba_bin" || ! -d "$env_prefix" ]]; then
      printf 'Micromamba or environment missing: %s. Check MAMBA_EXE and MAMBA_ROOT_PREFIX.\n' "$env_prefix" >&2
      exit 1
    fi
    export ERA5_ENV_READY=1
    exec "$mamba_bin" run -p "$env_prefix" bash "$ERA5_BUNDLE_ROOT/run_server.sh" "$@"
  fi
fi
rscript="$(command -v Rscript || true)"
if [[ -z "$rscript" ]]; then
  printf 'Rscript unavailable. Expected Micromamba environment: %s. See DEPLOY_MEMORY.txt.\n' "$env_prefix" >&2
  exit 1
fi
printf 'Bundle: %s\nRscript: %s\n' "$ERA5_BUNDLE_ROOT" "$rscript"
"$rscript" --version
if [[ "$mode" == env ]]; then exit 0; fi
mkdir -p logs results
case "$mode" in
  check)
    sha256sum -c SHA256SUMS.txt
    bash code/13_test_server_launcher.sh
    "$rscript" --vanilla code/06_test_regressions.R
    "$rscript" --vanilla code/10_test_parallel_runtime.R
    "$rscript" --vanilla code/12_test_memory_scheduler.R
    "$rscript" --vanilla code/04_verify_outputs.R "$config" pack
    "$rscript" --vanilla code/07_check_server.R "$config"
    ;;
  pilot)
    "$rscript" --vanilla code/03_run_ecs.R "${2:-config_pilot.R}"
    "$rscript" --vanilla code/04_verify_outputs.R "${2:-config_pilot.R}" run
    ;;
  full)
    "$rscript" --vanilla code/03_run_ecs.R "$config"
    "$rscript" --vanilla code/04_verify_outputs.R "$config" run
    ;;
  verify)
    "$rscript" --vanilla code/04_verify_outputs.R "$config" run
    ;;
esac
