#!/usr/bin/env bash
# Isolated launcher tests: no real R, Micromamba, models or data are modified.
set -euo pipefail
source_bundle=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
test_root=$(mktemp -d "${TMPDIR:-/tmp}/era5_launcher.XXXXXX")
trap 'rm -rf -- "$test_root"' EXIT
bundle="$test_root/project with spaces/03_server_bundle"
mkdir -p "$bundle" "$test_root/mamba/envs/era5/bin" "$test_root/system"
cp -- "$source_bundle/run_server.sh" "$bundle/"
export LAUNCHER_TEST_LOG="$test_root/calls"
cat > "$test_root/mamba/envs/era5/bin/Rscript" <<'SCRIPT'
#!/usr/bin/env bash
printf 'env_R:%s\n' "$*" >> "$LAUNCHER_TEST_LOG"
if [[ "$*" == --version ]]; then printf 'Rscript mock 4.5.3\n'; fi
if [[ "${TEST_FAIL_COMPUTE:-0}" == 1 && "$*" == *03_run_ecs.R* ]]; then exit 17; fi
SCRIPT
cat > "$test_root/system/Rscript" <<'SCRIPT'
#!/usr/bin/env bash
printf 'system_R:%s\n' "$*" >> "$LAUNCHER_TEST_LOG"
SCRIPT
cat > "$test_root/micromamba" <<'SCRIPT'
#!/usr/bin/env bash
set -euo pipefail
[[ "$1" == run && "$2" == -p ]]
export PATH="$3/bin:$PATH"
printf 'mamba:%s\n' "$3" >> "$LAUNCHER_TEST_LOG"
shift 3
exec "$@"
SCRIPT
chmod +x "$test_root/micromamba" "$test_root/system/Rscript" "$test_root/mamba/envs/era5/bin/Rscript"
export MAMBA_ROOT_PREFIX="$test_root/mamba" MAMBA_EXE="$test_root/micromamba"
export PATH="$test_root/system:$PATH"
unset ERA5_ENV_READY ERA5_MAMBA_ENV TEST_FAIL_COMPUTE
export ERA5_R_ENV=auto
checks=0
pass() { checks=$((checks+1)); printf 'PASS: %s\n' "$1"; }
cd -- "$test_root"
bash "$bundle/run_server.sh" env > "$test_root/env.log" 2>&1
[[ $(grep -c '^mamba:' "$LAUNCHER_TEST_LOG") == 1 ]]
grep -q '^env_R:--version$' "$LAUNCHER_TEST_LOG"
! grep -q '^system_R:' "$LAUNCHER_TEST_LOG"
grep -Fq "Bundle: $bundle" "$test_root/env.log"
pass 'fresh shell selects the existing Micromamba R over system R and resolves a spaced bundle path'
: > "$LAUNCHER_TEST_LOG"
bash "$bundle/run_server.sh" full > "$test_root/full.log" 2>&1
grep -q '^env_R:--vanilla code/03_run_ecs.R config_memory.R$' "$LAUNCHER_TEST_LOG"
grep -q '^env_R:--vanilla code/04_verify_outputs.R config_memory.R run$' "$LAUNCHER_TEST_LOG"
pass 'full uses the same environment for compute and verification'
: > "$LAUNCHER_TEST_LOG"
bash "$bundle/run_server.sh" pilot > "$test_root/pilot.log" 2>&1
grep -q '^env_R:--vanilla code/03_run_ecs.R config_pilot.R$' "$LAUNCHER_TEST_LOG"
bash "$bundle/run_server.sh" verify custom.R > "$test_root/verify.log" 2>&1
grep -q '^env_R:--vanilla code/04_verify_outputs.R custom.R run$' "$LAUNCHER_TEST_LOG"
pass 'pilot default and explicit configuration are forwarded intact'
: > "$LAUNCHER_TEST_LOG"
ERA5_R_ENV=current bash "$bundle/run_server.sh" env > "$test_root/current.log" 2>&1
grep -q '^system_R:--version$' "$LAUNCHER_TEST_LOG"
! grep -q '^mamba:' "$LAUNCHER_TEST_LOG"
pass 'explicit current environment selection avoids Micromamba'
if ERA5_R_ENV=micromamba MAMBA_ROOT_PREFIX="$test_root/missing" bash "$bundle/run_server.sh" env > "$test_root/missing.log" 2>&1; then exit 1; fi
grep -q 'Micromamba or environment missing' "$test_root/missing.log"
pass 'explicit missing Micromamba environment fails clearly'
: > "$LAUNCHER_TEST_LOG"
status=0
TEST_FAIL_COMPUTE=1 bash "$bundle/run_server.sh" full > "$test_root/fail.log" 2>&1 || status=$?
[[ "$status" == 17 ]]
! grep -q '04_verify_outputs.R' "$LAUNCHER_TEST_LOG"
pass 'compute failure keeps its exit status and skips verification'
if bash "$bundle/run_server.sh" invalid > "$test_root/invalid.log" 2>&1; then exit 1; fi
grep -q 'Usage:' "$test_root/invalid.log"
pass 'invalid mode fails before launching R'
printf 'SERVER LAUNCHER TESTS PASS: %s checks; mocked R/Micromamba.\n' "$checks"
