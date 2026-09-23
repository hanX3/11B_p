#!/usr/bin/env bash
# Runs the cross-section / reaction-channel validation macros and validates
# the resulting ROOT output with validation/cs_channel/validate_*.py.
#
# Usage: bash validation/cs_channel/run_cs_channel_validation.sh
# Must be run from (or discoverable relative to) the project root, since it
# shells out to ./batch.sh.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
cd "${PROJECT_ROOT}"

VALIDATION_DIR="validation/cs_channel"
REPORTS_DIR="${VALIDATION_DIR}/reports"
PLOTS_DIR="${VALIDATION_DIR}/plots"
SUMMARY_FILE="${REPORTS_DIR}/summary.txt"

mkdir -p "${REPORTS_DIR}" "${PLOTS_DIR}"
: > "${SUMMARY_FILE}"

THREADS="${THREADS:-4}"

require_macro() {
  local macro="$1"
  if [[ ! -f "${macro}" ]]; then
    echo "ERROR: macro not found: ${macro}" >&2
    exit 1
  fi
}

latest_reaction_root() {
  ls -t data/reaction_merged_*.root 2>/dev/null | head -1
}

run_scripts() {
  local tag="$1"
  local root_file="$2"
  shift 2
  local extra_args=("$@")

  echo "=== Validating ${root_file} (${tag}) ==="

  # A [FAIL] in one check must not abort the rest of the suite, otherwise a
  # single regression would hide the results of every later macro. Collect
  # exit codes and only fail the overall script at the very end.
  python3 "${VALIDATION_DIR}/validate_cs_components.py" "${root_file}" \
    | tee "${REPORTS_DIR}/${tag}_cs_components.txt" || true

  python3 "${VALIDATION_DIR}/validate_channel_sampling.py" "${root_file}" \
    | tee "${REPORTS_DIR}/${tag}_channel_sampling.txt" || true

  python3 "${VALIDATION_DIR}/validate_event_weight.py" "${root_file}" \
    | tee "${REPORTS_DIR}/${tag}_event_weight.txt" || true

  if [[ "${#extra_args[@]}" -gt 0 ]]; then
    python3 "${VALIDATION_DIR}/validate_gamma_branching.py" "${root_file}" "${extra_args[@]}" \
      | tee "${REPORTS_DIR}/${tag}_gamma_branching.txt" || true
  fi
}

run_one() {
  local macro="$1"
  local tag="$2"
  shift 2
  local extra_args=("$@")

  require_macro "${macro}"
  echo "=== Running ${macro} ==="
  ./batch.sh "${THREADS}" "${macro}"

  local root_file
  root_file="$(latest_reaction_root)"
  if [[ -z "${root_file}" ]]; then
    echo "ERROR: no reaction_merged_*.root file produced by ${macro}" >&2
    exit 1
  fi

  run_scripts "${tag}" "${root_file}" "${extra_args[@]}"
}

run_one macros/validation_directdecay_fraction0.mac directdecay_fraction0
run_one macros/validation_directdecay_default.mac directdecay_default
run_one macros/validation_directdecay_fraction10pct.mac directdecay_fraction10pct
run_one macros/validation_675scale_half.mac scale675_half
run_one macros/validation_675scale_double.mac scale675_double
run_one macros/validation_165gamma_only.mac gamma165 --macro-hint 165gamma_only
run_one macros/validation_675gamma_only.mac gamma675 --macro-hint 675gamma_only

echo "============================================================"
echo "Collecting summary"
echo "============================================================"

{
  echo "cs_channel validation summary"
  echo "generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo
  for report in "${REPORTS_DIR}"/*.txt; do
    [[ "$(basename "${report}")" == "summary.txt" ]] && continue
    echo "--- $(basename "${report}") ---"
    grep -E '^\[(PASS|WARN|FAIL)\]' "${report}" || echo "  (no PASS/WARN/FAIL lines found)"
    echo
  done
} | tee "${SUMMARY_FILE}"

echo "Summary written to ${SUMMARY_FILE}"

if grep -q '^\[FAIL\]' "${SUMMARY_FILE}"; then
  echo "One or more checks FAILED; see ${SUMMARY_FILE}" >&2
  exit 1
fi
