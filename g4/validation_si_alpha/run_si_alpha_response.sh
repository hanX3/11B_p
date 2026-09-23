#!/usr/bin/env bash
set -euo pipefail

# Run Si alpha-response validation for 165-keV and 675-keV model-only cases.
# Usage from project root:
#   ./validation_si_alpha/run_si_alpha_response.sh
#
# Fixed defaults for this validation:
#   - Geant4 threads: 8
#   - Si selection for plots: all Si rings
#   - spectra: unweighted hit counts

THREADS=8
RING="all"
WEIGHTED=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
PLOT_SCRIPT="${SCRIPT_DIR}/plot_si_alpha_response.py"

mkdir -p "${DATA_DIR}"

if [[ ! -x "${PROJECT_ROOT}/batch.sh" ]]; then
  echo "ERROR: ${PROJECT_ROOT}/batch.sh not found or not executable." >&2
  exit 1
fi

if [[ ! -x "${PLOT_SCRIPT}" ]]; then
  echo "ERROR: ${PLOT_SCRIPT} not found or not executable." >&2
  exit 1
fi

snapshot_root_files() {
  find "${DATA_DIR}" -maxdepth 1 -type f -name "*.root" -print | sort -u
}

find_newest_matching() {
  local list_file="$1"
  local prefix="$2"
  awk -v p="/${prefix}_" '$0 ~ p && $0 ~ /\.root$/ {print}' "${list_file}" | sort | tail -n 1
}

run_case() {
  local label="$1"
  local macro="$2"
  local channel="$3"
  local before after new_files reaction_root event_root plot_args

  echo "============================================================"
  echo "Running ${label}"
  echo "Macro   : ${macro}"
  echo "Channel : ${channel}"
  echo "Threads : ${THREADS}"
  echo "Si ring : ${RING}"
  echo "============================================================"

  before="$(mktemp)"
  after="$(mktemp)"
  new_files="$(mktemp)"
  snapshot_root_files > "${before}"

  (
    cd "${PROJECT_ROOT}"
    ./batch.sh "${THREADS}" "${macro}"
  )

  snapshot_root_files > "${after}"
  comm -13 "${before}" "${after}" > "${new_files}"

  echo "New ROOT files for ${label}:"
  sed 's/^/  /' "${new_files}"

  reaction_root="$(find_newest_matching "${new_files}" reaction || true)"
  event_root="$(find_newest_matching "${new_files}" event || true)"

  if [[ -z "${reaction_root}" ]]; then
    echo "ERROR: no new reaction ROOT file detected for ${label}." >&2
    exit 1
  fi
  if [[ -z "${event_root}" ]]; then
    echo "ERROR: no new event ROOT file detected for ${label}. Check /output/saveEvent true." >&2
    exit 1
  fi

  echo "Selected reaction ROOT: ${reaction_root}"
  echo "Selected event ROOT   : ${event_root}"

  plot_args=(
    python3 "${PLOT_SCRIPT}"
    --reaction-root "${reaction_root}"
    --event-root "${event_root}"
    --channel "${channel}"
    --label "${label}"
    --ring "${RING}"
    --outdir "${SCRIPT_DIR}"
  )
  if [[ "${WEIGHTED}" == "1" ]]; then
    plot_args+=(--weighted)
  fi

  "${plot_args[@]}"

  rm -f "${before}" "${after}" "${new_files}"
}

run_case "si_alpha_165_model_angular" "validation_si_alpha/si_alpha_165_model_angular.mac" 0
run_case "si_alpha_675_model_angular" "validation_si_alpha/si_alpha_675_model_angular.mac" 1

cat <<EOF2
============================================================
Si alpha-response validation complete.
Output files are written directly in:
  ${SCRIPT_DIR}

Most important files:
  *_summary.txt
  *_si_total.png
  *_si_by_branch_counts.png
  *_si_by_component_counts.png
  *_si_sector_energy_map.png
EOF2
