#!/usr/bin/env bash
set -euo pipefail

THREADS=8

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
BATCH_SH="${PROJECT_ROOT}/batch.sh"
PLOT_SCRIPT="${SCRIPT_DIR}/plot_gamma_response.py"

if [[ ! -x "${BATCH_SH}" ]]; then
  echo "ERROR: batch.sh not found or not executable: ${BATCH_SH}" >&2
  exit 1
fi

if [[ ! -f "${PLOT_SCRIPT}" ]]; then
  echo "ERROR: plot script not found: ${PLOT_SCRIPT}" >&2
  exit 1
fi

mkdir -p "${DATA_DIR}"

snapshot_files() {
  local pattern="$1"
  find "${DATA_DIR}" -maxdepth 1 -type f -name "${pattern}" -print | sort
}

pick_new_file() {
  local before_list="$1"
  local after_list="$2"
  local kind="$3"
  local new_files
  new_files="$(comm -13 \
    <(printf "%s\n" "${before_list}" | sed '/^$/d') \
    <(printf "%s\n" "${after_list}" | sed '/^$/d'))"

  if [[ -z "${new_files}" ]]; then
    echo "ERROR: no new ${kind} ROOT file detected." >&2
    exit 1
  fi

  printf "%s\n" "${new_files}" | tail -n 1
}

run_case() {
  local resonance="$1"
  local macro="${SCRIPT_DIR}/gamma_${resonance}_response.mac"
  local prefix="${SCRIPT_DIR}/gamma_${resonance}"

  if [[ ! -f "${macro}" ]]; then
    echo "ERROR: macro not found: ${macro}" >&2
    exit 1
  fi

  echo "============================================================"
  echo "Gamma-response validation: ${resonance} keV"
  echo "Macro   : ${macro}"
  echo "Threads : ${THREADS}"
  echo "Output  : ${prefix}_*.png / ${prefix}_summary.txt"
  echo "============================================================"

  local before_events before_reactions after_events after_reactions event_root reaction_root
  before_events="$(snapshot_files 'event_merged_*.root')"
  before_reactions="$(snapshot_files 'reaction_merged_*.root')"

  cd "${PROJECT_ROOT}"
  "${BATCH_SH}" "${THREADS}" "${macro}"

  after_events="$(snapshot_files 'event_merged_*.root')"
  after_reactions="$(snapshot_files 'reaction_merged_*.root')"

  event_root="$(pick_new_file "${before_events}" "${after_events}" event)"
  reaction_root="$(pick_new_file "${before_reactions}" "${after_reactions}" reaction)"

  echo "Event ROOT    : ${event_root}"
  echo "Reaction ROOT : ${reaction_root}"

  python3 "${PLOT_SCRIPT}" \
    --event-root "${event_root}" \
    --reaction-root "${reaction_root}" \
    --resonance "${resonance}" \
    --out-prefix "${prefix}"
}

run_case 165
sleep 1
run_case 675

echo "============================================================"
echo "Gamma-response validation complete."
echo "Results are in: ${SCRIPT_DIR}"
echo "============================================================"
