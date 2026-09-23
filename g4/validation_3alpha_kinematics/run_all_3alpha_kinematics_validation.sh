#!/usr/bin/env bash
set -euo pipefail

# Run generator-level three-alpha kinematics validation.
# Usage from project root:
#   ./validation_3alpha_kinematics/run_all_3alpha_kinematics_validation.sh

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${THIS_DIR}/.." && pwd)"

case "${PROJECT_ROOT}" in
  /tmp|/tmp/*)
    echo "ERROR: refusing to run from /tmp. Run inside the project directory."
    exit 1
    ;;
esac

cd "${PROJECT_ROOT}"

if [[ ! -x "./batch.sh" ]]; then
  echo "ERROR: ./batch.sh not found or not executable."
  echo "Run: chmod +x batch.sh"
  exit 1
fi

if ! python3 - <<'PY' >/dev/null 2>&1
import uproot, numpy, matplotlib
PY
then
  echo "ERROR: python3 packages are missing. Required: uproot numpy matplotlib"
  exit 1
fi

DATA_DIR="${PROJECT_ROOT}/data"
mkdir -p "${DATA_DIR}"
THREADS="${THREADS:-8}"

declare -A ROOTS

newest_merged_after_marker() {
  local marker="$1"
  local merged
  merged="$(find "${DATA_DIR}" -maxdepth 1 -type f -name 'reaction_merged_*.root' -newer "${marker}" -printf '%T@ %p\n' \
    | sort -n \
    | tail -n 1 \
    | cut -d' ' -f2-)"
  echo "${merged}"
}

run_one() {
  local key="$1"
  local macro="$2"
  local marker="${THIS_DIR}/.marker_${key}"

  echo "============================================================"
  echo "Running ${key}"
  echo "Macro: validation_3alpha_kinematics/${macro}"
  echo "Threads: ${THREADS}"
  echo "============================================================"

  touch "${marker}"
  ./batch.sh "${THREADS}" "validation_3alpha_kinematics/${macro}"

  local merged
  merged="$(newest_merged_after_marker "${marker}")"
  rm -f "${marker}"

  if [[ -z "${merged}" ]]; then
    echo "ERROR: no new merged reaction ROOT file found in ${DATA_DIR} after ${key}."
    exit 1
  fi

  ROOTS["${key}"]="${merged}"
  echo "ROOT for ${key}: ${merged}"
}

run_one 165_model validation_kinematics_165_model.mac
run_one 675_model validation_kinematics_675_model.mac
run_one background_phase_space validation_kinematics_background_phase_space.mac

echo "============================================================"
echo "Plotting and summarizing three-alpha kinematics"
echo "============================================================"

python3 "${THIS_DIR}/plot_3alpha_kinematics.py" \
  --sample "165_model:${ROOTS[165_model]}:0" \
  --sample "675_model:${ROOTS[675_model]}:1" \
  --sample "background_phase_space:${ROOTS[background_phase_space]}:2" \
  --output-dir "${THIS_DIR}"

echo "============================================================"
echo "Three-alpha kinematics validation runs finished."
echo "ROOT files remain in: ${DATA_DIR}"
echo "PNG/TXT outputs are in: validation_3alpha_kinematics/"
echo "Generated outputs:"
ls -1 "${THIS_DIR}"/*.png "${THIS_DIR}"/*.txt 2>/dev/null || true
echo "============================================================"
