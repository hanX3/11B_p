#!/usr/bin/env bash
set -euo pipefail

# Run all angular-distribution validation jobs.
# Usage from project root:
#   ./validation_angle_distribution/run_all_angle_validation.sh

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

DATA_DIR="${PROJECT_ROOT}/data"
mkdir -p "${DATA_DIR}"
THREADS=8

# Store only the merged ROOT files in data/. PNG files are written directly in validation_angle_distribution/.
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
  echo "Macro: validation_angle_distribution/${macro}"
  echo "Threads: ${THREADS}"
  echo "============================================================"

  touch "${marker}"
  ./batch.sh "${THREADS}" "validation_angle_distribution/${macro}"

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

plot_pair() {
  local label="$1"
  local off_key="$2"
  local on_key="$3"
  local channel="$4"
  local branch="$5"
  local output_png="${THIS_DIR}/${6}"
  local gamma_resonance="${7:-}"
  local gamma_branch="${8:-}"

  echo "============================================================"
  echo "Plotting ${label}"
  echo "============================================================"

  local cmd=(
    python3 "${THIS_DIR}/plot_compare_angular_cos.py"
    --off "${ROOTS[${off_key}]}"
    --on "${ROOTS[${on_key}]}"
    --channel "${channel}"
    --branch "${branch}"
    --output "${output_png}"
  )

  if [[ -n "${gamma_resonance}" ]]; then
    cmd+=(--gamma-resonance "${gamma_resonance}")
  fi
  if [[ -n "${gamma_branch}" ]]; then
    cmd+=(--gamma-branch "${gamma_branch}")
  fi

  "${cmd[@]}"
}

run_one 165_primary_off   validation_angle_165_primary_off.mac
run_one 165_primary_on    validation_angle_165_primary_on.mac
plot_pair 165_primary 165_primary_off 165_primary_on 0 cos_theta_primary_cm compare_165_primary_cos.png

run_one 675_primary_off   validation_angle_675_primary_off.mac
run_one 675_primary_on    validation_angle_675_primary_on.mac
plot_pair 675_primary 675_primary_off 675_primary_on 1 cos_theta_primary_cm compare_675_primary_cos.png

run_one 165_secondary_off validation_angle_165_secondary_off.mac
run_one 165_secondary_on  validation_angle_165_secondary_on.mac
plot_pair 165_secondary 165_secondary_off 165_secondary_on 0 cos_theta_secondary_correlation compare_165_secondary_cos.png

run_one 675_secondary_off validation_angle_675_secondary_off.mac
run_one 675_secondary_on  validation_angle_675_secondary_on.mac
plot_pair 675_secondary 675_secondary_off 675_secondary_on 1 cos_theta_secondary_correlation compare_675_secondary_cos.png

run_one 165_gamma_off     validation_angle_165_gamma_off.mac
run_one 165_gamma_on      validation_angle_165_gamma_on.mac
plot_pair 165_gamma 165_gamma_off 165_gamma_on 3 cos_theta_gamma_cm compare_165_gamma_cos.png 165 1

run_one 675_gamma_off     validation_angle_675_gamma_off.mac
run_one 675_gamma_on      validation_angle_675_gamma_on.mac
plot_pair 675_gamma 675_gamma_off 675_gamma_on 3 cos_theta_gamma_cm compare_675_gamma_cos.png 675

echo "============================================================"
echo "All angular-distribution validation runs finished."
echo "ROOT files remain in: ${DATA_DIR}"
echo "PNG files are in: validation_angle_distribution/"
echo "Generated PNGs:"
ls -1 "${THIS_DIR}"/*.png 2>/dev/null || true
echo "============================================================"
