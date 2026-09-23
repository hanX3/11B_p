#!/usr/bin/env bash
set -euo pipefail

events=5000
strict_events=5000
jobs=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --events)
      events="$2"
      shift 2
      ;;
    --strict-events)
      strict_events="$2"
      shift 2
      ;;
    --jobs)
      jobs="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: bash validation/angular_decay/run_angular_decay_validation.sh [--events N] [--strict-events N] [--jobs J]"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"
plots_dir="$script_dir/plots"
output_dir="$script_dir/output"
work_dir="$script_dir/work"
report="$script_dir/angular_decay_validation_report.md"
summary_json="$script_dir/angular_decay_validation_summary.json"

mkdir -p "$plots_dir" "$output_dir" "$work_dir"
export MPLCONFIGDIR="$script_dir/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

if [[ ! -x "$repo_root/build/HB" ]]; then
  echo "Missing executable: $repo_root/build/HB" >&2
  echo "Build first, for example: cmake --build build" >&2
  exit 1
fi

latest_reaction_root() {
  find "$repo_root/data" -maxdepth 1 -type f -name 'reaction_merged_*.root' -printf '%T@ %p\n' 2>/dev/null \
    | sort -n \
    | tail -n 1 \
    | cut -d' ' -f2-
}

write_temp_macro() {
  local src="$1"
  local dst="$2"
  local n_events="$3"
  local injection="$4"
  awk -v events="$n_events" -v injection="$injection" '
    /^\/run\/beamOn[[:space:]]+/ {
      if (injection != "") print injection
      print "/run/beamOn " events
      next
    }
    { print }
  ' "$src" > "$dst"
}

run_macro() {
  local label="$1"
  local macro="$2"
  local n_events="$3"
  local injection="$4"

  local before after tmp link_name
  before="$(latest_reaction_root || true)"
  tmp="$(mktemp "/tmp/angular_${label}.XXXXXX.mac")"
  write_temp_macro "$repo_root/$macro" "$tmp" "$n_events" "$injection"
  echo "Running $label from $macro with $n_events events" >&2
  (cd "$repo_root" && bash ./batch.sh "$jobs" "$tmp" > "$work_dir/${label}_geant4.log" 2>&1)
  after="$(latest_reaction_root || true)"
  if [[ -z "$after" || "$after" == "$before" ]]; then
    echo "Could not identify new merged reaction ROOT for $label" >&2
    exit 1
  fi
  link_name="$output_dir/${label}_$(basename "$after")"
  ln -sfn "$after" "$link_name"
  printf '%s\n' "$after"
}

alpha_injection_common=$'/h11b/enableDirectDecay false\n/h11b/165SequentialDecayFraction 1.0\n/h11b/675SequentialDecayFraction 1.0\n/h11b/enable165GammaCapture false\n/h11b/enable675GammaCapture false'
legacy_675_injection=$'/h11b/enableDirectDecay false\n/h11b/165SequentialDecayFraction 1.0\n/h11b/675SequentialDecayFraction 1.0\n/h11b/enable165GammaCapture false\n/h11b/enable675GammaCapture false\n/h11b/675AlphaDecayModel legacyLegendreA2A4'

declare -A roots

roots[165_primary_off]="$(run_macro 165_primary_off macros/validation_angle_165_primary_off.mac "$events" "$alpha_injection_common")"
roots[165_primary_on]="$(run_macro 165_primary_on macros/validation_angle_165_primary_on.mac "$events" "$alpha_injection_common")"
roots[165_secondary_off]="$(run_macro 165_secondary_off macros/validation_angle_165_secondary_off.mac "$events" "$alpha_injection_common")"
roots[165_secondary_on]="$(run_macro 165_secondary_on macros/validation_angle_165_secondary_on.mac "$events" "$alpha_injection_common")"
roots[675_primary_off]="$(run_macro 675_primary_off macros/validation_angle_675_primary_off.mac "$events" "$legacy_675_injection")"
roots[675_primary_on]="$(run_macro 675_primary_on macros/validation_angle_675_primary_on.mac "$events" "$legacy_675_injection")"
roots[675_secondary_off]="$(run_macro 675_secondary_off macros/validation_angle_675_secondary_off.mac "$events" "$legacy_675_injection")"
roots[675_secondary_on]="$(run_macro 675_secondary_on macros/validation_angle_675_secondary_on.mac "$events" "$legacy_675_injection")"

strict_labels=(default l1only l3only incoherent no_sym phase0)
strict_macros=(
  macros/validation_675alpha_strict_default.mac
  macros/validation_675alpha_strict_l1only.mac
  macros/validation_675alpha_strict_l3only.mac
  macros/validation_675alpha_strict_incoherent.mac
  macros/validation_675alpha_strict_no_sym.mac
  macros/validation_675alpha_strict_phase0.mac
)

for i in "${!strict_labels[@]}"; do
  label="strict_${strict_labels[$i]}"
  roots[$label]="$(run_macro "$label" "${strict_macros[$i]}" "$strict_events" "")"
done

legendre_jsons=()
run_legendre() {
  local name="$1"
  local off="$2"
  local on="$3"
  local observable="$4"
  local kind="$5"
  local select_args="$6"
  local mode_branch="$7"
  local on_mode="$8"
  local coeff_args="$9"
  local out="$plots_dir/legendre_${name}_on_off.png"
  local json="$work_dir/legendre_${name}.json"

  # shellcheck disable=SC2086
  set +e
  python3 "$script_dir/validate_legendre_distributions.py" \
    --off "$off" --on "$on" --observable "$observable" --kind "$kind" \
    --test-name "$name Legendre on/off validation" \
    $select_args $coeff_args \
    --mode-branch "$mode_branch" --expected-on-mode "$on_mode" \
    --out "$out" --json "$json" > "$work_dir/legendre_${name}.txt"
  set -e
  legendre_jsons+=("$json")
}

run_legendre 165_primary "${roots[165_primary_off]}" "${roots[165_primary_on]}" cos_theta_primary_cm a1a2 "--select reaction_channel=0" primary_angular_mode 1 "--a1-branch primary_a1_used --a2-branch primary_a2_used"
run_legendre 165_secondary "${roots[165_secondary_off]}" "${roots[165_secondary_on]}" cos_theta_secondary_correlation a2a4 "--select reaction_channel=0 --select branch_id=1" secondary_angular_model 1 "--a2-branch secondary_a2_used --a4-branch secondary_a4_used"
run_legendre 675_primary "${roots[675_primary_off]}" "${roots[675_primary_on]}" cos_theta_primary_cm a1a2 "--select reaction_channel=1" primary_angular_mode 1 "--a1-branch primary_a1_used --a2-branch primary_a2_used"
run_legendre 675_secondary "${roots[675_secondary_off]}" "${roots[675_secondary_on]}" cos_chi_secondary_8be a2a4 "--select reaction_channel=1 --select branch_id=1" secondary_angular_model 1 "--a2-branch secondary_a2_used --a4-branch secondary_a4_used"

strict_input_args=()
for label in "${strict_labels[@]}"; do
  strict_input_args+=("$label:${roots[strict_$label]}")
done
strict_json="$work_dir/strict_summary.json"
set +e
python3 "$script_dir/validate_675_strict_model.py" \
  --inputs "${strict_input_args[@]}" \
  --outdir "$plots_dir" \
  --json "$strict_json" \
  --max-attempts 10000 > "$work_dir/strict_summary.txt"
set -e

python3 - "$summary_json" "$report" "$events" "$strict_events" "$jobs" "${legendre_jsons[@]}" "$strict_json" <<'PY'
import json
import subprocess
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
report_path = Path(sys.argv[2])
events = sys.argv[3]
strict_events = sys.argv[4]
jobs = sys.argv[5]
json_paths = [Path(p) for p in sys.argv[6:]]

items = []
for path in json_paths:
    with open(path) as f:
        data = json.load(f)
    if "tests" in data:
        items.extend(data["tests"])
        strict = data
    else:
        items.append(data)

rank = {"PASS": 0, "WARN": 1, "FAIL": 2}
overall = "PASS"
for item in items:
    if rank[item["status"]] > rank[overall]:
        overall = item["status"]

summary = {
    "overall_status": overall,
    "events_per_legendre_macro": int(events),
    "events_per_strict_macro": int(strict_events),
    "jobs": int(jobs),
    "tests": items,
}
with open(summary_path, "w") as f:
    json.dump(summary, f, indent=2, sort_keys=True)
    f.write("\n")

commit = subprocess.check_output(["git", "log", "--oneline", "-1"], text=True).strip()
lines = [
    "# Angular Distribution and Decay Model Validation Report",
    "",
    "## Build / Run Information",
    f"- Git commit: `{commit}`",
    f"- Legendre macro events: `{events}`",
    f"- Strict macro events: `{strict_events}`",
    f"- Jobs: `{jobs}`",
    "",
    "## Legendre Angular Distribution Validation",
]
for item in items[:4]:
    lines.extend([
        f"### {item['test_name']}",
        f"- Status: `{item['status']}`",
        f"- Observable: `{item.get('observable', '')}`",
        f"- n_off / n_on: `{item.get('n_off', '')}` / `{item.get('n_on', '')}`",
        f"- coefficients: `{item.get('coefficients', '')}`",
        f"- off chi2/ndf: `{item.get('off_chi2_ndf', float('nan')):.4g}`",
        f"- on chi2/ndf: `{item.get('on_chi2_ndf', float('nan')):.4g}`",
        f"- on/off KS: `{item.get('on_off_ks_statistic', float('nan')):.4g}`",
        f"- plot: `{item.get('plot', '')}`",
        "",
    ])

lines.append("## 675 Strict Three-Alpha Model Validation")
for item in items[4:]:
    lines.append(f"- [{item['status']}] {item['test_name']}: " + ", ".join(f"{k}={v}" for k, v in item.items() if k not in ("test_name", "status")))

lines.extend(["", "## Final Summary"])
for item in items:
    lines.append(f"[{item['status']}] {item['test_name']}")
lines.append(f"\nOverall status: `{overall}`")
report_path.write_text("\n".join(lines) + "\n")
print(f"Wrote {summary_path}")
print(f"Wrote {report_path}")
sys.exit(0 if overall != "FAIL" else 1)
PY

echo "Wrote $summary_json"
echo "Wrote $report"
