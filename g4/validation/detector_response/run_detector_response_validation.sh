#!/usr/bin/env bash
set -euo pipefail

step_events=1000
stats_events=20000
jobs=4
skip_build=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step-events)
      step_events="$2"
      shift 2
      ;;
    --stats-events)
      stats_events="$2"
      shift 2
      ;;
    --jobs)
      jobs="$2"
      shift 2
      ;;
    --skip-build)
      skip_build=1
      shift
      ;;
    -h|--help)
      echo "Usage: bash validation/detector_response/run_detector_response_validation.sh [--step-events N] [--stats-events N] [--jobs J] [--skip-build]"
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
logs_dir="$script_dir/logs"
results_dir="$script_dir/results"
plots_dir="$script_dir/plots"
work_dir="$script_dir/work"
report="$script_dir/detector_response_validation_report.md"
summary_json="$results_dir/detector_response_summary.json"
summary_md="$results_dir/detector_response_summary.md"
checks_csv="$results_dir/detector_response_checks.csv"

mkdir -p "$logs_dir" "$results_dir" "$plots_dir" "$work_dir"
export MPLCONFIGDIR="$script_dir/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

if [[ "$skip_build" -eq 0 ]]; then
  echo "Building project" >&2
  cmake --build "$repo_root/build" > "$logs_dir/build.log" 2>&1
fi

if [[ ! -x "$repo_root/build/HB" ]]; then
  echo "Missing executable: $repo_root/build/HB" >&2
  echo "Build first, for example: cmake --build build" >&2
  exit 1
fi

snapshot_roots() {
  find "$repo_root/data" -maxdepth 1 -type f -name '*.root' -print 2>/dev/null | sort -u
}

write_temp_macro() {
  local src="$1"
  local dst="$2"
  local n_events="$3"
  awk -v events="$n_events" '
    /^\/run\/beamOn[[:space:]]+/ {
      print "/run/beamOn " events
      next
    }
    { print }
  ' "$src" > "$dst"
}

newest_new_class_file() {
  local before_file="$1"
  local after_file="$2"
  local class_name="$3"
  comm -13 "$before_file" "$after_file" \
    | awk -v class_name="$class_name" '
        {
          n=split($0, parts, "/");
          base=parts[n];
          if (base ~ "^" class_name "_merged_.*\\.root$" || base ~ "^" class_name "_.*_t[0-9]+\\.root$") {
            print $0;
          }
        }' \
    | while IFS= read -r path; do
        printf '%s %s\n' "$(stat -c '%Y' "$path")" "$path"
      done \
    | sort -n \
    | tail -n 1 \
    | cut -d' ' -f2-
}

run_macro() {
  local label="$1"
  local macro="$2"
  local n_events="$3"
  local tmp before after

  tmp="$(mktemp "/tmp/detector_${label}.XXXXXX.mac")"
  before="$work_dir/${label}_before.txt"
  after="$work_dir/${label}_after.txt"
  snapshot_roots > "$before"
  write_temp_macro "$macro" "$tmp" "$n_events"
  echo "Running $label with $n_events events" >&2
  (cd "$repo_root" && bash ./batch.sh "$jobs" "$tmp" > "$logs_dir/${label}.log" 2>&1)
  snapshot_roots > "$after"
  {
    printf 'reaction=%s\n' "$(newest_new_class_file "$before" "$after" reaction || true)"
    printf 'event=%s\n' "$(newest_new_class_file "$before" "$after" event || true)"
    printf 'track=%s\n' "$(newest_new_class_file "$before" "$after" track || true)"
    printf 'step=%s\n' "$(newest_new_class_file "$before" "$after" step || true)"
  } > "$results_dir/${label}_root_files.env"
}

validate_one() {
  local label="$1"
  local scenario="$2"
  local expect_track="$3"
  local expect_step="$4"
  # shellcheck disable=SC1090
  source "$results_dir/${label}_root_files.env"
  if [[ -z "${reaction:-}" || -z "${event:-}" ]]; then
    echo "$label did not produce required reaction/event ROOT files" >&2
    exit 1
  fi
  local args=(
    "$script_dir/validate_detector_response.py"
    --scenario "$scenario"
    --reaction "$reaction"
    --event "$event"
    --json "$results_dir/${label}_summary.json"
    --summary-md "$results_dir/${label}_summary.md"
    --checks-csv "$results_dir/${label}_checks.csv"
  )
  if [[ -n "${track:-}" ]]; then
    args+=(--track "$track")
  fi
  if [[ -n "${step:-}" ]]; then
    args+=(--step "$step")
  fi
  if [[ "$expect_track" == "1" ]]; then
    args+=(--expect-track)
  fi
  if [[ "$expect_step" == "1" ]]; then
    args+=(--expect-step)
  fi

  set +e
  python3 "${args[@]}" > "$logs_dir/${label}_validate.log" 2>&1
  local status=$?
  set -e

  local plot_args=(
    "$script_dir/plot_detector_response.py"
    --label "$label"
    --reaction "$reaction"
    --event "$event"
    --outdir "$plots_dir"
  )
  if [[ -n "${track:-}" ]]; then
    plot_args+=(--track "$track")
  fi
  if [[ -n "${step:-}" ]]; then
    plot_args+=(--step "$step")
  fi
  python3 "${plot_args[@]}" > "$logs_dir/${label}_plots.log" 2>&1 || true
  return "$status"
}

run_macro alpha_step "$script_dir/macros/validation_alpha_step_debug.mac" "$step_events"
run_macro gamma_step "$script_dir/macros/validation_gamma_step_debug.mac" "$step_events"
run_macro event_stats "$script_dir/macros/validation_event_response_stats.mac" "$stats_events"

validation_failed=0
validate_one alpha_step alpha_step 1 1 || validation_failed=1
validate_one gamma_step gamma_step 1 1 || validation_failed=1
validate_one event_stats event_stats 0 0 || validation_failed=1

python3 - "$summary_json" "$summary_md" "$checks_csv" "$report" "$step_events" "$stats_events" "$jobs" "$results_dir"/alpha_step_summary.json "$results_dir"/gamma_step_summary.json "$results_dir"/event_stats_summary.json <<'PY'
import csv
import json
import subprocess
import sys
from pathlib import Path

summary_json = Path(sys.argv[1])
summary_md = Path(sys.argv[2])
checks_csv = Path(sys.argv[3])
report = Path(sys.argv[4])
step_events = int(sys.argv[5])
stats_events = int(sys.argv[6])
jobs = int(sys.argv[7])
inputs = [Path(p) for p in sys.argv[8:]]

rank = {"PASS": 0, "WARN": 1, "FAIL": 2}
scenarios = []
overall = "PASS"
tests = []
for path in inputs:
    data = json.loads(path.read_text())
    scenarios.append(data)
    if rank[data["overall_status"]] > rank[overall]:
        overall = data["overall_status"]
    for item in data["tests"]:
        row = dict(item)
        row["scenario"] = data["scenario"]
        tests.append(row)

try:
    commit = subprocess.check_output(["git", "log", "--oneline", "-1"], text=True).strip()
except Exception:
    commit = "unknown"

combined = {
    "overall_status": overall,
    "git_commit": commit,
    "step_events": step_events,
    "stats_events": stats_events,
    "jobs": jobs,
    "scenarios": scenarios,
    "tests": tests,
}
summary_json.write_text(json.dumps(combined, indent=2, sort_keys=True) + "\n")

with checks_csv.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["scenario", "status", "name", "detail", "affected_entries", "related_branch"])
    writer.writeheader()
    for item in tests:
        writer.writerow({
            "scenario": item.get("scenario", ""),
            "status": item.get("status", ""),
            "name": item.get("name", ""),
            "detail": item.get("detail", ""),
            "affected_entries": item.get("affected_entries", ""),
            "related_branch": item.get("related_branch", ""),
        })

lines = [
    "# Detector response validation report",
    "",
    "## 1. Code version and run configuration",
    f"- Git commit: `{commit}`",
    f"- Step-debug events per macro: `{step_events}`",
    f"- Event-stat events: `{stats_events}`",
    f"- Jobs: `{jobs}`",
    "- Macros:",
    "  - `validation/detector_response/macros/validation_alpha_step_debug.mac`",
    "  - `validation/detector_response/macros/validation_gamma_step_debug.mac`",
    "  - `validation/detector_response/macros/validation_event_response_stats.mac`",
    "- Optional no-reaction macro: not included; no runtime command was found that cleanly disables only the custom H11B reaction process for an otherwise identical transport run.",
    "",
    "## 2. ROOT tree and branch inventory",
]
for scenario in scenarios:
    lines.append(f"### {scenario['scenario']}")
    for label, info in scenario["files"].items():
        if not info.get("path"):
            lines.append(f"- {label}: not produced/requested")
            continue
        branches = ", ".join(info.get("branches", []))
        lines.append(f"- {label}: `{info['path']}` entries=`{info['entries']}` branches: {branches}")

lines.extend(["", "## 3. Event / track / step consistency"])
for scenario in scenarios:
    lines.append(f"### {scenario['scenario']}")
    for item in scenario["tests"]:
        if any(key in item["name"] for key in ["event ids", "step-track", "track id", "track event", "track length", "track global", "parent", "step event", "step track"]):
            lines.append(f"- [{item['status']}] {item['name']}: {item.get('detail', '')}")

lines.extend(["", "## 4. Energy deposition checks"])
for scenario in scenarios:
    lines.append(f"### {scenario['scenario']}")
    for item in scenario["tests"]:
        if any(key in item["name"] for key in ["energy", "edep", "hit records", "numeric finite"]):
            lines.append(f"- [{item['status']}] {item['name']}: {item.get('detail', '')}")

lines.extend(["", "## 5. Detector hit checks"])
for scenario in scenarios:
    lines.append(f"### {scenario['scenario']}")
    for item in scenario["tests"]:
        if any(key in item["name"] for key in ["detector", "hit", "alpha", "gamma"]):
            lines.append(f"- [{item['status']}] {item['name']}: {item.get('detail', '')}")
    det_summary = scenario.get("detector_summary", {})
    if det_summary:
        lines.append("")
        lines.append("| detector_type | name | hits | events | total_edep | copy_min | copy_max |")
        lines.append("| --- | --- | ---: | ---: | ---: | ---: | ---: |")
        for key, row in sorted(det_summary.items(), key=lambda kv: int(kv[0])):
            lines.append(f"| {key} | {row['name']} | {row['hits']} | {row['events']} | {row['total_edep']:.8g} | {row['copy_min']} | {row['copy_max']} |")

lines.extend(["", "## 6. Transport path checks"])
for scenario in scenarios:
    lines.append(f"### {scenario['scenario']}")
    ss = scenario.get("string_summary", {})
    for key in ["track_particles", "step_particles", "step_volumes", "step_processes"]:
        if key in ss:
            top = ", ".join(f"{name}:{count}" for name, count in list(ss[key].items())[:12])
            lines.append(f"- {key}: {top}")

lines.extend(["", "## 7. Plots"])
for path in sorted(Path("validation/detector_response/plots").glob("*.png")):
    lines.append(f"- `{path}`")

lines.extend(["", "## 8. PASS / WARN / FAIL summary"])
for scenario in scenarios:
    lines.append(f"- {scenario['scenario']}: `{scenario['overall_status']}`")
lines.append(f"- Detector response validation: `{overall}`")

lines.extend(["", "## 9. Issues requiring source-code inspection"])
issues = [item for item in tests if item.get("status") in {"WARN", "FAIL"}]
if not issues:
    lines.append("- None.")
else:
    for item in issues:
        lines.append(f"- [{item['status']}] {item['scenario']} / {item['name']}: {item.get('detail', '')}")
lines.extend([
    "",
    "## Conclusion",
    f"Detector response validation: {overall}",
])
if overall == "PASS":
    lines.append("The detector-response output is internally consistent for the tested configurations.")
elif overall == "WARN":
    lines.append("The detector-response output is usable for follow-up analysis with the documented limitations.")
else:
    lines.append("At least one detector-response check failed; inspect the FAIL entries before using these outputs for physics analysis.")

text = "\n".join(lines) + "\n"
summary_md.write_text(text)
report.write_text(text)
PY

echo "Detector response validation summary: $summary_json"
echo "Detector response validation report: $report"
exit "$validation_failed"
