#!/usr/bin/env python3
"""Validate event_weight / event_sampling_weight against bias factors.

See gpt/codex_cs_channel_validation.md section 5.4.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from _common import Report, load_branches, open_tree, max_abs_diff, DEFAULT_TOL

BRANCHES = [
    "reaction_channel",
    "event_weight",
    "event_sampling_weight",
    "cross_section_bias_factor",
    "gamma_bias_factor",
]


def validate_file(path: Path) -> Report:
    report = Report(path)
    tree = open_tree(path)
    data, missing = load_branches(tree, BRANCHES)
    if missing:
        report.emit("FAIL", "required branches present (see MISSING BRANCH lines above)")
        return report
    if tree.num_entries == 0:
        report.emit("WARN", "no reaction entries; event weight checks skipped")
        return report

    channels = data["reaction_channel"]
    bias = data["cross_section_bias_factor"]
    gamma_bias = data["gamma_bias_factor"]

    three_alpha_mask = np.isin(channels, [0, 1, 2])
    if np.any(three_alpha_mask):
        expected = 1.0 / bias[three_alpha_mask]
        diff = max_abs_diff(data["event_weight"][three_alpha_mask], expected)
        if np.isfinite(diff) and diff < DEFAULT_TOL:
            report.emit("PASS", f"3-alpha event weight consistency (max abs diff = {diff:.3g})")
        else:
            report.emit("FAIL", f"3-alpha event weight consistency (max abs diff = {diff:.3g})")
            bad_idx = int(np.nanargmax(np.abs(data["event_weight"][three_alpha_mask] - expected)))
            report.fail_detail(
                "event_weight", float(expected[bad_idx]), float(data["event_weight"][three_alpha_mask][bad_idx]),
                entry_index=bad_idx, related_setting="/h11b/crossSectionBiasFactor",
            )
    else:
        report.emit("WARN", "no 3-alpha (reaction_channel in {0,1,2}) events in this file")

    gamma_mask = channels == 3
    if np.any(gamma_mask):
        expected_gamma = 1.0 / (bias[gamma_mask] * gamma_bias[gamma_mask])
        diff_gamma = max_abs_diff(data["event_weight"][gamma_mask], expected_gamma)
        if np.isfinite(diff_gamma) and diff_gamma < DEFAULT_TOL:
            report.emit("PASS", f"gamma event weight consistency (max abs diff = {diff_gamma:.3g})")
        else:
            report.emit("FAIL", f"gamma event weight consistency (max abs diff = {diff_gamma:.3g})")
            bad_idx = int(np.nanargmax(np.abs(data["event_weight"][gamma_mask] - expected_gamma)))
            report.fail_detail(
                "event_weight", float(expected_gamma[bad_idx]), float(data["event_weight"][gamma_mask][bad_idx]),
                entry_index=bad_idx, related_setting="/h11b/crossSectionBiasFactor, /h11b/gammaBiasFactor",
            )
    else:
        report.emit("WARN", "no gamma-capture (reaction_channel=3) events in this file")

    diff_sampling = max_abs_diff(data["event_sampling_weight"], data["event_weight"])
    if np.isfinite(diff_sampling) and diff_sampling < DEFAULT_TOL:
        report.emit("PASS", "event_sampling_weight consistency")
        print("    NOTE: event_sampling_weight is identical to event_weight in current implementation")
    else:
        report.emit("FAIL", "event_sampling_weight consistency")
        report.fail_detail("event_sampling_weight", "== event_weight", f"max abs diff = {diff_sampling:.3g}")

    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_files", nargs="+", type=Path)
    args = parser.parse_args()

    overall = 0
    for path in args.root_files:
        if not path.exists():
            print(f"{path}: file not found")
            overall = max(overall, 1)
            continue
        report = validate_file(path)
        print(report.summary_line())
        overall = max(overall, {"PASS": 0, "WARN": 0, "FAIL": 1}[report.status])
    return overall


if __name__ == "__main__":
    sys.exit(main())
