#!/usr/bin/env python3
"""Validate cross-section component algebra against H11BCrossSection::CalculateComponents().

See gpt/codex_cs_channel_validation.md section 3 and 5.1 for the formulas
checked here.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from _common import Report, load_branches, open_tree, DEFAULT_TOL

BRANCHES = [
    "sigma_165_model_b",
    "sigma_675_model_b",
    "sigma_165_total_b",
    "sigma_675_total_b",
    "sigma_165_used_b",
    "sigma_675_used_b",
    "sigma_165_directdecay_b",
    "sigma_675_directdecay_b",
    "sigma_directdecay_b",
    "sigma_total_used_b",
    "sigma_165_sampling_b",
    "sigma_675_sampling_b",
    "sigma_165_directdecay_sampling_b",
    "sigma_675_directdecay_sampling_b",
    "sigma_directdecay_sampling_b",
    "sigma_3alpha_sampling_total_b",
    "sigma_gamma_165_0_physical_b",
    "sigma_gamma_165_1_physical_b",
    "sigma_gamma_675_physical_b",
    "sigma_gamma_165_0_sampling_b",
    "sigma_gamma_165_1_sampling_b",
    "sigma_gamma_675_sampling_b",
    "sigma_gamma_sampling_total_b",
    "sigma_total_physical_all_b",
    "sigma_total_sampling_all_b",
    "cross_section_bias_factor",
    "gamma_bias_factor",
    "scale_factor_165",
    "scale_factor_675",
    "enable_direct_decay",
    "sequential_decay_fraction_165",
    "sequential_decay_fraction_675",
    "direct_decay_fraction_165",
    "direct_decay_fraction_675",
]


def validate_file(path: Path) -> Report:
    report = Report(path)
    tree = open_tree(path)
    data, missing = load_branches(tree, BRANCHES)
    if missing:
        report.emit("FAIL", "required branches present (see MISSING BRANCH lines above)")
        return report
    if tree.num_entries == 0:
        report.emit("WARN", "no reaction entries; numerical checks skipped")
        return report

    seq165 = data["sequential_decay_fraction_165"]
    seq675 = data["sequential_decay_fraction_675"]
    direct165 = data["direct_decay_fraction_165"]
    direct675 = data["direct_decay_fraction_675"]
    bias = data["cross_section_bias_factor"]
    gamma_bias = data["gamma_bias_factor"]

    # 3.1 total = scale * model
    sigma_165_total_expected = data["scale_factor_165"] * data["sigma_165_model_b"]
    sigma_675_total_expected = data["scale_factor_675"] * data["sigma_675_model_b"]
    ok_total = True
    ok_total &= report.check_close(
        "sigma_165_total = scale_factor_165 * sigma_165_model",
        data["sigma_165_total_b"], sigma_165_total_expected, "sigma_165_total_b",
    )
    ok_total &= report.check_close(
        "sigma_675_total = scale_factor_675 * sigma_675_model",
        data["sigma_675_total_b"], sigma_675_total_expected, "sigma_675_total_b",
    )
    if ok_total:
        report.emit("PASS", "cross-section component algebra")
    else:
        report.emit("FAIL", "cross-section component algebra")

    # 3.2 sequential / direct split
    sigma_165_used_expected = seq165 * sigma_165_total_expected
    sigma_675_used_expected = seq675 * sigma_675_total_expected
    sigma_165_directdecay_expected = direct165 * sigma_165_total_expected
    sigma_675_directdecay_expected = direct675 * sigma_675_total_expected
    sigma_directdecay_expected = sigma_165_directdecay_expected + sigma_675_directdecay_expected
    sigma_total_used_expected = sigma_165_used_expected + sigma_675_used_expected + sigma_directdecay_expected

    ok_split = True
    ok_split &= report.check_close(
        "sigma_165_used = seq165 * sigma_165_total", data["sigma_165_used_b"], sigma_165_used_expected,
        "sigma_165_used_b", related_setting="/h11b/165SequentialDecayFraction",
    )
    ok_split &= report.check_close(
        "sigma_675_used = seq675 * sigma_675_total", data["sigma_675_used_b"], sigma_675_used_expected,
        "sigma_675_used_b", related_setting="/h11b/675SequentialDecayFraction",
    )
    ok_split &= report.check_close(
        "sigma_165_directdecay = (1-seq165) * sigma_165_total", data["sigma_165_directdecay_b"],
        sigma_165_directdecay_expected, "sigma_165_directdecay_b",
        related_setting="/h11b/enableDirectDecay, /h11b/165SequentialDecayFraction",
    )
    ok_split &= report.check_close(
        "sigma_675_directdecay = (1-seq675) * sigma_675_total", data["sigma_675_directdecay_b"],
        sigma_675_directdecay_expected, "sigma_675_directdecay_b",
        related_setting="/h11b/enableDirectDecay, /h11b/675SequentialDecayFraction",
    )
    ok_split &= report.check_close(
        "sigma_directdecay = sigma_165_directdecay + sigma_675_directdecay", data["sigma_directdecay_b"],
        sigma_directdecay_expected, "sigma_directdecay_b",
    )
    ok_split &= report.check_close(
        "sigma_total_used = sigma_165_used + sigma_675_used + sigma_directdecay",
        data["sigma_total_used_b"], sigma_total_used_expected, "sigma_total_used_b",
    )
    # When enableDirectDecay, sequential + direct must reconstitute the un-split total.
    sigma_total_expected_from_total = sigma_165_total_expected + sigma_675_total_expected
    ok_split &= report.check_close(
        "sigma_total_used == sigma_165_total + sigma_675_total (direct decay is a branch, not extra background)",
        data["sigma_total_used_b"], sigma_total_expected_from_total, "sigma_total_used_b",
    )
    if not (bias >= 1.0 - DEFAULT_TOL).all():
        report.emit("WARN", "cross_section_bias_factor < 1 encountered")
    if ok_split:
        report.emit("PASS", "sequential/direct split")
    else:
        report.emit("FAIL", "sequential/direct split")

    # Degenerate case: enableDirectDecay = false must force seq=1, direct sigma=0.
    disabled_mask = data["enable_direct_decay"] == 0
    if np.any(disabled_mask):
        ok_disabled = True
        ok_disabled &= report.check_close(
            "seq165 == 1 when enableDirectDecay=false", seq165[disabled_mask],
            np.ones_like(seq165[disabled_mask]), "sequential_decay_fraction_165",
        )
        ok_disabled &= report.check_close(
            "seq675 == 1 when enableDirectDecay=false", seq675[disabled_mask],
            np.ones_like(seq675[disabled_mask]), "sequential_decay_fraction_675",
        )
        ok_disabled &= report.check_close(
            "sigma_directdecay == 0 when enableDirectDecay=false", data["sigma_directdecay_b"][disabled_mask],
            np.zeros_like(data["sigma_directdecay_b"][disabled_mask]), "sigma_directdecay_b",
        )
        if ok_disabled:
            report.emit("PASS", "enableDirectDecay=false degenerate case")
        else:
            report.emit("FAIL", "enableDirectDecay=false degenerate case")

    # 3.3 sampling cross sections
    sigma_165_sampling_expected = data["sigma_165_used_b"] * bias
    sigma_675_sampling_expected = data["sigma_675_used_b"] * bias
    sigma_165_directdecay_sampling_expected = data["sigma_165_directdecay_b"] * bias
    sigma_675_directdecay_sampling_expected = data["sigma_675_directdecay_b"] * bias
    sigma_directdecay_sampling_expected = (
        sigma_165_directdecay_sampling_expected + sigma_675_directdecay_sampling_expected
    )
    sigma_3alpha_sampling_total_expected = (
        sigma_165_sampling_expected + sigma_675_sampling_expected + sigma_directdecay_sampling_expected
    )

    ok_sampling = True
    ok_sampling &= report.check_close(
        "sigma_165_sampling = sigma_165_used * bias", data["sigma_165_sampling_b"],
        sigma_165_sampling_expected, "sigma_165_sampling_b", related_setting="/h11b/crossSectionBiasFactor",
    )
    ok_sampling &= report.check_close(
        "sigma_675_sampling = sigma_675_used * bias", data["sigma_675_sampling_b"],
        sigma_675_sampling_expected, "sigma_675_sampling_b", related_setting="/h11b/crossSectionBiasFactor",
    )
    ok_sampling &= report.check_close(
        "sigma_165_directdecay_sampling = sigma_165_directdecay * bias",
        data["sigma_165_directdecay_sampling_b"], sigma_165_directdecay_sampling_expected,
        "sigma_165_directdecay_sampling_b",
    )
    ok_sampling &= report.check_close(
        "sigma_675_directdecay_sampling = sigma_675_directdecay * bias",
        data["sigma_675_directdecay_sampling_b"], sigma_675_directdecay_sampling_expected,
        "sigma_675_directdecay_sampling_b",
    )
    ok_sampling &= report.check_close(
        "sigma_directdecay_sampling = sum of directdecay sampling parts",
        data["sigma_directdecay_sampling_b"], sigma_directdecay_sampling_expected,
        "sigma_directdecay_sampling_b",
    )
    ok_sampling &= report.check_close(
        "sigma_3alpha_sampling_total = sum of 165/675/directdecay sampling",
        data["sigma_3alpha_sampling_total_b"], sigma_3alpha_sampling_total_expected,
        "sigma_3alpha_sampling_total_b",
    )
    if ok_sampling:
        report.emit("PASS", "sampling cross-section algebra")
    else:
        report.emit("FAIL", "sampling cross-section algebra")

    # Gamma sampling algebra
    total_bias = bias * gamma_bias
    sigma_gamma_165_0_sampling_expected = data["sigma_gamma_165_0_physical_b"] * total_bias
    sigma_gamma_165_1_sampling_expected = data["sigma_gamma_165_1_physical_b"] * total_bias
    sigma_gamma_675_sampling_expected = data["sigma_gamma_675_physical_b"] * total_bias
    sigma_gamma_sampling_total_expected = (
        sigma_gamma_165_0_sampling_expected + sigma_gamma_165_1_sampling_expected + sigma_gamma_675_sampling_expected
    )

    ok_gamma = True
    ok_gamma &= report.check_close(
        "sigma_gamma_165_0_sampling = physical * bias * gammaBias", data["sigma_gamma_165_0_sampling_b"],
        sigma_gamma_165_0_sampling_expected, "sigma_gamma_165_0_sampling_b",
        related_setting="/h11b/gammaBiasFactor",
    )
    ok_gamma &= report.check_close(
        "sigma_gamma_165_1_sampling = physical * bias * gammaBias", data["sigma_gamma_165_1_sampling_b"],
        sigma_gamma_165_1_sampling_expected, "sigma_gamma_165_1_sampling_b",
        related_setting="/h11b/gammaBiasFactor",
    )
    ok_gamma &= report.check_close(
        "sigma_gamma_675_sampling = physical * bias * gammaBias", data["sigma_gamma_675_sampling_b"],
        sigma_gamma_675_sampling_expected, "sigma_gamma_675_sampling_b",
        related_setting="/h11b/gammaBiasFactor",
    )
    ok_gamma &= report.check_close(
        "sigma_gamma_sampling_total = sum of gamma sampling parts", data["sigma_gamma_sampling_total_b"],
        sigma_gamma_sampling_total_expected, "sigma_gamma_sampling_total_b",
    )
    if ok_gamma:
        report.emit("PASS", "gamma sampling algebra")
    else:
        report.emit("FAIL", "gamma sampling algebra")

    # Total sampling cross section
    sigma_total_sampling_all_expected = data["sigma_3alpha_sampling_total_b"] + data["sigma_gamma_sampling_total_b"]
    ok_grand_total = report.check_close(
        "sigma_total_sampling_all = sigma_3alpha_sampling_total + sigma_gamma_sampling_total",
        data["sigma_total_sampling_all_b"], sigma_total_sampling_all_expected, "sigma_total_sampling_all_b",
    )
    if ok_grand_total:
        report.emit("PASS", "total sampling cross section")
    else:
        report.emit("FAIL", "total sampling cross section")

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
