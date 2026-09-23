#!/usr/bin/env python3
"""Validate gamma-capture resonance/branch selection against sampling cross sections.

See gpt/codex_cs_channel_validation.md section 5.3.  gamma_resonance and
gamma_branch numeric codes come from include/Constants.hh
(H11BGammaResonance, H11BGammaBranch).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from _common import Report, load_branches, open_tree

BRANCHES = [
    "reaction_channel",
    "gamma_resonance",
    "gamma_branch",
    "gamma_branch_fraction_used",
    "gamma_relative_intensity_used",
    "gamma_primary_energy_MeV",
    "gamma_final_state_energy_MeV",
    "n_prompt_gammas",
    "sigma_gamma_165_0_sampling_b",
    "sigma_gamma_165_1_sampling_b",
    "sigma_gamma_675_sampling_b",
    "sigma_gamma_sampling_total_b",
]

GAMMA_RESONANCE_165 = 165
GAMMA_RESONANCE_675 = 675
GAMMA_BRANCH_165_GROUND = 1
GAMMA_BRANCH_165_FIRST_EXCITED = 2
PULL_FAIL_THRESHOLD = 5.0


def binomial_pull(observed_count: int, expected_p: float, n_total: int) -> float:
    if n_total == 0:
        return float("nan")
    expected_count = expected_p * n_total
    sigma = np.sqrt(max(expected_p * (1.0 - expected_p), 0.0) * n_total)
    if sigma <= 0.0:
        return 0.0 if observed_count == expected_count else float("inf")
    return (observed_count - expected_count) / sigma


def validate_file(path: Path, macro_hint: str | None) -> Report:
    report = Report(path)
    tree = open_tree(path)
    data, missing = load_branches(tree, BRANCHES)
    if missing:
        report.emit("FAIL", "required branches present (see MISSING BRANCH lines above)")
        return report

    gamma_mask = data["reaction_channel"] == 3
    n_gamma = int(np.count_nonzero(gamma_mask))
    if n_gamma == 0:
        if macro_hint is not None:
            # *_gamma_only.mac macros exist specifically to isolate gamma-capture
            # events; producing none contradicts their documented purpose
            # (gpt/codex_cs_channel_validation.md section 7.5), so this is a
            # real failure, not just an empty sample.
            report.emit("FAIL", f"no gamma-capture (reaction_channel=3) events in this file (macro_hint={macro_hint})")
        else:
            report.emit("WARN", "no gamma-capture (reaction_channel=3) events in this file")
        return report

    gamma_resonance = data["gamma_resonance"][gamma_mask]
    gamma_branch = data["gamma_branch"][gamma_mask]

    sigma_g165_0 = np.mean(data["sigma_gamma_165_0_sampling_b"])
    sigma_g165_1 = np.mean(data["sigma_gamma_165_1_sampling_b"])
    sigma_g675 = np.mean(data["sigma_gamma_675_sampling_b"])
    sigma_g_total = np.mean(data["sigma_gamma_sampling_total_b"])

    if sigma_g_total <= 0.0:
        report.emit("FAIL", "sigma_gamma_sampling_total_b <= 0; cannot form gamma branch probabilities")
        return report

    p_165_g0 = sigma_g165_0 / sigma_g_total
    p_165_g1 = sigma_g165_1 / sigma_g_total
    p_675 = sigma_g675 / sigma_g_total

    n_165_g0 = int(np.count_nonzero((gamma_resonance == GAMMA_RESONANCE_165) & (gamma_branch == GAMMA_BRANCH_165_GROUND)))
    n_165_g1 = int(
        np.count_nonzero((gamma_resonance == GAMMA_RESONANCE_165) & (gamma_branch == GAMMA_BRANCH_165_FIRST_EXCITED))
    )
    n_675 = int(np.count_nonzero(gamma_resonance == GAMMA_RESONANCE_675))

    print("branch          expected_probability  observed_probability  count  pull  status")
    ok_resonance = True
    for label, count, expected_p in (
        ("165_gamma0", n_165_g0, p_165_g0),
        ("165_gamma1", n_165_g1, p_165_g1),
        ("675", n_675, p_675),
    ):
        observed_p = count / n_gamma
        pull = binomial_pull(count, expected_p, n_gamma)
        status = "FAIL" if abs(pull) > PULL_FAIL_THRESHOLD else "PASS"
        if status == "FAIL":
            ok_resonance = False
        print(f"{label:<15} {expected_p:.6g}  {observed_p:.6g}  {count}  {pull:.3g}  {status}")

    if ok_resonance:
        report.emit("PASS", "gamma resonance selection")
    else:
        report.emit("FAIL", "gamma resonance selection")

    # Branch selection: every gamma event should carry a non-None branch id,
    # and 675 events should only use 675 branch codes (>= Gamma675ToGroundState).
    invalid_165_branch = np.count_nonzero(
        (gamma_resonance == GAMMA_RESONANCE_165)
        & ~np.isin(gamma_branch, [GAMMA_BRANCH_165_GROUND, GAMMA_BRANCH_165_FIRST_EXCITED])
    )
    invalid_675_branch = np.count_nonzero((gamma_resonance == GAMMA_RESONANCE_675) & (gamma_branch < 3))
    if invalid_165_branch == 0 and invalid_675_branch == 0:
        report.emit("PASS", "gamma branch selection")
    else:
        report.emit("FAIL", "gamma branch selection")
        if invalid_165_branch:
            print(f"    {invalid_165_branch} events with gamma_resonance=165 have an unexpected gamma_branch code")
        if invalid_675_branch:
            print(f"    {invalid_675_branch} events with gamma_resonance=675 have an unexpected gamma_branch code")

    # Macro-specific expectations for the *_gamma_only.mac validation macros.
    if macro_hint == "165gamma_only":
        bad = n_675
        if bad == 0:
            report.emit("PASS", "gamma-only macro behavior (validation_165gamma_only.mac: no 675 gamma branch)")
        else:
            report.emit("FAIL", "gamma-only macro behavior (validation_165gamma_only.mac: no 675 gamma branch)")
            print(f"    found {bad} gamma_resonance=675 events; expected 0")
    elif macro_hint == "675gamma_only":
        bad = n_165_g0 + n_165_g1
        if bad == 0:
            report.emit("PASS", "gamma-only macro behavior (validation_675gamma_only.mac: no 165 gamma branch)")
        else:
            report.emit("FAIL", "gamma-only macro behavior (validation_675gamma_only.mac: no 165 gamma branch)")
            print(f"    found {bad} gamma_resonance=165 events; expected 0")
    else:
        report.emit("WARN", "no --macro-hint given; skipped gamma-only macro behavior check")

    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_files", nargs="+", type=Path)
    parser.add_argument(
        "--macro-hint",
        choices=("165gamma_only", "675gamma_only"),
        default=None,
        help="Which gamma-only validation macro produced this file, if any.",
    )
    args = parser.parse_args()

    overall = 0
    for path in args.root_files:
        if not path.exists():
            print(f"{path}: file not found")
            overall = max(overall, 1)
            continue
        report = validate_file(path, args.macro_hint)
        print(report.summary_line())
        overall = max(overall, {"PASS": 0, "WARN": 0, "FAIL": 1}[report.status])
    return overall


if __name__ == "__main__":
    sys.exit(main())
