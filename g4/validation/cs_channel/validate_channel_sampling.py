#!/usr/bin/env python3
"""Validate sampled reaction_channel event counts against sampling cross sections.

See gpt/codex_cs_channel_validation.md section 4 and 5.2.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from _common import Report, load_branches, open_tree

BRANCHES = [
    "reaction_channel",
    "resonance_id",
    "sigma_165_sampling_b",
    "sigma_675_sampling_b",
    "sigma_directdecay_sampling_b",
    "sigma_165_directdecay_sampling_b",
    "sigma_675_directdecay_sampling_b",
    "sigma_3alpha_sampling_total_b",
    "sigma_gamma_sampling_total_b",
    "sigma_total_sampling_all_b",
]

PULL_FAIL_THRESHOLD = 5.0


def binomial_pull(observed_count: int, expected_p: float, n_total: int) -> float:
    if n_total == 0:
        return float("nan")
    expected_count = expected_p * n_total
    sigma = np.sqrt(max(expected_p * (1.0 - expected_p), 0.0) * n_total)
    if sigma <= 0.0:
        return 0.0 if observed_count == expected_count else float("inf")
    return (observed_count - expected_count) / sigma


def validate_file(path: Path) -> Report:
    report = Report(path)
    tree = open_tree(path)
    data, missing = load_branches(tree, BRANCHES)
    if missing:
        report.emit("FAIL", "required branches present (see MISSING BRANCH lines above)")
        return report
    n_total = tree.num_entries
    if n_total == 0:
        report.emit("WARN", "no reaction entries; sampling checks skipped")
        return report

    channels = data["reaction_channel"]
    resonance = data["resonance_id"]

    # sigma_*_sampling_b vary event-to-event: the proton loses energy inside
    # the target before reacting, and near a narrow resonance that straggling
    # changes the cross-section components substantially. Each event's actual
    # sampling probability is sigma_channel_sampling_b[i] / sigma_total_sampling_all_b[i],
    # so the expected aggregate probability is the mean of that per-event
    # ratio -- mean(branch)/mean(total_branch) is a biased (Jensen's-gap)
    # estimator here and was found to disagree with observed counts by >20
    # sigma even though the underlying simulation is correct.
    sigma_total_all = data["sigma_total_sampling_all_b"]
    if np.any(sigma_total_all <= 0.0):
        report.emit("FAIL", "sigma_total_sampling_all_b <= 0 for some events; cannot form probabilities")
        return report

    p_165 = float(np.mean(data["sigma_165_sampling_b"] / sigma_total_all))
    p_675 = float(np.mean(data["sigma_675_sampling_b"] / sigma_total_all))
    p_directdecay = float(np.mean(data["sigma_directdecay_sampling_b"] / sigma_total_all))
    p_gamma = float(np.mean(data["sigma_gamma_sampling_total_b"] / sigma_total_all))

    counts = {c: int(np.count_nonzero(channels == c)) for c in (0, 1, 2, 3)}

    print("channel  expected_probability  observed_probability  count  pull  status")
    rows = [
        (0, "165", p_165),
        (1, "675", p_675),
        (2, "directdecay", p_directdecay),
        (3, "gamma", p_gamma),
    ]
    ok_counts = True
    for channel, label, expected_p in rows:
        observed_p = counts[channel] / n_total
        pull = binomial_pull(counts[channel], expected_p, n_total)
        status = "FAIL" if abs(pull) > PULL_FAIL_THRESHOLD else "PASS"
        if status == "FAIL":
            ok_counts = False
        print(f"{channel} ({label})  {expected_p:.6g}  {observed_p:.6g}  {counts[channel]}  {pull:.3g}  {status}")

    if ok_counts:
        report.emit("PASS", "reaction channel count vs sampling probability")
    else:
        report.emit("FAIL", "reaction channel count vs sampling probability")

    # Probability normalization
    prob_sum = p_165 + p_675 + p_directdecay + p_gamma
    if abs(prob_sum - 1.0) < 1e-6:
        report.emit("PASS", f"sampling probability normalization (sum = {prob_sum:.9g})")
    else:
        report.emit("FAIL", f"sampling probability normalization (sum = {prob_sum:.9g})")
        report.fail_detail("p_165+p_675+p_directdecay+p_gamma", 1.0, prob_sum)

    # Direct decay resonance split, using resonance_id for channel==2 events.
    direct_mask = channels == 2
    n_direct = int(np.count_nonzero(direct_mask))
    if n_direct > 0:
        n_direct_165 = int(np.count_nonzero(resonance[direct_mask] == 165))
        n_direct_675 = int(np.count_nonzero(resonance[direct_mask] == 675))
        unresolved = n_direct - n_direct_165 - n_direct_675
        if unresolved != 0:
            report.emit(
                "WARN",
                f"{unresolved} direct-decay events have resonance_id outside {{165, 675}}",
            )

        # As above, use the per-event ratio conditioned on the events that
        # actually landed in the direct-decay channel (the same subset whose
        # 165-vs-675 coin flip this ratio drove), not a global mean/mean.
        direct_total = data["sigma_directdecay_sampling_b"][direct_mask]
        if np.all(direct_total > 0.0):
            expected_frac_165 = float(np.mean(data["sigma_165_directdecay_sampling_b"][direct_mask] / direct_total))
            expected_frac_675 = float(np.mean(data["sigma_675_directdecay_sampling_b"][direct_mask] / direct_total))
            observed_frac_165 = n_direct_165 / n_direct
            observed_frac_675 = n_direct_675 / n_direct
            pull_165 = binomial_pull(n_direct_165, expected_frac_165, n_direct)
            pull_675 = binomial_pull(n_direct_675, expected_frac_675, n_direct)
            print(
                "direct-from-165: expected="
                f"{expected_frac_165:.6g} observed={observed_frac_165:.6g} pull={pull_165:.3g}"
            )
            print(
                "direct-from-675: expected="
                f"{expected_frac_675:.6g} observed={observed_frac_675:.6g} pull={pull_675:.3g}"
            )
            if abs(pull_165) <= PULL_FAIL_THRESHOLD and abs(pull_675) <= PULL_FAIL_THRESHOLD and unresolved == 0:
                report.emit("PASS", "direct decay resonance split")
            else:
                report.emit("FAIL", "direct decay resonance split")
        else:
            report.emit("WARN", "direct decay sampling total is zero; cannot check resonance split")
    else:
        report.emit("WARN", "no direct-decay (channel=2) events in this file; resonance split not tested")

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
