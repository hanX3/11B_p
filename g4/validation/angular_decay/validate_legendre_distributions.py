#!/usr/bin/env python3
"""Validate an on/off Legendre angular-distribution pair."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import _common as c


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--off", required=True)
    parser.add_argument("--on", required=True)
    parser.add_argument("--observable", required=True)
    parser.add_argument("--kind", choices=["a1a2", "a2a4"], required=True)
    parser.add_argument("--test-name", required=True)
    parser.add_argument("--select", action="append", default=[], help="Selection branch=value, repeatable")
    parser.add_argument("--a1-branch", default=None)
    parser.add_argument("--a2-branch", default=None)
    parser.add_argument("--a4-branch", default=None)
    parser.add_argument("--mode-branch", default=None)
    parser.add_argument("--expected-on-mode", type=int, default=None)
    parser.add_argument("--expected-off-mode", type=int, default=0)
    parser.add_argument("--out", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--bins", type=int, default=30)
    return parser.parse_args()


def coefficient_values(data, mask, args):
    if args.kind == "a1a2":
        a1 = c.safe_median(data, args.a1_branch, mask, 0.0)
        a2 = c.safe_median(data, args.a2_branch, mask, 0.0)
        return (a1, a2)
    a2 = c.safe_median(data, args.a2_branch, mask, 0.0)
    a4 = c.safe_median(data, args.a4_branch, mask, 0.0)
    return (a2, a4)


def main() -> int:
    args = parse_args()
    branches = [args.observable] + [s.split("=", 1)[0] for s in args.select]
    for branch in (args.a1_branch, args.a2_branch, args.a4_branch, args.mode_branch):
        if branch:
            branches.append(branch)

    off_data, off_missing, _ = c.load_arrays(args.off, branches)
    on_data, on_missing, _ = c.load_arrays(args.on, branches)
    missing = sorted(set(off_missing + on_missing))
    if missing:
        result = {"test_name": args.test_name, "status": "FAIL", "missing_branches": missing}
        c.write_json(args.json, result)
        print(f"[FAIL] {args.test_name}: missing branches {missing}")
        return 1

    off_mask = c.apply_selection(off_data, args.select)
    on_mask = c.apply_selection(on_data, args.select)
    off_values = c.finite(off_data[args.observable][off_mask])
    on_values = c.finite(on_data[args.observable][on_mask])
    off_values = off_values[(off_values >= -1.0) & (off_values <= 1.0)]
    on_values = on_values[(on_values >= -1.0) & (on_values <= 1.0)]

    coeffs = coefficient_values(on_data, on_mask, args)
    edges = np.linspace(-1.0, 1.0, args.bins + 1)
    uniform_probs, _, _ = c.expected_probabilities(edges, "uniform", ())
    theory_probs, min_raw_weight, max_clipped_weight = c.expected_probabilities(edges, args.kind, coeffs)

    off_chi2, off_ndf = c.chi2_ndf(off_values, edges, uniform_probs)
    on_chi2, on_ndf = c.chi2_ndf(on_values, edges, theory_probs, n_parameters=2)
    off_ks = c.ks_one_sample(off_values, "uniform", ())
    on_ks = c.ks_one_sample(on_values, args.kind, coeffs)
    on_off_ks = c.ks_two_sample(on_values, off_values)

    mode_checks = []
    if args.mode_branch:
        off_mode = c.safe_median(off_data, args.mode_branch, off_mask, float("nan"))
        on_mode = c.safe_median(on_data, args.mode_branch, on_mask, float("nan"))
        mode_checks.append(abs(off_mode - args.expected_off_mode) < 1e-9)
        if args.expected_on_mode is not None:
            mode_checks.append(abs(on_mode - args.expected_on_mode) < 1e-9)
    else:
        off_mode = on_mode = float("nan")

    statuses = []
    if off_values.size < 50 or on_values.size < 50:
        statuses.append("WARN")
    if off_chi2 < 3.0 and off_ks < 0.15:
        statuses.append("PASS")
    else:
        statuses.append("FAIL")
    if on_chi2 < 5.0 and on_ks < 0.15:
        statuses.append("PASS")
    else:
        statuses.append("FAIL")
    if on_off_ks > 0.03:
        statuses.append("PASS")
    else:
        statuses.append("WARN")
    if mode_checks and not all(mode_checks):
        statuses.append("FAIL")
    status = c.combine_status(statuses)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        centers = 0.5 * (edges[:-1] + edges[1:])
        width = edges[1] - edges[0]
        off_hist, _ = np.histogram(off_values, bins=edges, density=True)
        on_hist, _ = np.histogram(on_values, bins=edges, density=True)
        theory_density = theory_probs / width

        fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.0, 6.0), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
        ax.step(centers, off_hist, where="mid", label=f"off n={off_values.size}")
        ax.step(centers, on_hist, where="mid", label=f"on n={on_values.size}")
        ax.plot(centers, theory_density, "k--", label=f"theory {args.kind} coeffs={tuple(round(v, 4) for v in coeffs)}")
        ax.set_ylabel("density")
        ax.set_title(args.test_name)
        ax.legend(fontsize=8)
        ratio = np.divide(on_hist, theory_density, out=np.full_like(on_hist, np.nan), where=theory_density > 0)
        rx.axhline(1.0, color="k", linewidth=0.8)
        rx.plot(centers, ratio, "o", ms=3)
        rx.set_xlabel(args.observable)
        rx.set_ylabel("on/theory")
        fig.tight_layout()
        fig.savefig(args.out, dpi=140)
        plt.close(fig)
    except ImportError:
        pass

    result = {
        "test_name": args.test_name,
        "status": status,
        "observable": args.observable,
        "n_off": int(off_values.size),
        "n_on": int(on_values.size),
        "coefficients": list(coeffs),
        "off_chi2_ndf": off_chi2,
        "on_chi2_ndf": on_chi2,
        "off_ks_statistic": off_ks,
        "on_ks_statistic": on_ks,
        "on_off_ks_statistic": on_off_ks,
        "off_mode_median": off_mode,
        "on_mode_median": on_mode,
        "minimum_raw_weight": min_raw_weight,
        "maximum_clipped_weight": max_clipped_weight,
        "off_moments": c.moment_summary(off_values),
        "on_moments": c.moment_summary(on_values),
        "plot": args.out,
        "missing_branches": [],
    }
    c.write_json(args.json, result)
    print(f"[{status}] {args.test_name}: n_off={off_values.size} n_on={on_values.size} off_chi2/ndf={off_chi2:.3g} on_chi2/ndf={on_chi2:.3g} on/off_KS={on_off_ks:.3g}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())

