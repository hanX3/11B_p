#!/usr/bin/env python3
"""Validate 675-keV strict three-alpha model diagnostics and responses."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import _common as c


BRANCHES = [
    "reaction_channel", "resonance_id",
    "h11b675_alpha_decay_model",
    "h11b675_strict_coherent_l13",
    "h11b675_strict_permutation_symmetrized",
    "h11b675_strict_l1_fraction",
    "h11b675_strict_l13_phase",
    "h11b675_strict_weight",
    "h11b675_strict_weight_max",
    "h11b675_strict_sampling_attempts",
    "h11b675_strict_max_sampling_attempts",
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "eaa_8Be",
    "opening_angle_alpha12_cm", "opening_angle_alpha13_cm", "opening_angle_alpha23_cm",
    "cos_chi_secondary_8be",
]

EXPECTED = {
    "default": {"l1": 0.76, "coherent": 1, "sym": 1, "phase": 4.209734155810323},
    "l1only": {"l1": 1.0, "coherent": 1, "sym": 1, "phase": 4.209734155810323},
    "l3only": {"l1": 0.0, "coherent": 1, "sym": 1, "phase": 4.209734155810323},
    "incoherent": {"l1": 0.76, "coherent": 0, "sym": 1, "phase": 4.209734155810323},
    "no_sym": {"l1": 0.76, "coherent": 1, "sym": 0, "phase": 4.209734155810323},
    "phase0": {"l1": 0.76, "coherent": 1, "sym": 1, "phase": 0.0},
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, help="label:path")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--max-attempts", type=int, default=10000)
    return parser.parse_args()


def parse_inputs(items):
    result = {}
    for item in items:
        if ":" not in item:
            raise SystemExit(f"--inputs entries must be label:path: {item}")
        label, path = item.split(":", 1)
        result[label] = path
    return result


def filtered_sample(path):
    data, missing, _ = c.load_arrays(path, BRANCHES)
    hard_missing = [m for m in missing if m != "h11b675_strict_max_sampling_attempts"]
    if hard_missing:
        raise SystemExit(f"{path}: missing required branches {hard_missing}")
    mask = (data["reaction_channel"] == 1)
    if "resonance_id" in data:
        mask &= data["resonance_id"] == 675
    return data, mask, missing


def dalitz(data, mask):
    e1 = data["e_3alpha_cm_alpha1"][mask]
    e2 = data["e_3alpha_cm_alpha2"][mask]
    e3 = data["e_3alpha_cm_alpha3"][mask]
    total = e1 + e2 + e3
    x = np.sqrt(3.0) * (e2 - e3) / total
    y = (2.0 * e1 - e2 - e3) / total
    finite = np.isfinite(x) & np.isfinite(y)
    return x[finite], y[finite]


def hist2(x, y):
    h, _, _ = np.histogram2d(x, y, bins=50, range=[[-1.8, 1.8], [-1.8, 1.8]], density=True)
    return h


def hist_distance(a, b):
    denom = a + b
    mask = denom > 0
    return float(0.5 * np.mean(((a[mask] - b[mask]) ** 2) / denom[mask])) if np.any(mask) else float("nan")


def plot_dalitz(outdir, label, x, y):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    ax.hist2d(x, y, bins=60, range=[[-1.8, 1.8], [-1.8, 1.8]], cmap="viridis")
    ax.set_xlabel("sqrt(3) * (E2 - E3) / sum(E)")
    ax.set_ylabel("(2E1 - E2 - E3) / sum(E)")
    ax.set_title(f"675 strict {label} (n={x.size})")
    fig.tight_layout()
    path = Path(outdir) / f"dalitz_675_strict_{label}.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return str(path)


def plot_diff(outdir, label, h_a, h_b):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    image = ax.imshow((h_a - h_b).T, origin="lower", extent=[-1.8, 1.8, -1.8, 1.8], cmap="coolwarm", aspect="auto")
    fig.colorbar(image, ax=ax, label="density difference")
    ax.set_xlabel("Dalitz x")
    ax.set_ylabel("Dalitz y")
    ax.set_title(label.replace("_", " "))
    fig.tight_layout()
    path = Path(outdir) / f"dalitz_diff_{label}.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return str(path)


def main() -> int:
    args = parse_args()
    inputs = parse_inputs(args.inputs)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = {}
    tests = []
    plots = []
    missing_notes = []

    for label, path in inputs.items():
        data, mask, missing = filtered_sample(path)
        if "h11b675_strict_max_sampling_attempts" in missing:
            missing_notes.append(f"{label}: missing h11b675_strict_max_sampling_attempts; used CLI --max-attempts={args.max_attempts}")
        n = int(np.sum(mask))
        x, y = dalitz(data, mask)
        h = hist2(x, y)
        plot = plot_dalitz(outdir, label, x, y)
        plots.append(plot)

        weight = data["h11b675_strict_weight"][mask]
        weight_max = data["h11b675_strict_weight_max"][mask]
        attempts = data["h11b675_strict_sampling_attempts"][mask]
        param_status = "PASS"
        details = {}
        expected = EXPECTED.get(label, {})
        if expected:
            l1 = c.safe_median(data, "h11b675_strict_l1_fraction", mask)
            coherent = c.safe_median(data, "h11b675_strict_coherent_l13", mask)
            sym = c.safe_median(data, "h11b675_strict_permutation_symmetrized", mask)
            phase = c.safe_median(data, "h11b675_strict_l13_phase", mask)
            details = {"l1_fraction": l1, "coherent": coherent, "symmetrized": sym, "phase_rad": phase}
            checks = [
                abs(l1 - expected["l1"]) < 1e-9,
                abs(coherent - expected["coherent"]) < 1e-9,
                abs(sym - expected["sym"]) < 1e-9,
                abs(phase - expected["phase"]) < 1e-9,
            ]
            if not all(checks):
                param_status = "FAIL"

        weight_status = "PASS"
        if n == 0 or np.any(weight < -1e-12) or np.any(weight - weight_max > 1e-9):
            weight_status = "FAIL"
        elif np.mean(attempts == args.max_attempts) >= 1e-3 or np.mean(attempts) >= 100:
            weight_status = "WARN"

        tests.append({
            "test_name": f"675_strict_{label}_parameters",
            "status": param_status,
            "n_events": n,
            **details,
        })
        tests.append({
            "test_name": f"675_strict_{label}_rejection_sampling",
            "status": weight_status,
            "n_events": n,
            "weight_min": float(np.min(weight)) if n else float("nan"),
            "weight_max_observed": float(np.max(weight)) if n else float("nan"),
            "weight_max_limit_min": float(np.min(weight_max)) if n else float("nan"),
            "mean_sampling_attempts": float(np.mean(attempts)) if n else float("nan"),
            "fraction_at_max_attempts": float(np.mean(attempts == args.max_attempts)) if n else float("nan"),
        })

        samples[label] = {
            "data": data,
            "mask": mask,
            "n": n,
            "x": x,
            "y": y,
            "hist": h,
            "plot": plot,
        }

    def comparison(name, a, b, observable, threshold=0.05):
        stat = c.ks_two_sample(samples[a]["data"][observable][samples[a]["mask"]], samples[b]["data"][observable][samples[b]["mask"]])
        status = "PASS" if stat > threshold else "FAIL"
        tests.append({"test_name": name, "status": status, "ks_statistic": stat, "threshold": threshold, "observable": observable})

    if "l1only" in samples and "l3only" in samples:
        d = hist_distance(samples["l1only"]["hist"], samples["l3only"]["hist"])
        tests.append({"test_name": "675 strict L1-only vs L3-only Dalitz response", "status": "PASS" if d > 0.01 else "FAIL", "hist_chi_distance": d, "threshold": 0.01})
        comparison("675 strict L1-only vs L3-only cos_chi response", "l1only", "l3only", "cos_chi_secondary_8be")
    if "default" in samples and "incoherent" in samples:
        d = hist_distance(samples["default"]["hist"], samples["incoherent"]["hist"])
        tests.append({"test_name": "675 strict coherent vs incoherent Dalitz response", "status": "PASS" if d > 0.002 else "FAIL", "hist_chi_distance": d, "threshold": 0.002})
        comparison("675 strict coherent vs incoherent eaa response", "default", "incoherent", "eaa_8Be", threshold=0.03)
        plots.append(plot_diff(outdir, "default_minus_incoherent", samples["default"]["hist"], samples["incoherent"]["hist"]))
    if "default" in samples and "phase0" in samples:
        d = hist_distance(samples["default"]["hist"], samples["phase0"]["hist"])
        tests.append({"test_name": "675 strict phase response Dalitz", "status": "PASS" if d > 0.002 else "FAIL", "hist_chi_distance": d, "threshold": 0.002})
        comparison("675 strict phase response cos_chi", "default", "phase0", "cos_chi_secondary_8be", threshold=0.03)
        plots.append(plot_diff(outdir, "default_minus_phase0", samples["default"]["hist"], samples["phase0"]["hist"]))
    if "default" in samples and "no_sym" in samples:
        plots.append(plot_diff(outdir, "default_minus_no_sym", samples["default"]["hist"], samples["no_sym"]["hist"]))

    if "default" in samples:
        d = samples["default"]["data"]
        m = samples["default"]["mask"]
        e1, e2, e3 = d["e_3alpha_cm_alpha1"][m], d["e_3alpha_cm_alpha2"][m], d["e_3alpha_cm_alpha3"][m]
        a12, a13, a23 = d["opening_angle_alpha12_cm"][m], d["opening_angle_alpha13_cm"][m], d["opening_angle_alpha23_cm"][m]
        kse = max(c.ks_two_sample(e1, e2), c.ks_two_sample(e1, e3), c.ks_two_sample(e2, e3))
        ksa = max(c.ks_two_sample(a12, a13), c.ks_two_sample(a12, a23), c.ks_two_sample(a13, a23))
        # Mean-gap check as a significance (z-score), not a fixed absolute MeV
        # tolerance.  The CM alpha energy is broad (sigma ~ 1.5 MeV at 675 keV),
        # so the standard error of each mean is sigma/sqrt(n); a fixed 0.02 MeV
        # tolerance sits *below* that noise floor at n>~20k and produces a
        # ~25-30% spurious-FAIL rate on a perfectly symmetric generator.  A
        # z-score is self-calibrating: a genuine primary/secondary asymmetry
        # (e.g. permutation_symmetrized=false) shows up at >100 sigma, while
        # statistical fluctuations of a symmetric sample stay within a few sigma.
        def _mean_gap_sigma(x, y):
            gap = float(np.mean(x) - np.mean(y))
            se = float(np.sqrt(np.var(x) / x.size + np.var(y) / y.size)) if x.size and y.size else float("inf")
            return gap, se, (abs(gap) / se if se > 0 else float("inf"))

        gaps = [_mean_gap_sigma(e1, e2), _mean_gap_sigma(e1, e3), _mean_gap_sigma(e2, e3)]
        mean_gap = max(abs(g) for g, _, _ in gaps)
        mean_gap_sigma_max = max(z for _, _, z in gaps)
        threshold_ks = 0.08 if samples["default"]["n"] < 2000 else 0.03
        # 5-sigma keeps the false-FAIL probability negligible while still
        # catching any real (>~10 keV at this n) permutation asymmetry.
        threshold_mean_sigma = 5.0
        tests.append({
            "test_name": "675 strict permutation symmetrization",
            "status": "PASS" if kse < threshold_ks and ksa < threshold_ks and mean_gap_sigma_max < threshold_mean_sigma else "FAIL",
            "energy_ks_max": kse,
            "opening_angle_ks_max": ksa,
            "energy_mean_gap_max_MeV": mean_gap,
            "energy_mean_gap_max_sigma": mean_gap_sigma_max,
            "ks_threshold": threshold_ks,
            "mean_gap_threshold_sigma": threshold_mean_sigma,
        })

    status = c.combine_status(t["status"] for t in tests)
    result = {"status": status, "tests": tests, "plots": plots, "missing_branch_notes": missing_notes}
    c.write_json(args.json, result)
    for test in tests:
        print(f"[{test['status']}] {test['test_name']}")
    if missing_notes:
        for note in missing_notes:
            print(f"[INFO] {note}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())

