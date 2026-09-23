#!/usr/bin/env python3
"""Direct 3-alpha phase-space validation (task doc section 7).

Applies to reaction_channel==2.  The script rebuilds the three alpha
four-vectors from lab kinetic energies and directions, boosts them to the
3-alpha CM frame, checks energy/momentum closure, checks pair invariant-mass
limits, and reconstructs simple Dalitz coordinates from CM kinetic-energy
fractions.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

DEFAULT_TOLERANCE_MEV = 1.0e-6
KINETIC_FLOOR_MEV = -1.0e-9

BRANCHES = [
    "reaction_channel", "resonance_id",
    "e_alpha1", "e_alpha2", "e_alpha3",
    "theta_lab_alpha1", "theta_lab_alpha2", "theta_lab_alpha3",
    "phi_lab_alpha1", "phi_lab_alpha2", "phi_lab_alpha3",
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "projectile_kinetic_lab", "projectile_px_lab", "projectile_py_lab", "projectile_pz_lab", "projectile_p_lab",
    "e_cm_p11B",
]


def maybe_plot(outdir: Path | None, x, y, energy_residual, p_abs) -> None:
    if outdir is None:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plots")
        return

    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots()
    ax.scatter(x, y, s=2, alpha=0.35)
    ax.set_xlabel("Dalitz x = sqrt(3) * (epsilon2 - epsilon3)")
    ax.set_ylabel("Dalitz y = 2*epsilon1 - epsilon2 - epsilon3")
    ax.set_title("direct 3-alpha Dalitz coordinates")
    fig.savefig(outdir / "direct_dalitz.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(energy_residual, bins=60)
    ax.set_xlabel("sum_T_cm - T_available_cm (MeV)")
    ax.set_ylabel("events")
    ax.set_title("direct 3-alpha CM energy residual")
    fig.savefig(outdir / "direct_energy_sum_residual.png", dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(p_abs, bins=60)
    ax.set_xlabel("|Delta_p_lab| (MeV/c)")
    ax.set_ylabel("events")
    ax.set_title("direct 3-alpha lab momentum closure")
    fig.savefig(outdir / "direct_momentum_closure.png", dpi=120)
    plt.close(fig)


def main() -> int:
    parser = c.common_argparser(__doc__)
    args = parser.parse_args()

    tree = c.open_tree(Path(args.root))
    data, missing = c.load_branches(tree, BRANCHES, args.max_events)
    if missing:
        raise SystemExit(f"required branches missing: {missing}")

    masses_pm = c.derive_projectile_masses(data)
    m_p, m_11b = masses_pm["m_p"][0], masses_pm["m_11B"][0]

    channel_arg = c.parse_channel_arg(args.channel)
    channels = [2] if channel_arg is None else [channel_arg]
    tol = args.tolerance if args.tolerance is not None else DEFAULT_TOLERANCE_MEV
    report = c.Report(Path(args.report) if args.report else Path(args.root))
    outdir = Path(args.outdir) if args.outdir else None

    for channel in channels:
        if channel != 2:
            report.emit("WARN", f"channel {channel} is not direct 3-alpha; skipped")
            continue

        mask = data["reaction_channel"] == 2
        n = int(np.sum(mask))
        if n == 0:
            report.emit("WARN", "channel 2: no events in this file")
            continue

        d = {k: v[mask] for k, v in data.items()}
        sqrt_s = c.derive_sqrt_s(d, m_p, m_11b)
        m_alpha = c.derive_alpha_mass(
            sqrt_s, d["e_3alpha_cm_alpha1"], d["e_3alpha_cm_alpha2"], d["e_3alpha_cm_alpha3"]
        )[0]
        masses = {"m_p": masses_pm["m_p"], "m_11B": masses_pm["m_11B"], "m_alpha": (m_alpha, 0.0)}
        c.print_mass_table(masses)
        report.note(f"-- reaction_channel=2: {n} events --")

        e1, px1, py1, pz1 = c.spherical_lab_vector(d["e_alpha1"], d["theta_lab_alpha1"], d["phi_lab_alpha1"], m_alpha)
        e2, px2, py2, pz2 = c.spherical_lab_vector(d["e_alpha2"], d["theta_lab_alpha2"], d["phi_lab_alpha2"], m_alpha)
        e3, px3, py3, pz3 = c.spherical_lab_vector(d["e_alpha3"], d["theta_lab_alpha3"], d["phi_lab_alpha3"], m_alpha)

        e_initial = d["projectile_kinetic_lab"] + m_p + m_11b
        px_i, py_i, pz_i = d["projectile_px_lab"], d["projectile_py_lab"], d["projectile_pz_lab"]
        delta_e = (e1 + e2 + e3) - e_initial
        delta_px = px1 + px2 + px3 - px_i
        delta_py = py1 + py2 + py3 - py_i
        delta_pz = pz1 + pz2 + pz3 - pz_i
        delta_p_abs = np.sqrt(delta_px ** 2 + delta_py ** 2 + delta_pz ** 2)

        report.check_residual("direct Delta_E_lab", delta_e, tol, "e_alpha1/2/3 + theta/phi")
        report.check_residual("direct Delta_px_lab", delta_px, tol, "e_alpha1/2/3 + theta/phi")
        report.check_residual("direct Delta_py_lab", delta_py, tol, "e_alpha1/2/3 + theta/phi")
        report.check_residual("direct Delta_pz_lab", delta_pz, tol, "e_alpha1/2/3 + theta/phi")
        report.check_residual("direct |Delta_p_lab|", delta_p_abs, tol, "e_alpha1/2/3 + theta/phi")

        beta_x, beta_y, beta_z = c.boost_vector(e_initial, px_i, py_i, pz_i)
        e1c, px1c, py1c, pz1c = c.boost(e1, px1, py1, pz1, -beta_x, -beta_y, -beta_z)
        e2c, px2c, py2c, pz2c = c.boost(e2, px2, py2, pz2, -beta_x, -beta_y, -beta_z)
        e3c, px3c, py3c, pz3c = c.boost(e3, px3, py3, pz3, -beta_x, -beta_y, -beta_z)
        e_cm_sum = e1c + e2c + e3c
        px_cm_sum = px1c + px2c + px3c
        py_cm_sum = py1c + py2c + py3c
        pz_cm_sum = pz1c + pz2c + pz3c
        p_cm_abs = np.sqrt(px_cm_sum ** 2 + py_cm_sum ** 2 + pz_cm_sum ** 2)

        report.check_residual("direct Delta_E_cm", e_cm_sum - sqrt_s, tol, "lab vectors boosted to CM")
        report.check_residual("direct |p_cm_total|", p_cm_abs, tol, "lab vectors boosted to CM")

        t1, t2, t3 = e1c - m_alpha, e2c - m_alpha, e3c - m_alpha
        sum_t = t1 + t2 + t3
        t_available = sqrt_s - 3.0 * m_alpha
        energy_residual = sum_t - t_available
        report.check_residual("direct sum_T_alpha_cm_minus_available", energy_residual, tol, "boosted lab vectors")
        report.check_residual("direct stored Tcm alpha1 vs reconstructed", d["e_3alpha_cm_alpha1"] - t1, tol, "e_3alpha_cm_alpha1")
        report.check_residual("direct stored Tcm alpha2 vs reconstructed", d["e_3alpha_cm_alpha2"] - t2, tol, "e_3alpha_cm_alpha2")
        report.check_residual("direct stored Tcm alpha3 vs reconstructed", d["e_3alpha_cm_alpha3"] - t3, tol, "e_3alpha_cm_alpha3")

        min_t = np.minimum(np.minimum(t1, t2), t3)
        report.check_residual("direct negative CM kinetic energy floor", np.minimum(min_t - KINETIC_FLOOR_MEV, 0.0), tol, "reconstructed CM alpha kinetic energies")

        m12 = c.invariant_mass(e1c + e2c, px1c + px2c, py1c + py2c, pz1c + pz2c)
        m13 = c.invariant_mass(e1c + e3c, px1c + px3c, py1c + py3c, pz1c + pz3c)
        m23 = c.invariant_mass(e2c + e3c, px2c + px3c, py2c + py3c, pz2c + pz3c)
        lower = 2.0 * m_alpha
        upper = sqrt_s - m_alpha
        pair_low_violation = np.maximum.reduce([lower - m12, lower - m13, lower - m23, np.zeros(n)])
        pair_high_violation = np.maximum.reduce([m12 - upper, m13 - upper, m23 - upper, np.zeros(n)])
        report.check_residual("direct pair invariant mass lower-bound violation", pair_low_violation, tol, "alpha-pair invariant masses")
        report.check_residual("direct pair invariant mass upper-bound violation", pair_high_violation, tol, "alpha-pair invariant masses")

        eps1 = t1 / sum_t
        eps2 = t2 / sum_t
        eps3 = t3 / sum_t
        eps_sum = eps1 + eps2 + eps3
        x = np.sqrt(3.0) * (eps2 - eps3)
        y = 2.0 * eps1 - eps2 - eps3
        eps_violation = np.maximum.reduce([-eps1, -eps2, -eps3, eps1 - 1.0, eps2 - 1.0, eps3 - 1.0, np.zeros(n)])
        report.check_residual("direct Dalitz epsilon sum minus 1", eps_sum - 1.0, tol, "reconstructed CM kinetic energies")
        report.check_residual("direct Dalitz epsilon physical-bound violation", eps_violation, tol, "reconstructed Dalitz epsilons")
        report.note("  direct Dalitz variables are reconstructed from epsilon_i = T_i_cm/sum(T_cm); no dalitz_x/dalitz_y ROOT branches are present.")

        if "resonance_id" in d:
            unique, counts = np.unique(d["resonance_id"], return_counts=True)
            report.note("  direct resonance_id counts: " + ", ".join(f"{int(u)}:{int(cn)}" for u, cn in zip(unique, counts)))

        maybe_plot(outdir, x, y, energy_residual, delta_p_abs)

    print(report.summary_line())
    if args.report:
        report.write()
    return 0 if report.status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
