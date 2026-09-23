#!/usr/bin/env python3
"""Core kinematics check (task doc section 9): for every 3-alpha-channel
event (reaction_channel 0/1/2), verify

    P_alpha1_lab + P_alpha2_lab + P_alpha3_lab == P_projectile_lab + P_target_lab

and that boosting the final state into the frame where the initial state is
at rest gives CM total energy == sqrt_s and CM total momentum == 0.

Usage:
    python3 validation/kinematics/scripts/validate_four_momentum_closure.py \
        --root data/reaction_merged_xxx.root --channel all \
        --outdir validation/kinematics/plots --report validation/kinematics/reports/x.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

DEFAULT_TOLERANCE_MEV = 1.0e-6

BRANCHES = [
    "reaction_channel",
    "e_alpha1", "e_alpha2", "e_alpha3",
    "theta_lab_alpha1", "theta_lab_alpha2", "theta_lab_alpha3",
    "phi_lab_alpha1", "phi_lab_alpha2", "phi_lab_alpha3",
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "projectile_kinetic_lab", "projectile_px_lab", "projectile_py_lab", "projectile_pz_lab", "projectile_p_lab",
    "e_cm_p11B",
]


def maybe_plot(outdir: Path | None, channel_label: str, delta_e, delta_p) -> None:
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
    ax.hist(delta_e, bins=60)
    ax.set_xlabel("Delta_E_lab (MeV)")
    ax.set_ylabel("events")
    ax.set_title(f"channel {channel_label}: lab energy closure residual")
    fig.savefig(outdir / f"channel{channel_label}_delta_E_lab.png", dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(delta_p, bins=60)
    ax.set_xlabel("|Delta_p_lab| (MeV/c)")
    ax.set_ylabel("events")
    ax.set_title(f"channel {channel_label}: lab momentum closure residual")
    fig.savefig(outdir / f"channel{channel_label}_delta_p_lab.png", dpi=120)
    plt.close(fig)


def run_for_channel(all_data: dict[str, np.ndarray], channel: int, masses: dict, report: c.Report, tol: float, outdir: Path | None) -> None:
    data = c.select_channel(all_data, channel)
    n = data["reaction_channel"].size
    if n == 0:
        report.emit("WARN", f"channel {channel}: no events in this file")
        return

    m_p, m_11b, m_alpha = masses["m_p"][0], masses["m_11B"][0], masses["m_alpha"][0]

    e1, px1, py1, pz1 = c.spherical_lab_vector(data["e_alpha1"], data["theta_lab_alpha1"], data["phi_lab_alpha1"], m_alpha)
    e2, px2, py2, pz2 = c.spherical_lab_vector(data["e_alpha2"], data["theta_lab_alpha2"], data["phi_lab_alpha2"], m_alpha)
    e3, px3, py3, pz3 = c.spherical_lab_vector(data["e_alpha3"], data["theta_lab_alpha3"], data["phi_lab_alpha3"], m_alpha)

    e_final = e1 + e2 + e3
    px_final = px1 + px2 + px3
    py_final = py1 + py2 + py3
    pz_final = pz1 + pz2 + pz3

    e_initial = data["projectile_kinetic_lab"] + m_p + m_11b
    px_initial = data["projectile_px_lab"]
    py_initial = data["projectile_py_lab"]
    pz_initial = data["projectile_pz_lab"]

    delta_e = e_final - e_initial
    delta_px = px_final - px_initial
    delta_py = py_final - py_initial
    delta_pz = pz_final - pz_initial
    delta_p_abs = np.sqrt(delta_px ** 2 + delta_py ** 2 + delta_pz ** 2)

    report.note(f"-- reaction_channel={channel}: {n} events --")
    report.check_residual(f"channel{channel} Delta_E_lab", delta_e, tol, "e_alpha1/2/3 + theta/phi")
    report.check_residual(f"channel{channel} Delta_px_lab", delta_px, tol, "e_alpha1/2/3 + theta/phi")
    report.check_residual(f"channel{channel} Delta_py_lab", delta_py, tol, "e_alpha1/2/3 + theta/phi")
    report.check_residual(f"channel{channel} Delta_pz_lab", delta_pz, tol, "e_alpha1/2/3 + theta/phi")
    report.check_residual(f"channel{channel} Delta_p_abs_lab", delta_p_abs, tol, "e_alpha1/2/3 + theta/phi")

    # CM frame: boost the final-state alphas by minus the initial-state boost
    # vector (CLHEP convention: boost(v) moves v INTO the frame comoving with
    # v's argument, so boosting by -beta_initial brings lab 4-vectors into
    # the frame where the initial state is at rest).
    beta_x, beta_y, beta_z = c.boost_vector(e_initial, px_initial, py_initial, pz_initial)
    e1c, px1c, py1c, pz1c = c.boost(e1, px1, py1, pz1, -beta_x, -beta_y, -beta_z)
    e2c, px2c, py2c, pz2c = c.boost(e2, px2, py2, pz2, -beta_x, -beta_y, -beta_z)
    e3c, px3c, py3c, pz3c = c.boost(e3, px3, py3, pz3, -beta_x, -beta_y, -beta_z)

    sqrt_s = c.derive_sqrt_s(data, m_p, m_11b)
    e_cm_final = e1c + e2c + e3c
    px_cm_final = px1c + px2c + px3c
    py_cm_final = py1c + py2c + py3c
    pz_cm_final = pz1c + pz2c + pz3c
    p_cm_abs = np.sqrt(px_cm_final ** 2 + py_cm_final ** 2 + pz_cm_final ** 2)

    delta_e_cm = e_cm_final - sqrt_s
    report.check_residual(f"channel{channel} Delta_E_cm (vs sqrt_s)", delta_e_cm, tol, "e_cm_p11B + derived m_p/m_11B")
    report.check_residual(f"channel{channel} |p_cm_total| (should be 0)", p_cm_abs, tol, "e_cm_p11B + derived m_p/m_11B")

    t_available_cm = sqrt_s - 3.0 * m_alpha
    sum_t_cm = (data["e_3alpha_cm_alpha1"] + data["e_3alpha_cm_alpha2"] + data["e_3alpha_cm_alpha3"])
    residual_energy_budget = sum_t_cm - t_available_cm
    report.check_residual(
        f"channel{channel} sum_T_alpha_cm_minus_available",
        residual_energy_budget, tol, "e_3alpha_cm_alpha1/2/3 vs sqrt_s - 3*m_alpha",
    )

    maybe_plot(outdir, str(channel), delta_e, delta_p_abs)


def main() -> int:
    parser = c.common_argparser(__doc__)
    args = parser.parse_args()

    tree = c.open_tree(Path(args.root))
    data, missing = c.load_branches(tree, BRANCHES, args.max_events)
    if missing:
        raise SystemExit(f"required branches missing: {missing}")

    m_p, m_11b = c.derive_projectile_masses(data)["m_p"], c.derive_projectile_masses(data)["m_11B"]
    three_alpha = c.select_channel(data, None)
    alpha_mask = np.isin(data["reaction_channel"], [0, 1, 2])
    three_alpha_data = {k: v[alpha_mask] for k, v in data.items()}
    sqrt_s_for_mass = c.derive_sqrt_s(three_alpha_data, m_p[0], m_11b[0])
    m_alpha = c.derive_alpha_mass(
        sqrt_s_for_mass,
        three_alpha_data["e_3alpha_cm_alpha1"],
        three_alpha_data["e_3alpha_cm_alpha2"],
        three_alpha_data["e_3alpha_cm_alpha3"],
    )
    masses = {"m_p": m_p, "m_11B": m_11b, "m_alpha": m_alpha}

    report_path = Path(args.report) if args.report else None
    report = c.Report(report_path if report_path else Path(args.root))
    c.print_mass_table(masses)
    for line_key, (value, spread) in masses.items():
        report.note(f"mass {line_key} = {value:.9f} MeV (mad={spread:.3e} MeV)")

    tol = args.tolerance if args.tolerance is not None else DEFAULT_TOLERANCE_MEV
    outdir = Path(args.outdir) if args.outdir else None
    channel_arg = c.parse_channel_arg(args.channel)
    channels = [0, 1, 2] if channel_arg is None else [channel_arg]

    for channel in channels:
        if channel not in (0, 1, 2):
            report.emit("WARN", f"channel {channel} is not a 3-alpha channel; four-momentum closure not applicable here (see validate_gamma_capture_kinematics.py)")
            continue
        run_for_channel(data, channel, masses, report, tol, outdir)

    print(report.summary_line())
    if report_path:
        report.write()
    return 0 if report.status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
