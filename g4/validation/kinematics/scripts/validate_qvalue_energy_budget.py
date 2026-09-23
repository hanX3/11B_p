#!/usr/bin/env python3
"""Q-value and total-energy-budget check (task doc section 5).

Verifies, for 3-alpha channel events:
    Q_3alpha = m_p + m_11B - 3*m_alpha
    T_available_cm = sqrt_s - 3*m_alpha == e_cm_p11B + Q_3alpha
    T_alpha1_cm + T_alpha2_cm + T_alpha3_cm ~= T_available_cm

Usage:
    python3 validation/kinematics/scripts/validate_qvalue_energy_budget.py \
        --root data/reaction_merged_xxx.root --channel all
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
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "projectile_kinetic_lab", "projectile_p_lab", "e_cm_p11B",
]


def main() -> int:
    parser = c.common_argparser(__doc__)
    args = parser.parse_args()

    tree = c.open_tree(Path(args.root))
    data, missing = c.load_branches(tree, BRANCHES, args.max_events)
    if missing:
        raise SystemExit(f"required branches missing: {missing}")

    masses_pm = c.derive_projectile_masses(data)
    m_p, m_11b = masses_pm["m_p"][0], masses_pm["m_11B"][0]

    alpha_mask = np.isin(data["reaction_channel"], [0, 1, 2])
    three_alpha_data = {k: v[alpha_mask] for k, v in data.items()}
    if three_alpha_data["reaction_channel"].size == 0:
        print("WARN: no 3-alpha-channel events in this file; cannot derive m_alpha or run Q-value checks")
        report = c.Report(Path(args.report) if args.report else Path(args.root))
        report.emit("WARN", "no 3-alpha channel events in this file")
        print(report.summary_line())
        if args.report:
            report.write()
        return 0

    sqrt_s = c.derive_sqrt_s(three_alpha_data, m_p, m_11b)
    m_alpha, m_alpha_spread = c.derive_alpha_mass(
        sqrt_s,
        three_alpha_data["e_3alpha_cm_alpha1"],
        three_alpha_data["e_3alpha_cm_alpha2"],
        three_alpha_data["e_3alpha_cm_alpha3"],
    )

    masses = {"m_p": (m_p, masses_pm["m_p"][1]), "m_11B": (m_11b, masses_pm["m_11B"][1]), "m_alpha": (m_alpha, m_alpha_spread)}
    c.print_mass_table(masses)

    q_3alpha = m_p + m_11b - 3.0 * m_alpha
    print(f"Q_3alpha = m_p + m_11B - 3*m_alpha = {q_3alpha:.6f} MeV")

    tol = args.tolerance if args.tolerance is not None else DEFAULT_TOLERANCE_MEV
    report = c.Report(Path(args.report) if args.report else Path(args.root))
    report.note(f"mass m_p = {m_p:.9f} MeV")
    report.note(f"mass m_11B = {m_11b:.9f} MeV")
    report.note(f"mass m_alpha = {m_alpha:.9f} MeV")
    report.note(f"Q_3alpha = {q_3alpha:.9f} MeV")

    channel_arg = c.parse_channel_arg(args.channel)
    channels = [0, 1, 2] if channel_arg is None else [channel_arg]

    for channel in channels:
        if channel not in (0, 1, 2):
            report.emit("WARN", f"channel {channel} is not a 3-alpha channel; Q-value budget not applicable")
            continue
        mask = three_alpha_data["reaction_channel"] == channel
        n = int(np.sum(mask))
        if n == 0:
            report.emit("WARN", f"channel {channel}: no events in this file")
            continue

        sqrt_s_ch = sqrt_s[mask]
        e_cm_ch = three_alpha_data["e_cm_p11B"][mask]
        t_available_direct = sqrt_s_ch - 3.0 * m_alpha
        t_available_from_q = e_cm_ch + q_3alpha
        report.note(f"-- reaction_channel={channel}: {n} events --")
        report.check_residual(
            f"channel{channel} T_available_cm identity (sqrt_s - 3m_alpha == e_cm_p11B + Q_3alpha)",
            t_available_direct - t_available_from_q, tol, "e_cm_p11B, derived masses",
        )

        sum_t_cm = (
            three_alpha_data["e_3alpha_cm_alpha1"][mask]
            + three_alpha_data["e_3alpha_cm_alpha2"][mask]
            + three_alpha_data["e_3alpha_cm_alpha3"][mask]
        )
        report.check_residual(
            f"channel{channel} sum(T_alpha_cm) vs T_available_cm",
            sum_t_cm - t_available_direct, tol, "e_3alpha_cm_alpha1/2/3",
        )

    print(report.summary_line())
    if args.report:
        report.write()
    return 0 if report.status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
