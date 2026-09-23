#!/usr/bin/env python3
"""Gamma-capture kinematics validation (task doc section 8).

The C++ generator uses recoil-aware two-body gamma kinematics for the primary
transition and explicitly generates the recoil 12C ion.  For cascade branches,
the secondary gamma is generated in the recoiling intermediate 12C rest frame.
This script checks those assumptions against ROOT diagnostics.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

DEFAULT_TOLERANCE_MEV = 1.0e-6
# Stored theta values are recomputed from boosted vectors through acos; the
# 675-keV cascade check can accumulate O(1e-6 rad) roundoff while energy closure
# remains at O(1e-13 MeV).  Keep the angle tolerance tight but above that floor.
DEFAULT_ANGLE_TOLERANCE = 2.0e-6

BRANCHES = [
    "reaction_channel", "resonance_id",
    "gamma_resonance", "gamma_branch", "gamma_angular_mode",
    "n_prompt_gammas", "gamma1_energy", "gamma2_energy",
    "gamma1_theta_lab", "gamma2_theta_lab", "gamma1_phi_lab", "gamma2_phi_lab",
    "gamma1_theta_cm", "gamma2_theta_cm", "cos_theta_gamma_cm",
    "gamma_final_state_energy_MeV", "gamma_primary_energy_MeV",
    "gamma_relative_intensity_used", "gamma_branch_fraction_used",
    "gamma_branch_is_upper_limit", "gamma_branch_from_relative_table",
    "gamma_cascade_generated",
    "projectile_kinetic_lab", "projectile_px_lab", "projectile_py_lab", "projectile_pz_lab", "projectile_p_lab",
    "e_cm_p11B",
]


def photon_vector(energy, theta, phi):
    px = energy * np.sin(theta) * np.cos(phi)
    py = energy * np.sin(theta) * np.sin(phi)
    pz = energy * np.cos(theta)
    return energy, px, py, pz


def angle_residual(theta_a, theta_b):
    return np.arctan2(np.sin(theta_a - theta_b), np.cos(theta_a - theta_b))


def main() -> int:
    parser = c.common_argparser(__doc__)
    args = parser.parse_args()

    tree = c.open_tree(Path(args.root))
    data, missing = c.load_branches(tree, BRANCHES, args.max_events)
    if missing:
        raise SystemExit(f"required branches missing: {missing}")

    masses_pm = c.derive_projectile_masses(data)
    m_p, m_11b = masses_pm["m_p"][0], masses_pm["m_11B"][0]
    tol = args.tolerance if args.tolerance is not None else DEFAULT_TOLERANCE_MEV
    report = c.Report(Path(args.report) if args.report else Path(args.root))

    channel_arg = c.parse_channel_arg(args.channel)
    channels = [3] if channel_arg is None else [channel_arg]

    for channel in channels:
        if channel != 3:
            report.emit("WARN", f"channel {channel} is not gamma capture; skipped")
            continue

        mask = data["reaction_channel"] == 3
        n = int(np.sum(mask))
        if n == 0:
            report.emit("WARN", "channel 3: no events in this file")
            continue

        d = {k: v[mask] for k, v in data.items()}
        report.note(f"-- reaction_channel=3: {n} events --")
        report.note("  source implementation: recoil-aware two-body primary gamma; explicit recoil 12C; cascade gamma generated in intermediate 12C rest frame.")

        sqrt_s = c.derive_sqrt_s(d, m_p, m_11b)
        e_initial = d["projectile_kinetic_lab"] + m_p + m_11b
        px_i, py_i, pz_i = d["projectile_px_lab"], d["projectile_py_lab"], d["projectile_pz_lab"]
        beta_x, beta_y, beta_z = c.boost_vector(e_initial, px_i, py_i, pz_i)

        # Primary gamma: reconstruct CM vector by boosting the stored lab vector.
        eg1, g1x, g1y, g1z = photon_vector(d["gamma1_energy"], d["gamma1_theta_lab"], d["gamma1_phi_lab"])
        eg1_cm, g1x_cm, g1y_cm, g1z_cm = c.boost(eg1, g1x, g1y, g1z, -beta_x, -beta_y, -beta_z)
        g1p_cm = np.sqrt(g1x_cm ** 2 + g1y_cm ** 2 + g1z_cm ** 2)
        theta1_cm = np.arccos(np.clip(g1z_cm / np.where(g1p_cm > 0, g1p_cm, np.nan), -1.0, 1.0))

        report.check_residual("gamma primary E_cm reconstructed vs stored", eg1_cm - d["gamma_primary_energy_MeV"], tol, "gamma1_energy/theta_lab/phi_lab")
        report.check_residual("gamma primary theta_cm reconstructed vs stored", angle_residual(theta1_cm, d["gamma1_theta_cm"]), DEFAULT_ANGLE_TOLERANCE, "gamma1_theta_cm")

        # The stored primary energy is recoil-aware:
        #   E_gamma = (M_i^2 - M_f^2) / (2*M_i)
        # so M_f can be reconstructed per event.  M_f - excitation should be
        # a single effective 12C ground mass across branches.
        m_final_reco = np.sqrt(np.clip(sqrt_s ** 2 - 2.0 * sqrt_s * d["gamma_primary_energy_MeV"], 0.0, None))
        m_c12_ground_reco = m_final_reco - d["gamma_final_state_energy_MeV"]
        ground_median = float(np.median(m_c12_ground_reco))
        primary_expected = (sqrt_s ** 2 - (ground_median + d["gamma_final_state_energy_MeV"]) ** 2) / (2.0 * sqrt_s)
        report.note(f"  reconstructed m_12C_ground median = {ground_median:.9f} MeV")
        report.check_residual("gamma recoil-aware primary energy formula", d["gamma_primary_energy_MeV"] - primary_expected, tol, "gamma_primary_energy_MeV")
        report.check_residual("gamma M_final - excitation constant", m_c12_ground_reco - ground_median, tol, "gamma_final_state_energy_MeV")

        cos_mask = np.isfinite(d["cos_theta_gamma_cm"])
        if np.any(cos_mask):
            beam_p = np.sqrt(px_i ** 2 + py_i ** 2 + pz_i ** 2)
            cos_reco = (g1x_cm * px_i + g1y_cm * py_i + g1z_cm * pz_i) / np.where(g1p_cm * beam_p > 0, g1p_cm * beam_p, np.nan)
            report.check_residual("gamma cos_theta_gamma_cm reconstructed vs stored", cos_reco[cos_mask] - d["cos_theta_gamma_cm"][cos_mask], DEFAULT_ANGLE_TOLERANCE, "cos_theta_gamma_cm")

        cascade_mask = d["gamma_cascade_generated"] == 1
        if np.any(cascade_mask):
            m = cascade_mask
            eg2, g2x, g2y, g2z = photon_vector(d["gamma2_energy"][m], d["gamma2_theta_lab"][m], d["gamma2_phi_lab"][m])
            eg2_cm, g2x_cm, g2y_cm, g2z_cm = c.boost(eg2, g2x, g2y, g2z, -beta_x[m], -beta_y[m], -beta_z[m])
            g2p_cm = np.sqrt(g2x_cm ** 2 + g2y_cm ** 2 + g2z_cm ** 2)
            theta2_cm = np.arccos(np.clip(g2z_cm / np.where(g2p_cm > 0, g2p_cm, np.nan), -1.0, 1.0))
            report.check_residual("gamma cascade theta2_cm reconstructed vs stored", angle_residual(theta2_cm, d["gamma2_theta_cm"][m]), DEFAULT_ANGLE_TOLERANCE, "gamma2_theta_cm")

            # Intermediate recoil after the first primary gamma in the initial CM.
            m_intermediate = m_final_reco[m]
            e_int_cm = sqrt_s[m] - eg1_cm[m]
            px_int_cm = -g1x_cm[m]
            py_int_cm = -g1y_cm[m]
            pz_int_cm = -g1z_cm[m]
            beta_ix, beta_iy, beta_iz = c.boost_vector(e_int_cm, px_int_cm, py_int_cm, pz_int_cm)
            eg2_int, _, _, _ = c.boost(eg2_cm, g2x_cm, g2y_cm, g2z_cm, -beta_ix, -beta_iy, -beta_iz)
            m_ground = m_intermediate - d["gamma_final_state_energy_MeV"][m]
            eg2_expected = (m_intermediate ** 2 - m_ground ** 2) / (2.0 * m_intermediate)
            report.check_residual("gamma cascade E2 in intermediate rest frame", eg2_int - eg2_expected, tol, "gamma2_energy/theta_lab/phi_lab")

        noncascade_mask = d["gamma_cascade_generated"] == 0
        if np.any(noncascade_mask):
            report.check_residual("gamma non-cascade gamma2_energy is zero", d["gamma2_energy"][noncascade_mask], tol, "gamma2_energy")

        for name in ("gamma_resonance", "gamma_branch", "n_prompt_gammas"):
            unique, counts = np.unique(d[name], return_counts=True)
            report.note(f"  {name} counts: " + ", ".join(f"{int(u)}:{int(cn)}" for u, cn in zip(unique, counts)))

    print(report.summary_line())
    if args.report:
        report.write()
    return 0 if report.status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
