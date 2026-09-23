#!/usr/bin/env python3
"""Sequential 165/675-keV alpha-decay kinematics check (task doc section 6).

Applies to reaction_channel 0 (165 sequential) and 1 (675 sequential,
including the SymmetrizedCoherentL1L3 "strict" model, which src/H11BReaction.cc
also tags as reaction_channel==1 -- see h11b675_alpha_decay_model). By
construction in both code paths, alpha1 is the primary breakup alpha
(12C* -> alpha1 + 8Be*) and alpha2/alpha3 are the pair from 8Be* -> alpha+alpha,
so:

    invariant_mass(alpha2, alpha3) == 2*m_alpha + eaa_8Be

is checked directly against the stored eaa_8Be branch, along with the
cos_theta_primary_cm and cos_chi_secondary_8be helicity-angle branches
reconstructed independently from the lab-frame alpha kinematics.

Usage:
    python3 validation/kinematics/scripts/validate_sequential_decay_kinematics.py \
        --root data/reaction_merged_xxx.root --channel all
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

DEFAULT_TOLERANCE_MEV = 1.0e-6
DEFAULT_COS_TOLERANCE = 1.0e-6
EX8BE_GROUND_ABOVE_2ALPHA_KEV = 91.84  # include/Constants.hh Ex8BeGroundAbove2Alpha

BRANCHES = [
    "reaction_channel", "resonance_id", "branch_id", "h11b675_alpha_decay_model",
    "e_alpha1", "e_alpha2", "e_alpha3",
    "theta_lab_alpha1", "theta_lab_alpha2", "theta_lab_alpha3",
    "phi_lab_alpha1", "phi_lab_alpha2", "phi_lab_alpha3",
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "eaa_8Be", "ex_8Be",
    "cos_theta_primary_cm", "cos_chi_secondary_8be",
    "projectile_kinetic_lab", "projectile_px_lab", "projectile_py_lab", "projectile_pz_lab", "projectile_p_lab",
    "e_cm_p11B",
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

    tol = args.tolerance if args.tolerance is not None else DEFAULT_TOLERANCE_MEV
    report = c.Report(Path(args.report) if args.report else Path(args.root))

    channel_arg = c.parse_channel_arg(args.channel)
    channels = [0, 1] if channel_arg is None else [channel_arg]

    for channel in channels:
        if channel not in (0, 1):
            report.emit("WARN", f"channel {channel} is not a sequential-decay channel (0/1); skipped")
            continue

        mask = data["reaction_channel"] == channel
        n = int(np.sum(mask))
        if n == 0:
            report.emit("WARN", f"channel {channel}: no events in this file")
            continue

        d = {k: v[mask] for k, v in data.items()}
        report.note(f"-- reaction_channel={channel}: {n} events --")

        strict_frac = float(np.mean(d["h11b675_alpha_decay_model"] == 1)) if channel == 1 else 0.0
        branch0_frac = float(np.mean(d["branch_id"] == 0))
        report.note(
            f"  branch_id==0 (alpha0/ground-8Be) fraction = {branch0_frac:.4f}; "
            f"h11b675_alpha_decay_model==1 (SymmetrizedCoherentL1L3) fraction = {strict_frac:.4f}"
        )

        m_alpha = c.derive_alpha_mass(
            c.derive_sqrt_s(d, m_p, m_11b), d["e_3alpha_cm_alpha1"], d["e_3alpha_cm_alpha2"], d["e_3alpha_cm_alpha3"]
        )[0]

        e1, px1, py1, pz1 = c.spherical_lab_vector(d["e_alpha1"], d["theta_lab_alpha1"], d["phi_lab_alpha1"], m_alpha)
        e2, px2, py2, pz2 = c.spherical_lab_vector(d["e_alpha2"], d["theta_lab_alpha2"], d["phi_lab_alpha2"], m_alpha)
        e3, px3, py3, pz3 = c.spherical_lab_vector(d["e_alpha3"], d["theta_lab_alpha3"], d["phi_lab_alpha3"], m_alpha)

        # 1-3: four-momentum closure and CM totals (already covered in depth
        # by validate_four_momentum_closure.py; repeated compactly here since
        # section 6.3 lists it as part of this script's core judgment too).
        e_initial = d["projectile_kinetic_lab"] + m_p + m_11b
        px_i, py_i, pz_i = d["projectile_px_lab"], d["projectile_py_lab"], d["projectile_pz_lab"]
        delta_e = (e1 + e2 + e3) - e_initial
        delta_p = np.sqrt((px1 + px2 + px3 - px_i) ** 2 + (py1 + py2 + py3 - py_i) ** 2 + (pz1 + pz2 + pz3 - pz_i) ** 2)
        report.check_residual(f"channel{channel} lab four-momentum closure |Delta_p|", delta_p, tol, "e_alpha1/2/3+theta+phi")
        report.check_residual(f"channel{channel} lab energy closure Delta_E", delta_e, tol, "e_alpha1/2/3")

        # 4-5: eaa_8Be from the alpha2/alpha3 pair invariant mass (frame-independent).
        # eaa_8Be is *defined* in src/H11BReaction.cc as
        # ex_8Be + Ex8BeGroundAbove2Alpha, a hardcoded literature constant
        # (91.84 keV), not as invariant_mass(pair) - 2*m_alpha using
        # G4IonTable's own 8He/8Be masses. The two need not agree exactly:
        # any gap between them is the actual difference between G4's ion-table
        # 8Be mass and "2*m_alpha_iontable + 91.84 keV". Check that this gap
        # is a *pure additive constant* (i.e. the eaa_8Be formula still tracks
        # ex_8Be correctly and no random/model-dependent error is present),
        # and separately report its size.
        pair_mass = c.invariant_mass(e2 + e3, px2 + px3, py2 + py3, pz2 + pz3)
        eaa_reconstructed = pair_mass - 2.0 * m_alpha
        raw_offset = eaa_reconstructed - d["eaa_8Be"]
        offset_median = float(np.median(raw_offset))
        report.note(
            f"  channel{channel} eaa_8Be vs invariant-mass offset: median={offset_median * 1000:.4f} keV "
            f"(constant gap between G4IonTable's 8Be mass and the hardcoded Ex8BeGroundAbove2Alpha=91.84 keV "
            f"constant in src/Constants.hh; physically negligible next to the 1.513 MeV 8Be(2+) width, "
            f"but see gpt/codex_reaction_kinematics_validation.md section 15 point 7)"
        )
        report.check_residual(
            f"channel{channel} eaa_8Be formula self-consistency (offset-subtracted)",
            raw_offset - offset_median, tol, "eaa_8Be vs alpha2+alpha3 invariant mass",
        )
        if abs(offset_median) >= 1.0e-3:
            report.emit("WARN", f"channel{channel} eaa_8Be absolute offset from invariant mass = {offset_median * 1000:.4f} keV (>= 1 keV)")
        # Ground-state (branch_id==0) events must sit at the 8Be(g.s.) pole
        # (this is tautological from the eaa_8Be formula itself; the physically
        # meaningful version of this check is the invariant-mass reconstruction above).
        ground_mask = d["branch_id"] == 0
        if np.any(ground_mask):
            report.check_residual(
                f"channel{channel} branch_id==0 eaa_8Be == Ex8BeGroundAbove2Alpha ({EX8BE_GROUND_ABOVE_2ALPHA_KEV} keV, by formula)",
                (d["eaa_8Be"][ground_mask] - EX8BE_GROUND_ABOVE_2ALPHA_KEV / 1000.0), tol, "eaa_8Be, branch_id",
            )

        # 6: primary/secondary role is defined structurally by construction
        # (alpha1 = 12C*->alpha+8Be* primary; alpha2/alpha3 = 8Be* pair); the
        # eaa_8Be check above is exactly the test that alpha2/alpha3 form the
        # 8Be-consistent pair, so no separate check is needed here.

        # 7: cos_theta_primary_cm = cos(angle between alpha1's CM direction and beam axis).
        # The beam is along a fixed lab direction and the initial-state boost is
        # purely along that axis, so the beam unit vector is identical in lab
        # and CM frames -- no extra boost of the axis itself is needed.
        e_total_i = e1 + e2 + e3
        beta_x, beta_y, beta_z = c.boost_vector(e_total_i, px1 + px2 + px3, py1 + py2 + py3, pz1 + pz2 + pz3)
        e1c, px1c, py1c, pz1c = c.boost(e1, px1, py1, pz1, -beta_x, -beta_y, -beta_z)
        p1c_mag = np.sqrt(px1c ** 2 + py1c ** 2 + pz1c ** 2)
        beam_p = np.sqrt(px_i ** 2 + py_i ** 2 + pz_i ** 2)
        cos_primary_reconstructed = (px1c * px_i + py1c * py_i + pz1c * pz_i) / np.where((p1c_mag * beam_p) > 0, p1c_mag * beam_p, np.nan)
        report.check_residual(
            f"channel{channel} cos_theta_primary_cm reconstructed vs stored",
            cos_primary_reconstructed - d["cos_theta_primary_cm"], DEFAULT_COS_TOLERANCE, "cos_theta_primary_cm",
        )

        # 8: cos_chi_secondary_8be = cos(angle between the alpha2/alpha3 pair's
        # CM-frame direction and alpha2's direction in the pair's own rest frame).
        #
        # cos_chi is a rotation-sensitive quantity, not a boost-path-independent
        # scalar like mass or energy, so the reconstruction must mirror the
        # *exact* sequence of boosts src/H11BReaction.cc applies, or it picks
        # up a spurious Thomas-Wigner rotation relative to the stored value:
        #   - regular sequential path (h11b675_alpha_decay_model != 1, i.e. all
        #     165 events and legacy 675 events): the code boosts alpha2/alpha3
        #     directly between the "8Be rest frame" and the lab with ONE boost
        #     vector each way (lv_lab_alphaN.boost(lv_lab_8Be.boostVector())),
        #     so the mirror reconstruction is a single direct lab<->pair-rest
        #     boost.
        #   - SymmetrizedCoherentL1L3 strict model (h11b675_alpha_decay_model
        #     == 1): the code boosts alpha2/alpha3 pair-rest -> overall CM ->
        #     lab as two separate .boost() calls, and cos_chi_secondary_8be is
        #     computed by undoing only the pair-rest -> overall-CM leg, so the
        #     mirror reconstruction must also go through the overall CM frame
        #     as an explicit intermediate step (lab -> overall CM -> pair rest).
        #
        # h11b675_alpha_decay_model is filled by FillRuntimeConfigDiagnostics()
        # from the *global* /h11b/675AlphaDecayModel setting on every event
        # regardless of channel/branch, not per-event routing info -- so it
        # must be combined with channel==1 (675) and branch_id==1
        # (is_alpha1_branch) to identify which events actually went through
        # Generate675SymmetrizedCoherent3Alpha().
        strict_mask = (channel == 1) & (d["branch_id"] == 1) & (d["h11b675_alpha_decay_model"] == 1)
        cos_chi_reconstructed = np.full(n, np.nan)

        if np.any(~strict_mask):
            m = ~strict_mask
            # dir_cm_alpha2: single direct lab -> pair-rest-frame boost.
            beta_px, beta_py, beta_pz = c.boost_vector(e2[m] + e3[m], px2[m] + px3[m], py2[m] + py3[m], pz2[m] + pz3[m])
            e2p, px2p, py2p, pz2p = c.boost(e2[m], px2[m], py2[m], pz2[m], -beta_px, -beta_py, -beta_pz)
            p2p_mag = np.sqrt(px2p ** 2 + py2p ** 2 + pz2p ** 2)
            # axis_8be_recoil_cm: the pair's direction in the *overall* CM
            # frame (single direct lab -> overall-CM boost), not its raw lab
            # direction.
            e_pair_cm, px_pair_cm, py_pair_cm, pz_pair_cm = c.boost(
                e2[m] + e3[m], px2[m] + px3[m], py2[m] + py3[m], pz2[m] + pz3[m], -beta_x[m], -beta_y[m], -beta_z[m]
            )
            pair_cm_mag = np.sqrt(px_pair_cm ** 2 + py_pair_cm ** 2 + pz_pair_cm ** 2)
            cos_chi_reconstructed[m] = (px2p * px_pair_cm + py2p * py_pair_cm + pz2p * pz_pair_cm) / np.where(
                (p2p_mag * pair_cm_mag) > 0, p2p_mag * pair_cm_mag, np.nan
            )

        if np.any(strict_mask):
            m = strict_mask
            e2c, px2c, py2c, pz2c = c.boost(e2[m], px2[m], py2[m], pz2[m], -beta_x[m], -beta_y[m], -beta_z[m])
            e3c, px3c, py3c, pz3c = c.boost(e3[m], px3[m], py3[m], pz3[m], -beta_x[m], -beta_y[m], -beta_z[m])
            e_pair_c, px_pair_c, py_pair_c, pz_pair_c = e2c + e3c, px2c + px3c, py2c + py3c, pz2c + pz3c
            pair_c_mag = np.sqrt(px_pair_c ** 2 + py_pair_c ** 2 + pz_pair_c ** 2)
            beta_px, beta_py, beta_pz = c.boost_vector(e_pair_c, px_pair_c, py_pair_c, pz_pair_c)
            e2p, px2p, py2p, pz2p = c.boost(e2c, px2c, py2c, pz2c, -beta_px, -beta_py, -beta_pz)
            p2p_mag = np.sqrt(px2p ** 2 + py2p ** 2 + pz2p ** 2)
            cos_chi_reconstructed[m] = (px2p * px_pair_c + py2p * py_pair_c + pz2p * pz_pair_c) / np.where(
                (p2p_mag * pair_c_mag) > 0, p2p_mag * pair_c_mag, np.nan
            )

        report.check_residual(
            f"channel{channel} cos_chi_secondary_8be reconstructed vs stored",
            cos_chi_reconstructed - d["cos_chi_secondary_8be"], DEFAULT_COS_TOLERANCE, "cos_chi_secondary_8be",
        )

    print(report.summary_line())
    if args.report:
        report.write()
    return 0 if report.status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
