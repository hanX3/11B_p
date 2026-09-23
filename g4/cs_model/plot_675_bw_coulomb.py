#!/usr/bin/env python3
"""Plot the 11B(p,3alpha) cross section around the 675-keV resonance.

Every literature point is reconstructed from numerical quantities printed in
the original papers:

1. Becker, Rolfs and Trautvetter, Z. Phys. A 327, 341-355 (1987), Table 2
   - published quantities: E_cm, S_alpha0(E), S_alpha1(E)
   - conversion:
         sigma(E) = [S_alpha0(E) + S_alpha1(E)]
                    / [E * exp(4.7528 / sqrt(E))]
     where E is in MeV, S is in MeV barn, and sigma is in barn.
     (4.7528/sqrt(E) is Becker's own Gamow factor 2*pi*eta for p+11B.)
   - proton laboratory energy:
         E_p,lab = (12/11) E_cm  (mass-number ratio; the exact nuclear-mass
         ratio used elsewhere in this script differs by ~0.05 %).

2. Taskaev et al., Nucl. Instrum. Methods Phys. Res. B 555 (2024)
   165490, Tables 5 and 8
   - published quantities: sigma_alpha0 and sigma_alpha1 in mb
   - plotted total:
         sigma_total = sigma_alpha0 + sigma_alpha1.

3. Sikora and Weller, J. Fusion Energ. 35, 538-543 (2016), Table 1
   - published quantity: Legendre coefficient A0 in mb/sr
   - plotted total:
         sigma_total = (4*pi/3) A0.
     See the SIKORA_MULTIPLICITY note below for why the factor is 4*pi/3
     and not the textbook 4*pi.
   - these are model-corrected experimental results because the unmeasured
     low-energy part of the alpha spectrum was completed with a reaction model.

The comparison curve, also used by the runtime C++ cross-section model, is a
single-level Breit-Wigner expression with an energy-dependent s-wave proton
width evaluated using a low-energy Coulomb penetrability approximation.

Usage
-----
Terminal:
    python3 plot_675_bw_coulomb.py --show
    python3 plot_675_bw_coulomb.py --output fig.png

Jupyter / IPython:
    from plot_675_bw_coulomb import plot_675_bw_coulomb
    fig = plot_675_bw_coulomb()                       # inline, default parameters
    fig = plot_675_bw_coulomb(gamma_p_resonance_keV=200, exclude_sikora=True)
    # or simply:  %run plot_675_bw_coulomb.py

Dependencies:
    python3 -m pip install --user numpy matplotlib
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib


# ---------------------------------------------------------------------------
# Backend selection (works in a head-less terminal AND in Jupyter)
# ---------------------------------------------------------------------------
#   * Jupyter / IPython  -> leave the inline (or user) backend untouched.
#   * Terminal with a display -> keep the interactive backend so --show works.
#   * Head-less terminal (no DISPLAY) -> force Agg so savefig still works and
#     plt.show() is a harmless no-op.
# savefig() works on every backend, so figures are always written to disk.

def _in_notebook() -> bool:
    """True only inside a Jupyter/IPython ZMQ kernel."""
    try:
        from IPython import get_ipython

        shell = get_ipython()
        return shell is not None and shell.__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


IN_NOTEBOOK = _in_notebook()

if (
    not IN_NOTEBOOK
    and os.name != "nt"
    and not os.environ.get("DISPLAY")
    and not os.environ.get("MPLBACKEND")
):
    matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (import after backend choice)
import numpy as np  # noqa: E402


# ---------------------------------------------------------------------------
# Physical constants and frame conversion
# ---------------------------------------------------------------------------

HBARC_MEV_FM = 197.3269804
AMU_MEV = 931.49410242
ALPHA_FINE_STRUCTURE = 1.0 / 137.035999084

# Real masses (AME2020).  The reduced mass that enters the wave number k and
# the Sommerfeld parameter eta is the *nuclear* reduced mass, so the atomic
# masses are corrected for the electrons.
ELECTRON_MASS_U = 5.48579909065e-4
PROTON_MASS_U = 1.007276466879                       # bare proton
B11_ATOMIC_MASS_U = 11.009305166                     # neutral 11B atom
B11_NUCLEAR_MASS_U = B11_ATOMIC_MASS_U - 5.0 * ELECTRON_MASS_U

PROJECTILE_MASS_MEV = PROTON_MASS_U * AMU_MEV
TARGET_MASS_MEV = B11_NUCLEAR_MASS_U * AMU_MEV

REDUCED_MASS_MEV = (
    PROJECTILE_MASS_MEV
    * TARGET_MASS_MEV
    / (PROJECTILE_MASS_MEV + TARGET_MASS_MEV)
)

Z_PROJECTILE = 1.0
Z_TARGET = 5.0

# Non-relativistic lab <-> c.m. energy conversion using the exact nuclear
# masses.  Numerically this is 0.91616, i.e. within 0.05 % of 11/12; the
# Becker table (which quotes E_cm) is converted with the conventional 12/11.
LAB_TO_CM = TARGET_MASS_MEV / (PROJECTILE_MASS_MEV + TARGET_MASS_MEV)
CM_TO_LAB = 1.0 / LAB_TO_CM

# Sikora-Weller multiplicity factor -----------------------------------------
#   The textbook total cross section from the 0th-order Legendre coefficient
#   A0 of dsigma/dOmega is  sigma = 4*pi*A0.  Sikora & Weller, however, report
#   their absolute cross section in terms of the TOTAL number of alpha
#   particles detected (to avoid ambiguity in how many alphas a single
#   reaction emits).  Because 11B(p,3alpha) yields THREE alphas per reaction,
#   the alpha-yield-integrated 4*pi*A0 overcounts the reaction cross section by
#   the 3-alpha multiplicity.  Dividing it out gives the reaction cross
#   section:
#         sigma_reaction = (4*pi / 3) * A0 .
#   Cross-check: with 4*pi/3 the Sikora points agree with Becker (peak
#   ~1.3-1.4 barn, the accepted absolute scale); using the bare 4*pi would
#   overshoot Becker and Taskaev by a factor of 3.
SIKORA_MULTIPLICITY = 3.0


# ---------------------------------------------------------------------------
# Becker et al. (1987), Table 2
# ---------------------------------------------------------------------------

BECKER_E_CM_KEV = np.array(
    [
        292.0,
        300.0,
        350.0,
        400.0,
        450.0,
        500.0,
        550.0,
        600.0,
        642.0,
        700.0,
        750.0,
        800.0,
        850.0,
        900.0,
    ],
    dtype=float,
)

BECKER_S_ALPHA0_MEV_B = np.array(
    [
        1.91,
        1.82,
        1.78,
        1.64,
        1.56,
        1.38,
        1.11,
        1.35,
        1.14,
        1.12,
        1.29,
        1.00,
        0.94,
        0.91,
    ],
    dtype=float,
)

BECKER_S_ALPHA0_ERR_MEV_B = np.array(
    [
        0.10,
        0.09,
        0.09,
        0.09,
        0.08,
        0.08,
        0.12,
        0.08,
        0.06,
        0.06,
        0.08,
        0.06,
        0.05,
        0.05,
    ],
    dtype=float,
)

BECKER_S_ALPHA1_MEV_B = np.array(
    [
        271.0,
        290.0,
        293.0,
        332.0,
        361.0,
        377.0,
        357.0,
        340.0,
        239.0,
        135.0,
        68.0,
        55.0,
        42.0,
        35.0,
    ],
    dtype=float,
)

BECKER_S_ALPHA1_ERR_MEV_B = np.array(
    [
        14.0,
        15.0,
        15.0,
        17.0,
        18.0,
        19.0,
        18.0,
        17.0,
        12.0,
        7.0,
        4.0,
        3.0,
        2.0,
        2.0,
    ],
    dtype=float,
)


# ---------------------------------------------------------------------------
# Taskaev et al. (2024), Tables 5 and 8
# ---------------------------------------------------------------------------

TASKAEV_E_LAB_KEV = np.array(
    [
        355.0,
        461.0,
        565.0,
        668.0,
        771.0,
        873.0,
        975.0,
    ],
    dtype=float,
)

TASKAEV_E_LAB_ERR_KEV = np.array(
    [
        45.0,
        39.0,
        35.0,
        32.0,
        29.0,
        27.0,
        25.0,
    ],
    dtype=float,
)

TASKAEV_SIGMA_ALPHA0_MB = np.array(
    [
        1.40,
        2.02,
        2.53,
        3.31,
        3.45,
        3.04,
        2.90,
    ],
    dtype=float,
)

TASKAEV_SIGMA_ALPHA0_ERR_MB = np.array(
    [
        0.24,
        0.34,
        0.42,
        0.55,
        0.57,
        0.50,
        0.49,
    ],
    dtype=float,
)

TASKAEV_SIGMA_ALPHA1_MB = np.array(
    [
        148.0,
        357.0,
        598.0,
        668.0,
        386.0,
        234.0,
        171.0,
    ],
    dtype=float,
)

TASKAEV_SIGMA_ALPHA1_ERR_MB = np.array(
    [
        19.0,
        42.0,
        88.0,
        89.0,
        58.0,
        34.0,
        22.0,
    ],
    dtype=float,
)


# ---------------------------------------------------------------------------
# Sikora and Weller (2016), Table 1
# ---------------------------------------------------------------------------

SIKORA_E_LAB_MEV = np.array(
    [
        0.30,
        0.40,
        0.49,
        0.57,
        0.65,
        0.73,
        0.80,
        0.88,
        0.94,
        1.00,
    ],
    dtype=float,
)

SIKORA_A0_MB_PER_SR = np.array(
    [
        31.4329,
        93.5936,
        173.7688,
        285.1283,
        333.7993,
        273.8339,
        172.2051,
        110.7989,
        79.4042,
        75.4242,
    ],
    dtype=float,
)

SIKORA_A0_ERR_MB_PER_SR = np.array(
    [
        0.0683,
        0.5210,
        1.0916,
        0.8436,
        0.6866,
        0.6153,
        0.4377,
        0.3744,
        0.2940,
        0.4428,
    ],
    dtype=float,
)


# ---------------------------------------------------------------------------
# Literature-data reconstruction
# ---------------------------------------------------------------------------

def becker_cross_section() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Becker proton-lab energy, total cross section, and uncertainty."""
    energy_cm_mev = BECKER_E_CM_KEV / 1000.0
    # Becker's own convention for the lab energy is the mass-number ratio 12/11.
    energy_lab_kev = BECKER_E_CM_KEV * (12.0 / 11.0)

    total_s = (
        BECKER_S_ALPHA0_MEV_B
        + BECKER_S_ALPHA1_MEV_B
    )
    total_s_error = np.hypot(
        BECKER_S_ALPHA0_ERR_MEV_B,
        BECKER_S_ALPHA1_ERR_MEV_B,
    )

    denominator = (
        energy_cm_mev
        * np.exp(
            4.7528 / np.sqrt(energy_cm_mev)
        )
    )
    sigma_barn = total_s / denominator
    sigma_error_barn = total_s_error / denominator

    return (
        energy_lab_kev,
        sigma_barn,
        sigma_error_barn,
    )


def taskaev_cross_section() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return Taskaev lab energy, x error, total cross section, and y error."""
    sigma_mb = (
        TASKAEV_SIGMA_ALPHA0_MB
        + TASKAEV_SIGMA_ALPHA1_MB
    )
    sigma_error_mb = np.hypot(
        TASKAEV_SIGMA_ALPHA0_ERR_MB,
        TASKAEV_SIGMA_ALPHA1_ERR_MB,
    )

    return (
        TASKAEV_E_LAB_KEV,
        TASKAEV_E_LAB_ERR_KEV,
        sigma_mb / 1000.0,
        sigma_error_mb / 1000.0,
    )


def sikora_cross_section() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Sikora lab energy, total cross section, and statistical error.

    sigma = (4*pi / 3) * A0 : the 4*pi integrates the isotropic (0th-order)
    Legendre term, the 1/3 removes the 3-alpha multiplicity (see the
    SIKORA_MULTIPLICITY note above).
    """
    factor = 4.0 * np.pi / SIKORA_MULTIPLICITY
    energy_lab_kev = 1000.0 * SIKORA_E_LAB_MEV
    sigma_barn = (
        factor * SIKORA_A0_MB_PER_SR / 1000.0
    )
    sigma_error_barn = (
        factor * SIKORA_A0_ERR_MB_PER_SR / 1000.0
    )

    return (
        energy_lab_kev,
        sigma_barn,
        sigma_error_barn,
    )


# ---------------------------------------------------------------------------
# Coulomb-suppressed Breit-Wigner curve
# ---------------------------------------------------------------------------

def wave_number_fm_inverse(
    energy_cm_mev: np.ndarray,
) -> np.ndarray:
    return np.sqrt(
        2.0
        * REDUCED_MASS_MEV
        * energy_cm_mev
    ) / HBARC_MEV_FM


def sommerfeld_parameter(
    energy_cm_mev: np.ndarray,
) -> np.ndarray:
    velocity_over_c = np.sqrt(
        2.0
        * energy_cm_mev
        / REDUCED_MASS_MEV
    )
    return (
        Z_PROJECTILE
        * Z_TARGET
        * ALPHA_FINE_STRUCTURE
        / velocity_over_c
    )


def s_wave_penetrability_approx(
    energy_cm_mev: np.ndarray,
    channel_radius_fm: float,
) -> np.ndarray:
    """Low-energy approximation P0 proportional to rho*C0(eta)^2.

    NOTE: this hard-wires the s-wave (l=0) proton channel appropriate to the
    675-keV (2-) resonance.  Re-using this routine for a different partial wave
    would require the corresponding P_l.
    """
    k = wave_number_fm_inverse(
        energy_cm_mev
    )
    rho = k * float(channel_radius_fm)

    eta = sommerfeld_parameter(
        energy_cm_mev
    )
    two_pi_eta = 2.0 * np.pi * eta
    c0_squared = (
        two_pi_eta
        / np.expm1(two_pi_eta)
    )
    return rho * c0_squared


def breit_wigner_cross_section_barn(
    energy_lab_kev: np.ndarray,
    resonance_energy_lab_kev: float,
    gamma_p_resonance_kev: float,
    gamma_out_kev: float,
    spin_factor: float,
    channel_radius_fm: float,
) -> np.ndarray:
    energy_lab_kev = np.asarray(
        energy_lab_kev,
        dtype=float,
    )
    energy_cm_mev = (
        energy_lab_kev
        * LAB_TO_CM
        / 1000.0
    )
    resonance_cm_mev = (
        float(resonance_energy_lab_kev)
        * LAB_TO_CM
        / 1000.0
    )

    penetrability = s_wave_penetrability_approx(
        energy_cm_mev,
        channel_radius_fm,
    )
    penetrability_resonance = float(
        s_wave_penetrability_approx(
            np.array(
                [resonance_cm_mev],
                dtype=float,
            ),
            channel_radius_fm,
        )[0]
    )

    gamma_p_kev = (
        float(gamma_p_resonance_kev)
        * penetrability
        / penetrability_resonance
    )
    gamma_total_kev = (
        gamma_p_kev
        + float(gamma_out_kev)
    )

    energy_offset_kev = (
        energy_cm_mev
        - resonance_cm_mev
    ) * 1000.0

    denominator = (
        energy_offset_kev**2
        + 0.25 * gamma_total_kev**2
    )

    k = wave_number_fm_inverse(
        energy_cm_mev
    )
    cross_section_fm2 = (
        np.pi
        / k**2
        * float(spin_factor)
        * gamma_p_kev
        * float(gamma_out_kev)
        / denominator
    )

    # 1 barn = 100 fm^2.
    return cross_section_fm2 / 100.0


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 14,
            "legend.fontsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 1.1,
            "lines.linewidth": 1.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot exact literature-derived p+11B cross-section points "
            "around the 675-keV resonance."
        )
    )
    parser.add_argument("--output", type=Path, default=Path("p11b_675keV_exact_literature.png"))
    parser.add_argument("--resonance-energy-lab-keV", type=float, default=675.0)
    parser.add_argument("--gamma-p-resonance-keV", type=float, default=150.0)
    parser.add_argument("--gamma-out-keV", type=float, default=150.0)
    parser.add_argument("--spin-factor", type=float, default=5.0 / 8.0)
    parser.add_argument("--channel-radius-fm", type=float, default=4.5)
    parser.add_argument("--energy-min-keV", type=float, default=300.0)
    parser.add_argument("--energy-max-keV", type=float, default=1000.0)
    parser.add_argument("--ymax-barn", type=float, default=1.48)
    parser.add_argument("--no-title", action="store_true")
    parser.add_argument(
        "--exclude-sikora",
        action="store_true",
        help="omit the model-corrected Sikora-Weller points",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="skip writing the figure/CSV to disk (useful in notebooks)",
    )
    parser.add_argument("--show", action="store_true")

    # In a notebook, sys.argv carries the kernel's own flags (e.g. -f ...).
    # Fall back to defaults there and ignore any stray/unknown flags.
    if argv is None:
        argv = [] if IN_NOTEBOOK else sys.argv[1:]
    args, _unknown = parser.parse_known_args(argv)
    return args


def validate_args(
    args: argparse.Namespace,
) -> None:
    if args.energy_max_keV <= args.energy_min_keV:
        raise SystemExit("--energy-max-keV must exceed --energy-min-keV")
    if args.ymax_barn <= 0.0:
        raise SystemExit("--ymax-barn must be positive")
    if args.resonance_energy_lab_keV <= 0.0:
        raise SystemExit("--resonance-energy-lab-keV must be positive")
    if args.gamma_p_resonance_keV <= 0.0 or args.gamma_out_keV <= 0.0:
        raise SystemExit("Breit-Wigner widths must be positive")
    if args.spin_factor <= 0.0:
        raise SystemExit("--spin-factor must be positive")
    if args.channel_radius_fm <= 0.0:
        raise SystemExit("--channel-radius-fm must be positive")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def plot_675_bw_coulomb(
    argv: list[str] | None = None,
    **overrides,
) -> "plt.Figure":
    """Build the 675-keV cross-section figure and return it.

    Terminal:
        python3 plot_675_bw_coulomb.py --show --output fig.png

    Jupyter / IPython (returns the Figure, shown inline):
        fig = plot_675_bw_coulomb()
        fig = plot_675_bw_coulomb(gamma_p_resonance_keV=200, exclude_sikora=True)

    Any command-line option is available as a keyword (dashes -> underscores),
    e.g. --gamma-p-resonance-keV becomes gamma_p_resonance_keV. By default the
    figure is written to disk; pass no_save=True to skip that.
    """
    args = parse_args(argv)
    for key, value in overrides.items():
        if not hasattr(args, key):
            raise TypeError(f"unknown option: {key!r}")
        setattr(args, key, Path(value) if key == "output" else value)
    validate_args(args)

    configure_plot_style()

    energy_lab_kev = np.linspace(
        args.energy_min_keV,
        args.energy_max_keV,
        1401,
    )
    bw_sigma_barn = breit_wigner_cross_section_barn(
        energy_lab_kev,
        args.resonance_energy_lab_keV,
        args.gamma_p_resonance_keV,
        args.gamma_out_keV,
        args.spin_factor,
        args.channel_radius_fm,
    )

    becker_e, becker_sigma, becker_error = becker_cross_section()
    (
        taskaev_e,
        taskaev_xerr,
        taskaev_sigma,
        taskaev_error,
    ) = taskaev_cross_section()
    sikora_e, sikora_sigma, sikora_error = sikora_cross_section()

    figure, axis = plt.subplots(figsize=(7, 4.6))

    axis.plot(
        energy_lab_kev,
        bw_sigma_barn,
        label="TUNL-parameter Breit-Wigner",
    )
    axis.errorbar(
        becker_e,
        becker_sigma,
        yerr=becker_error,
        fmt="s",
        linestyle="none",
        capsize=3,
        markersize=4.8,
        label="Becker et al. (1987)",
    )
    axis.errorbar(
        taskaev_e,
        taskaev_sigma,
        xerr=taskaev_xerr,
        yerr=taskaev_error,
        fmt="o",
        linestyle="none",
        capsize=3,
        markersize=5.2,
        label="Taskaev et al. (2024)",
    )
    if not args.exclude_sikora:
        axis.errorbar(
            sikora_e,
            sikora_sigma,
            yerr=sikora_error,
            fmt="^",
            linestyle="none",
            capsize=2.5,
            markersize=5.0,
            markerfacecolor="none",
            label="Sikora-Weller (model-corrected, 2016)",
        )

    axis.set_xlim(args.energy_min_keV, args.energy_max_keV)
    axis.set_ylim(0.0, args.ymax_barn)
    axis.set_xlabel(
        r"Proton laboratory energy, $E_{p,\mathrm{lab}}$ (keV)"
    )
    axis.set_ylabel(
        r"Total cross section, $\sigma$ (barn)"
    )
    axis.grid(alpha=0.25)
    axis.legend(loc="upper right", frameon=False)
    figure.tight_layout()

    if not args.no_save:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.output, dpi=300, bbox_inches="tight")

    peak_index = int(np.argmax(bw_sigma_barn))
    print("Breit-Wigner peak")
    print(f"  E_lab = {energy_lab_kev[peak_index]:.3f} keV")
    print(f"  sigma = {bw_sigma_barn[peak_index]:.6f} barn")
    if not args.no_save:
        print(f"Wrote figure: {args.output}")

    # --show behaves correctly on every backend: interactive backends open a
    # window; on Agg it is a harmless no-op (the figure was already saved).
    interactive = matplotlib.get_backend().lower() != "agg"
    if args.show and interactive:
        plt.show()
    elif args.show and not interactive:
        print(
            "[info] --show requested but a non-interactive (Agg) backend is "
            "active; the figure was saved to disk instead."
        )

    return figure


if __name__ == "__main__":
    if IN_NOTEBOOK:
        # `%run` inside Jupyter: show inline, keep the kernel alive.
        plot_675_bw_coulomb()
        plt.show()
    else:
        plot_675_bw_coulomb()
        plt.close("all")
        raise SystemExit(0)
