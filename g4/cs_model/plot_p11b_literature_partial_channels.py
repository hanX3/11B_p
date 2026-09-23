#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_p11b_literature_partial_channels.py

Plot literature partial-channel cross sections for 11B(p,alpha) up to Ep_lab = 1 MeV.

Channels:
  alpha0: 11B(p,alpha0)8Be(g.s.)
  alpha1: 11B(p,alpha1)8Be*(2+)

Included:
  - Becker 1987 alpha0 and alpha1 curves from published low-energy empirical fits
  - Taskaev 2024 alpha0 and alpha1 total cross-section points
  - Munch 2020 alpha0 points
  - Munch & Fynbo 2018 alpha0 and alpha1 resonance points near Ep_lab = 0.162 MeV

Run:
  python3 plot_p11b_literature_partial_channels.py

Output:
  p11b_literature_partial_channels.png
"""

import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import matplotlib.pyplot as plt


EG_MEV = 22.589
LAB_TO_CM = 11.0 / 12.0
CM_TO_LAB = 12.0 / 11.0


def sigma_from_s(E_cm_MeV, S_MeV_b):
    """
    Convert S(E) to cross section.

    sigma(E) = S(E) / E * exp[-sqrt(EG/E)]

    E_cm_MeV : center-of-mass energy in MeV
    S_MeV_b  : S(E) in MeV barn
    return   : cross section in barn
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    S = np.asarray(S_MeV_b, dtype=float)

    sigma = np.full_like(E, np.nan, dtype=float)
    mask = E > 0.0
    sigma[mask] = S[mask] / E[mask] * np.exp(-np.sqrt(EG_MEV / E[mask]))
    return sigma


# ----------------------------------------------------------------------
# Becker 1987 partial-channel empirical fits
# ----------------------------------------------------------------------
def becker1987_s_alpha0(E_cm_MeV):
    """
    Becker 1987 alpha0 S(E), E in MeV, S in MeV barn.
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    return (
        2.1
        - 1.26 * E
        - 0.14 * E**2
        + 0.69e-3 / ((E - 0.148) ** 2 + 7.13e-6)
    )


def becker1987_s_alpha1(E_cm_MeV):
    """
    Becker 1987 alpha1 S(E), E in MeV, S in MeV barn.

    The published alpha1 expression is used only below E_cm = 0.5 MeV.
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    return (
        195.0
        + 241.0 * E
        + 231.0 * E**2
        + 1.76e-2 / ((E - 0.148) ** 2 + 5.52e-6)
    )


# ----------------------------------------------------------------------
# Taskaev 2024 partial-channel points
# Tables 5 and 8, values in mb.
# Energies are the average proton energies in the boron layer.
# ----------------------------------------------------------------------
TASKAEV_Ep_keV = np.array([
    75, 134, 247, 355, 461, 565, 668, 771, 873, 975
], dtype=float)

TASKAEV_ALPHA0_mb = np.array([
    0.017, 0.27, 0.64, 1.40, 2.02, 2.53, 3.31, 3.45, 3.04, 2.90
], dtype=float)

TASKAEV_ALPHA1_mb = np.array([
    0.56, 6.4, 40.0, 148.0, 357.0, 598.0, 668.0, 386.0, 234.0, 171.0
], dtype=float)


# ----------------------------------------------------------------------
# Munch 2020 alpha0 points
# Table 1, values in mb.
# ----------------------------------------------------------------------
MUNCH2020_ALPHA0_Ep_keV = np.array([
    500, 600, 683, 800, 850, 900, 950, 1000
], dtype=float)

MUNCH2020_ALPHA0_mb = np.array([
    2.6, 3.9, 4.6, 5.9, 5.7, 6.2, 6.5, 6.7
], dtype=float)


# ----------------------------------------------------------------------
# Munch & Fynbo 2018 resonance points near Ep_lab = 162 keV
# Recommended peak cross sections:
#   sigma(p,alpha0) = 2.03 mb
#   sigma(p,alpha1) = 39 mb
# ----------------------------------------------------------------------
MUNCH2018_RESONANCE_Ep_MeV = np.array([0.162], dtype=float)
MUNCH2018_ALPHA0_b = np.array([2.03e-3], dtype=float)
MUNCH2018_ALPHA1_b = np.array([39.0e-3], dtype=float)


def main():
    ep_max_lab = 1.0

    # Becker curves: alpha1 expression is limited to E_cm < 0.5 MeV.
    e_becker_cm = np.linspace(0.022, 0.499, 2500)
    ep_becker_lab = CM_TO_LAB * e_becker_cm
    sigma_becker_alpha0 = sigma_from_s(e_becker_cm, becker1987_s_alpha0(e_becker_cm))
    sigma_becker_alpha1 = sigma_from_s(e_becker_cm, becker1987_s_alpha1(e_becker_cm))

    # Taskaev points.
    ep_taskaev_lab = TASKAEV_Ep_keV / 1000.0
    sigma_taskaev_alpha0 = TASKAEV_ALPHA0_mb / 1000.0
    sigma_taskaev_alpha1 = TASKAEV_ALPHA1_mb / 1000.0

    # Munch 2020 alpha0 points.
    ep_munch2020_lab = MUNCH2020_ALPHA0_Ep_keV / 1000.0
    sigma_munch2020_alpha0 = MUNCH2020_ALPHA0_mb / 1000.0

    plt.figure(figsize=(9.4, 6.0))

    # alpha0
    plt.plot(
        ep_becker_lab,
        sigma_becker_alpha0,
        linestyle="-",
        linewidth=2.0,
        label=r"Becker 1987 $\alpha_0$"
    )
    plt.scatter(
        ep_taskaev_lab,
        sigma_taskaev_alpha0,
        marker="o",
        s=34,
        label=r"Taskaev 2024 $\alpha_0$"
    )
    plt.scatter(
        ep_munch2020_lab,
        sigma_munch2020_alpha0,
        marker="^",
        s=42,
        label=r"Munch 2020 $\alpha_0$"
    )
    plt.scatter(
        MUNCH2018_RESONANCE_Ep_MeV,
        MUNCH2018_ALPHA0_b,
        marker="D",
        s=44,
        label=r"Munch & Fynbo 2018 $\alpha_0$"
    )

    # alpha1
    plt.plot(
        ep_becker_lab,
        sigma_becker_alpha1,
        linestyle="--",
        linewidth=2.0,
        label=r"Becker 1987 $\alpha_1$"
    )
    plt.scatter(
        ep_taskaev_lab,
        sigma_taskaev_alpha1,
        marker="s",
        s=34,
        label=r"Taskaev 2024 $\alpha_1$"
    )
    plt.scatter(
        MUNCH2018_RESONANCE_Ep_MeV,
        MUNCH2018_ALPHA1_b,
        marker="P",
        s=54,
        label=r"Munch & Fynbo 2018 $\alpha_1$"
    )

    plt.yscale("log")
    plt.xlim(0.0, ep_max_lab)
    plt.ylim(1e-5, 3.0)

    plt.xlabel(r"proton lab energy $E_p$ (MeV)")
    plt.ylabel(r"partial cross section $\sigma$ (barn)")
    plt.title(r"$^{11}$B$(p,\alpha)^{8}$Be partial channels, $E_p \leq 1$ MeV")

    plt.grid(True, which="both", alpha=0.3)
    plt.legend(ncol=2, fontsize=9)
    plt.tight_layout()

    output = "p11b_literature_partial_channels.png"
    plt.savefig(output, dpi=300)
    print(f"Saved: {output}")
    plt.close()


if __name__ == "__main__":
    main()
