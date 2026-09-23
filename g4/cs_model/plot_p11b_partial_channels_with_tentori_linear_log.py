#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_p11b_partial_channels_with_tentori_linear_log.py

Plot 11B(p,alpha)8Be partial-channel literature data together with the
Tentori & Belloni 2023 total reference cross-section curve.

Energy axis:
  proton laboratory energy Ep_lab, 0--1 MeV

Outputs:
  p11b_partial_channels_with_tentori_linear.png
  p11b_partial_channels_with_tentori_log.png

Run:
  python3 plot_p11b_partial_channels_with_tentori_linear_log.py

Notes:
  - Becker 1987 alpha0 and alpha1 are partial-channel curves.
  - Taskaev 2024 alpha0 and alpha1 are partial-channel points.
  - Munch 2020 is alpha0 only.
  - Munch & Fynbo 2018 gives the 165-keV resonance alpha0 and alpha1 points.
  - Tentori & Belloni 2023 is a total p-11B fusion cross-section curve.
    It is not an alpha0/alpha1 decomposition.
"""

import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import matplotlib.pyplot as plt


EG_MEV = 22.589
LAB_TO_CM = 11.0 / 12.0
CM_TO_LAB = 12.0 / 11.0


def sigma_from_sfactor(E_cm_MeV, S_MeV_b):
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
# Becker 1987: partial-channel empirical fits
# ----------------------------------------------------------------------
def becker1987_s_alpha0(E_cm_MeV):
    """
    Becker 1987 alpha0 S(E).

    E_cm_MeV : center-of-mass energy in MeV
    return   : S_alpha0(E) in MeV barn
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
    Becker 1987 alpha1 S(E).

    E_cm_MeV : center-of-mass energy in MeV
    return   : S_alpha1(E) in MeV barn

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
# Tentori & Belloni 2023: total reference curve
# ----------------------------------------------------------------------
def tentori2023_s_total(E_cm_MeV):
    """
    Tentori & Belloni 2023 total reference S(E).

    This is the total p-11B fusion cross-section parameterization.
    It is not an alpha0 or alpha1 partial-channel curve.

    The low-energy 148-keV resonance term uses the Nevins-Swain parameters.
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    E_keV = 1000.0 * E
    S = np.full_like(E, np.nan, dtype=float)

    # E <= 0.400 MeV
    C0 = 197.0
    C1 = 0.269
    C2 = 2.54e-4

    # Nevins-Swain 148 keV resonance term
    AL = 1.82e4
    EL = 148.0
    dEL = 2.35

    # 0.400 < E <= 0.668 MeV
    E1_break = 0.400
    E2_break = 0.668
    D0 = 346.0
    D1 = 150.0
    D2 = -59.9
    D5 = -0.460

    # E > 0.668 MeV
    B = 0.381
    A = np.array([1.98e6, 3.89e6, 1.36e6, 3.71e6])
    ER = np.array([640.9, 1211.0, 2340.0, 3294.0])
    dE = np.array([85.5, 414.0, 221.0, 351.0])

    m1 = E <= E1_break
    if np.any(m1):
        x = E_keV[m1]
        S[m1] = C0 + C1 * x + C2 * x**2
        S[m1] += AL / ((x - EL) ** 2 + dEL**2)

    m2 = (E > E1_break) & (E <= E2_break)
    if np.any(m2):
        x = (E_keV[m2] - 400.0) / 100.0
        S[m2] = D0 + D1 * x + D2 * x**2 + D5 * x**5

    m3 = E > E2_break
    if np.any(m3):
        x = E_keV[m3]
        s3 = np.full_like(x, B, dtype=float)
        for Ak, Erk, dek in zip(A, ER, dE):
            s3 += Ak / ((x - Erk) ** 2 + dek**2)
        S[m3] = s3

    return S


def tentori2023_total_cross_section(E_cm_MeV):
    return sigma_from_sfactor(E_cm_MeV, tentori2023_s_total(E_cm_MeV))


# ----------------------------------------------------------------------
# Taskaev 2024: partial-channel points
# Tables 5 and 8, cross sections in mb.
# Energies are average proton energies in the boron layer.
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
# Munch 2020: alpha0 points
# Table 1, cross sections in mb.
# ----------------------------------------------------------------------
MUNCH2020_ALPHA0_Ep_keV = np.array([
    500, 600, 683, 800, 850, 900, 950, 1000
], dtype=float)

MUNCH2020_ALPHA0_mb = np.array([
    2.6, 3.9, 4.6, 5.9, 5.7, 6.2, 6.5, 6.7
], dtype=float)


# ----------------------------------------------------------------------
# Munch & Fynbo 2018: 165-keV resonance points
# Recommended cross sections near Ep_lab = 0.162 MeV:
#   sigma(p,alpha0) = 2.03 mb
#   sigma(p,alpha1) = 39 mb
# ----------------------------------------------------------------------
MUNCH2018_Ep_MeV = np.array([0.162], dtype=float)
MUNCH2018_ALPHA0_b = np.array([2.03e-3], dtype=float)
MUNCH2018_ALPHA1_b = np.array([39.0e-3], dtype=float)


def prepare_data():
    ep_max_lab = 1.0

    # Becker partial-channel curves.
    # Becker alpha1 formula is restricted to E_cm < 0.5 MeV.
    e_becker_cm = np.linspace(0.022, 0.499, 2500)
    ep_becker_lab = CM_TO_LAB * e_becker_cm
    sigma_becker_alpha0 = sigma_from_sfactor(e_becker_cm, becker1987_s_alpha0(e_becker_cm))
    sigma_becker_alpha1 = sigma_from_sfactor(e_becker_cm, becker1987_s_alpha1(e_becker_cm))

    # Tentori total curve.
    e_tentori_cm = np.linspace(0.02, LAB_TO_CM * ep_max_lab, 5000)
    ep_tentori_lab = CM_TO_LAB * e_tentori_cm
    sigma_tentori_total = tentori2023_total_cross_section(e_tentori_cm)

    # Taskaev points.
    ep_taskaev_lab = TASKAEV_Ep_keV / 1000.0
    sigma_taskaev_alpha0 = TASKAEV_ALPHA0_mb / 1000.0
    sigma_taskaev_alpha1 = TASKAEV_ALPHA1_mb / 1000.0

    # Munch 2020 alpha0 points.
    ep_munch2020_lab = MUNCH2020_ALPHA0_Ep_keV / 1000.0
    sigma_munch2020_alpha0 = MUNCH2020_ALPHA0_mb / 1000.0

    return {
        "ep_max_lab": ep_max_lab,
        "ep_becker_lab": ep_becker_lab,
        "sigma_becker_alpha0": sigma_becker_alpha0,
        "sigma_becker_alpha1": sigma_becker_alpha1,
        "ep_tentori_lab": ep_tentori_lab,
        "sigma_tentori_total": sigma_tentori_total,
        "ep_taskaev_lab": ep_taskaev_lab,
        "sigma_taskaev_alpha0": sigma_taskaev_alpha0,
        "sigma_taskaev_alpha1": sigma_taskaev_alpha1,
        "ep_munch2020_lab": ep_munch2020_lab,
        "sigma_munch2020_alpha0": sigma_munch2020_alpha0,
    }


def draw_common(data, yscale, output):
    """
    Draw and save one figure.

    yscale:
      "linear" or "log"
    output:
      output PNG filename
    """
    plt.figure(figsize=(10.0, 6.3))

    # Total reference.
    plt.plot(
        data["ep_tentori_lab"],
        data["sigma_tentori_total"],
        linestyle="-.",
        linewidth=2.4,
        label=r"Tentori & Belloni 2023 total"
    )

    # alpha0 data and curve.
    plt.plot(
        data["ep_becker_lab"],
        data["sigma_becker_alpha0"],
        linestyle="-",
        linewidth=2.0,
        label=r"Becker 1987 $\alpha_0$"
    )
    plt.scatter(
        data["ep_taskaev_lab"],
        data["sigma_taskaev_alpha0"],
        marker="o",
        s=34,
        label=r"Taskaev 2024 $\alpha_0$"
    )
    plt.scatter(
        data["ep_munch2020_lab"],
        data["sigma_munch2020_alpha0"],
        marker="^",
        s=42,
        label=r"Munch 2020 $\alpha_0$"
    )
    plt.scatter(
        MUNCH2018_Ep_MeV,
        MUNCH2018_ALPHA0_b,
        marker="D",
        s=46,
        label=r"Munch & Fynbo 2018 $\alpha_0$"
    )

    # alpha1 data and curve.
    plt.plot(
        data["ep_becker_lab"],
        data["sigma_becker_alpha1"],
        linestyle="--",
        linewidth=2.0,
        label=r"Becker 1987 $\alpha_1$"
    )
    plt.scatter(
        data["ep_taskaev_lab"],
        data["sigma_taskaev_alpha1"],
        marker="s",
        s=34,
        label=r"Taskaev 2024 $\alpha_1$"
    )
    plt.scatter(
        MUNCH2018_Ep_MeV,
        MUNCH2018_ALPHA1_b,
        marker="P",
        s=56,
        label=r"Munch & Fynbo 2018 $\alpha_1$"
    )

    if yscale == "log":
        plt.yscale("log")
        plt.ylim(1e-5, 3.0)
        title_suffix = "log scale"
    elif yscale == "linear":
        plt.yscale("linear")
        plt.ylim(0.0, 1.6)
        title_suffix = "linear scale"
    else:
        raise ValueError(f"Unknown yscale: {yscale}")

    plt.xlim(0.0, data["ep_max_lab"])

    plt.xlabel(r"proton lab energy $E_p$ (MeV)")
    plt.ylabel(r"cross section $\sigma$ (barn)")
    plt.title(
        r"$^{11}$B$(p,\alpha)^{8}$Be partial channels "
        rf"with total reference, $E_p \leq 1$ MeV, {title_suffix}"
    )

    plt.grid(True, which="both", alpha=0.3)
    plt.legend(ncol=2, fontsize=8.5)
    plt.tight_layout()

    plt.savefig(output, dpi=300)
    print(f"Saved: {output}")
    plt.close()


def main():
    data = prepare_data()

    draw_common(
        data,
        yscale="linear",
        output="p11b_partial_channels_with_tentori_linear.png"
    )

    draw_common(
        data,
        yscale="log",
        output="p11b_partial_channels_with_tentori_log.png"
    )


if __name__ == "__main__":
    main()
