#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_p11b_literature_total_xs_1MeV.py

Plot literature total cross sections of 11B(p,3alpha) up to Ep_lab = 1 MeV.

Included:
  - Becker 1987 original low-energy total cross section
  - Sikora & Weller 2016 evaluated total cross-section points
  - Tentori & Belloni 2023 reference total cross section

Run:
  python3 plot_p11b_literature_total.py

Output:
  p11b_literature_total_cross_sections.png
"""

import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import matplotlib.pyplot as plt


EG_MEV = 22.589
LAB_TO_CM = 11.0 / 12.0
CM_TO_LAB = 12.0 / 11.0


def sigma_from_literature_parameterization(E_cm_MeV, parameter_MeV_b):
    """
    Convert the Becker/Tentori published parameterization to cross section.

    E_cm_MeV       : center-of-mass energy in MeV
    parameter_MeV_b: published parameterization value in MeV barn
    return         : cross section in barn
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    P = np.asarray(parameter_MeV_b, dtype=float)

    sigma = np.full_like(E, np.nan, dtype=float)
    mask = E > 0.0
    sigma[mask] = P[mask] / E[mask] * np.exp(-np.sqrt(EG_MEV / E[mask]))
    return sigma


# ----------------------------------------------------------------------
# Becker 1987 original low-energy total cross section
# ----------------------------------------------------------------------
def becker1987_total_cross_section(E_cm_MeV):
    """
    Becker 1987 original low-energy empirical total cross-section fit.

    E_cm_MeV : center-of-mass energy in MeV
    return   : total cross section in barn
    """
    E = np.asarray(E_cm_MeV, dtype=float)

    alpha0_parameter = (
        2.1
        - 1.26 * E
        - 0.14 * E**2
        + 0.69e-3 / ((E - 0.148) ** 2 + 7.13e-6)
    )

    alpha1_parameter = (
        195.0
        + 241.0 * E
        + 231.0 * E**2
        + 1.76e-2 / ((E - 0.148) ** 2 + 5.52e-6)
    )

    total_parameter = alpha0_parameter + alpha1_parameter
    return sigma_from_literature_parameterization(E, total_parameter)


# ----------------------------------------------------------------------
# Sikora & Weller 2016 total cross-section points
# ----------------------------------------------------------------------
SIKORA_Ep_lab_MeV = np.array([
    0.15, 0.22, 0.25, 0.30, 0.40, 0.49, 0.57,
    0.65, 0.73, 0.80, 0.88, 0.94, 1.00
])

SIKORA_A0_mbsr = np.array([
    1.3532, 8.9100, 14.0372, 31.4329, 93.5936, 173.7688,
    285.1283, 333.7993, 273.8339, 172.2051, 110.7989,
    79.4042, 75.4242
])


def sikora2016_total_cross_section_points():
    """
    Sikora & Weller 2016, Ep_lab <= 1 MeV.

    A0 is in mb/sr. The total cross section is:
        sigma = 4*pi/3 * A0

    return:
        Ep_lab in MeV
        sigma in barn
    """
    sigma_mb = 4.0 * np.pi / 3.0 * SIKORA_A0_mbsr
    sigma_barn = sigma_mb / 1000.0
    return SIKORA_Ep_lab_MeV, sigma_barn


# ----------------------------------------------------------------------
# Tentori & Belloni 2023 reference total cross section
# ----------------------------------------------------------------------
def tentori2023_total_cross_section(E_cm_MeV):
    """
    Tentori & Belloni 2023 reference total cross-section parameterization.

    E_cm_MeV : center-of-mass energy in MeV
    return   : total cross section in barn
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    E_keV = 1000.0 * E

    parameter = np.full_like(E, np.nan, dtype=float)

    # E <= 0.400 MeV
    C0 = 197.0
    C1 = 0.269
    C2 = 2.54e-4
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

    mask1 = E <= E1_break
    if np.any(mask1):
        x = E_keV[mask1]
        parameter[mask1] = C0 + C1 * x + C2 * x**2
        parameter[mask1] += AL / ((x - EL) ** 2 + dEL**2)

    mask2 = (E > E1_break) & (E <= E2_break)
    if np.any(mask2):
        x = (E_keV[mask2] - 400.0) / 100.0
        parameter[mask2] = D0 + D1 * x + D2 * x**2 + D5 * x**5

    mask3 = E > E2_break
    if np.any(mask3):
        x = E_keV[mask3]
        p3 = np.full_like(x, B, dtype=float)
        for Ak, Erk, dek in zip(A, ER, dE):
            p3 += Ak / ((x - Erk) ** 2 + dek**2)
        parameter[mask3] = p3

    return sigma_from_literature_parameterization(E, parameter)


def main():
    ep_max_lab = 1.0

    # Becker curve. The Becker low-energy formula is only used below E_cm = 0.5 MeV.
    e_becker_cm = np.linspace(0.022, 0.499, 2500)
    ep_becker_lab = CM_TO_LAB * e_becker_cm
    sigma_becker = becker1987_total_cross_section(e_becker_cm)

    # Sikora points.
    ep_sikora_lab, sigma_sikora = sikora2016_total_cross_section_points()

    # Tentori curve up to Ep_lab = 1 MeV.
    e_tentori_cm = np.linspace(0.02, LAB_TO_CM * ep_max_lab, 5000)
    ep_tentori_lab = CM_TO_LAB * e_tentori_cm
    sigma_tentori = tentori2023_total_cross_section(e_tentori_cm)

    plt.figure(figsize=(9, 5.8))

    plt.plot(
        ep_tentori_lab,
        sigma_tentori,
        linewidth=2.2,
        label="Tentori & Belloni 2023"
    )

    plt.scatter(
        ep_sikora_lab,
        sigma_sikora,
        s=34,
        marker="o",
        label="Sikora & Weller 2016"
    )

    plt.plot(
        ep_becker_lab,
        sigma_becker,
        linestyle="--",
        linewidth=2.0,
        label="Becker 1987 original"
    )

    plt.yscale("log")
    plt.xlim(0.0, ep_max_lab)
    plt.ylim(1e-5, 3.0)

    plt.xlabel(r"proton lab energy $E_p$ (MeV)")
    plt.ylabel(r"total cross section $\sigma$ (barn)")
    plt.title(r"$^{11}$B$(p,3\alpha)$ total cross section, $E_p \leq 1$ MeV")

    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output = "p11b_literature_total_cross_sections.png"
    plt.savefig(output, dpi=300)
    print(f"Saved: {output}")
    plt.close()


if __name__ == "__main__":
    main()
