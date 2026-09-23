#!/usr/bin/env python3
"""Plot the 162/165 keV single-level Breit-Wigner cross section for the slide.

Reproduces H11BCrossSection::GetSigma165 exactly:
  sigma(E) = pi/k^2 * omega * Gamma_p(E) * Gamma_exit
             / [ (E-Er)^2 + Gamma_tot(E)^2/4 ]
with the energy-dependent entrance width
  Gamma_p(E) = Gamma_p(Er) * P1(E)/P1(Er)
from the tabulated full Coulomb penetrability, and a constant-width
Breit-Wigner drawn for comparison.

Outputs (written next to this script):
  p11b_cs_162bw_slide.log.png
"""

from __future__ import annotations

import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
SCRIPT_DIR = THIS_FILE.parent
PROJECT_ROOT = THIS_FILE.parents[1]
COULOMB_DATA = PROJECT_ROOT / "src" / "CoulombPenetrabilityData.cc"

LAB_TO_CM = 11.0 / 12.0

AMU_C2_KEV = 931_494.10242
HBARC_KEV_FM = 197_326.9804
FM2_TO_CM2 = 1.0e-26
BARN_CM2 = 1.0e-24

# Constants.hh values (H11B165*)
RESONANCE_ECM_KEV = 147.95
TOTAL_WIDTH_KEV = 5.3
PROTON_WIDTH_KEV = 0.0215
ALPHA0_WIDTH_KEV = 0.26
ALPHA1_WIDTH_KEV = 5.0
SPIN_STAT_FACTOR = 5.0 / 8.0
ENTRANCE_ORBITAL_L = 1

# plot range (c.m. energy)
ECM_MIN_KEV = 50.0
ECM_MAX_KEV = 400.0
ECM_STEP_KEV = 0.05

OUTPUT_PLOT = "p11b_cs_162bw_slide.log.png"


def extract_cpp_array(text: str, name: str) -> np.ndarray:
    pattern = rf"(?:static\s+)?const double {re.escape(name)}\[\] = \{{(.*?)\}};"
    match = re.search(pattern, text, flags=re.S)
    if not match:
        raise RuntimeError(f"Could not find array {name} in {COULOMB_DATA}")
    values = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", match.group(1))
    return np.array([float(v) for v in values], dtype=float)


def load_p11b_tables() -> tuple[np.ndarray, np.ndarray]:
    text = COULOMB_DATA.read_text()
    energy = extract_cpp_array(text, "P11B_energy_keV")
    log_p = extract_cpp_array(text, f"P11B_L{ENTRANCE_ORBITAL_L}_logP")
    if energy.size != log_p.size:
        raise RuntimeError("P11B table size mismatch")
    return energy, log_p


def interp_log_p(ecm_kev: np.ndarray, table_e: np.ndarray, table_logp: np.ndarray) -> np.ndarray:
    """log-log interpolation, clamped at the table ends (matches the C++ code)."""
    log_e = np.log(np.clip(ecm_kev, table_e[0], table_e[-1]))
    return np.interp(log_e, np.log(table_e), table_logp)


def penetrability_ratio(ecm_kev: np.ndarray, er_kev: float,
                        table_e: np.ndarray, table_logp: np.ndarray) -> np.ndarray:
    log_p_e = interp_log_p(np.asarray(ecm_kev, dtype=float), table_e, table_logp)
    log_p_r = interp_log_p(np.array([er_kev]), table_e, table_logp)[0]
    return np.exp(log_p_e - log_p_r)


def sigma_bw_barn(ecm_kev: np.ndarray, ratio: np.ndarray | float) -> np.ndarray:
    """Single-level BW in barn; ratio = P1(E)/P1(Er) (1.0 -> constant width)."""
    ecm = np.asarray(ecm_kev, dtype=float)
    reduced_mass_kev = 1.0 * 11.0 / 12.0 * AMU_C2_KEV
    pi_over_k2_fm2 = math.pi * HBARC_KEV_FM**2 / (2.0 * reduced_mass_kev) / ecm

    gamma_exit = ALPHA0_WIDTH_KEV + ALPHA1_WIDTH_KEV
    gamma_p_e = PROTON_WIDTH_KEV * np.asarray(ratio, dtype=float)
    missing = max(0.0, TOTAL_WIDTH_KEV - PROTON_WIDTH_KEV - gamma_exit)
    gamma_tot_e = gamma_p_e + gamma_exit + missing

    denom = (ecm - RESONANCE_ECM_KEV) ** 2 + gamma_tot_e**2 / 4.0
    bw = SPIN_STAT_FACTOR * gamma_p_e * gamma_exit / denom
    return pi_over_k2_fm2 * bw * FM2_TO_CM2 / BARN_CM2


def main() -> None:
    table_e, table_logp = load_p11b_tables()

    ecm = np.arange(ECM_MIN_KEV, ECM_MAX_KEV + ECM_STEP_KEV, ECM_STEP_KEV)
    ratio = penetrability_ratio(ecm, RESONANCE_ECM_KEV, table_e, table_logp)

    sigma_full = sigma_bw_barn(ecm, ratio)          # energy-dependent Gamma_p
    sigma_const = sigma_bw_barn(ecm, 1.0)           # constant-width reference

    ep_lab = ecm / LAB_TO_CM

    # ---- summary ----
    i_peak = int(np.argmax(sigma_full))
    lines = [
        "162/165 keV single-level Breit-Wigner (slide figure)",
        "",
        f"Er            = {RESONANCE_ECM_KEV:.2f} keV (c.m.)"
        f"  ->  Ep(lab) = {RESONANCE_ECM_KEV / LAB_TO_CM:.2f} keV",
        f"Gamma_tot(Er) = {TOTAL_WIDTH_KEV:.2f} keV",
        f"Gamma_p(Er)   = {PROTON_WIDTH_KEV * 1e3:.1f} eV   (L = {ENTRANCE_ORBITAL_L}, full Coulomb)",
        f"Gamma_exit    = {ALPHA0_WIDTH_KEV + ALPHA1_WIDTH_KEV:.2f} keV (alpha0 + alpha1, constant)",
        f"omega         = {SPIN_STAT_FACTOR:.4f}",
        "",
        f"peak: sigma = {sigma_full[i_peak]:.4f} b at Ecm = {ecm[i_peak]:.2f} keV"
        f" (Ep lab = {ep_lab[i_peak]:.2f} keV)",
        "",
        "suppression of the full-Coulomb curve relative to constant width:",
    ]
    for e_probe in (60.0, 80.0, 100.0, 120.0):
        s_f = sigma_bw_barn(np.array([e_probe]),
                            penetrability_ratio(np.array([e_probe]), RESONANCE_ECM_KEV,
                                                table_e, table_logp))[0]
        s_c = sigma_bw_barn(np.array([e_probe]), 1.0)[0]
        lines.append(f"  Ecm = {e_probe:5.1f} keV : full/const = {s_f / s_c:.3e}")
    summary = "\n".join(lines) + "\n"
    print(summary)

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.plot(ecm, sigma_const, "--", color="#C77D2E", lw=1.8,
            label=r"constant width $\Gamma_p(E_r)$")
    ax.plot(ecm, sigma_full, "-", color="#0E7C6B", lw=2.6,
            label=r"$\Gamma_p(E)=\Gamma_p(E_r)\,P_1(E)/P_1(E_r)$")
    ax.axvline(RESONANCE_ECM_KEV, color="#1A2E35", lw=1.0, ls=":", alpha=0.7)
    ax.annotate(rf"$E_r$ = {RESONANCE_ECM_KEV:.1f} keV",
                xy=(RESONANCE_ECM_KEV, 0.0), xycoords=("data", "axes fraction"),
                xytext=(6, 8), textcoords="offset points",
                fontsize=10, color="#1A2E35")

    ax.set_yscale("log")
    ax.set_xlim(ECM_MIN_KEV, ECM_MAX_KEV)
    ax.set_xlabel(r"$E_{\mathrm{c.m.}}$ (keV)", fontsize=12)
    ax.set_ylabel(r"$\sigma(p,\alpha)$ (barn)", fontsize=12)
    ax.set_title(r"$^{11}$B(p,$\alpha$) 162 keV resonance: single-level BW",
                 fontsize=12)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=10, loc="lower right")
    fig.tight_layout()
    fig.savefig(SCRIPT_DIR / OUTPUT_PLOT, dpi=200)
    print(f"wrote {SCRIPT_DIR / OUTPUT_PLOT}")


if __name__ == "__main__":
    main()
