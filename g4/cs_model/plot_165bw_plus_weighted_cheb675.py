#!/usr/bin/env python3
"""Plot the final 165 BW + weighted_chebE16 analytic fit675 cross section."""

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
EG_MEV = 22.589

AMU_C2_KEV = 931_494.10242
HBARC_KEV_FM = 197_326.9804
FM2_TO_CM2 = 1.0e-26
BARN_CM2 = 1.0e-24

ENERGY_MIN_MEV = 0.020
ENERGY_MAX_MEV = 1.000
ENERGY_STEP_MEV = 0.0002

H11B_165_RESONANCE_ECM_KEV = 147.95
H11B_165_TOTAL_WIDTH_KEV = 5.3
H11B_165_PROTON_WIDTH_KEV = 0.0215
H11B_165_ALPHA0_WIDTH_KEV = 0.26
H11B_165_ALPHA1_WIDTH_KEV = 5.0
H11B_165_SPIN_STAT_FACTOR = 5.0 / 8.0
H11B_165_ENTRANCE_ORBITAL_L = 1

WEIGHTED_CHEB_E16_COEFFICIENTS = np.array(
    [
        -4.240991093218323,
        6.5335030639422,
        -5.348623593275032,
        2.836731408213122,
        -1.88972601027126,
        1.556619927270825,
        -0.9501800114612392,
        0.6222767937946161,
        -0.4865823078410484,
        0.3336851443104574,
        -0.2357424108239984,
        0.1936547242371347,
        -0.1385260429805155,
        0.1008831241333581,
        -0.05830063682899136,
        0.02350578311364693,
        0.009138850080015581,
    ],
    dtype=float,
)

OUTPUT_COMPARE_LOG = "p11b_cs_165bw_plus_weighted_cheb675_vs_tentori.log.png"
OUTPUT_MODEL_LOG = "p11b_cs_165bw_plus_weighted_cheb675_model.log.png"
OUTPUT_SUMMARY = "p11b_cs_165bw_plus_weighted_cheb675_summary.txt"


def make_energy_grid() -> np.ndarray:
    n_steps = int(round((ENERGY_MAX_MEV - ENERGY_MIN_MEV) / ENERGY_STEP_MEV))
    return ENERGY_MIN_MEV + ENERGY_STEP_MEV * np.arange(n_steps + 1, dtype=float)


def extract_cpp_array(text: str, name: str) -> np.ndarray:
    pattern = rf"(?:static\s+)?const double {re.escape(name)}\[\] = \{{(.*?)\}};"
    match = re.search(pattern, text, flags=re.S)
    if not match:
        raise RuntimeError(f"Could not find array {name} in {COULOMB_DATA}")

    number_pattern = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
    return np.array([float(value) for value in re.findall(number_pattern, match.group(1))], dtype=float)


def load_p11b_penetrability_tables() -> tuple[np.ndarray, dict[int, np.ndarray]]:
    text = COULOMB_DATA.read_text()
    return (
        extract_cpp_array(text, "P11B_energy_keV"),
        {
            0: extract_cpp_array(text, "P11B_L0_logP"),
            1: extract_cpp_array(text, "P11B_L1_logP"),
        },
    )


def interpolate_logp(energy_keV: np.ndarray, table_energy_keV: np.ndarray, table_logp: np.ndarray) -> np.ndarray:
    clipped_energy = np.clip(energy_keV, table_energy_keV[0], table_energy_keV[-1])
    return np.interp(np.log(clipped_energy), np.log(table_energy_keV), table_logp)


def penetrability_ratio(
    energy_keV: np.ndarray,
    reference_energy_keV: float,
    orbital_l: int,
    table_energy_keV: np.ndarray,
    table_logp_by_l: dict[int, np.ndarray],
) -> np.ndarray:
    logp = table_logp_by_l[orbital_l]
    logp_e = interpolate_logp(energy_keV, table_energy_keV, logp)
    logp_ref = interpolate_logp(np.array([reference_energy_keV]), table_energy_keV, logp)[0]
    return np.where(energy_keV > 0.0, np.exp(logp_e - logp_ref), 0.0)


def sigma_from_sfactor(e_cm_mev: np.ndarray, s_factor_mev_b: np.ndarray) -> np.ndarray:
    sigma = np.zeros_like(e_cm_mev, dtype=float)
    valid = (e_cm_mev > 0.0) & (s_factor_mev_b > 0.0)
    sigma[valid] = (
        s_factor_mev_b[valid] / e_cm_mev[valid] * np.exp(-np.sqrt(EG_MEV / e_cm_mev[valid]))
    )
    return sigma


def tentori2023_s_total(e_cm_mev: np.ndarray) -> np.ndarray:
    e = np.asarray(e_cm_mev, dtype=float)
    x_keV = 1000.0 * e
    s = np.zeros_like(e, dtype=float)

    mask1 = e <= 0.400
    if np.any(mask1):
        x = x_keV[mask1]
        s[mask1] = 197.0 + 0.269 * x + 2.54e-4 * x * x
        s[mask1] += 1.82e4 / ((x - 148.0) ** 2 + 2.35**2)

    mask2 = (e > 0.400) & (e <= 0.668)
    if np.any(mask2):
        x = (x_keV[mask2] - 400.0) / 100.0
        s[mask2] = 346.0 + 150.0 * x - 59.9 * x * x - 0.460 * x**5

    mask3 = e > 0.668
    if np.any(mask3):
        x = x_keV[mask3]
        s3 = np.full_like(x, 0.381, dtype=float)
        for amplitude, center, width in (
            (1.98e6, 640.9, 85.5),
            (3.89e6, 1211.0, 414.0),
            (1.36e6, 2340.0, 221.0),
            (3.71e6, 3294.0, 351.0),
        ):
            s3 += amplitude / ((x - center) ** 2 + width**2)
        s[mask3] = s3

    return s


def tentori2023_total_cross_section_b(e_cm_mev: np.ndarray) -> np.ndarray:
    return sigma_from_sfactor(e_cm_mev, tentori2023_s_total(e_cm_mev))


def breit_wigner_sigma_cm2(
    ecm_keV: np.ndarray,
    resonance_energy_keV: float,
    total_width_keV: float,
    entrance_width_at_resonance_keV: float,
    exit_width_at_resonance_keV: float,
    spin_stat_factor: float,
    entrance_orbital_l: int,
    table_energy_keV: np.ndarray,
    table_logp_by_l: dict[int, np.ndarray],
) -> np.ndarray:
    sigma = np.zeros_like(ecm_keV, dtype=float)
    valid = ecm_keV > 0.0
    if not np.any(valid):
        return sigma

    reduced_mass_keV = 1.0 * 11.0 / 12.0 * AMU_C2_KEV
    gamma_entrance = entrance_width_at_resonance_keV * penetrability_ratio(
        ecm_keV,
        resonance_energy_keV,
        entrance_orbital_l,
        table_energy_keV,
        table_logp_by_l,
    )
    missing_width = max(0.0, total_width_keV - entrance_width_at_resonance_keV - exit_width_at_resonance_keV)
    gamma_total = gamma_entrance + exit_width_at_resonance_keV + missing_width
    pi_over_k2_fm2 = math.pi * HBARC_KEV_FM * HBARC_KEV_FM / (2.0 * reduced_mass_keV * ecm_keV)
    denominator = (ecm_keV - resonance_energy_keV) ** 2 + gamma_total * gamma_total / 4.0
    bw_factor = spin_stat_factor * gamma_entrance * exit_width_at_resonance_keV / denominator
    sigma[valid] = pi_over_k2_fm2[valid] * bw_factor[valid] * FM2_TO_CM2
    return sigma


def sigma_165_bw_b(ep_lab_mev: np.ndarray, table_energy_keV: np.ndarray, table_logp_by_l: dict[int, np.ndarray]) -> np.ndarray:
    ecm_keV = 1000.0 * LAB_TO_CM * ep_lab_mev
    sigma_cm2 = breit_wigner_sigma_cm2(
        ecm_keV,
        H11B_165_RESONANCE_ECM_KEV,
        H11B_165_TOTAL_WIDTH_KEV,
        H11B_165_PROTON_WIDTH_KEV,
        H11B_165_ALPHA0_WIDTH_KEV + H11B_165_ALPHA1_WIDTH_KEV,
        H11B_165_SPIN_STAT_FACTOR,
        H11B_165_ENTRANCE_ORBITAL_L,
        table_energy_keV,
        table_logp_by_l,
    )
    return np.where(ecm_keV >= 1.0, sigma_cm2 / BARN_CM2, 0.0)


def chebyshev_t_series(z: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    result = np.full_like(z, coefficients[0], dtype=float)
    if len(coefficients) == 1:
        return result

    t_nm2 = np.ones_like(z)
    t_nm1 = z
    result += coefficients[1] * t_nm1
    for n in range(2, len(coefficients)):
        t_n = 2.0 * z * t_nm1 - t_nm2
        result += coefficients[n] * t_n
        t_nm2 = t_nm1
        t_nm1 = t_n
    return result


def sigma_weighted_cheb675_b(ep_lab_mev: np.ndarray) -> np.ndarray:
    z = 2.0 * (ep_lab_mev - ENERGY_MIN_MEV) / (ENERGY_MAX_MEV - ENERGY_MIN_MEV) - 1.0
    return np.exp(chebyshev_t_series(z, WEIGHTED_CHEB_E16_COEFFICIENTS))


def calculate_curves() -> dict[str, np.ndarray]:
    ep_lab_mev = make_energy_grid()
    table_energy_keV, table_logp_by_l = load_p11b_penetrability_tables()
    sigma_165 = sigma_165_bw_b(ep_lab_mev, table_energy_keV, table_logp_by_l)
    sigma_675 = sigma_weighted_cheb675_b(ep_lab_mev)
    return {
        "ep_lab_mev": ep_lab_mev,
        "sigma_tentori_b": tentori2023_total_cross_section_b(LAB_TO_CM * ep_lab_mev),
        "sigma_165_bw_b": sigma_165,
        "sigma_675_weighted_cheb_b": sigma_675,
        "sigma_165_plus_weighted_cheb675_b": sigma_165 + sigma_675,
    }


def write_summary(curves: dict[str, np.ndarray]) -> None:
    ep = curves["ep_lab_mev"]
    main675 = (ep >= 0.550) & (ep <= 0.750)
    peak_idx = np.flatnonzero(main675)[int(np.argmax(curves["sigma_675_weighted_cheb_b"][main675]))]
    lines = [
        "Current 165 BW + weighted_chebE16 analytic fit675",
        "",
        "Function:",
        "  sigma_fit675(Ep_lab) = exp(sum c_n T_n(z)) barn",
        "  z = 2*(Ep_lab - 0.020)/(1.000 - 0.020) - 1",
        "  degree = 16",
        "",
        "Main 675-region fit675 component:",
        f"  weighted_chebE16 peak at {ep[peak_idx]:.6f} MeV lab",
        f"  weighted_chebE16 peak = {curves['sigma_675_weighted_cheb_b'][peak_idx]:.6f} b",
        f"  165 BW + weighted_chebE16 at peak = {curves['sigma_165_plus_weighted_cheb675_b'][peak_idx]:.6f} b",
        f"  Tentori2023 total at peak = {curves['sigma_tentori_b'][peak_idx]:.6f} b",
        "",
        "Generated plots:",
        f"  {OUTPUT_COMPARE_LOG}",
        f"  {OUTPUT_MODEL_LOG}",
    ]
    (SCRIPT_DIR / OUTPUT_SUMMARY).write_text("\n".join(lines) + "\n")


def draw_compare_log(curves: dict[str, np.ndarray]) -> None:
    ep = curves["ep_lab_mev"]
    plt.figure(figsize=(10.8, 6.6))
    plt.plot(ep, curves["sigma_tentori_b"], linewidth=2.7, label="Tentori 2023 total")
    plt.plot(ep, curves["sigma_165_bw_b"], linewidth=1.9, label="current 165 BW")
    plt.plot(ep, curves["sigma_675_weighted_cheb_b"], linewidth=2.1, label="weighted_chebE16 fit675")
    plt.plot(
        ep,
        curves["sigma_165_plus_weighted_cheb675_b"],
        linewidth=2.0,
        linestyle="--",
        label="current 165 BW + weighted_chebE16 fit675",
    )
    plt.yscale("log")
    positive = curves["sigma_tentori_b"][curves["sigma_tentori_b"] > 0.0]
    plt.ylim(max(1.0e-8, 0.2 * float(np.min(positive))), 2.0 * float(np.max(curves["sigma_tentori_b"])))
    plt.xlabel("proton lab energy (MeV)")
    plt.ylabel("cross section (barn)")
    plt.title("Current 165 BW plus weighted_chebE16 fit675 vs Tentori total")
    plt.xlim(ENERGY_MIN_MEV, ENERGY_MAX_MEV)
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(frameon=True, fontsize=8.2)
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / OUTPUT_COMPARE_LOG, dpi=220)
    plt.close()


def draw_model_log(curves: dict[str, np.ndarray]) -> None:
    ep = curves["ep_lab_mev"]
    plt.figure(figsize=(10.8, 6.6))
    plt.plot(ep, curves["sigma_165_bw_b"], linewidth=2.1, label="current 165 BW")
    plt.plot(ep, curves["sigma_675_weighted_cheb_b"], linewidth=2.3, label="weighted_chebE16 fit675")
    plt.plot(
        ep,
        curves["sigma_165_plus_weighted_cheb675_b"],
        linewidth=2.3,
        linestyle="--",
        label="current 165 BW + weighted_chebE16 fit675",
    )
    plt.yscale("log")
    positive = curves["sigma_165_plus_weighted_cheb675_b"][curves["sigma_165_plus_weighted_cheb675_b"] > 0.0]
    plt.ylim(max(1.0e-8, 0.2 * float(np.min(positive))), 2.0 * float(np.max(curves["sigma_165_plus_weighted_cheb675_b"])))
    plt.xlabel("proton lab energy (MeV)")
    plt.ylabel("cross section (barn)")
    plt.title("Current model: 165 BW plus weighted_chebE16 fit675")
    plt.xlim(ENERGY_MIN_MEV, ENERGY_MAX_MEV)
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(frameon=True, fontsize=8.2)
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / OUTPUT_MODEL_LOG, dpi=220)
    plt.close()


def main() -> None:
    curves = calculate_curves()
    write_summary(curves)
    draw_compare_log(curves)
    draw_model_log(curves)

    for line in (SCRIPT_DIR / OUTPUT_SUMMARY).read_text().splitlines():
        print(line)


if __name__ == "__main__":
    main()
