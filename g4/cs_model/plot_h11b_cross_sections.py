#!/usr/bin/env python3
"""Plot the p + 11B cross-section decomposition used by this Geant4 app.

The C++ model is implemented in src/H11BCrossSection.cc.  It is a sum of two
single-level Breit-Wigner components.  The "with Coulomb" version uses the same
precomputed p+11B penetrability tables as the simulation to make the entrance
proton width energy dependent:

    Gamma_p(E) = Gamma_p(Er) * P_L(E) / P_L(Er)

The "without Coulomb" version keeps Gamma_p fixed at Gamma_p(Er).

The current simulation normalizes the total reaction probability to the
Tentori & Belloni 2023 evaluated total cross section and treats the residual
between that total and the explicit 165/675 components as a phenomenological
background 3alpha channel.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
COULOMB_DATA = PROJECT_ROOT / "src" / "CoulombPenetrabilityData.cc"

H11B_165_RESONANCE_ECM_KEV = 148.3
H11B_165_TOTAL_WIDTH_KEV = 5.3
H11B_165_PROTON_WIDTH_KEV = 0.0215
H11B_165_ALPHA0_WIDTH_KEV = 0.26
H11B_165_ALPHA1_WIDTH_KEV = 5.0
H11B_165_SPIN_STAT_FACTOR = 5.0 / 8.0
H11B_165_ENTRANCE_ORBITAL_L = 1

H11B_675_RESONANCE_ECM_KEV = 618.75
H11B_675_TOTAL_WIDTH_KEV = 300.0
H11B_675_PROTON_WIDTH_KEV = 150.0
H11B_675_ALPHA0_WIDTH_KEV = 0.0
H11B_675_ALPHA1_WIDTH_KEV = 150.0
H11B_675_SPIN_STAT_FACTOR = 5.0 / 8.0
H11B_675_ENTRANCE_ORBITAL_L = 0

AMU_C2_KEV = 931_494.10242
HBARC_KEV_FM = 197_326.9804
FM2_TO_CM2 = 1.0e-26
BARN_CM2 = 1.0e-24
LAB_TO_CM = 11.0 / 12.0
CM_TO_LAB = 12.0 / 11.0
EG_MEV = 22.589


def extract_cpp_array(text: str, name: str) -> np.ndarray:
    pattern = rf"(?:static\s+)?const double {re.escape(name)}\[\] = \{{(.*?)\}};"
    match = re.search(pattern, text, flags=re.S)
    if not match:
        raise ValueError(f"Could not find array {name} in {COULOMB_DATA}")

    number_pattern = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
    return np.array([float(value) for value in re.findall(number_pattern, match.group(1))], dtype=float)


def load_p11b_penetrability_tables() -> tuple[np.ndarray, dict[int, np.ndarray]]:
    text = COULOMB_DATA.read_text()
    energy = extract_cpp_array(text, "P11B_energy_keV")
    logp_l0 = extract_cpp_array(text, "P11B_L0_logP")
    logp_l1 = extract_cpp_array(text, "P11B_L1_logP")
    return energy, {0: logp_l0, 1: logp_l1}


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
    if orbital_l not in table_logp_by_l:
        raise ValueError(f"No p+11B penetrability table for L={orbital_l}")

    logp = table_logp_by_l[orbital_l]
    logp_e = interpolate_logp(energy_keV, table_energy_keV, logp)
    logp_ref = interpolate_logp(np.array([reference_energy_keV]), table_energy_keV, logp)[0]
    ratio = np.exp(logp_e - logp_ref)
    return np.where(energy_keV > 0.0, ratio, 0.0)


def breit_wigner_sigma_cm2(
    ecm_keV: np.ndarray,
    resonance_energy_keV: float,
    total_width_keV: float,
    entrance_width_at_resonance_keV: float,
    exit_width_at_resonance_keV: float,
    spin_stat_factor: float,
    entrance_orbital_l: int,
    use_coulomb: bool,
    table_energy_keV: np.ndarray,
    table_logp_by_l: dict[int, np.ndarray],
) -> np.ndarray:
    sigma = np.zeros_like(ecm_keV, dtype=float)
    valid = (
        (ecm_keV > 0.0)
        & (resonance_energy_keV > 0.0)
        & (total_width_keV > 0.0)
        & (entrance_width_at_resonance_keV > 0.0)
        & (exit_width_at_resonance_keV > 0.0)
    )
    if not np.any(valid):
        return sigma

    project_a = 1.0
    target_a = 11.0
    reduced_mass_keV = project_a * target_a / (project_a + target_a) * AMU_C2_KEV

    gamma_entrance = np.full_like(ecm_keV, entrance_width_at_resonance_keV, dtype=float)
    if use_coulomb:
        gamma_entrance *= penetrability_ratio(
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


def components(ep_lab_keV: np.ndarray, use_coulomb: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table_energy_keV, table_logp_by_l = load_p11b_penetrability_tables()
    ecm_keV = 11.0 / 12.0 * ep_lab_keV

    sigma_165 = breit_wigner_sigma_cm2(
        ecm_keV,
        H11B_165_RESONANCE_ECM_KEV,
        H11B_165_TOTAL_WIDTH_KEV,
        H11B_165_PROTON_WIDTH_KEV,
        H11B_165_ALPHA0_WIDTH_KEV + H11B_165_ALPHA1_WIDTH_KEV,
        H11B_165_SPIN_STAT_FACTOR,
        H11B_165_ENTRANCE_ORBITAL_L,
        use_coulomb,
        table_energy_keV,
        table_logp_by_l,
    )

    sigma_675 = breit_wigner_sigma_cm2(
        ecm_keV,
        H11B_675_RESONANCE_ECM_KEV,
        H11B_675_TOTAL_WIDTH_KEV,
        H11B_675_PROTON_WIDTH_KEV,
        H11B_675_ALPHA0_WIDTH_KEV + H11B_675_ALPHA1_WIDTH_KEV,
        H11B_675_SPIN_STAT_FACTOR,
        H11B_675_ENTRANCE_ORBITAL_L,
        use_coulomb,
        table_energy_keV,
        table_logp_by_l,
    )

    sigma_165 = np.where((ecm_keV >= 1.0) & (ecm_keV <= 400.0), sigma_165, 0.0)
    sigma_675 = np.where((ecm_keV >= 1.0) & (ecm_keV <= 3500.0), sigma_675, 0.0)

    return sigma_165, sigma_675, sigma_165 + sigma_675


def cm2_to_barn(sigma_cm2: np.ndarray) -> np.ndarray:
    return sigma_cm2 / BARN_CM2


def sigma_from_sfactor_mev(E_cm_MeV: np.ndarray, S_MeV_b: np.ndarray) -> np.ndarray:
    E = np.asarray(E_cm_MeV, dtype=float)
    S = np.asarray(S_MeV_b, dtype=float)

    sigma = np.full_like(E, np.nan, dtype=float)
    valid = E > 0.0
    sigma[valid] = S[valid] / E[valid] * np.exp(-np.sqrt(EG_MEV / E[valid]))
    return sigma


def tentori2023_s_total(E_cm_MeV: np.ndarray) -> np.ndarray:
    E = np.asarray(E_cm_MeV, dtype=float)
    E_keV = 1000.0 * E
    S = np.full_like(E, np.nan, dtype=float)

    C0 = 197.0
    C1 = 0.269
    C2 = 2.54e-4
    AL = 1.82e4
    EL = 148.0
    dEL = 2.35

    E1_break = 0.400
    E2_break = 0.668
    D0 = 346.0
    D1 = 150.0
    D2 = -59.9
    D5 = -0.460

    B = 0.381
    A = np.array([1.98e6, 3.89e6, 1.36e6, 3.71e6])
    ER = np.array([640.9, 1211.0, 2340.0, 3294.0])
    dE = np.array([85.5, 414.0, 221.0, 351.0])

    mask1 = E <= E1_break
    if np.any(mask1):
        x = E_keV[mask1]
        S[mask1] = C0 + C1 * x + C2 * x**2
        S[mask1] += AL / ((x - EL) ** 2 + dEL**2)

    mask2 = (E > E1_break) & (E <= E2_break)
    if np.any(mask2):
        x = (E_keV[mask2] - 400.0) / 100.0
        S[mask2] = D0 + D1 * x + D2 * x**2 + D5 * x**5

    mask3 = E > E2_break
    if np.any(mask3):
        x = E_keV[mask3]
        s3 = np.full_like(x, B, dtype=float)
        for Ak, Erk, dek in zip(A, ER, dE):
            s3 += Ak / ((x - Erk) ** 2 + dek**2)
        S[mask3] = s3

    return S


def tentori2023_total_cross_section_b(ep_lab_keV: np.ndarray) -> np.ndarray:
    ecm_MeV = LAB_TO_CM * ep_lab_keV / 1000.0
    return sigma_from_sfactor_mev(ecm_MeV, tentori2023_s_total(ecm_MeV))


def evaluated_channel_decomposition(
    sigma_165_model_b: np.ndarray,
    sigma_675_model_b: np.ndarray,
    sigma_eval_b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sigma_model_sum = sigma_165_model_b + sigma_675_model_b
    sigma_165_used = sigma_165_model_b.copy()
    sigma_675_used = sigma_675_model_b.copy()
    sigma_background = np.zeros_like(sigma_eval_b, dtype=float)

    no_total = sigma_eval_b <= 0.0
    has_background = (~no_total) & (sigma_model_sum <= sigma_eval_b)
    needs_scaling = (~no_total) & (sigma_model_sum > sigma_eval_b)

    sigma_165_used[no_total] = 0.0
    sigma_675_used[no_total] = 0.0
    sigma_background[has_background] = sigma_eval_b[has_background] - sigma_model_sum[has_background]

    scale = np.ones_like(sigma_eval_b, dtype=float)
    np.divide(sigma_eval_b, sigma_model_sum, out=scale, where=needs_scaling)
    sigma_165_used[needs_scaling] *= scale[needs_scaling]
    sigma_675_used[needs_scaling] *= scale[needs_scaling]

    sigma_used_sum = sigma_165_used + sigma_675_used + sigma_background
    return sigma_165_used, sigma_675_used, sigma_background, sigma_used_sum


def print_peak_summary(ep_lab_keV: np.ndarray) -> None:
    tentori_b = tentori2023_total_cross_section_b(ep_lab_keV)
    for use_coulomb in (False, True):
        s165, s675, stotal = components(ep_lab_keV, use_coulomb)
        s165_b = cm2_to_barn(s165)
        s675_b = cm2_to_barn(s675)
        stotal_b = cm2_to_barn(stotal)
        s165_used_b, s675_used_b, sbg_b, sused_sum_b = evaluated_channel_decomposition(s165_b, s675_b, tentori_b)
        tag = "with Coulomb" if use_coulomb else "constant widths"
        print(f"{tag} (physical, barns):")
        for name, curve in (
            ("raw 165", s165_b),
            ("raw 675", s675_b),
            ("raw total", stotal_b),
            ("Tentori", tentori_b),
            ("used 165", s165_used_b),
            ("used 675", s675_used_b),
            ("background", sbg_b),
            ("used sum", sused_sum_b),
        ):
            i = int(np.argmax(curve))
            print(f"  {name:>5s}: peak {curve[i]:.6g} b at Ep_lab = {ep_lab_keV[i]:.3f} keV")
        closure = np.nanmax(np.abs(sused_sum_b - tentori_b))
        print(f"  closure: max |used sum - Tentori| = {closure:.6e} b")


def plot(ep_lab_keV: np.ndarray, output: Path, yscale: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), constrained_layout=True)
    tentori_b = tentori2023_total_cross_section_b(ep_lab_keV)

    for ax, use_coulomb, title in (
        (axes[0], False, "Constant entrance width"),
        (axes[1], True, "Energy-dependent Coulomb penetrability"),
    ):
        s165, s675, stotal = components(ep_lab_keV, use_coulomb)
        s165_model_b = cm2_to_barn(s165)
        s675_model_b = cm2_to_barn(s675)
        stotal_model_b = cm2_to_barn(stotal)
        s165_used_b, s675_used_b, sbg_b, sused_sum_b = evaluated_channel_decomposition(
            s165_model_b,
            s675_model_b,
            tentori_b,
        )

        ax.plot(ep_lab_keV, tentori_b, color="black", lw=2.2, label="Tentori 2023 total")
        ax.plot(ep_lab_keV, sused_sum_b, color="#777777", lw=1.5, ls=":", label="used sum")
        ax.plot(ep_lab_keV, s165_used_b, color="#2764b8", lw=1.8, ls="--", label="used 165")
        ax.plot(ep_lab_keV, s675_used_b, color="#c43c39", lw=1.8, ls="-.", label="used 675")
        ax.plot(ep_lab_keV, sbg_b, color="#2f8f46", lw=1.9, label="background 3alpha")
        ax.plot(ep_lab_keV, stotal_model_b, color="#999999", lw=1.2, ls=(0, (4, 2)), label="raw model total")
        ax.axvline(CM_TO_LAB * H11B_165_RESONANCE_ECM_KEV, color="#2764b8", ls=":", lw=1.0)
        ax.axvline(CM_TO_LAB * H11B_675_RESONANCE_ECM_KEV, color="#c43c39", ls=":", lw=1.0)
        ax.set_title(title)
        ax.set_xlabel(r"proton lab energy $E_p$ (keV)")
        ax.set_yscale(yscale)
        if yscale == "linear":
            ax.set_ylim(0.0, max(1.6, 1.15 * np.nanmax(tentori_b)))
        else:
            ax.set_ylim(1e-5, max(30.0, 1.15 * np.nanmax(stotal_model_b)))
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, fontsize=8.2)

    axes[0].set_ylabel("physical cross section (barn)")
    fig.suptitle(r"$p + ^{11}$B evaluated total and phenomenological background decomposition")
    fig.savefig(output, dpi=220)
    print(f"Wrote {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--emin", type=float, default=1.0, help="minimum proton lab energy in keV")
    parser.add_argument("--emax", type=float, default=1000.0, help="maximum proton lab energy in keV")
    parser.add_argument("--npoints", type=int, default=5000, help="number of energy grid points")
    parser.add_argument("--output", type=Path, default=Path("cs_model/h11b_cross_section_model.png"))
    parser.add_argument("--yscale", choices=("linear", "log"), default="linear", help="y-axis scale")
    args = parser.parse_args()

    ep_lab_keV = np.linspace(args.emin, args.emax, args.npoints)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print_peak_summary(ep_lab_keV)
    plot(ep_lab_keV, args.output, args.yscale)


if __name__ == "__main__":
    main()
