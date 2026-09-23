#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_p11b_model_vs_tentori.py

Compare the p + 11B cross-section model used in the Geant4 simulation
with the Tentori & Belloni 2023 total reference cross-section curve.

Model/decomposition curves:
  - raw model 165 keV resonance
  - raw model 675 keV resonance
  - raw model total = 165 + 675
  - used 165 and 675 components after normalization to the evaluated total
  - phenomenological background 3alpha = Tentori total - used 165 - used 675

Literature curve, 1 line:
  - Tentori & Belloni 2023 total p-11B fusion cross section

Default:
  - energy axis: proton lab energy Ep_lab in MeV
  - energy range: 0--1 MeV
  - simulation model: energy-dependent Coulomb penetrability entrance width
  - output two figures: linear and log y scales

Run inside the Geant4 project directory, or from cs_model/:
  python3 plot_p11b_model_vs_tentori.py

Optional:
  python3 plot_p11b_model_vs_tentori.py --constant-width
  python3 plot_p11b_model_vs_tentori.py --emax 1000
  python3 plot_p11b_model_vs_tentori.py --coulomb-data ../src/CoulombPenetrabilityData.cc

Outputs:
  p11b_model_vs_tentori_linear.png
  p11b_model_vs_tentori_log.png
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


# ----------------------------------------------------------------------
# Simulation model constants.
# Same constants as plot_h11b_cross_sections.py / H11BCrossSection.cc model.
# ----------------------------------------------------------------------
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

# Gamow energy used by Tentori & Belloni 2023 for p + 11B.
EG_MEV = 22.589


# ----------------------------------------------------------------------
# Coulomb penetrability table utilities.
# ----------------------------------------------------------------------
def find_coulomb_data(user_path: str | None) -> Path:
    """
    Find src/CoulombPenetrabilityData.cc robustly.

    This allows the script to be run either from the project root or from
    cs_model/.
    """
    if user_path:
        path = Path(user_path).expanduser().resolve()
        if path.exists():
            return path
        raise FileNotFoundError(f"--coulomb-data does not exist: {path}")

    script_dir = Path(__file__).resolve().parent
    cwd = Path.cwd().resolve()

    candidates = [
        cwd / "src" / "CoulombPenetrabilityData.cc",
        cwd.parent / "src" / "CoulombPenetrabilityData.cc",
        script_dir / "src" / "CoulombPenetrabilityData.cc",
        script_dir.parent / "src" / "CoulombPenetrabilityData.cc",
        script_dir.parent.parent / "src" / "CoulombPenetrabilityData.cc",
    ]

    for path in candidates:
        if path.exists():
            return path

    searched = "\n  ".join(str(p) for p in candidates)
    raise FileNotFoundError(
        "Could not find CoulombPenetrabilityData.cc. Searched:\n"
        f"  {searched}\n"
        "Run from the Geant4 project root/cs_model directory, or pass:\n"
        "  --coulomb-data /path/to/src/CoulombPenetrabilityData.cc"
    )


def extract_cpp_array(text: str, name: str, source: Path) -> np.ndarray:
    pattern = rf"(?:static\s+)?const double {re.escape(name)}\[\] = \{{(.*?)\}};"
    match = re.search(pattern, text, flags=re.S)
    if not match:
        raise ValueError(f"Could not find array {name} in {source}")

    number_pattern = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
    return np.array([float(value) for value in re.findall(number_pattern, match.group(1))], dtype=float)


def load_p11b_penetrability_tables(coulomb_data: Path) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    text = coulomb_data.read_text()
    energy = extract_cpp_array(text, "P11B_energy_keV", coulomb_data)
    logp_l0 = extract_cpp_array(text, "P11B_L0_logP", coulomb_data)
    logp_l1 = extract_cpp_array(text, "P11B_L1_logP", coulomb_data)
    return energy, {0: logp_l0, 1: logp_l1}


def interpolate_logp(
    energy_keV: np.ndarray,
    table_energy_keV: np.ndarray,
    table_logp: np.ndarray,
) -> np.ndarray:
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

    table_logp = table_logp_by_l[orbital_l]
    logp_e = interpolate_logp(energy_keV, table_energy_keV, table_logp)
    logp_ref = interpolate_logp(np.array([reference_energy_keV]), table_energy_keV, table_logp)[0]
    ratio = np.exp(logp_e - logp_ref)
    return np.where(energy_keV > 0.0, ratio, 0.0)


# ----------------------------------------------------------------------
# Simulation Breit-Wigner model.
# ----------------------------------------------------------------------
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

    projectile_a = 1.0
    target_a = 11.0
    reduced_mass_keV = projectile_a * target_a / (projectile_a + target_a) * AMU_C2_KEV

    gamma_entrance = np.full_like(ecm_keV, entrance_width_at_resonance_keV, dtype=float)
    if use_coulomb:
        gamma_entrance *= penetrability_ratio(
            ecm_keV,
            resonance_energy_keV,
            entrance_orbital_l,
            table_energy_keV,
            table_logp_by_l,
        )

    missing_width = max(
        0.0,
        total_width_keV - entrance_width_at_resonance_keV - exit_width_at_resonance_keV,
    )
    gamma_total = gamma_entrance + exit_width_at_resonance_keV + missing_width

    pi_over_k2_fm2 = math.pi * HBARC_KEV_FM * HBARC_KEV_FM / (
        2.0 * reduced_mass_keV * ecm_keV
    )
    denominator = (ecm_keV - resonance_energy_keV) ** 2 + gamma_total * gamma_total / 4.0
    bw_factor = spin_stat_factor * gamma_entrance * exit_width_at_resonance_keV / denominator

    sigma[valid] = pi_over_k2_fm2[valid] * bw_factor[valid] * FM2_TO_CM2
    return sigma


def model_components(
    ep_lab_keV: np.ndarray,
    use_coulomb: bool,
    table_energy_keV: np.ndarray,
    table_logp_by_l: dict[int, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return model 165, 675, and total cross sections in barn.
    """
    ecm_keV = LAB_TO_CM * ep_lab_keV

    sigma_165_cm2 = breit_wigner_sigma_cm2(
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

    sigma_675_cm2 = breit_wigner_sigma_cm2(
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

    # Same cutoffs as the original simulation plotting script.
    sigma_165_cm2 = np.where((ecm_keV >= 1.0) & (ecm_keV <= 400.0), sigma_165_cm2, 0.0)
    sigma_675_cm2 = np.where((ecm_keV >= 1.0) & (ecm_keV <= 3500.0), sigma_675_cm2, 0.0)

    sigma_165_b = sigma_165_cm2 / BARN_CM2
    sigma_675_b = sigma_675_cm2 / BARN_CM2
    sigma_total_b = sigma_165_b + sigma_675_b

    return sigma_165_b, sigma_675_b, sigma_total_b


# ----------------------------------------------------------------------
# Tentori & Belloni 2023 total reference curve.
# ----------------------------------------------------------------------
def sigma_from_sfactor_mev(E_cm_MeV: np.ndarray, S_MeV_b: np.ndarray) -> np.ndarray:
    """
    Convert S(E) in MeV b to sigma in barn.
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    S = np.asarray(S_MeV_b, dtype=float)

    sigma = np.full_like(E, np.nan, dtype=float)
    valid = E > 0.0
    sigma[valid] = S[valid] / E[valid] * np.exp(-np.sqrt(EG_MEV / E[valid]))
    return sigma


def tentori2023_s_total(E_cm_MeV: np.ndarray) -> np.ndarray:
    """
    Tentori & Belloni 2023 total reference S(E) parameterization.

    The low-energy 148-keV resonance term uses the Nevins-Swain parameters.
    """
    E = np.asarray(E_cm_MeV, dtype=float)
    E_keV = 1000.0 * E
    S = np.full_like(E, np.nan, dtype=float)

    # E <= 0.400 MeV.
    C0 = 197.0
    C1 = 0.269
    C2 = 2.54e-4

    # Nevins-Swain 148 keV resonance term.
    AL = 1.82e4
    EL = 148.0
    dEL = 2.35

    # 0.400 < E <= 0.668 MeV.
    E1_break = 0.400
    E2_break = 0.668
    D0 = 346.0
    D1 = 150.0
    D2 = -59.9
    D5 = -0.460

    # E > 0.668 MeV.
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


def tentori2023_total_cross_section(ep_lab_keV: np.ndarray) -> np.ndarray:
    """
    Return Tentori 2023 total cross section in barn on an Ep_lab grid.
    """
    ecm_MeV = LAB_TO_CM * ep_lab_keV / 1000.0
    return sigma_from_sfactor_mev(ecm_MeV, tentori2023_s_total(ecm_MeV))


def evaluated_channel_decomposition(
    sigma_165_model_b: np.ndarray,
    sigma_675_model_b: np.ndarray,
    sigma_eval_b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Apply the same channel-normalization logic as H11BCrossSection.cc.

    Returns:
        sigma_165_used, sigma_675_used, sigma_background, sigma_used_sum
        all in barn.
    """
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


# ----------------------------------------------------------------------
# Plotting.
# ----------------------------------------------------------------------
def print_peak_summary(
    ep_lab_keV: np.ndarray,
    s165_model_b,
    s675_model_b,
    stotal_model_b,
    tentori_b,
    s165_used_b,
    s675_used_b,
    sbg_b,
    sused_sum_b,
    tag: str,
) -> None:
    print(f"{tag}:")
    for name, curve in (
        ("raw model 165", s165_model_b),
        ("raw model 675", s675_model_b),
        ("raw model total", stotal_model_b),
        ("Tentori total", tentori_b),
        ("used 165", s165_used_b),
        ("used 675", s675_used_b),
        ("background", sbg_b),
        ("used sum", sused_sum_b),
    ):
        i = int(np.nanargmax(curve))
        print(f"  {name:>13s}: peak {curve[i]:.6g} b at Ep_lab = {ep_lab_keV[i] / 1000.0:.6f} MeV")

    closure = np.nanmax(np.abs(sused_sum_b - tentori_b))
    print(f"  {'closure':>13s}: max |used sum - Tentori| = {closure:.6e} b")


def draw_figure(
    ep_lab_MeV: np.ndarray,
    s165_model_b: np.ndarray,
    s675_model_b: np.ndarray,
    stotal_model_b: np.ndarray,
    tentori_b: np.ndarray,
    s165_used_b: np.ndarray,
    s675_used_b: np.ndarray,
    sbg_b: np.ndarray,
    sused_sum_b: np.ndarray,
    yscale: str,
    output: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 6.2))

    ax.plot(ep_lab_MeV, tentori_b, color="black", linewidth=2.5, label="Tentori 2023 evaluated total")
    ax.plot(ep_lab_MeV, sused_sum_b, color="#777777", linewidth=1.7, linestyle=":", label="used sum")
    ax.plot(ep_lab_MeV, s165_used_b, color="#2764b8", linewidth=2.0, linestyle="--", label="used 165 component")
    ax.plot(ep_lab_MeV, s675_used_b, color="#c43c39", linewidth=2.0, linestyle="-.", label="used 675 component")
    ax.plot(ep_lab_MeV, sbg_b, color="#2f8f46", linewidth=2.1, label="phenom. background 3alpha")
    ax.plot(ep_lab_MeV, stotal_model_b, color="#999999", linewidth=1.4, linestyle=(0, (4, 2)), label="raw model total")

    # Resonance positions in lab energy.
    ax.axvline(CM_TO_LAB * H11B_165_RESONANCE_ECM_KEV / 1000.0, linestyle=":", linewidth=1.0)
    ax.axvline(CM_TO_LAB * H11B_675_RESONANCE_ECM_KEV / 1000.0, linestyle=":", linewidth=1.0)

    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"proton lab energy $E_p$ (MeV)")
    ax.set_ylabel(r"cross section $\sigma$ (barn)")

    ymax = max(np.nanmax(stotal_model_b), np.nanmax(tentori_b), np.nanmax(sbg_b))

    if yscale == "linear":
        ax.set_yscale("linear")
        ax.set_ylim(0.0, max(1.6, 1.15 * ymax))
        suffix = "linear scale"
    elif yscale == "log":
        ax.set_yscale("log")
        ax.set_ylim(1e-5, max(3.0, 1.15 * ymax))
        suffix = "log scale"
    else:
        raise ValueError(f"Unknown yscale: {yscale}")

    ax.set_title(r"$p+^{11}$B evaluated total and phenomenological background, $E_p \leq 1$ MeV, " + suffix)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8.7, ncol=2)

    fig.tight_layout()
    fig.savefig(output, dpi=300)
    print(f"Saved: {output}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--emin", type=float, default=1.0, help="minimum proton lab energy in keV")
    parser.add_argument("--emax", type=float, default=1000.0, help="maximum proton lab energy in keV")
    parser.add_argument("--npoints", type=int, default=5000, help="number of energy grid points")
    parser.add_argument(
        "--constant-width",
        action="store_true",
        help="use constant entrance widths instead of energy-dependent Coulomb penetrability",
    )
    parser.add_argument(
        "--coulomb-data",
        type=str,
        default=None,
        help="path to src/CoulombPenetrabilityData.cc",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="p11b_model_vs_tentori",
        help="prefix for output PNG files",
    )
    args = parser.parse_args()

    coulomb_data = find_coulomb_data(args.coulomb_data)
    table_energy_keV, table_logp_by_l = load_p11b_penetrability_tables(coulomb_data)

    use_coulomb = not args.constant_width
    model_tag = "energy-dependent Coulomb penetrability" if use_coulomb else "constant entrance widths"

    ep_lab_keV = np.linspace(args.emin, args.emax, args.npoints)
    ep_lab_MeV = ep_lab_keV / 1000.0

    s165_model_b, s675_model_b, stotal_model_b = model_components(
        ep_lab_keV,
        use_coulomb,
        table_energy_keV,
        table_logp_by_l,
    )

    tentori_b = tentori2023_total_cross_section(ep_lab_keV)
    s165_used_b, s675_used_b, sbg_b, sused_sum_b = evaluated_channel_decomposition(
        s165_model_b,
        s675_model_b,
        tentori_b,
    )

    print(f"Using Coulomb data: {coulomb_data}")
    print_peak_summary(
        ep_lab_keV,
        s165_model_b,
        s675_model_b,
        stotal_model_b,
        tentori_b,
        s165_used_b,
        s675_used_b,
        sbg_b,
        sused_sum_b,
        f"Model ({model_tag})",
    )

    draw_figure(
        ep_lab_MeV,
        s165_model_b,
        s675_model_b,
        stotal_model_b,
        tentori_b,
        s165_used_b,
        s675_used_b,
        sbg_b,
        sused_sum_b,
        yscale="linear",
        output=Path(f"{args.output_prefix}_linear.png"),
    )

    draw_figure(
        ep_lab_MeV,
        s165_model_b,
        s675_model_b,
        stotal_model_b,
        tentori_b,
        s165_used_b,
        s675_used_b,
        sbg_b,
        sused_sum_b,
        yscale="log",
        output=Path(f"{args.output_prefix}_log.png"),
    )


if __name__ == "__main__":
    main()
