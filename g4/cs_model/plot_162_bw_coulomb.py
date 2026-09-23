#!/usr/bin/env python3
from __future__ import annotations

import argparse
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

# fallback channel radius
CHANNEL_RADIUS_FM = 4.5
Z_PROJECTILE = 1.0
Z_TARGET = 5.0
ALPHA_FINE_STRUCTURE = 1.0 / 137.035999084

# plot range (c.m. energy)
ECM_MIN_KEV = 130.0
ECM_MAX_KEV = 165.0
ECM_STEP_KEV = 0.02

OUTPUT_PLOT = "p11b_cs_162_bw_coulomb.png"

# Near-resonance gas-target points only
BECKER_NEAR_RESONANCE = np.array(
    [
        (143.0, 23.0, 1.1, 940.0, 50.0),
        (145.0, 42.0, 2.0, 1620.0, 80.0),
        (146.5, 64.0, 3.0, 2390.0, 120.0),
        (148.0, 98.0, 5.0, 3290.0, 160.0),
        (150.0, 76.0, 4.0, 2600.0, 130.0),
        (155.0, 16.7, 0.8, 463.0, 23.0),
    ],
    dtype=float,
)


def becker_to_sigma_barn(table: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    e_cm_kev = table[:, 0]
    s_total = table[:, 1] + table[:, 3]
    s_total_err = np.hypot(table[:, 2], table[:, 4])

    e_mev = e_cm_kev / 1000.0
    denom = e_mev * np.exp(4.7528 / np.sqrt(e_mev))
    return e_cm_kev, s_total / denom, s_total_err / denom


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
    log_e = np.log(np.clip(ecm_kev, table_e[0], table_e[-1]))
    return np.interp(log_e, np.log(table_e), table_logp)


def _penetrability_l1_mpmath(ecm_kev: np.ndarray) -> np.ndarray:
    import mpmath as mp

    reduced_mass_kev = 1.0 * 11.0 / 12.0 * AMU_C2_KEV
    arr = np.atleast_1d(np.asarray(ecm_kev, dtype=float))
    out = np.empty_like(arr)

    for i, e in enumerate(arr):
        k = math.sqrt(2.0 * reduced_mass_kev * e) / HBARC_KEV_FM
        rho = k * CHANNEL_RADIUS_FM
        eta = (
            Z_PROJECTILE * Z_TARGET * ALPHA_FINE_STRUCTURE
            * math.sqrt(reduced_mass_kev / (2.0 * e))
        )
        f = float(mp.coulombf(ENTRANCE_ORBITAL_L, eta, rho))
        g = float(mp.coulombg(ENTRANCE_ORBITAL_L, eta, rho))
        out[i] = rho / (f * f + g * g)

    return out.reshape(np.asarray(ecm_kev, dtype=float).shape)


def penetrability_ratio(ecm_kev: np.ndarray, er_kev: float) -> np.ndarray:
    ecm_kev = np.asarray(ecm_kev, dtype=float)
    try:
        table_e, table_logp = load_p11b_tables()
        log_p_e = interp_log_p(ecm_kev, table_e, table_logp)
        log_p_r = interp_log_p(np.array([er_kev]), table_e, table_logp)[0]
        return np.exp(log_p_e - log_p_r)
    except (FileNotFoundError, RuntimeError, OSError):
        print(
            f"[warn] {COULOMB_DATA.name} not found - using mpmath fallback "
            f"(R = {CHANNEL_RADIUS_FM} fm)."
        )
        p_e = _penetrability_l1_mpmath(ecm_kev)
        p_r = _penetrability_l1_mpmath(np.array([er_kev]))[0]
        return p_e / p_r


def sigma_bw_barn(ecm_kev: np.ndarray, ratio: np.ndarray | float) -> np.ndarray:
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


def plot_162_bw_coulomb(output=None):
    ecm = np.arange(ECM_MIN_KEV, ECM_MAX_KEV + ECM_STEP_KEV, ECM_STEP_KEV)
    ratio = penetrability_ratio(ecm, RESONANCE_ECM_KEV)
    sigma_full = sigma_bw_barn(ecm, ratio)

    e_res, sig_res, err_res = becker_to_sigma_barn(BECKER_NEAR_RESONANCE)

    i_peak = int(np.argmax(sigma_full))
    print("162 keV BW near-resonance plot")
    print(f"Er = {RESONANCE_ECM_KEV:.2f} keV (c.m.)")
    print(f"peak sigma = {sigma_full[i_peak]:.4f} b at Ecm = {ecm[i_peak]:.2f} keV")
    print(f"near-resonance data points = {e_res.size}")

    fig, ax = plt.subplots(figsize=(7.2, 5.4))

    # only the Coulomb-suppressed BW curve
    ax.plot(
        ecm,
        sigma_full,
        "-",
        color="#0E7C6B",
        lw=2.6,
        label=r"Breit-Wigner with $\Gamma_p(E)=\Gamma_p(E_r)\,P_1(E)/P_1(E_r)$",
    )

    # only near-resonance points
    ax.errorbar(
        e_res,
        sig_res,
        yerr=err_res,
        fmt="o",
        ms=5.2,
        color="#1A2E35",
        linestyle="none",
        capsize=2.5,
        zorder=5,
        label="Becker et al. (1987), near resonance",
    )

    ax.axvline(
        RESONANCE_ECM_KEV,
        color="#1A2E35",
        lw=1.0,
        ls=":",
        alpha=0.7,
    )
    ax.annotate(
        rf"$E_r$ = {RESONANCE_ECM_KEV:.1f} keV",
        xy=(RESONANCE_ECM_KEV, 0.0),
        xycoords=("data", "axes fraction"),
        xytext=(6, 8),
        textcoords="offset points",
        fontsize=10,
        color="#1A2E35",
    )

    ax.set_xlim(ECM_MIN_KEV, ECM_MAX_KEV)
    ax.set_xlabel(r"$E_{\mathrm{c.m.}}$ (keV)", fontsize=12)
    ax.set_ylabel(r"$\sigma(p,\alpha)$ (barn)", fontsize=12)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9, loc="upper right")

    fig.tight_layout()

    output_path = Path(output) if output is not None else SCRIPT_DIR / OUTPUT_PLOT
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"wrote {output_path}")

    return fig


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot the 162-keV Breit-Wigner cross section near resonance."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SCRIPT_DIR / OUTPUT_PLOT,
        help="Output image filename",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively",
    )
    args = parser.parse_args()

    fig = plot_162_bw_coulomb(output=args.output)

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
