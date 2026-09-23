#!/usr/bin/env python3
"""Plot the current p+11B cross-section model (162-keV penetrability-corrected
Breit-Wigner + Wang-2026 S-factor 675-region component, as implemented in
H11BCrossSection::GetSigma162 / GetWang675CrossSection) against the two
experimental datasets used to anchor the Wang parameterization:

    Becker 1987     : H. W. Becker, et al., Z. Phys. A 327, 341 (1987).
                      EXFOR entry A0413. Total 11B(p,3a) cross section,
                      E_cm = 22 - 1100 keV.
    Mazzucconi 2025 : D. Mazzucconi, et al., Eur. Phys. J. A 61, 114 (2025).
                      Table 2, dsigma/dOmega at 60 deg (lab); total
                      approximated here as 4*pi * dsigma/dOmega(60deg)
                      (60 deg is near the P2 node at 54.7 deg, so the P2
                      anisotropy largely cancels; odd/higher terms survive).

Data files (same directory as this script):
    becker1987.csv         : Table 2 of the paper (per-channel S factors,
                             c.m. energies); loader sums a0+a1 and converts
                             to sigma via the paper's Gamow convention
    mazzucconi2025_table2.csv : Ep_MeV, dEp_MeV, a0_mbsr, da0, a1_mbsr, da1,
                             tot_mbsr, dtot  (already provided)

Outputs (cs_model conventions):
    p11b_cs_162bw_plus_675wang.log.png
    p11b_cs_162bw_plus_675wang.lin.png
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
EG_MEV = 22.589

AMU_C2_KEV = 931_494.10242
HBARC_KEV_FM = 197_326.9804
FM2_TO_CM2 = 1.0e-26
BARN_CM2 = 1.0e-24

ENERGY_MIN_MEV = 0.0
ENERGY_MAX_MEV = 1.000
ENERGY_STEP_MEV = 0.0002

H11B_162_RESONANCE_ECM_KEV = 147.95
H11B_162_TOTAL_WIDTH_KEV = 5.3
H11B_162_PROTON_WIDTH_KEV = 0.0215
H11B_162_ALPHA0_WIDTH_KEV = 0.26
H11B_162_ALPHA1_WIDTH_KEV = 5.0
H11B_162_SPIN_STAT_FACTOR = 5.0 / 8.0
H11B_162_ENTRANCE_ORBITAL_L = 1

OUTPUT_COMPARE_LOG = "p11b_cs_162bw_plus_675wang.log.png"
OUTPUT_COMPARE_LIN = "p11b_cs_162bw_plus_675wang.lin.png"

C_TOTAL = "#1A2E35"
C_BW = "#4472C4"
C_WANG = "#0E7C6B"
C_BECKER = "#C77D2E"
C_MAZZ = "#990000"


# --------------------------------------------------------------------------
# penetrability tables (parsed from the C++ source, cs_model convention)
# --------------------------------------------------------------------------
def extract_cpp_array(text: str, name: str) -> np.ndarray:
    pattern = rf"(?:static\s+)?const double {re.escape(name)}\[\] = \{{(.*?)\}};"
    match = re.search(pattern, text, flags=re.S)
    if not match:
        raise RuntimeError(f"Could not find array {name} in {COULOMB_DATA}")
    number_pattern = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
    return np.array([float(v) for v in re.findall(number_pattern, match.group(1))], dtype=float)


def load_p11b_penetrability_tables() -> tuple[np.ndarray, dict[int, np.ndarray]]:
    text = COULOMB_DATA.read_text()
    return (
        extract_cpp_array(text, "P11B_energy_keV"),
        {
            0: extract_cpp_array(text, "P11B_L0_logP"),
            1: extract_cpp_array(text, "P11B_L1_logP"),
        },
    )


def interpolate_logp(energy_keV, table_energy_keV, table_logp):
    clipped = np.clip(energy_keV, table_energy_keV[0], table_energy_keV[-1])
    return np.interp(np.log(clipped), np.log(table_energy_keV), table_logp)


def penetrability_ratio(energy_keV, reference_energy_keV, orbital_l,
                        table_energy_keV, table_logp_by_l):
    logp = table_logp_by_l[orbital_l]
    logp_e = interpolate_logp(energy_keV, table_energy_keV, logp)
    logp_ref = interpolate_logp(np.array([reference_energy_keV]),
                                table_energy_keV, logp)[0]
    return np.where(energy_keV > 0.0, np.exp(logp_e - logp_ref), 0.0)


# --------------------------------------------------------------------------
# model: 162 BW (mirrors H11BCrossSection::GetSigma162)
# --------------------------------------------------------------------------
def breit_wigner_sigma_cm2(ecm_keV, resonance_energy_keV, total_width_keV,
                           entrance_width_keV, exit_width_keV,
                           spin_stat_factor, orbital_l,
                           table_energy_keV, table_logp_by_l):
    sigma = np.zeros_like(ecm_keV, dtype=float)
    valid = ecm_keV > 0.0
    if not np.any(valid):
        return sigma
    reduced_mass_keV = 1.0 * 11.0 / 12.0 * AMU_C2_KEV
    gamma_entrance = entrance_width_keV * penetrability_ratio(
        ecm_keV, resonance_energy_keV, orbital_l,
        table_energy_keV, table_logp_by_l)
    missing = max(0.0, total_width_keV - entrance_width_keV - exit_width_keV)
    gamma_total = gamma_entrance + exit_width_keV + missing
    with np.errstate(divide="ignore", invalid="ignore"):
        pi_over_k2_fm2 = math.pi * HBARC_KEV_FM ** 2 / (2.0 * reduced_mass_keV * ecm_keV)
    denom = (ecm_keV - resonance_energy_keV) ** 2 + gamma_total ** 2 / 4.0
    bw = spin_stat_factor * gamma_entrance * exit_width_keV / denom
    sigma[valid] = pi_over_k2_fm2[valid] * bw[valid] * FM2_TO_CM2
    return sigma


def sigma_162_bw_b(ep_lab_mev, table_energy_keV, table_logp_by_l):
    ecm_keV = 1000.0 * LAB_TO_CM * ep_lab_mev
    sigma_cm2 = breit_wigner_sigma_cm2(
        ecm_keV,
        H11B_162_RESONANCE_ECM_KEV,
        H11B_162_TOTAL_WIDTH_KEV,
        H11B_162_PROTON_WIDTH_KEV,
        H11B_162_ALPHA0_WIDTH_KEV + H11B_162_ALPHA1_WIDTH_KEV,
        H11B_162_SPIN_STAT_FACTOR,
        H11B_162_ENTRANCE_ORBITAL_L,
        table_energy_keV,
        table_logp_by_l,
    )
    return np.where(ecm_keV >= 1.0, sigma_cm2 / BARN_CM2, 0.0)


# --------------------------------------------------------------------------
# model: Wang-2026 675-region component
# (mirrors H11BCrossSection::GetWang675CrossSection: S1 quadratic background
#  WITHOUT the 148-keV Lorentzian, S2 polynomial, S3 constant + Lorentzians)
# --------------------------------------------------------------------------
def sigma_wang675_b(ep_lab_mev: np.ndarray) -> np.ndarray:
    e_cm = LAB_TO_CM * np.asarray(ep_lab_mev, dtype=float)
    x_keV = 1000.0 * e_cm
    s = np.zeros_like(e_cm)

    m1 = e_cm <= 0.400
    s[m1] = 197.0 + 0.240 * x_keV[m1] + 2.31e-4 * x_keV[m1] ** 2

    m2 = (e_cm > 0.400) & (e_cm <= 0.700)
    x = (x_keV[m2] - 400.0) / 100.0
    s[m2] = 330.2 + 102.436 * x - 58.481 * x ** 2 + 0.0933 * x ** 5

    m3 = e_cm > 0.700
    x3 = x_keV[m3]
    s3 = np.full_like(x3, 0.209689)
    for amp, center, width in ((2.0235e6, 622.2, 99.6),
                               (4.0102e6, 1388.4, 449.9),
                               (1.3220e6, 2492.4, 238.6),
                               (4.9451e6, 3528.6, 398.5),
                               (4.3430e5, 4703.6, 152.5)):
        s3 += amp / ((x3 - center) ** 2 + width ** 2)
    s[m3] = s3

    sigma = np.zeros_like(e_cm)
    ok = e_cm > 0.0
    sigma[ok] = s[ok] / e_cm[ok] * np.exp(-np.sqrt(EG_MEV / e_cm[ok]))
    return sigma


# --------------------------------------------------------------------------
# experimental datasets
# --------------------------------------------------------------------------
def load_csv(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        rows.append([float(v) for v in line.split(",")])
    return np.array(rows) if rows else None


def load_becker() -> dict | None:
    """becker1987.csv: Table 2 S factors (per channel), c.m. energies.

    Columns: e_cm_keV, s_a0_MeVb, ds_a0, s_a1_MeVb, ds_a1, flag.
    Total S = S(a0) + S(a1); sigma = S/E * exp(-sqrt(EG/E)) with the
    same Gamow constant as the model (2*pi*eta = 4.7528 E^-1/2 in the
    paper, i.e. EG = 4.7528^2 = 22.589 MeV).
    """
    arr = load_csv(SCRIPT_DIR / "becker1987.csv")
    if arr is None:
        return None
    e_cm_mev = arr[:, 0] / 1000.0
    s_tot = arr[:, 1] + arr[:, 3]
    ds_tot = np.sqrt(arr[:, 2] ** 2 + arr[:, 4] ** 2)
    gamow = np.exp(-np.sqrt(EG_MEV / e_cm_mev))
    sigma_b = s_tot / e_cm_mev * gamow
    dsigma_b = ds_tot / e_cm_mev * gamow
    return {"ep": e_cm_mev / LAB_TO_CM,
            "sigma_b": sigma_b,
            "dsigma_b": dsigma_b}


def load_mazzucconi() -> dict | None:
    """mazzucconi2025_table2.csv: total ~ 4*pi * dsigma/dOmega(60deg)."""
    arr = load_csv(SCRIPT_DIR / "mazzucconi2025_table2.csv")
    if arr is None:
        return None
    return {"ep": arr[:, 0],
            "dep": arr[:, 1],
            "sigma_b": 4.0 * math.pi * arr[:, 6] * 1e-3,
            "dsigma_b": 4.0 * math.pi * arr[:, 7] * 1e-3}


# --------------------------------------------------------------------------
def main() -> None:
    ep = ENERGY_MIN_MEV + ENERGY_STEP_MEV * np.arange(
        int(round((ENERGY_MAX_MEV - ENERGY_MIN_MEV) / ENERGY_STEP_MEV)) + 1)

    table_energy_keV, table_logp_by_l = load_p11b_penetrability_tables()
    sig_bw = sigma_162_bw_b(ep, table_energy_keV, table_logp_by_l)
    sig_wang = sigma_wang675_b(ep)
    sig_total = sig_bw + sig_wang

    becker = load_becker()
    mazz = load_mazzucconi()

    plt.figure(figsize=(8.4, 5.2))
    plt.semilogy(ep * 1e3, sig_total, lw=2.4, color=C_TOTAL,
                 label="model total (162 BW + Wang-2026 675)")

    if becker is not None:
        m = (becker["ep"] >= ENERGY_MIN_MEV) & (becker["ep"] <= ENERGY_MAX_MEV)
        plt.errorbar(becker["ep"][m] * 1e3, becker["sigma_b"][m],
                     yerr=becker["dsigma_b"][m], fmt="s", ms=3.5,
                     color=C_BECKER, lw=0.9, capsize=2,
                     label="Becker 1987 (EXFOR A0413)")
    if mazz is not None:
        m = mazz["ep"] <= ENERGY_MAX_MEV
        plt.errorbar(mazz["ep"][m] * 1e3, mazz["sigma_b"][m],
                     xerr=mazz["dep"][m] * 1e3, yerr=mazz["dsigma_b"][m],
                     fmt="o", ms=4, color=C_MAZZ, lw=0.9, capsize=2,
                     label=r"Mazzucconi 2025 ($4\pi\,d\sigma/d\Omega|_{60^\circ}$)")

    plt.xlabel(r"$E_p$ lab (keV)")
    plt.ylabel(r"$\sigma$ (b)")
    plt.xlim(ENERGY_MIN_MEV * 1e3, ENERGY_MAX_MEV * 1e3)
    plt.ylim(1e-13, 4)
    plt.legend(frameon=False, fontsize=9, loc="lower right")
    plt.title("p+11B model vs Becker 1987 / Mazzucconi 2025")
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / OUTPUT_COMPARE_LOG, dpi=200)

    # ---------------- linear-scale figure: 675 resonance region ----------
    lin_lo, lin_hi = 0.0, 1.00   # MeV lab
    plt.figure(figsize=(8.4, 5.2))
    zoom = (ep >= lin_lo) & (ep <= lin_hi)
    plt.plot(ep[zoom] * 1e3, sig_total[zoom], lw=2.4, color=C_TOTAL,
             label="model total (162 BW + Wang-2026 675)")
    if becker is not None:
        m = (becker["ep"] >= lin_lo) & (becker["ep"] <= lin_hi)
        plt.errorbar(becker["ep"][m] * 1e3, becker["sigma_b"][m],
                     yerr=becker["dsigma_b"][m], fmt="s", ms=4,
                     color=C_BECKER, lw=0.9, capsize=2,
                     label="Becker 1987")
    if mazz is not None:
        m = (mazz["ep"] >= lin_lo) & (mazz["ep"] <= lin_hi)
        plt.errorbar(mazz["ep"][m] * 1e3, mazz["sigma_b"][m],
                     xerr=mazz["dep"][m] * 1e3, yerr=mazz["dsigma_b"][m],
                     fmt="o", ms=4.5, color=C_MAZZ, lw=0.9, capsize=2,
                     label=r"Mazzucconi 2025 ($4\pi\,d\sigma/d\Omega|_{60^\circ}$)")
    plt.axvline(675, color="#999999", lw=0.9, ls=":")
    plt.xlabel(r"$E_p$ lab (keV)")
    plt.ylabel(r"$\sigma$ (b)")
    plt.xlim(lin_lo * 1e3, lin_hi * 1e3)
    plt.ylim(bottom=0)
    plt.legend(frameon=False, fontsize=9, loc="upper left")
    plt.title("675-keV resonance region (linear)")
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / OUTPUT_COMPARE_LIN, dpi=200)

    # ---------------- summary ----------------
    i_peak = np.argmax(sig_total[ep > 0.4])
    ep_peak = ep[ep > 0.4][i_peak]
    lines = [
        "162 BW + Wang-2026 675-region S-factor model",
        "",
        "Model components:",
        "  162 BW : penetrability-corrected single-level BW "
        f"(Er={H11B_162_RESONANCE_ECM_KEV} keV cm, Gtot={H11B_162_TOTAL_WIDTH_KEV} keV)",
        "  675    : Wang 2026 (arXiv:2601.00241) piecewise S-factor,",
        "           148-keV Lorentzian omitted (handled by the BW)",
        "",
        f"675-region peak: {ep_peak*1e3:.1f} keV lab, "
        f"{sig_total[ep > 0.4][i_peak]:.4f} b",
        f"model at Ep=675 keV lab: "
        f"{np.interp(0.675, ep, sig_total):.4f} b",
    ]
    for name, d in (("Becker 1987", becker), ("Mazzucconi 2025", mazz)):
        if d is None:
            lines.append(f"{name}: data file not found -- overlay skipped")
            continue
        m = (d["ep"] >= ENERGY_MIN_MEV) & (d["ep"] <= ENERGY_MAX_MEV)
        model_at = np.interp(d["ep"][m], ep, sig_total)
        resid = (d["sigma_b"][m] - model_at) / d["sigma_b"][m]
        lines.append(f"{name}: {int(m.sum())} points <= 1 MeV, "
                     f"mean (data-model)/data = {np.mean(resid)*100:+.1f} %, "
                     f"rms = {np.std(resid)*100:.1f} %")

    print("\n".join(lines))
    print(f"\nwrote {OUTPUT_COMPARE_LOG}")
    print(f"wrote {OUTPUT_COMPARE_LIN}")


if __name__ == "__main__":
    main()
