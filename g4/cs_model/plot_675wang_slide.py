#!/usr/bin/env python3
"""675-keV resonance region, linear scale, for the slide.

Single figure: the Wang-2026 675-region S-factor component (as implemented
in H11BCrossSection::GetWang675CrossSection; 148-keV Lorentzian omitted)
over Ep(lab) = 300 - 1000 keV, compared with:

    Becker 1987     : H. W. Becker, et al., Z. Phys. A 327, 341 (1987).
                      Table 2 S factors (a0 + a1), c.m. energies.
    Mazzucconi 2025 : D. Mazzucconi, et al., Eur. Phys. J. A 61, 114 (2025).
                      Table 2, total approximated as 4*pi*dsigma/dOmega(60deg).

Data files (same directory): becker1987.csv, mazzucconi2025_table2.csv.
Output: p11b_cs_675wang_slide.png
"""

from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent

LAB_TO_CM = 11.0 / 12.0
EG_MEV = 22.589

LIN_LO, LIN_HI = 0.300, 1.000   # MeV lab
OUTPUT = "p11b_cs_675wang_slide.png"

C_WANG = "#0E7C6B"
C_BECKER = "#C77D2E"
C_MAZZ = "#990000"

plt.rcParams.update({
    "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "figure.dpi": 100,
})


def sigma_wang675_b(ep_lab_mev: np.ndarray) -> np.ndarray:
    """Wang 2026 (arXiv:2601.00241) piecewise S-factor -> sigma (barn).

    Mirrors H11BCrossSection::GetWang675CrossSection (S1 quadratic
    background without the 148-keV Lorentzian; S2 polynomial; S3
    constant + Lorentzians).
    """
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
    """Table 2 S factors (a0 + a1, c.m. keV) -> sigma via the paper's
    Gamow convention (2*pi*eta = 4.7528 E^-1/2, i.e. EG = 22.589 MeV)."""
    arr = load_csv(SCRIPT_DIR / "becker1987.csv")
    if arr is None:
        return None
    e_cm_mev = arr[:, 0] / 1000.0
    s_tot = arr[:, 1] + arr[:, 3]
    ds_tot = np.sqrt(arr[:, 2] ** 2 + arr[:, 4] ** 2)
    gamow = np.exp(-np.sqrt(EG_MEV / e_cm_mev))
    return {"ep": e_cm_mev / LAB_TO_CM,
            "sigma_b": s_tot / e_cm_mev * gamow,
            "dsigma_b": ds_tot / e_cm_mev * gamow}


def load_mazzucconi() -> dict | None:
    """Total ~ 4*pi * dsigma/dOmega(60deg lab)."""
    arr = load_csv(SCRIPT_DIR / "mazzucconi2025_table2.csv")
    if arr is None:
        return None
    return {"ep": arr[:, 0],
            "dep": arr[:, 1],
            "sigma_b": 4.0 * math.pi * arr[:, 6] * 1e-3,
            "dsigma_b": 4.0 * math.pi * arr[:, 7] * 1e-3}


def main() -> None:
    ep = np.arange(LIN_LO, LIN_HI + 1e-9, 0.0002)
    sig = sigma_wang675_b(ep)

    becker = load_becker()
    mazz = load_mazzucconi()

    plt.figure(figsize=(8.4, 5.2))
    plt.plot(ep * 1e3, sig, lw=2.4, color=C_WANG,
             label="675-keV resonance region (Wang 2026)")

    if becker is not None:
        m = (becker["ep"] >= LIN_LO) & (becker["ep"] <= LIN_HI)
        plt.errorbar(becker["ep"][m] * 1e3, becker["sigma_b"][m],
                     yerr=becker["dsigma_b"][m], fmt="s", ms=4,
                     color=C_BECKER, lw=0.9, capsize=2,
                     label="Becker 1987")
    if mazz is not None:
        m = (mazz["ep"] >= LIN_LO) & (mazz["ep"] <= LIN_HI)
        plt.errorbar(mazz["ep"][m] * 1e3, mazz["sigma_b"][m],
                     xerr=mazz["dep"][m] * 1e3, yerr=mazz["dsigma_b"][m],
                     fmt="o", ms=4.5, color=C_MAZZ, lw=0.9, capsize=2,
                     label=r"Mazzucconi 2025 ($4\pi\,d\sigma/d\Omega|_{60^\circ}$)")

    plt.axvline(675, color="#999999", lw=0.9, ls=":")
    plt.xlabel(r"$E_p$ lab (keV)")
    plt.ylabel(r"$\sigma$ (b)")
    plt.xlim(LIN_LO * 1e3, LIN_HI * 1e3)
    plt.ylim(bottom=0)
    plt.legend(frameon=False, fontsize=10, loc="upper left")
    plt.title("675-keV resonance region")
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / OUTPUT, dpi=200)
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
