#!/usr/bin/env python3
"""162-keV gamma-channel validation figures, 11B(p,gamma)12C.

Produces THREE separate pngs, each sized for a 2-row x 3-column grid
on a 16:8 PowerPoint slide (panel aspect 4:3):

    <stem>_spectrum.png : LAB gamma energy spectrum, split by branch
                          (gamma0 16.11 MeV; cascade 11.66 + 4.44 MeV),
                          with rest-frame line positions as dashed
                          references and an inset zoom on the 16.11 MeV
                          peak showing the Doppler box vs the sharp
                          rest-frame line.
    <stem>_gamma0_angular.png : gamma0 cos(theta) in the c.m. vs beam
                          axis; OFF (isotropic) vs ON (A1/A2) overlay
                          with the analytic W(theta) curve.
    <stem>_cascade_isotropy.png : cascade gamma1/gamma2 cos(theta_cm),
                          both expected flat (current model has no
                          gamma-gamma correlation -- stated approximation).

Gamma events are selected with gamma_branch > 0:
    gamma_branch == 1 : gamma0          (n_prompt_gammas = 1)
    gamma_branch == 2 : gamma1 cascade  (n_prompt_gammas = 2)

Energies stored in the tree are LAB energies (Doppler-shifted);
gamma_primary_energy_MeV is the emitter-rest-frame energy of the
primary gamma (sharp). Spectra are unweighted: the gamma bias factor
is common to both gamma sub-branches, so their relative ratio is
unaffected.

Usage:
    python3 plot_gamma_162.py ON.root                       # no off overlay
    python3 plot_gamma_162.py ON.root OFF.root              # angular on/off
    python3 plot_gamma_162.py ON.root OFF.root --stem g162
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

# ---- styling: panels sit at 1/3 slide width ----
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.1,
    "legend.fontsize": 12,
    "figure.dpi": 100,
})

FIGSIZE = (16 / 3, 8 / 2)   # 4:3 panel for a 2x3 grid on a 16:8 slide

FIELDS = [
    "gamma_branch", "n_prompt_gammas",
    "gamma1_energy", "gamma2_energy",
    "gamma_primary_energy_MeV",
    "gamma1_theta_cm", "gamma2_theta_cm",
    "cos_theta_gamma_cm",
]

C_G0   = "#0E7C6B"   # teal   = gamma0
C_CASC = "#C77D2E"   # orange = cascade
C_OFF  = "#C77D2E"   # orange = isotropic baseline
C_ON   = "#0E7C6B"   # teal   = physics on
C_REF  = "#1A2E35"   # dark   = analytic curve / reference lines

# default gamma0 A1/A2 (run.mac, Craig-1956-equivalent)
A1_DEF = -0.17757009345794392
A2_DEF = 0.13084112149532712


def load(path):
    with uproot.open(path) as f:
        return f["tr"].arrays(FIELDS, library="np")


def clean_cos(c):
    c = np.asarray(c, dtype=float)
    return c[np.isfinite(c) & (c >= -1.0) & (c <= 1.0)]


def finite(x):
    x = np.asarray(x, dtype=float)
    return x[np.isfinite(x) & (x > 0)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("on_file", help="run with gamma0 W(theta) ON")
    ap.add_argument("off_file", nargs="?", default=None,
                    help="optional run with gamma0 W(theta) OFF (isotropic)")
    ap.add_argument("--stem", default="gamma_162")
    ap.add_argument("--a1", type=float, default=A1_DEF)
    ap.add_argument("--a2", type=float, default=A2_DEF)
    args = ap.parse_args()

    a_on = load(args.on_file)
    a_off = load(args.off_file) if args.off_file else None

    gb = np.asarray(a_on["gamma_branch"])
    m0 = gb == 1     # gamma0
    mc = gb == 2     # cascade
    if m0.sum() + mc.sum() == 0:
        sys.exit("no gamma events (gamma_branch > 0) in the first file.\n"
                 "  Did you run validation_162_gamma.mac?")

    # =====================================================================
    # Figure 1: LAB energy spectrum by branch, rest-frame refs, inset zoom
    # =====================================================================
    e_g0  = finite(np.asarray(a_on["gamma1_energy"])[m0])
    e_c1  = finite(np.asarray(a_on["gamma1_energy"])[mc])
    e_c2  = finite(np.asarray(a_on["gamma2_energy"])[mc])
    e_casc = np.concatenate([e_c1, e_c2])

    # rest-frame reference energies from the tree (sharp, no Doppler)
    ref_g0 = np.nanmedian(np.asarray(a_on["gamma_primary_energy_MeV"])[m0]) \
        if m0.sum() else np.nan
    ref_c1 = np.nanmedian(np.asarray(a_on["gamma_primary_energy_MeV"])[mc]) \
        if mc.sum() else np.nan
    # 4.44 MeV line rest-frame energy incl. recoil: (M_i^2 - M_g^2)/(2 M_i)
    M_GS, DE = 11177.93, 4.4389   # MeV
    ref_c2 = ((M_GS + DE) ** 2 - M_GS ** 2) / (2.0 * (M_GS + DE))

    emax = 17.0
    ebins = np.linspace(0, emax, 1700)   # 10 keV bins

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if e_g0.size:
        ax.hist(e_g0, bins=ebins, histtype="step", lw=1.8, color=C_G0)
    if e_casc.size:
        ax.hist(e_casc, bins=ebins, histtype="step", lw=1.8, color=C_CASC)
    from matplotlib.lines import Line2D
    handles = []
    if e_g0.size:
        handles.append(Line2D([], [], color=C_G0, lw=1.8,
                              label=r"$\gamma_0$ (16.11 $\to$ g.s.)"))
    if e_casc.size:
        handles.append(Line2D([], [], color=C_CASC, lw=1.8,
                              label=r"cascade (11.66 + 4.44)"))
    ax.legend(handles=handles, loc="upper left", frameon=False)
    ax.set_xlabel(r"$\gamma$ energy, LAB (MeV)")
    ax.set_ylabel("counts")
    ax.set_xlim(0, emax)
    ax.set_yscale("log")
    ymax = max(np.histogram(e_casc, bins=ebins)[0].max() if e_casc.size else 1,
               np.histogram(e_g0,   bins=ebins)[0].max() if e_g0.size else 1)
    ax.set_ylim(top=ymax * 30)   # headroom so the legend clears the peaks
    ax.set_title(r"$^{11}$B(p,$\gamma$)$^{12}$C - 162 keV")

    # inset: zoom on the 16.11 MeV peak (Doppler box vs sharp line)
    if e_g0.size and np.isfinite(ref_g0):
        axin = ax.inset_axes([0.75, 0.55, 0.22, 0.24])
        half = 0.06   # MeV, ~2x the Doppler half-width
        zbins = np.linspace(ref_g0 - half, ref_g0 + half, 60)
        axin.hist(e_g0, bins=zbins, histtype="step", lw=1.0, color=C_G0)
        axin.axvline(ref_g0, color=C_REF, lw=0.8, ls="--")
        axin.set_xlim(ref_g0 - half, ref_g0 + half)
        axin.set_xticks([round(ref_g0 - 0.04, 2), round(ref_g0 + 0.04, 2)])
        axin.tick_params(labelsize=7)
        axin.set_yticks([])

    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_spectrum.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_spectrum.png")

    # =====================================================================
    # Figure 2: gamma0 angular, ON vs OFF, analytic W(theta)
    # =====================================================================
    abins = np.linspace(-1, 1, 41)
    actr  = 0.5 * (abins[:-1] + abins[1:])

    def ang_hist(c):
        h, _ = np.histogram(c, bins=abins, density=True)
        return h

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if a_off is not None:
        gb_off = np.asarray(a_off["gamma_branch"])
        c_off = clean_cos(np.asarray(a_off["cos_theta_gamma_cm"])[gb_off == 1])
        ax.step(actr, ang_hist(c_off), where="mid", lw=2.0, color=C_OFF,
                label="isotropic (off)")
    c_on = clean_cos(np.asarray(a_on["cos_theta_gamma_cm"])[m0])
    ax.step(actr, ang_hist(c_on), where="mid", lw=2.0, color=C_ON,
            label=r"$W(\theta)$ on")

    x = np.linspace(-1, 1, 400)
    P1 = x
    P2 = 0.5 * (3 * x ** 2 - 1)
    W = 1 + args.a1 * P1 + args.a2 * P2
    W = W / (_trapz(W, x) / 2.0) * 0.5
    ax.plot(x, W, "--", lw=1.4, color=C_REF,
            label=r"$1+A_1P_1+A_2P_2$")

    ax.axhline(0.5, color="#999999", lw=0.8, ls=":")
    ax.set_xlabel(r"$\cos\theta_{\gamma}$  (c.m., vs beam axis)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1); ax.set_ylim(0, 1)
    ax.legend(loc="lower right", frameon=False)
    ax.set_title(r"$\gamma_0$ angular distribution")

    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_gamma0_angular.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_gamma0_angular.png")

    # =====================================================================
    # Figure 3: cascade isotropy check (gamma1 and gamma2)
    # =====================================================================
    fig, ax = plt.subplots(figsize=FIGSIZE)
    c1 = clean_cos(np.cos(np.asarray(a_on["gamma1_theta_cm"])[mc]))
    c2 = clean_cos(np.cos(np.asarray(a_on["gamma2_theta_cm"])[mc]))
    if c1.size:
        ax.step(actr, ang_hist(c1), where="mid", lw=2.0, color=C_G0,
                label=r"$\gamma_1$ (11.66 MeV)")
    if c2.size:
        ax.step(actr, ang_hist(c2), where="mid", lw=2.0, color=C_CASC,
                label=r"$\gamma_2$ (4.44 MeV)")
    ax.axhline(0.5, color="#999999", lw=0.8, ls=":")
    ax.set_xlabel(r"$\cos\theta_{\gamma}$  (c.m.)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1); ax.set_ylim(0, 1)
    ax.legend(loc="lower right", frameon=False)
    ax.set_title("cascade angular (model: isotropic)")

    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_cascade_isotropy.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_cascade_isotropy.png")

    print(f"  gamma0 events : {int(m0.sum())}")
    print(f"  cascade events: {int(mc.sum())}")
    if m0.sum() + mc.sum() > 0:
        print(f"  gamma0 fraction: "
              f"{m0.sum() / (m0.sum() + mc.sum()):.3f}")


if __name__ == "__main__":
    main()
