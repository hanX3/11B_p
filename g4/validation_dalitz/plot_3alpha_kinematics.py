#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot

# These constants are used only for diagnostics and plot reference lines.
# ROOT energies and angles come from Geant4; exact lab-momentum closure would
# require storing final-state four-vectors or the exact Geant4 ion masses.
PROTON_MASS_MEV = 938.2720813
ALPHA_MASS_MEV = 3727.3794066
REFERENCE_Q_P11B_TO_3ALPHA_MEV = 8.680
BE8_GS_RELATIVE_ENERGY_MEV = 0.09184
BE8_2PLUS_RELATIVE_ENERGY_MEV = 3.03

REQUIRED_BRANCHES = {
    "reaction_channel",
    "branch_id",
    "e_cm_p11B",
    "e_alpha1",
    "e_alpha2",
    "e_alpha3",
    "theta_lab_alpha1",
    "theta_lab_alpha2",
    "theta_lab_alpha3",
    "phi_lab_alpha1",
    "phi_lab_alpha2",
    "phi_lab_alpha3",
    "e_3alpha_cm_alpha1",
    "e_3alpha_cm_alpha2",
    "e_3alpha_cm_alpha3",
    "e_alpha1_cm",
    "e_alpha2_cm",
    "e_alpha3_cm",
    "opening_angle_alpha12_cm",
    "opening_angle_alpha13_cm",
    "opening_angle_alpha23_cm",
    "ex_8Be",
    "e_8be_excitation",
    "event_weight",
}


@dataclass
class Sample:
    label: str
    root_file: Path
    channel: int
    tree_name: str
    arrays: dict


# -----------------------------------------------------------------------------
# Loading and basic utilities
# -----------------------------------------------------------------------------

def sanitize_label(label: str) -> str:
    out = []
    for ch in label:
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out)


def resolve_tree(root_file: Path, requested_tree: str | None) -> str:
    with uproot.open(root_file) as f:
        keys = {key.split(";")[0]: key for key in f.keys()}
        if requested_tree is not None:
            if requested_tree in keys:
                return requested_tree
            raise KeyError(f"Requested tree '{requested_tree}' not found in {root_file}. Available keys: {list(keys)}")
        if "tr" in keys:
            return "tr"
        if "reaction" in keys:
            return "reaction"
        raise KeyError(f"No supported reaction tree found in {root_file}. Available keys: {list(keys)}")


def load_sample(spec: str, tree_name: str | None) -> Sample:
    parts = spec.split(":")
    if len(parts) < 3:
        raise ValueError("--sample must have format label:root_file:reaction_channel")
    label = parts[0]
    channel = int(parts[-1])
    root_file = Path(":".join(parts[1:-1])).expanduser()
    if not root_file.exists():
        raise FileNotFoundError(f"ROOT file not found: {root_file}")

    resolved_tree = resolve_tree(root_file, tree_name)
    with uproot.open(root_file) as f:
        tree = f[resolved_tree]
        available = set(tree.keys())
        missing = sorted(REQUIRED_BRANCHES - available)
        if missing:
            raise KeyError(
                f"Missing required branch(es) in {root_file}: {missing}\n"
                f"Available branches:\n" + "\n".join(sorted(available))
            )
        arrays = tree.arrays(sorted(REQUIRED_BRANCHES), library="np")

    mask = np.asarray(arrays["reaction_channel"]) == channel
    arrays = {key: np.asarray(value)[mask] for key, value in arrays.items()}
    n = len(arrays["reaction_channel"])
    if n == 0:
        raise RuntimeError(f"Sample {label} has zero entries after reaction_channel=={channel} selection.")
    return Sample(label=label, root_file=root_file, channel=channel, tree_name=resolved_tree, arrays=arrays)


def make_filtered_sample(sample: Sample, label: str, mask: np.ndarray) -> Sample:
    mask = np.asarray(mask, dtype=bool)
    arrays = {key: np.asarray(value)[mask] for key, value in sample.arrays.items()}
    if len(arrays["reaction_channel"]) == 0:
        raise RuntimeError(f"Filtered sample {label} has zero entries.")
    return Sample(label=label, root_file=sample.root_file, channel=sample.channel, tree_name=sample.tree_name, arrays=arrays)


def finite(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def finite_mask(*values: np.ndarray) -> np.ndarray:
    mask = np.ones(len(values[0]), dtype=bool)
    for value in values:
        mask &= np.isfinite(np.asarray(value, dtype=float))
    return mask


def branch_ids(sample: Sample) -> np.ndarray:
    return np.unique(np.asarray(sample.arrays["branch_id"], dtype=int))


# -----------------------------------------------------------------------------
# Energies, Dalitz coordinates, projections, and Eij
# -----------------------------------------------------------------------------

def cm_energies(sample: Sample) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a = sample.arrays
    return finite(a["e_3alpha_cm_alpha1"]), finite(a["e_3alpha_cm_alpha2"]), finite(a["e_3alpha_cm_alpha3"])


def energy_triplet(sample: Sample) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a = sample.arrays
    e1 = np.asarray(a["e_3alpha_cm_alpha1"], dtype=float)
    e2 = np.asarray(a["e_3alpha_cm_alpha2"], dtype=float)
    e3 = np.asarray(a["e_3alpha_cm_alpha3"], dtype=float)
    esum = e1 + e2 + e3
    mask = np.isfinite(e1) & np.isfinite(e2) & np.isfinite(e3) & np.isfinite(esum) & (esum > 0.0)
    return e1[mask], e2[mask], e3[mask]


def energy_sum_cm(sample: Sample) -> np.ndarray:
    e1, e2, e3 = energy_triplet(sample)
    return e1 + e2 + e3


def q_eff(sample: Sample) -> np.ndarray:
    a = sample.arrays
    e1 = np.asarray(a["e_3alpha_cm_alpha1"], dtype=float)
    e2 = np.asarray(a["e_3alpha_cm_alpha2"], dtype=float)
    e3 = np.asarray(a["e_3alpha_cm_alpha3"], dtype=float)
    ecm = np.asarray(a["e_cm_p11B"], dtype=float)
    esum = e1 + e2 + e3
    mask = np.isfinite(esum) & np.isfinite(ecm)
    values = esum[mask] - ecm[mask]
    return values[np.isfinite(values)]


def dalitz_xy_from_triplet(e1: np.ndarray, e2: np.ndarray, e3: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    esum = e1 + e2 + e3
    mask = np.isfinite(e1) & np.isfinite(e2) & np.isfinite(e3) & np.isfinite(esum) & (esum > 0.0)
    e1, e2, e3, esum = e1[mask], e2[mask], e3[mask], esum[mask]
    x = np.sqrt(3.0) * (e2 - e3) / esum
    y = (2.0 * e1 - e2 - e3) / esum
    return x, y


def dalitz_xy(sample: Sample, symmetrized: bool = False) -> tuple[np.ndarray, np.ndarray]:
    e1, e2, e3 = energy_triplet(sample)
    if not symmetrized:
        return dalitz_xy_from_triplet(e1, e2, e3)

    xs = []
    ys = []
    # Six permutations. This is a visualization symmetrization, not a coherent
    # amplitude-level Bose symmetrization.
    for a, b, c in ((e1, e2, e3), (e1, e3, e2), (e2, e1, e3), (e2, e3, e1), (e3, e1, e2), (e3, e2, e1)):
        x, y = dalitz_xy_from_triplet(a, b, c)
        xs.append(x)
        ys.append(y)
    return np.concatenate(xs), np.concatenate(ys)


def dalitz_rho_phi(sample: Sample, symmetrized: bool = True) -> tuple[np.ndarray, np.ndarray]:
    x, y = dalitz_xy(sample, symmetrized=symmetrized)
    rho = np.sqrt(x * x + y * y)
    phi = np.degrees(np.arctan2(y, x))
    # Fold to 0--60 deg to use the sixfold symmetry of three identical alphas.
    phi_folded = np.mod(phi, 60.0)
    mask = np.isfinite(rho) & np.isfinite(phi_folded)
    return rho[mask], phi_folded[mask]


def alpha_alpha_relative_energies(sample: Sample) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return nonrelativistic E12, E13, E23 from CM kinetic energies and opening angles.

    For equal masses, Eij = (Ti + Tj - 2 sqrt(Ti Tj) cos(theta_ij)) / 2.
    This uses generator-level CM angles stored in the ROOT tree.
    """
    a = sample.arrays
    e1 = np.asarray(a["e_3alpha_cm_alpha1"], dtype=float)
    e2 = np.asarray(a["e_3alpha_cm_alpha2"], dtype=float)
    e3 = np.asarray(a["e_3alpha_cm_alpha3"], dtype=float)
    th12 = np.asarray(a["opening_angle_alpha12_cm"], dtype=float)
    th13 = np.asarray(a["opening_angle_alpha13_cm"], dtype=float)
    th23 = np.asarray(a["opening_angle_alpha23_cm"], dtype=float)
    mask = finite_mask(e1, e2, e3, th12, th13, th23) & (e1 >= 0.0) & (e2 >= 0.0) & (e3 >= 0.0)
    e1, e2, e3, th12, th13, th23 = [v[mask] for v in (e1, e2, e3, th12, th13, th23)]

    def eij(ti: np.ndarray, tj: np.ndarray, theta: np.ndarray) -> np.ndarray:
        values = 0.5 * (ti + tj - 2.0 * np.sqrt(np.maximum(0.0, ti * tj)) * np.cos(theta))
        return values[np.isfinite(values)]

    return eij(e1, e2, th12), eij(e1, e3, th13), eij(e2, e3, th23)


# -----------------------------------------------------------------------------
# Lab momentum closure diagnostic
# -----------------------------------------------------------------------------

def momentum_diagnostics(sample: Sample) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a = sample.arrays
    t1 = np.asarray(a["e_alpha1"], dtype=float)
    t2 = np.asarray(a["e_alpha2"], dtype=float)
    t3 = np.asarray(a["e_alpha3"], dtype=float)
    th1 = np.asarray(a["theta_lab_alpha1"], dtype=float)
    th2 = np.asarray(a["theta_lab_alpha2"], dtype=float)
    th3 = np.asarray(a["theta_lab_alpha3"], dtype=float)
    ph1 = np.asarray(a["phi_lab_alpha1"], dtype=float)
    ph2 = np.asarray(a["phi_lab_alpha2"], dtype=float)
    ph3 = np.asarray(a["phi_lab_alpha3"], dtype=float)
    ecm = np.asarray(a["e_cm_p11B"], dtype=float)

    mask = finite_mask(t1, t2, t3, th1, th2, th3, ph1, ph2, ph3, ecm)
    t1, t2, t3, th1, th2, th3, ph1, ph2, ph3, ecm = [v[mask] for v in (t1, t2, t3, th1, th2, th3, ph1, ph2, ph3, ecm)]

    def p_from_t(t):
        return np.sqrt(np.maximum(0.0, t * t + 2.0 * ALPHA_MASS_MEV * t))

    def components(p, theta, phi):
        st = np.sin(theta)
        return p * st * np.cos(phi), p * st * np.sin(phi), p * np.cos(theta)

    p1, p2, p3 = p_from_t(t1), p_from_t(t2), p_from_t(t3)
    px1, py1, pz1 = components(p1, th1, ph1)
    px2, py2, pz2 = components(p2, th2, ph2)
    px3, py3, pz3 = components(p3, th3, ph3)

    px = px1 + px2 + px3
    py = py1 + py2 + py3
    pz = pz1 + pz2 + pz3
    pt = np.sqrt(px * px + py * py)

    # Approximate lab proton energy from E_cm using nonrelativistic two-body relation.
    # This is only for a convenient scale in the diagnostic plot.
    lab_energy_est = ecm * 12.0 / 11.0
    p_initial = np.sqrt(np.maximum(0.0, lab_energy_est * lab_energy_est + 2.0 * PROTON_MASS_MEV * lab_energy_est))
    delta_pz = pz - p_initial
    relative_delta_pz = np.divide(delta_pz, p_initial, out=np.full_like(delta_pz, np.nan), where=p_initial > 0.0)
    return pt[np.isfinite(pt)], delta_pz[np.isfinite(delta_pz)], relative_delta_pz[np.isfinite(relative_delta_pz)]


# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

def append_summary_for_sample(lines: list[str], sample: Sample, indent: str = "") -> None:
    a = sample.arrays
    esum = energy_sum_cm(sample)
    q = q_eff(sample)
    ecm = finite(a["e_cm_p11B"])
    branch = np.asarray(a["branch_id"], dtype=int)
    unique_branch, counts_branch = np.unique(branch, return_counts=True)
    pt, dpz, rdpz = momentum_diagnostics(sample)
    e12, e13, e23 = alpha_alpha_relative_energies(sample)
    emin = np.minimum(np.minimum(e12, e13), e23)

    lines.append(f"{indent}[{sample.label}]")
    lines.append(f"{indent}  ROOT file              : {sample.root_file}")
    lines.append(f"{indent}  tree                   : {sample.tree_name}")
    lines.append(f"{indent}  reaction_channel       : {sample.channel}")
    lines.append(f"{indent}  selected entries        : {len(a['reaction_channel'])}")
    lines.append(f"{indent}  branch_id counts        : " + ", ".join(f"{b}:{c}" for b, c in zip(unique_branch, counts_branch)))
    lines.append(f"{indent}  E_cm mean/std [MeV]     : {np.mean(ecm):.8f} / {np.std(ecm):.8e}")
    lines.append(f"{indent}  sum T_cm mean/std [MeV] : {np.mean(esum):.8f} / {np.std(esum):.8e}")
    lines.append(f"{indent}  Q_eff mean/std [MeV]    : {np.mean(q):.8f} / {np.std(q):.8e}")
    lines.append(f"{indent}  Q_eff median [MeV]      : {np.median(q):.8f}")
    lines.append(f"{indent}  Q_eff - Qref mean [MeV] : {np.mean(q - REFERENCE_Q_P11B_TO_3ALPHA_MEV):.8e}")
    lines.append(f"{indent}  min(Eij) median [MeV]   : {np.median(emin):.8e}")
    lines.append(f"{indent}  min(Eij) 5/95% [MeV]    : {np.percentile(emin, 5):.8e} / {np.percentile(emin, 95):.8e}")
    lines.append(f"{indent}  pT_final median [MeV/c] : {np.median(pt):.8e}")
    lines.append(f"{indent}  pT_final 95% [MeV/c]    : {np.percentile(pt, 95):.8e}")
    lines.append(f"{indent}  delta_pz median [MeV/c] : {np.median(dpz):.8e}")
    lines.append(f"{indent}  rel delta_pz median     : {np.median(rdpz):.8e}")
    lines.append("")


def write_summary(samples: list[Sample], output_dir: Path) -> None:
    lines = []
    lines.append("Three-alpha kinematics and Dalitz validation summary")
    lines.append("====================================================")
    lines.append("")
    lines.append(f"Reference Q(p+11B -> 3 alpha) used only for orientation: {REFERENCE_Q_P11B_TO_3ALPHA_MEV:.6f} MeV")
    lines.append(f"8Be(g.s.) E_alpha-alpha reference: {BE8_GS_RELATIVE_ENERGY_MEV:.5f} MeV")
    lines.append(f"8Be(2+) broad-reference centroid : {BE8_2PLUS_RELATIVE_ENERGY_MEV:.3f} MeV")
    lines.append("Energy-closure check mainly uses Q_eff = sum(T_alpha_cm) - E_cm(p11B).")
    lines.append("Symmetrized Dalitz plots fill each event six times as a visualization diagnostic; this is not a coherent amplitude-level Bose symmetrization.")
    lines.append("")

    for sample in samples:
        append_summary_for_sample(lines, sample)
        for b in branch_ids(sample):
            mask = np.asarray(sample.arrays["branch_id"], dtype=int) == b
            sub = make_filtered_sample(sample, f"{sample.label}_branch{b}", mask)
            append_summary_for_sample(lines, sub, indent="  ")

    path = output_dir / "three_alpha_kinematics_summary.txt"
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"saved: {path}")


# -----------------------------------------------------------------------------
# Plot functions
# -----------------------------------------------------------------------------

def savefig(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"saved: {path}")


def plot_energy_spectra(sample: Sample, output_dir: Path) -> None:
    e1, e2, e3 = cm_energies(sample)
    all_e = np.concatenate([e1, e2, e3])
    lo, hi = np.percentile(all_e, [0.1, 99.9])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = np.min(all_e), np.max(all_e)
    bins = np.linspace(lo, hi, 120)

    plt.figure(figsize=(8, 6))
    plt.hist(e1, bins=bins, histtype="step", density=True, linewidth=1.4, label="alpha1")
    plt.hist(e2, bins=bins, histtype="step", density=True, linewidth=1.4, label="alpha2")
    plt.hist(e3, bins=bins, histtype="step", density=True, linewidth=1.4, label="alpha3")
    plt.xlabel("Three-alpha CM kinetic energy [MeV]")
    plt.ylabel("Normalized counts")
    plt.title(f"CM alpha energy spectra: {sample.label}")
    plt.legend()
    savefig(output_dir / f"cm_energy_spectra_{sanitize_label(sample.label)}.png")


def plot_q_eff(sample: Sample, output_dir: Path) -> None:
    q = q_eff(sample)
    residual = q - np.median(q)
    width = np.percentile(np.abs(residual), 99.5)
    if not np.isfinite(width) or width <= 0.0:
        width = max(np.std(residual) * 5.0, 1e-9)
    bins = np.linspace(-width, width, 120)

    plt.figure(figsize=(8, 6))
    plt.hist(residual, bins=bins, histtype="step", linewidth=1.4)
    plt.xlabel("Q_eff - median(Q_eff) [MeV]")
    plt.ylabel("Counts")
    plt.title(f"Energy closure residual: {sample.label}")
    savefig(output_dir / f"energy_closure_residual_{sanitize_label(sample.label)}.png")


def plot_dalitz(sample: Sample, output_dir: Path, max_points: int, symmetrized: bool = False, prefix: str | None = None) -> None:
    x, y = dalitz_xy(sample, symmetrized=symmetrized)
    if len(x) > max_points:
        rng = np.random.default_rng(12345)
        idx = rng.choice(len(x), size=max_points, replace=False)
        x, y = x[idx], y[idx]

    mode = "symmetrized" if symmetrized else "labelled"
    filename_prefix = prefix if prefix is not None else ("dalitz_symmetrized" if symmetrized else "dalitz")

    plt.figure(figsize=(7, 6))
    plt.hist2d(x, y, bins=180, range=[[-1.1, 1.1], [-1.1, 1.1]])
    circle = plt.Circle((0.0, 0.0), 1.0, fill=False, linestyle="--", linewidth=0.8)
    plt.gca().add_patch(circle)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel(r"$\sqrt{3}(T_2-T_3)/(T_1+T_2+T_3)$")
    plt.ylabel(r"$(2T_1-T_2-T_3)/(T_1+T_2+T_3)$")
    plt.title(f"Dalitz plot ({mode}): {sample.label}")
    plt.colorbar(label="Counts")
    savefig(output_dir / f"{filename_prefix}_{sanitize_label(sample.label)}.png")


def plot_dalitz_projections(sample: Sample, output_dir: Path, symmetrized: bool = True) -> None:
    rho, phi = dalitz_rho_phi(sample, symmetrized=symmetrized)
    mode = "symmetrized" if symmetrized else "labelled"
    suffix = "symmetrized" if symmetrized else "labelled"

    plt.figure(figsize=(8, 6))
    plt.hist(rho, bins=np.linspace(0.0, 1.05, 106), histtype="step", density=True, linewidth=1.4)
    plt.xlabel(r"Dalitz radial coordinate $\rho$")
    plt.ylabel("Normalized counts")
    plt.title(f"Dalitz radial projection ({mode}): {sample.label}")
    savefig(output_dir / f"rho_projection_{suffix}_{sanitize_label(sample.label)}.png")

    plt.figure(figsize=(8, 6))
    plt.hist(phi, bins=np.linspace(0.0, 60.0, 61), histtype="step", density=True, linewidth=1.4)
    plt.xlabel(r"Folded Dalitz angle $\phi$ [deg]")
    plt.ylabel("Normalized counts")
    plt.title(f"Dalitz angular projection ({mode}): {sample.label}")
    savefig(output_dir / f"phi_projection_{suffix}_{sanitize_label(sample.label)}.png")


def plot_relative_energies(sample: Sample, output_dir: Path) -> None:
    e12, e13, e23 = alpha_alpha_relative_energies(sample)
    all_e = np.concatenate([e12, e13, e23])
    hi = np.percentile(all_e, 99.5)
    if not np.isfinite(hi) or hi <= 0.0:
        hi = max(np.max(all_e), 5.0)
    hi = max(hi * 1.05, 3.5)
    bins = np.linspace(0.0, hi, 160)

    plt.figure(figsize=(8, 6))
    plt.hist(e12, bins=bins, histtype="step", density=True, linewidth=1.4, label="E12")
    plt.hist(e13, bins=bins, histtype="step", density=True, linewidth=1.4, label="E13")
    plt.hist(e23, bins=bins, histtype="step", density=True, linewidth=1.4, label="E23")
    plt.axvline(BE8_GS_RELATIVE_ENERGY_MEV, linestyle="--", linewidth=1.0, label="8Be g.s. 92 keV")
    plt.axvline(BE8_2PLUS_RELATIVE_ENERGY_MEV, linestyle=":", linewidth=1.0, label="8Be 2+ ~3.03 MeV")
    plt.xlabel(r"Alpha-alpha relative energy $E_{ij}$ [MeV]")
    plt.ylabel("Normalized counts")
    plt.title(f"Alpha-alpha relative energies: {sample.label}")
    plt.legend()
    savefig(output_dir / f"relative_energy_eij_{sanitize_label(sample.label)}.png")

    emin = np.minimum(np.minimum(e12, e13), e23)
    plt.figure(figsize=(8, 6))
    plt.hist(emin, bins=bins, histtype="step", density=True, linewidth=1.4)
    plt.axvline(BE8_GS_RELATIVE_ENERGY_MEV, linestyle="--", linewidth=1.0, label="8Be g.s. 92 keV")
    plt.axvline(BE8_2PLUS_RELATIVE_ENERGY_MEV, linestyle=":", linewidth=1.0, label="8Be 2+ ~3.03 MeV")
    plt.xlabel(r"Minimum alpha-alpha relative energy $\min(E_{ij})$ [MeV]")
    plt.ylabel("Normalized counts")
    plt.title(f"Minimum alpha-alpha relative energy: {sample.label}")
    plt.legend()
    savefig(output_dir / f"relative_energy_min_eij_{sanitize_label(sample.label)}.png")


def plot_opening_angles(sample: Sample, output_dir: Path) -> None:
    a = sample.arrays
    ang12 = np.degrees(finite(a["opening_angle_alpha12_cm"]))
    ang13 = np.degrees(finite(a["opening_angle_alpha13_cm"]))
    ang23 = np.degrees(finite(a["opening_angle_alpha23_cm"]))
    bins = np.linspace(0.0, 180.0, 121)

    plt.figure(figsize=(8, 6))
    plt.hist(ang12, bins=bins, histtype="step", density=True, linewidth=1.4, label="alpha1-alpha2")
    plt.hist(ang13, bins=bins, histtype="step", density=True, linewidth=1.4, label="alpha1-alpha3")
    plt.hist(ang23, bins=bins, histtype="step", density=True, linewidth=1.4, label="alpha2-alpha3")
    plt.xlabel("Opening angle in three-alpha CM [deg]")
    plt.ylabel("Normalized counts")
    plt.title(f"Opening angles: {sample.label}")
    plt.legend()
    savefig(output_dir / f"opening_angles_{sanitize_label(sample.label)}.png")


def plot_ex8be(sample: Sample, output_dir: Path) -> None:
    ex = finite(sample.arrays["ex_8Be"])
    ex2 = finite(sample.arrays["e_8be_excitation"])
    if len(ex) == 0 and len(ex2) == 0:
        return
    if np.nanmax(np.concatenate([ex, ex2])) <= 1e-12 and sample.channel == 2:
        return

    hi = np.percentile(np.concatenate([ex, ex2]), 99.5)
    hi = max(hi, 0.2)
    bins = np.linspace(0.0, hi * 1.05, 120)
    plt.figure(figsize=(8, 6))
    plt.hist(ex, bins=bins, histtype="step", density=True, linewidth=1.4, label="ex_8Be")
    plt.hist(ex2, bins=bins, histtype="step", density=True, linewidth=1.4, label="e_8be_excitation")
    plt.xlabel("8Be excitation energy [MeV]")
    plt.ylabel("Normalized counts")
    plt.title(f"8Be excitation diagnostic: {sample.label}")
    plt.legend()
    savefig(output_dir / f"ex8be_{sanitize_label(sample.label)}.png")


def plot_momentum(sample: Sample, output_dir: Path) -> None:
    pt, dpz, rdpz = momentum_diagnostics(sample)

    hi = np.percentile(pt, 99.5)
    if not np.isfinite(hi) or hi <= 0.0:
        hi = max(np.max(pt), 1e-9)
    plt.figure(figsize=(8, 6))
    plt.hist(pt, bins=np.linspace(0.0, hi * 1.05, 120), histtype="step", linewidth=1.4)
    plt.xlabel(r"$|p_{x,final}+p_{y,final}|$ transverse residual [MeV/c]")
    plt.ylabel("Counts")
    plt.title(f"Lab transverse momentum closure: {sample.label}")
    savefig(output_dir / f"momentum_closure_transverse_{sanitize_label(sample.label)}.png")

    center = np.median(rdpz)
    residual = rdpz - center
    width = np.percentile(np.abs(residual), 99.5)
    if not np.isfinite(width) or width <= 0.0:
        width = max(np.std(residual) * 5.0, 1e-10)
    plt.figure(figsize=(8, 6))
    plt.hist(residual, bins=np.linspace(-width, width, 120), histtype="step", linewidth=1.4)
    plt.xlabel(r"$(p_{z,final}-p_{z,beam})/p_{z,beam}$ minus median")
    plt.ylabel("Counts")
    plt.title(f"Lab longitudinal momentum closure shape: {sample.label}")
    savefig(output_dir / f"momentum_closure_longitudinal_{sanitize_label(sample.label)}.png")


def plot_q_eff_all(samples: list[Sample], output_dir: Path) -> None:
    plt.figure(figsize=(8, 6))
    for sample in samples:
        q = q_eff(sample)
        plt.hist(q, bins=120, histtype="step", density=True, linewidth=1.4, label=sample.label)
    plt.axvline(REFERENCE_Q_P11B_TO_3ALPHA_MEV, linestyle="--", linewidth=1.0, label="Q reference")
    plt.xlabel("Q_eff = sum(T_alpha_cm) - E_cm(p11B) [MeV]")
    plt.ylabel("Normalized counts")
    plt.title("Effective Q-value consistency")
    plt.legend()
    savefig(output_dir / "q_eff_all_samples.png")


def plot_all_for_sample(sample: Sample, output_dir: Path, max_dalitz_points: int, include_heavy_plots: bool = True) -> None:
    plot_q_eff(sample, output_dir)
    plot_energy_spectra(sample, output_dir)
    plot_dalitz(sample, output_dir, max_dalitz_points, symmetrized=False, prefix="dalitz")
    plot_dalitz(sample, output_dir, max_dalitz_points, symmetrized=True, prefix="dalitz_symmetrized")
    plot_dalitz_projections(sample, output_dir, symmetrized=True)
    plot_relative_energies(sample, output_dir)
    plot_opening_angles(sample, output_dir)
    plot_ex8be(sample, output_dir)
    if include_heavy_plots:
        plot_momentum(sample, output_dir)


def plot_branch_breakdowns(sample: Sample, output_dir: Path, max_dalitz_points: int) -> None:
    branch = np.asarray(sample.arrays["branch_id"], dtype=int)
    for b in branch_ids(sample):
        mask = branch == b
        if np.count_nonzero(mask) == 0:
            continue
        sub = make_filtered_sample(sample, f"{sample.label}_branch{b}", mask)
        print(f"Branch breakdown: {sub.label}, entries={len(sub.arrays['reaction_channel'])}")
        plot_energy_spectra(sub, output_dir)
        plot_dalitz(sub, output_dir, max_dalitz_points, symmetrized=False, prefix="dalitz_labelled")
        plot_dalitz(sub, output_dir, max_dalitz_points, symmetrized=True, prefix="dalitz_symmetrized")
        plot_dalitz_projections(sub, output_dir, symmetrized=True)
        plot_relative_energies(sub, output_dir)
        plot_opening_angles(sub, output_dir)
        plot_ex8be(sub, output_dir)
        plot_momentum(sub, output_dir)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot and summarize generator-level three-alpha kinematics and Dalitz diagnostics.")
    parser.add_argument(
        "--sample",
        action="append",
        required=True,
        help="Sample specification: label:root_file:reaction_channel. Can be repeated.",
    )
    parser.add_argument("--output-dir", required=True, help="Directory for PNG and TXT outputs")
    parser.add_argument("--tree", default=None, help="TTree name. Default: auto-detect tr or reaction")
    parser.add_argument("--max-dalitz-points", type=int, default=120000, help="Maximum points used in each Dalitz hist2d")
    parser.add_argument("--no-branch-breakdown", action="store_true", help="Disable automatic branch_id-split diagnostic plots")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = [load_sample(spec, args.tree) for spec in args.sample]
    for sample in samples:
        print(f"Loaded {sample.label}: {sample.root_file}, tree={sample.tree_name}, channel={sample.channel}, entries={len(sample.arrays['reaction_channel'])}")

    write_summary(samples, output_dir)
    plot_q_eff_all(samples, output_dir)

    for sample in samples:
        plot_all_for_sample(sample, output_dir, args.max_dalitz_points)
        if not args.no_branch_breakdown:
            plot_branch_breakdowns(sample, output_dir, args.max_dalitz_points)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
