#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

ALPHA_MASS_MEV = 3727.3794066
PROTON_MASS_MEV = 938.2720813


REQUIRED = [
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
]


def resolve_tree(root_file, requested_tree=None):
    with uproot.open(root_file) as f:
        keys = {k.split(";")[0]: k for k in f.keys()}
        if requested_tree:
            if requested_tree not in keys:
                raise KeyError(f"Tree {requested_tree} not found. Available: {list(keys)}")
            return requested_tree
        if "tr" in keys:
            return "tr"
        if "reaction" in keys:
            return "reaction"
        raise KeyError(f"No tr/reaction tree found. Available: {list(keys)}")


def p_from_t_alpha(t):
    return np.sqrt(np.maximum(0.0, t * t + 2.0 * ALPHA_MASS_MEV * t))


def components_from_t_theta_phi(t, theta, phi):
    p = p_from_t_alpha(t)
    st = np.sin(theta)
    px = p * st * np.cos(phi)
    py = p * st * np.sin(phi)
    pz = p * np.cos(theta)
    return p, px, py, pz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root_file", help="reaction_merged_*.root for 165_model")
    ap.add_argument("--tree", default=None)
    ap.add_argument("--channel", type=int, default=0)
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--outdir", default="validation_dalitz/diagnose_165_momentum")
    args = ap.parse_args()

    root_file = Path(args.root_file).expanduser()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tree_name = resolve_tree(root_file, args.tree)

    with uproot.open(root_file) as f:
        tree = f[tree_name]
        available = set(tree.keys())

        missing = [b for b in REQUIRED if b not in available]
        if missing:
            print("Missing branches:", missing)
            print("Available branches:")
            for b in sorted(available):
                print("  ", b)
            raise SystemExit(1)

        branches = list(REQUIRED)
        if "event_id" in available:
            branches.append("event_id")
        if "run_id" in available:
            branches.append("run_id")

        arr = tree.arrays(branches, library="np", entry_start=0)

    ch = np.asarray(arr["reaction_channel"])
    mask = ch == args.channel

    # Keep original tree entry index for locating suspicious events.
    entry = np.nonzero(mask)[0]

    def get(name):
        return np.asarray(arr[name])[mask]

    branch_id = get("branch_id").astype(int)
    ecm = get("e_cm_p11B").astype(float)

    t1 = get("e_alpha1").astype(float)
    t2 = get("e_alpha2").astype(float)
    t3 = get("e_alpha3").astype(float)

    th1 = get("theta_lab_alpha1").astype(float)
    th2 = get("theta_lab_alpha2").astype(float)
    th3 = get("theta_lab_alpha3").astype(float)

    ph1 = get("phi_lab_alpha1").astype(float)
    ph2 = get("phi_lab_alpha2").astype(float)
    ph3 = get("phi_lab_alpha3").astype(float)

    finite = (
        np.isfinite(ecm)
        & np.isfinite(t1) & np.isfinite(t2) & np.isfinite(t3)
        & np.isfinite(th1) & np.isfinite(th2) & np.isfinite(th3)
        & np.isfinite(ph1) & np.isfinite(ph2) & np.isfinite(ph3)
    )

    entry = entry[finite]
    branch_id = branch_id[finite]
    ecm = ecm[finite]
    t1, t2, t3 = t1[finite], t2[finite], t3[finite]
    th1, th2, th3 = th1[finite], th2[finite], th3[finite]
    ph1, ph2, ph3 = ph1[finite], ph2[finite], ph3[finite]

    p1, px1, py1, pz1 = components_from_t_theta_phi(t1, th1, ph1)
    p2, px2, py2, pz2 = components_from_t_theta_phi(t2, th2, ph2)
    p3, px3, py3, pz3 = components_from_t_theta_phi(t3, th3, ph3)

    px = px1 + px2 + px3
    py = py1 + py2 + py3
    pz = pz1 + pz2 + pz3
    pt = np.sqrt(px * px + py * py)

    # Same rough scale used in the previous validation script.
    # This pz estimate is not the main diagnostic. pT is independent of beam pz.
    lab_energy_est = ecm * 12.0 / 11.0
    p_initial_est = np.sqrt(np.maximum(0.0, lab_energy_est * lab_energy_est + 2.0 * PROTON_MASS_MEV * lab_energy_est))
    delta_pz = pz - p_initial_est
    rel_delta_pz = delta_pz / p_initial_est

    # For diagnosing which alpha may be inconsistent in transverse momentum:
    # If alpha i were the bad one, its transverse vector should equal minus the other two.
    def transverse_angle(x, y):
        return np.arctan2(y, x)

    need1x, need1y = -(px2 + px3), -(py2 + py3)
    need2x, need2y = -(px1 + px3), -(py1 + py3)
    need3x, need3y = -(px1 + px2), -(py1 + py2)

    rec1_pt = np.sqrt(px1 * px1 + py1 * py1)
    rec2_pt = np.sqrt(px2 * px2 + py2 * py2)
    rec3_pt = np.sqrt(px3 * px3 + py3 * py3)

    need1_pt = np.sqrt(need1x * need1x + need1y * need1y)
    need2_pt = np.sqrt(need2x * need2x + need2y * need2y)
    need3_pt = np.sqrt(need3x * need3x + need3y * need3y)

    dpt1 = rec1_pt - need1_pt
    dpt2 = rec2_pt - need2_pt
    dpt3 = rec3_pt - need3_pt

    dphi1 = np.arctan2(np.sin(transverse_angle(px1, py1) - transverse_angle(need1x, need1y)),
                       np.cos(transverse_angle(px1, py1) - transverse_angle(need1x, need1y)))
    dphi2 = np.arctan2(np.sin(transverse_angle(px2, py2) - transverse_angle(need2x, need2y)),
                       np.cos(transverse_angle(px2, py2) - transverse_angle(need2x, need2y)))
    dphi3 = np.arctan2(np.sin(transverse_angle(px3, py3) - transverse_angle(need3x, need3y)),
                       np.cos(transverse_angle(px3, py3) - transverse_angle(need3x, need3y)))

    data = {
        "entry": entry,
        "branch_id": branch_id,
        "e_cm_p11B": ecm,
        "pT_residual": pt,
        "px_residual": px,
        "py_residual": py,
        "pz_sum": pz,
        "p_initial_est": p_initial_est,
        "delta_pz": delta_pz,
        "rel_delta_pz": rel_delta_pz,
        "e_alpha1": t1,
        "e_alpha2": t2,
        "e_alpha3": t3,
        "theta_lab_alpha1": th1,
        "theta_lab_alpha2": th2,
        "theta_lab_alpha3": th3,
        "phi_lab_alpha1": ph1,
        "phi_lab_alpha2": ph2,
        "phi_lab_alpha3": ph3,
        "p_alpha1": p1,
        "p_alpha2": p2,
        "p_alpha3": p3,
        "px_alpha1": px1,
        "py_alpha1": py1,
        "pz_alpha1": pz1,
        "px_alpha2": px2,
        "py_alpha2": py2,
        "pz_alpha2": pz2,
        "px_alpha3": px3,
        "py_alpha3": py3,
        "pz_alpha3": pz3,
        "dpt_if_alpha1_bad": dpt1,
        "dpt_if_alpha2_bad": dpt2,
        "dpt_if_alpha3_bad": dpt3,
        "dphi_if_alpha1_bad_rad": dphi1,
        "dphi_if_alpha2_bad_rad": dphi2,
        "dphi_if_alpha3_bad_rad": dphi3,
    }

    if "event_id" in arr:
        data["event_id"] = np.asarray(arr["event_id"])[mask][finite]
    if "run_id" in arr:
        data["run_id"] = np.asarray(arr["run_id"])[mask][finite]

    df = pd.DataFrame(data)

    df_sorted = df.sort_values("pT_residual", ascending=False)
    top = df_sorted.head(args.top)

    csv_all = outdir / "all_165_momentum_residual.csv"
    csv_top = outdir / "top_165_momentum_residual.csv"
    txt = outdir / "summary_165_momentum_residual.txt"

    df.to_csv(csv_all, index=False)
    top.to_csv(csv_top, index=False)

    with open(txt, "w") as w:
        w.write("165 lab momentum residual diagnostic\n")
        w.write("====================================\n\n")
        w.write(f"ROOT file: {root_file}\n")
        w.write(f"tree     : {tree_name}\n")
        w.write(f"channel  : {args.channel}\n")
        w.write(f"entries  : {len(df)}\n\n")
        for bid in sorted(df["branch_id"].unique()):
            sub = df[df["branch_id"] == bid]
            w.write(f"branch_id={bid}\n")
            w.write(f"  entries       = {len(sub)}\n")
            w.write(f"  pT median     = {sub['pT_residual'].median():.8e} MeV/c\n")
            w.write(f"  pT 90%        = {sub['pT_residual'].quantile(0.90):.8e} MeV/c\n")
            w.write(f"  pT 95%        = {sub['pT_residual'].quantile(0.95):.8e} MeV/c\n")
            w.write(f"  pT 99%        = {sub['pT_residual'].quantile(0.99):.8e} MeV/c\n")
            w.write(f"  pT max        = {sub['pT_residual'].max():.8e} MeV/c\n")
            w.write(f"  delta_pz med  = {sub['delta_pz'].median():.8e} MeV/c\n")
            w.write("\n")

        w.write("\nTop events by pT_residual\n")
        w.write("-------------------------\n")
        w.write(top.to_string(index=False))
        w.write("\n")

    # Plot 1: pT histogram by branch.
    plt.figure(figsize=(8, 5))
    bins = np.linspace(0, np.nanpercentile(df["pT_residual"], 99.8), 100)
    for bid in sorted(df["branch_id"].unique()):
        sub = df[df["branch_id"] == bid]
        plt.hist(sub["pT_residual"], bins=bins, histtype="step", density=True, label=f"branch {bid}")
    plt.xlabel("pT residual [MeV/c]")
    plt.ylabel("normalized counts")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "pt_residual_hist_by_branch.png", dpi=160)
    plt.close()

    # Plot 2: residual vector in px-py.
    plt.figure(figsize=(6, 6))
    take = df.sample(min(len(df), 50000), random_state=1)
    plt.scatter(take["px_residual"], take["py_residual"], s=1, alpha=0.25)
    plt.xlabel("px residual [MeV/c]")
    plt.ylabel("py residual [MeV/c]")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(outdir / "px_py_residual_scatter.png", dpi=160)
    plt.close()

    # Plot 3: pT versus branch / energies / angles.
    for name in [
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
    ]:
        plt.figure(figsize=(7, 5))
        take = df.sample(min(len(df), 60000), random_state=2)
        plt.scatter(take[name], take["pT_residual"], s=1, alpha=0.25)
        plt.xlabel(name)
        plt.ylabel("pT residual [MeV/c]")
        plt.tight_layout()
        plt.savefig(outdir / f"pt_vs_{name}.png", dpi=160)
        plt.close()

    # Plot 4: which alpha would be most inconsistent if treated as the bad vector.
    plt.figure(figsize=(8, 5))
    bins = np.linspace(-0.5, 0.5, 120)
    plt.hist(df["dpt_if_alpha1_bad"], bins=bins, histtype="step", density=True, label="alpha1")
    plt.hist(df["dpt_if_alpha2_bad"], bins=bins, histtype="step", density=True, label="alpha2")
    plt.hist(df["dpt_if_alpha3_bad"], bins=bins, histtype="step", density=True, label="alpha3")
    plt.xlabel("recorded pT_i - required pT_i [MeV/c]")
    plt.ylabel("normalized counts")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "which_alpha_transverse_magnitude_mismatch.png", dpi=160)
    plt.close()

    print(f"Wrote:")
    print(f"  {csv_all}")
    print(f"  {csv_top}")
    print(f"  {txt}")
    print(f"  plots in {outdir}")


if __name__ == "__main__":
    main()
