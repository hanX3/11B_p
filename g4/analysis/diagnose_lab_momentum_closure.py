#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

ALPHA_MASS_MEV = 3727.3794066

REQUIRED_ALPHA_BRANCHES = [
    "reaction_channel",
    "branch_id",
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

PROJECTILE_BRANCHES = [
    "projectile_px_lab",
    "projectile_py_lab",
    "projectile_pz_lab",
    "projectile_p_lab",
    "projectile_kinetic_lab",
    "projectile_theta_lab",
    "projectile_phi_lab",
]


def p_from_t_alpha(t_mev):
    return np.sqrt(np.maximum(0.0, t_mev * t_mev + 2.0 * ALPHA_MASS_MEV * t_mev))


def components_from_t_theta_phi(t_mev, theta, phi):
    p = p_from_t_alpha(t_mev)
    st = np.sin(theta)
    px = p * st * np.cos(phi)
    py = p * st * np.sin(phi)
    pz = p * np.cos(theta)
    return p, px, py, pz


def resolve_tree(root_path, requested=None):
    with uproot.open(root_path) as f:
        names = {k.split(";")[0] for k in f.keys()}
        if requested:
            if requested not in names:
                raise RuntimeError(f"Requested tree {requested} not found. Available: {sorted(names)}")
            return requested
        if "tr" in names:
            return "tr"
        if "reaction" in names:
            return "reaction"
        raise RuntimeError(f"No tr/reaction tree found. Available: {sorted(names)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root_file")
    ap.add_argument("--tree", default=None)
    ap.add_argument("--channel", type=int, default=0)
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--outdir", default="diagnose_lab_momentum_closure")
    args = ap.parse_args()

    root_file = Path(args.root_file)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tree_name = resolve_tree(root_file, args.tree)

    with uproot.open(root_file) as f:
        tree = f[tree_name]
        available = set(tree.keys())

        missing = [b for b in REQUIRED_ALPHA_BRANCHES if b not in available]
        if missing:
            raise RuntimeError(f"Missing alpha branches: {missing}")

        has_projectile = all(b in available for b in PROJECTILE_BRANCHES)
        if not has_projectile:
            missing_projectile = [b for b in PROJECTILE_BRANCHES if b not in available]
            print(f"WARNING: missing projectile branches: {missing_projectile}")

        branches = list(REQUIRED_ALPHA_BRANCHES)
        if has_projectile:
            branches += PROJECTILE_BRANCHES

        if "event_id" in available:
            branches.append("event_id")
        if "run_id" in available:
            branches.append("run_id")
        if "event" in available:
            branches.append("event")
        if "e_cm_p11B" in available:
            branches.append("e_cm_p11B")

        arr = tree.arrays(branches, library="np")

    mask = np.asarray(arr["reaction_channel"]) == args.channel

    def get(name):
        return np.asarray(arr[name])[mask]

    entry = np.nonzero(mask)[0]
    branch_id = get("branch_id").astype(int)

    e1 = get("e_alpha1").astype(float)
    e2 = get("e_alpha2").astype(float)
    e3 = get("e_alpha3").astype(float)

    th1 = get("theta_lab_alpha1").astype(float)
    th2 = get("theta_lab_alpha2").astype(float)
    th3 = get("theta_lab_alpha3").astype(float)

    ph1 = get("phi_lab_alpha1").astype(float)
    ph2 = get("phi_lab_alpha2").astype(float)
    ph3 = get("phi_lab_alpha3").astype(float)

    finite = (
        np.isfinite(e1) & np.isfinite(e2) & np.isfinite(e3)
        & np.isfinite(th1) & np.isfinite(th2) & np.isfinite(th3)
        & np.isfinite(ph1) & np.isfinite(ph2) & np.isfinite(ph3)
    )

    entry = entry[finite]
    branch_id = branch_id[finite]
    e1, e2, e3 = e1[finite], e2[finite], e3[finite]
    th1, th2, th3 = th1[finite], th2[finite], th3[finite]
    ph1, ph2, ph3 = ph1[finite], ph2[finite], ph3[finite]

    p1, px1, py1, pz1 = components_from_t_theta_phi(e1, th1, ph1)
    p2, px2, py2, pz2 = components_from_t_theta_phi(e2, th2, ph2)
    p3, px3, py3, pz3 = components_from_t_theta_phi(e3, th3, ph3)

    px_sum = px1 + px2 + px3
    py_sum = py1 + py2 + py3
    pz_sum = pz1 + pz2 + pz3

    if has_projectile:
        proj_px = get("projectile_px_lab").astype(float)[finite]
        proj_py = get("projectile_py_lab").astype(float)[finite]
        proj_pz = get("projectile_pz_lab").astype(float)[finite]
        mode = "subtract_projectile_momentum"
    else:
        proj_px = np.zeros_like(px_sum)
        proj_py = np.zeros_like(py_sum)
        proj_pz = np.zeros_like(pz_sum)
        mode = "zero_projectile_transverse_assumption"

    rx = px_sum - proj_px
    ry = py_sum - proj_py
    rz = pz_sum - proj_pz
    rt = np.sqrt(rx * rx + ry * ry)
    r3 = np.sqrt(rx * rx + ry * ry + rz * rz)

    data = {
        "entry": entry,
        "branch_id": branch_id,
        "pT_residual": rt,
        "p3_residual": r3,
        "rx": rx,
        "ry": ry,
        "rz": rz,
        "px_sum_alpha": px_sum,
        "py_sum_alpha": py_sum,
        "pz_sum_alpha": pz_sum,
        "projectile_px_lab": proj_px,
        "projectile_py_lab": proj_py,
        "projectile_pz_lab": proj_pz,
        "e_alpha1": e1,
        "e_alpha2": e2,
        "e_alpha3": e3,
        "theta_lab_alpha1": th1,
        "theta_lab_alpha2": th2,
        "theta_lab_alpha3": th3,
        "phi_lab_alpha1": ph1,
        "phi_lab_alpha2": ph2,
        "phi_lab_alpha3": ph3,
        "p_alpha1": p1,
        "p_alpha2": p2,
        "p_alpha3": p3,
    }

    if has_projectile:
        for name in PROJECTILE_BRANCHES:
            data[name] = get(name).astype(float)[finite]
    if "event_id" in arr:
        data["event_id"] = np.asarray(arr["event_id"])[mask][finite]
    if "run_id" in arr:
        data["run_id"] = np.asarray(arr["run_id"])[mask][finite]
    if "event" in arr:
        data["event"] = np.asarray(arr["event"])[mask][finite]
    if "e_cm_p11B" in arr:
        data["e_cm_p11B"] = np.asarray(arr["e_cm_p11B"])[mask][finite]

    df = pd.DataFrame(data)
    top = df.sort_values("pT_residual", ascending=False).head(args.top)

    df.to_csv(outdir / "all_lab_momentum_closure.csv", index=False)
    top.to_csv(outdir / "top_lab_momentum_closure.csv", index=False)

    with open(outdir / "summary_lab_momentum_closure.txt", "w") as w:
        w.write("Lab momentum closure diagnostic\n")
        w.write("===============================\n\n")
        w.write(f"ROOT file: {root_file}\n")
        w.write(f"tree     : {tree_name}\n")
        w.write(f"channel  : {args.channel}\n")
        w.write(f"mode     : {mode}\n")
        w.write(f"entries  : {len(df)}\n\n")

        for bid in sorted(df["branch_id"].unique()):
            sub = df[df["branch_id"] == bid]
            w.write(f"branch_id={bid}\n")
            w.write(f"  entries   = {len(sub)}\n")
            w.write(f"  pT median = {sub['pT_residual'].median():.8e} MeV/c\n")
            w.write(f"  pT 90%    = {sub['pT_residual'].quantile(0.90):.8e} MeV/c\n")
            w.write(f"  pT 95%    = {sub['pT_residual'].quantile(0.95):.8e} MeV/c\n")
            w.write(f"  pT 99%    = {sub['pT_residual'].quantile(0.99):.8e} MeV/c\n")
            w.write(f"  pT max    = {sub['pT_residual'].max():.8e} MeV/c\n")
            w.write(f"  p3 median = {sub['p3_residual'].median():.8e} MeV/c\n")
            w.write(f"  p3 95%    = {sub['p3_residual'].quantile(0.95):.8e} MeV/c\n")
            w.write(f"  p3 max    = {sub['p3_residual'].max():.8e} MeV/c\n\n")

        w.write("Top events by pT_residual\n")
        w.write("-------------------------\n")
        w.write(top.to_string(index=False))
        w.write("\n")

    upper = np.nanpercentile(df["pT_residual"], 99.9) if len(df) else 1.0
    if not np.isfinite(upper) or upper <= 0:
        upper = 1.0

    plt.figure(figsize=(8, 5))
    bins = np.linspace(0, upper, 120)
    for bid in sorted(df["branch_id"].unique()):
        sub = df[df["branch_id"] == bid]
        plt.hist(sub["pT_residual"], bins=bins, histtype="step", density=True, label=f"branch {bid}")
    plt.xlabel("pT residual after projectile subtraction [MeV/c]")
    plt.ylabel("normalized counts")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "pt_residual_hist_by_branch.png", dpi=160)
    plt.close()

    sample = df.sample(min(len(df), 60000), random_state=1) if len(df) else df

    plt.figure(figsize=(6, 6))
    plt.scatter(sample["rx"], sample["ry"], s=1, alpha=0.25)
    plt.xlabel("rx = sum(px_alpha) - px_projectile [MeV/c]")
    plt.ylabel("ry = sum(py_alpha) - py_projectile [MeV/c]")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(outdir / "rx_ry_residual_scatter.png", dpi=160)
    plt.close()

    if has_projectile:
        for xcol in ["projectile_p_lab", "projectile_theta_lab", "projectile_phi_lab", "projectile_kinetic_lab"]:
            plt.figure(figsize=(7, 5))
            plt.scatter(sample[xcol], sample["pT_residual"], s=1, alpha=0.25)
            plt.xlabel(xcol)
            plt.ylabel("pT residual [MeV/c]")
            plt.tight_layout()
            plt.savefig(outdir / f"pt_vs_{xcol}.png", dpi=160)
            plt.close()

    print("Wrote:")
    print(f"  {outdir / 'summary_lab_momentum_closure.txt'}")
    print(f"  {outdir / 'top_lab_momentum_closure.csv'}")
    print(f"  {outdir / 'pt_residual_hist_by_branch.png'}")
    print(f"  {outdir / 'rx_ry_residual_scatter.png'}")


if __name__ == "__main__":
    main()
