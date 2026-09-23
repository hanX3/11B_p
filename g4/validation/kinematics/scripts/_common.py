"""Shared helpers for validation/kinematics/scripts/*.py.

Loads branches from the "tr" TTree in a reaction_merged_*.root file with
uproot, reports MISSING BRANCH instead of crashing, derives the particle
masses actually used by the C++ generator from the stored kinematic
branches themselves (so Python can never silently disagree with
src/H11BReaction.cc about which mass convention -- atomic vs. nuclear -- is
in use, per gpt/codex_reaction_kinematics_validation.md section 13), and
prints PASS/WARN/FAIL blocks in the format that document requires.

All energies/momenta are MeV and MeV/c with c = 1, matching Geant4's
internal unit system (the ROOT branches are filled directly from
G4LorentzVector values without conversion, except the few branches whose
names already end in _MeV or _keV).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

TREE_NAME = "tr"


class PyRootBranch:
    def __init__(self, tree: "PyRootTree", name: str):
        self.tree = tree
        self.name = name

    def array(self, library: str = "np", entry_stop: Optional[int] = None) -> np.ndarray:
        if library != "np":
            raise ValueError("PyROOT backend only supports library='np'")
        return self.tree.arrays([self.name], entry_stop).get(self.name, np.asarray([]))


class PyRootTree:
    """Minimal uproot-like wrapper backed by ROOT.RDataFrame.

    Some ROOT 6.34 files produced by this Geant4 job hang when branch baskets
    are read with the installed uproot/fsspec stack.  PyROOT reads the same
    files promptly, so the validators use this wrapper first and keep uproot as
    a fallback for environments without PyROOT.
    """

    def __init__(self, path: Path, root_module):
        self.path = Path(path)
        self.ROOT = root_module
        self._file = root_module.TFile.Open(str(self.path))
        if not self._file or self._file.IsZombie():
            raise SystemExit(f"{path}: could not open ROOT file")
        self._tree = self._file.Get(TREE_NAME)
        if not self._tree:
            raise SystemExit(f"{path}: missing tree '{TREE_NAME}'")
        self.num_entries = int(self._tree.GetEntries())
        self._keys = [branch.GetName() for branch in self._tree.GetListOfBranches()]

    def keys(self) -> list[str]:
        return list(self._keys)

    def __getitem__(self, name: str) -> PyRootBranch:
        if name not in self._keys:
            raise KeyError(name)
        return PyRootBranch(self, name)

    def arrays(self, names: Iterable[str], entry_stop: Optional[int] = None) -> dict[str, np.ndarray]:
        names = [name for name in names if name in self._keys]
        if not names:
            return {}
        dataframe = self.ROOT.RDataFrame(TREE_NAME, str(self.path))
        if entry_stop is not None:
            dataframe = dataframe.Range(int(entry_stop))
        arrays = dataframe.AsNumpy(names)
        return {name: np.asarray(arrays[name]) for name in names}


def load_uproot():
    try:
        import uproot  # type: ignore
    except ImportError as exc:
        raise SystemExit("uproot is not installed; cannot validate ROOT files.") from exc
    return uproot


def open_tree(path: Path):
    try:
        import ROOT  # type: ignore
        ROOT.gROOT.SetBatch(True)
        return PyRootTree(path, ROOT)
    except ImportError:
        pass

    uproot = load_uproot()
    root_file = uproot.open(path)
    if TREE_NAME not in root_file:
        raise SystemExit(f"{path}: missing tree '{TREE_NAME}'")
    return root_file[TREE_NAME]


def load_branches(tree, names: Iterable[str], max_events: Optional[int] = None) -> tuple[dict[str, np.ndarray], list[str]]:
    """Return (data, missing) where data holds only the branches that exist."""
    available = set(tree.keys())
    names = list(dict.fromkeys(names))
    missing = [name for name in names if name not in available]
    for name in missing:
        print(f"MISSING BRANCH: {name}")
    existing = [name for name in names if name in available]
    if hasattr(tree, "arrays"):
        data = tree.arrays(existing, max_events)
    else:
        entry_stop = max_events
        data = {
            name: tree[name].array(library="np", entry_stop=entry_stop)
            for name in existing
        }
    return data, missing


def select_channel(data: dict[str, np.ndarray], channel: Optional[int]) -> dict[str, np.ndarray]:
    """Return a copy of data filtered to reaction_channel == channel (or unfiltered if None)."""
    if channel is None or "reaction_channel" not in data:
        return data
    mask = data["reaction_channel"] == channel
    return {name: values[mask] if values.shape[:1] == mask.shape else values for name, values in data.items()}


def common_argparser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--root", required=True, help="Input reaction_merged_*.root file")
    parser.add_argument("--channel", default="all", help="reaction_channel to select: 0/1/2/3/all")
    parser.add_argument("--max-events", type=int, default=None, help="Maximum number of events to read")
    parser.add_argument("--outdir", default=None, help="Directory for output plots")
    parser.add_argument("--report", default=None, help="Path to write the text report")
    parser.add_argument("--tolerance", type=float, default=None, help="Override the default numeric tolerance (MeV)")
    return parser


def parse_channel_arg(value: str) -> Optional[int]:
    if value is None or value.lower() == "all":
        return None
    return int(value)


# ---------------------------------------------------------------------------
# Lorentz-vector helpers.  All functions are vectorized over numpy arrays and
# follow CLHEP::HepLorentzVector::boost() sign conventions: boost(v, beta)
# moves v from the frame in which `beta` is the velocity of the origin's rest
# frame into the frame where that origin is at rest (i.e. the same operation
# as G4LorentzVector::boost(G4ThreeVector) in the C++ generator).
# ---------------------------------------------------------------------------


def spherical_lab_vector(kinetic_energy: np.ndarray, theta: np.ndarray, phi: np.ndarray, mass: float):
    """Rebuild (E, px, py, pz) from stored kinetic energy + direction + a mass.

    Matches src/H11BReaction.cc: e_alphaN branches are lv.e() - m_4He, and
    theta_lab_alphaN/phi_lab_alphaN are G4LorentzVector::theta()/phi() of the
    same lv, so E_total = T + m and p = sqrt(T^2 + 2*T*m) recovers it exactly.
    """
    e_total = kinetic_energy + mass
    p = np.sqrt(np.clip(kinetic_energy * kinetic_energy + 2.0 * kinetic_energy * mass, 0.0, None))
    px = p * np.sin(theta) * np.cos(phi)
    py = p * np.sin(theta) * np.sin(phi)
    pz = p * np.cos(theta)
    return e_total, px, py, pz


def invariant_mass(e: np.ndarray, px: np.ndarray, py: np.ndarray, pz: np.ndarray) -> np.ndarray:
    m2 = e * e - (px * px + py * py + pz * pz)
    return np.sqrt(np.clip(m2, 0.0, None))


def boost(e, px, py, pz, beta_x, beta_y, beta_z):
    """Boost (E, p) by velocity vector (beta_x, beta_y, beta_z), CLHEP convention."""
    beta2 = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z
    gamma = 1.0 / np.sqrt(np.clip(1.0 - beta2, 1e-300, None))
    bp = beta_x * px + beta_y * py + beta_z * pz
    # gamma2 = (gamma - 1) / beta2, guarding beta2 -> 0 (no motion).
    gamma2 = np.where(beta2 > 0.0, (gamma - 1.0) / np.where(beta2 > 0.0, beta2, 1.0), 0.0)
    px2 = px + gamma2 * bp * beta_x + gamma * beta_x * e
    py2 = py + gamma2 * bp * beta_y + gamma * beta_y * e
    pz2 = pz + gamma2 * bp * beta_z + gamma * beta_z * e
    e2 = gamma * (e + bp)
    return e2, px2, py2, pz2


def boost_vector(e, px, py, pz):
    """G4LorentzVector::boostVector(): velocity of this 4-vector's rest frame."""
    return px / e, py / e, pz / e


# ---------------------------------------------------------------------------
# Mass derivation.  Every formula below is the algebraic inverse of a
# formula actually used in src/H11BReaction.cc (see README.md "Mass
# constants used" for the derivation), so a mismatch here means the Python
# and C++ side disagree about physical masses, not a genuine kinematics bug.
# ---------------------------------------------------------------------------


def _robust_center_spread(values: np.ndarray) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan"), float("nan")
    median = float(np.median(values))
    spread = float(np.median(np.abs(values - median)))
    return median, spread


def derive_projectile_masses(data: dict[str, np.ndarray]) -> dict[str, tuple[float, float]]:
    """Derive m_p and m_11B from projectile_kinetic_lab / projectile_p_lab / e_cm_p11B.

    From src/H11BReaction.cc: p_lab = sqrt(T^2 + 2*T*m_p), so
        m_p = (p_lab^2 - T^2) / (2*T).
    And e_cm_p11B = sqrt_s - m_p - m_11B with
        s = m_p^2 + m_11B^2 + 2*m_11B*(m_p + T),
    which solves to
        m_11B = (C^2 + 2*C*m_p) / (2*(T - C))   where C = e_cm_p11B.
    """
    result: dict[str, tuple[float, float]] = {}
    T = data.get("projectile_kinetic_lab")
    p_lab = data.get("projectile_p_lab")
    C = data.get("e_cm_p11B")
    if T is None or p_lab is None:
        return result

    mask = T > 0
    m_p_samples = (p_lab[mask] ** 2 - T[mask] ** 2) / (2.0 * T[mask])
    result["m_p"] = _robust_center_spread(m_p_samples)

    if C is not None:
        m_p_center = result["m_p"][0]
        denom = T[mask] - C[mask]
        valid = np.abs(denom) > 1e-12
        m_11b_samples = (C[mask][valid] ** 2 + 2.0 * C[mask][valid] * m_p_center) / (2.0 * denom[valid])
        result["m_11B"] = _robust_center_spread(m_11b_samples)

    return result


def derive_sqrt_s(data: dict[str, np.ndarray], m_p: float, m_11b: float) -> np.ndarray:
    """sqrt_s = e_cm_p11B + m_p + m_11B, matching reaction_data.e_cm_p11B fills."""
    return data["e_cm_p11B"] + m_p + m_11b


def derive_alpha_mass(sqrt_s: np.ndarray, e_cm_alpha1: np.ndarray, e_cm_alpha2: np.ndarray, e_cm_alpha3: np.ndarray) -> tuple[float, float]:
    """m_alpha = (sqrt_s - sum(T_cm_i)) / 3 for any 3-alpha-channel event.

    Energy conservation in ReactionKinematic()/GenerateThreeBodyPhaseSpace()/
    Generate675SymmetrizedCoherent3Alpha() guarantees
        sum(E_cm_i) = sqrt_s  =>  sum(T_cm_i) + 3*m_alpha = sqrt_s
    exactly for the true alpha mass, so this recovers it from data alone.
    """
    samples = (sqrt_s - (e_cm_alpha1 + e_cm_alpha2 + e_cm_alpha3)) / 3.0
    return _robust_center_spread(samples)


def print_mass_table(masses: dict[str, tuple[float, float]]) -> None:
    print("Mass constants used (derived from stored kinematic branches, MeV, c=1):")
    for name, (value, spread) in masses.items():
        print(f"  {name:10s} = {value:.9f} MeV  (median-abs-deviation across events = {spread:.3e} MeV)")


# ---------------------------------------------------------------------------
# PASS/WARN/FAIL residual reporting, per section 14 of the task doc:
#   PASS: |residual| < tolerance for essentially all events (fail fraction < 1e-4)
#   WARN: fail fraction in [1e-4, 1e-2), or a variable is missing but derivable
#   FAIL: fail fraction >= 1e-2, or a systematic (mean-level) bias is evident
# ---------------------------------------------------------------------------


class Report:
    def __init__(self, path: Path):
        self.path = path
        self.status = "PASS"
        self.lines: list[str] = []

    def _bump(self, level: str) -> None:
        order = {"PASS": 0, "WARN": 1, "FAIL": 2}
        if order[level] > order[self.status]:
            self.status = level

    def emit(self, level: str, label: str) -> None:
        self._bump(level)
        line = f"[{level}] {label}"
        print(line)
        self.lines.append(line)

    def note(self, text: str) -> None:
        print(text)
        self.lines.append(text)

    def check_residual(
        self,
        label: str,
        residual: np.ndarray,
        tolerance: float,
        branch: str,
        related_setting: Optional[str] = None,
        systematic_tolerance: Optional[float] = None,
    ) -> bool:
        residual = np.asarray(residual, dtype=float)
        finite = residual[np.isfinite(residual)]
        if finite.size == 0:
            self.emit("WARN", f"{label}: no finite residuals to check")
            return False

        abs_res = np.abs(finite)
        n_total = finite.size
        n_fail = int(np.sum(abs_res >= tolerance))
        fail_frac = n_fail / n_total
        mean = float(np.mean(finite))
        std = float(np.std(finite))
        rms = float(np.sqrt(np.mean(finite ** 2)))
        max_abs = float(np.max(abs_res))
        p95 = float(np.percentile(abs_res, 95))
        p99 = float(np.percentile(abs_res, 99))

        sys_tol = systematic_tolerance if systematic_tolerance is not None else tolerance
        systematic_bias = abs(mean) >= sys_tol

        self.note(
            f"  {label}: mean={mean:.3e} std={std:.3e} rms={rms:.3e} max_abs={max_abs:.3e} "
            f"p95={p95:.3e} p99={p99:.3e} n_fail={n_fail}/{n_total} ({fail_frac:.3e}) tol={tolerance:.3e}"
        )

        if fail_frac >= 1e-2 or systematic_bias:
            level = "FAIL"
        elif fail_frac >= 1e-4:
            level = "WARN"
        else:
            level = "PASS"

        self.emit(level, f"{label} (fail_frac={fail_frac:.3e}, mean={mean:.3e}, tol={tolerance:.3e})")
        if level == "FAIL":
            bad_idx = int(np.argmax(abs_res))
            print(f"    branch name        : {branch}")
            print(f"    expected value     : 0")
            print(f"    observed value     : {finite[bad_idx]}")
            print(f"    absolute difference: {abs_res[bad_idx]}")
            print(f"    entry index        : {bad_idx}")
            if related_setting is not None:
                print(f"    related macro setting: {related_setting}")
        return level == "PASS"

    def summary_line(self) -> str:
        return f"[{self.status}] {self.path}"

    def write(self) -> None:
        if self.path is None:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "w") as f:
            f.write("\n".join(self.lines) + "\n")
