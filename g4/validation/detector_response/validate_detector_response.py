#!/usr/bin/env python3
"""Validate detector response ROOT outputs.

The project writes separate reaction/event/track/step ROOT files.  Each file
contains a TTree named "tr".  This validator uses PyROOT first because the ROOT
files produced by this Geant4 application have been more reliable with PyROOT
than with uproot in this environment; it falls back to uproot when PyROOT is not
available.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np

TREE_NAME = "tr"
RANK = {"PASS": 0, "WARN": 1, "FAIL": 2}
DETECTOR_NAMES = {1: "Si", 2: "LaBr3", 3: "HPGe", 4: "BGO", 5: "NaI"}
_VOL_JIT_DECLARED = False


class CheckSet:
    def __init__(self) -> None:
        self.tests: list[dict[str, object]] = []

    def add(self, status: str, name: str, detail: str = "", **extra: object) -> None:
        item: dict[str, object] = {"status": status, "name": name, "detail": detail}
        item.update(extra)
        self.tests.append(item)

    @property
    def overall(self) -> str:
        status = "PASS"
        for item in self.tests:
            current = str(item["status"])
            if RANK[current] > RANK[status]:
                status = current
        return status


class RootTree:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.backend = "pyroot"
        try:
            import ROOT  # type: ignore

            ROOT.gROOT.SetBatch(True)
            self.ROOT = ROOT
            self.file = ROOT.TFile.Open(str(self.path))
            if not self.file or self.file.IsZombie():
                raise RuntimeError(f"could not open {self.path}")
            self.tree = self.file.Get(TREE_NAME)
            if not self.tree:
                raise RuntimeError(f"{self.path}: missing tree '{TREE_NAME}'")
            self.branches = [b.GetName() for b in self.tree.GetListOfBranches()]
            self.entries = int(self.tree.GetEntries())
            return
        except ImportError:
            pass

        try:
            import uproot  # type: ignore
        except ImportError as exc:
            raise RuntimeError("neither PyROOT nor uproot is available") from exc

        self.backend = "uproot"
        self.uproot_file = uproot.open(self.path)
        if TREE_NAME not in self.uproot_file:
            raise RuntimeError(f"{self.path}: missing tree '{TREE_NAME}'")
        self.tree = self.uproot_file[TREE_NAME]
        self.branches = [str(k).split(";")[0] for k in self.tree.keys()]
        self.entries = int(self.tree.num_entries)

    def arrays(self, names: Iterable[str], limit: int | None = None) -> dict[str, np.ndarray]:
        names = [name for name in dict.fromkeys(names) if name in self.branches]
        if not names:
            return {}
        if self.backend == "pyroot":
            dataframe = self.ROOT.RDataFrame(TREE_NAME, str(self.path))
            if limit is not None:
                dataframe = dataframe.Range(int(limit))
            out = dataframe.AsNumpy(names)
            return {name: np.asarray(out[name]) for name in names}
        return {
            name: self.tree[name].array(library="np", entry_stop=limit)
            for name in names
        }

    def strings(self, name: str, limit: int | None = None) -> np.ndarray:
        if name not in self.branches:
            return np.asarray([], dtype=object)
        if self.backend == "uproot":
            try:
                return np.asarray(self.tree[name].array(library="np", entry_stop=limit), dtype=str)
            except Exception:
                return np.asarray([], dtype=object)
        n = self.entries if limit is None else min(self.entries, int(limit))
        values = []
        for i in range(n):
            self.tree.GetEntry(i)
            raw = getattr(self.tree, name)
            if isinstance(raw, bytes):
                values.append(raw.decode(errors="replace").rstrip("\x00"))
            elif hasattr(raw, "__len__") and raw.__class__.__name__ == "LowLevelView":
                values.append(bytes(raw).split(b"\x00", 1)[0].decode(errors="replace"))
            else:
                values.append(str(raw).rstrip("\x00"))
        return np.asarray(values, dtype=object)

    def sensitive_step_edep_by_event(self) -> dict[int, float] | None:
        """Sum step `de` per event over sensitive crystal volumes (Si/LaBr3/HPGe
        `*_phy`, excluding the `*_al_shell_phy` housings).  Returns None if the
        required branches are missing.  Uses a JIT-compiled volume-name filter on
        the PyROOT backend so it stays fast on multi-million-row step files (the
        per-entry Python string decode would be far too slow here); falls back to
        a vectorized bytes decode on the uproot backend.
        """
        if not {"event", "de", "volume"}.issubset(self.branches):
            return None
        if self.backend == "pyroot":
            global _VOL_JIT_DECLARED
            if not _VOL_JIT_DECLARED:
                self.ROOT.gInterpreter.Declare(
                    r"""
bool h11b_vol_is_sensitive(const ROOT::VecOps::RVec<char>& v){
  std::string s(v.begin(), v.end());
  auto z = s.find('\0'); if (z != std::string::npos) s = s.substr(0, z);
  if (s.size() < 4 || s.compare(s.size()-4, 4, "_phy") != 0) return false;
  if (s.find("al_shell") != std::string::npos) return false;
  return s.rfind("Si_",0)==0 || s.rfind("LaBr3_",0)==0 || s.rfind("HPGe_",0)==0;
}
"""
                )
                _VOL_JIT_DECLARED = True
            df = self.ROOT.RDataFrame(TREE_NAME, str(self.path)).Define("is_sens", "h11b_vol_is_sensitive(volume)")
            out = df.AsNumpy(["event", "de", "is_sens"])
            ev, de, sens = out["event"].astype(np.int64), out["de"].astype(float), out["is_sens"].astype(bool)
        else:
            ev = np.asarray(self.tree["event"].array(library="np"), dtype=np.int64)
            de = np.asarray(self.tree["de"].array(library="np"), dtype=float)
            raw = self.tree["volume"].array(library="np")
            names = np.asarray([bytes(r).split(b"\x00", 1)[0].decode(errors="replace") if not isinstance(r, str) else r for r in raw])
            sens = np.asarray([
                nm.endswith("_phy") and "al_shell" not in nm and (nm.startswith("Si_") or nm.startswith("LaBr3_") or nm.startswith("HPGe_"))
                for nm in names
            ])
        ev, de = ev[sens], de[sens]
        totals: dict[int, float] = {}
        for e, d in zip(ev.tolist(), de.tolist()):
            totals[e] = totals.get(e, 0.0) + d
        return totals


def load_tree(path: str | None, checks: CheckSet, label: str, required: bool) -> RootTree | None:
    if not path:
        if required:
            checks.add("FAIL", f"{label} file present", "required path was not provided")
        else:
            checks.add("PASS", f"{label} file absent", "not requested for this scenario")
        return None
    p = Path(path)
    if not p.exists():
        checks.add("FAIL", f"{label} file readable", f"missing file: {p}")
        return None
    try:
        tree = RootTree(p)
    except Exception as exc:
        checks.add("FAIL", f"{label} file readable", str(exc))
        return None
    checks.add("PASS", f"{label} tree readable", f"{p.name}: {tree.entries} entries, backend={tree.backend}")
    return tree


def require_branches(checks: CheckSet, label: str, tree: RootTree | None, branches: Iterable[str]) -> None:
    if tree is None:
        return
    missing = [b for b in branches if b not in tree.branches]
    if missing:
        checks.add("FAIL", f"{label} required branches", "missing " + ", ".join(missing))
    else:
        checks.add("PASS", f"{label} required branches", f"{len(list(branches))} required branches present")


def numeric_sanity(checks: CheckSet, label: str, data: dict[str, np.ndarray], branches: Iterable[str]) -> None:
    bad = []
    for branch in branches:
        values = data.get(branch)
        if values is None or values.size == 0:
            continue
        if np.issubdtype(values.dtype, np.number):
            n_bad = int((~np.isfinite(values)).sum())
            if n_bad:
                bad.append(f"{branch}:{n_bad}")
    if bad:
        checks.add("FAIL", f"{label} numeric finite", "; ".join(bad))
    else:
        checks.add("PASS", f"{label} numeric finite", "no NaN/Inf in loaded numeric branches")


def unique_pairs(a: np.ndarray, b: np.ndarray) -> set[tuple[int, int]]:
    return set(zip(a.astype(int).tolist(), b.astype(int).tolist()))


def detector_summary(event: dict[str, np.ndarray]) -> dict[str, object]:
    if not event or "detector_type" not in event:
        return {}
    rows: dict[int, dict[str, object]] = {}
    det_types = event["detector_type"].astype(int)
    energies = event.get("e", np.zeros_like(det_types, dtype=float))
    copy_no = event.get("copy_no", np.full_like(det_types, -1))
    for det_type in sorted(set(det_types.tolist())):
        mask = det_types == det_type
        copies = copy_no[mask]
        rows[int(det_type)] = {
            "name": DETECTOR_NAMES.get(int(det_type), "unknown"),
            "hits": int(mask.sum()),
            "events": int(len(set(event["event"][mask].astype(int).tolist()))) if "event" in event else 0,
            "total_edep": float(np.sum(energies[mask])),
            "copy_min": int(np.min(copies)) if copies.size else -1,
            "copy_max": int(np.max(copies)) if copies.size else -1,
        }
    return rows


def check_ids(checks: CheckSet, event: dict[str, np.ndarray]) -> None:
    if not event:
        return
    e = event
    det = e.get("detector_type")
    copy_no = e.get("copy_no")
    ring = e.get("ring_id")
    module = e.get("module_id")
    segment = e.get("segment_id")
    if det is None or copy_no is None:
        checks.add("FAIL", "event detector identifiers", "detector_type or copy_no branch is missing")
        return
    # Widened decimal copy-number layout (see include/DetectorChannel.hh):
    #   copy_no = type*1e8 + array*1e7 + ring*1e5 + module*1e3 + segment
    # segment now spans 3 decimal digits (0..999) to hold DSSD strip ids, so the
    # per-field decode below must use the wider place values.
    TYPE_UNIT = 100_000_000
    ARRAY_UNIT = 10_000_000
    RING_UNIT = 100_000
    MODULE_UNIT = 1_000
    cn = copy_no.astype(np.int64)
    invalid_types = ~np.isin(det.astype(int), [1, 2, 3])
    invalid_copy = cn < TYPE_UNIT
    decoded_type = cn // TYPE_UNIT
    mismatch_type = decoded_type != det.astype(int)
    bad_geo = np.zeros(det.shape, dtype=bool)
    if ring is not None:
        bad_geo |= ring.astype(int) < 0
        # ring branch must agree with the ring field encoded in copy_no
        bad_geo |= (cn % TYPE_UNIT % ARRAY_UNIT) // RING_UNIT != ring.astype(np.int64)
    if module is not None:
        bad_geo |= module.astype(int) < 0
        bad_geo |= (cn % RING_UNIT) // MODULE_UNIT != module.astype(np.int64)
    if segment is not None:
        bad_geo |= segment.astype(int) < 0
        # segment branch must agree with the segment field encoded in copy_no
        bad_geo |= cn % MODULE_UNIT != segment.astype(np.int64)
    n_bad = int((invalid_types | invalid_copy | mismatch_type | bad_geo).sum())
    if n_bad:
        checks.add(
            "FAIL",
            "event detector identifiers",
            f"{n_bad} entries have invalid detector_type/copy_no/ring/module/segment encoding",
            affected_entries=n_bad,
            related_branch="detector_type,copy_no,ring_id,module_id,segment_id",
        )
    else:
        checks.add("PASS", "event detector identifiers", "detector_type and encoded copy_no are consistent")


def check_event_hits(checks: CheckSet, event: dict[str, np.ndarray], scenario: str) -> None:
    if not event or "event" not in event or "e" not in event:
        checks.add("FAIL", "event hit records", "event/e branches are unavailable")
        return
    n = len(event["event"])
    if n == 0:
        checks.add("FAIL", "event hit records", "event tree has zero detector hits")
        return
    energy = event["e"]
    negative = int((energy < -1e-12).sum())
    if negative:
        checks.add("FAIL", "event hit energy nonnegative", f"{negative} negative hit energies")
    else:
        checks.add("PASS", "event hit energy nonnegative", f"{n} detector hit records")
    high = int((energy > 50.0).sum())
    if high:
        checks.add("FAIL", "event hit energy budget", f"{high} hits exceed 50 MeV")
    elif np.max(energy) > 30.0:
        checks.add("WARN", "event hit energy budget", f"max hit energy is {float(np.max(energy)):.6g} MeV")
    else:
        checks.add("PASS", "event hit energy budget", f"max hit energy is {float(np.max(energy)):.6g} MeV")

    if "detector_type" in event:
        det = event["detector_type"].astype(int)
        if scenario == "alpha_step":
            n_si = int((det == 1).sum())
            checks.add("PASS" if n_si else "FAIL", "alpha charged-particle detector hits", f"Si hits: {n_si}")
        if scenario == "gamma_step":
            n_gamma_det = int(np.isin(det, [2, 3]).sum())
            checks.add("PASS" if n_gamma_det else "FAIL", "gamma detector hits", f"LaBr3/HPGe hits: {n_gamma_det}")

    if {"event", "detector_type", "copy_no"}.issubset(event):
        keys = list(zip(event["event"].astype(int), event["detector_type"].astype(int), event["copy_no"].astype(int)))
        dup = len(keys) - len(set(keys))
        if dup:
            checks.add("FAIL", "event hit uniqueness", f"{dup} duplicate (event, detector_type, copy_no) hit rows")
        else:
            checks.add("PASS", "event hit uniqueness", "one event-level hit row per detector copy")

    # Parent/track-id sanity on the hit tree (task section 5.2): the track tree
    # omits parent_id, but each hit records the depositing track's track_id and
    # parent_id.  A primary proton hit has parent_id==0; secondaries have
    # parent_id>0 referencing their creator.  Uninitialized/negative ids or a
    # nonpositive track_id would indicate a bookkeeping bug.
    if {"track_id", "parent_id"}.issubset(event):
        tid = event["track_id"].astype(int)
        pid = event["parent_id"].astype(int)
        bad_tid = int((tid <= 0).sum())
        bad_pid = int((pid < 0).sum())
        if bad_tid or bad_pid:
            checks.add(
                "FAIL", "hit track/parent ids",
                f"track_id<=0: {bad_tid}, parent_id<0: {bad_pid}",
                affected_entries=bad_tid + bad_pid, related_branch="track_id,parent_id",
            )
        else:
            n_primary = int((pid == 0).sum())
            checks.add(
                "PASS", "hit track/parent ids",
                f"all {tid.size} hits have track_id>0 and parent_id>=0 ({n_primary} primary-parent hits)",
            )


def check_track_step(checks: CheckSet, track: dict[str, np.ndarray], step: dict[str, np.ndarray]) -> None:
    if track:
        if "event" in track:
            neg = int((track["event"] < 0).sum())
            checks.add("FAIL" if neg else "PASS", "track event ids nonnegative", f"negative entries: {neg}")
        if "track" in track:
            bad = int((track["track"] <= 0).sum())
            checks.add("FAIL" if bad else "PASS", "track ids positive", f"nonpositive entries: {bad}")
        if {"event", "track"}.issubset(track):
            pairs = list(zip(track["event"].astype(int), track["track"].astype(int)))
            dup = len(pairs) - len(set(pairs))
            checks.add("FAIL" if dup else "PASS", "track id uniqueness", f"duplicate (event, track) rows: {dup}")
        if "length" in track:
            bad_len = int((track["length"] < -1e-12).sum())
            checks.add("FAIL" if bad_len else "PASS", "track length nonnegative", f"negative entries: {bad_len}")
        if "ts" in track:
            bad_t = int((track["ts"] < -1e-12).sum())
            checks.add("FAIL" if bad_t else "PASS", "track global time nonnegative", f"negative entries: {bad_t}")
        # The track tree does not store parent_id (schema choice, not a defect);
        # the parent/track-id sanity of the reaction products is validated on the
        # hit tree instead (check_hit_parentage), which does carry both fields.
        checks.add("PASS", "track parent branch", "track tree stores no parent_id; parent-id sanity checked on the hit tree")

    if step:
        if "event" in step:
            neg = int((step["event"] < 0).sum())
            checks.add("FAIL" if neg else "PASS", "step event ids nonnegative", f"negative entries: {neg}")
        if "track" in step:
            bad = int((step["track"] <= 0).sum())
            checks.add("FAIL" if bad else "PASS", "step track ids positive", f"nonpositive entries: {bad}")
        if "de" in step:
            bad_de = int((step["de"] < -1e-12).sum())
            checks.add("FAIL" if bad_de else "PASS", "step edep nonnegative", f"negative entries: {bad_de}")
        if "length" in step:
            bad_len = int((step["length"] < -1e-12).sum())
            checks.add("FAIL" if bad_len else "PASS", "step length nonnegative", f"negative entries: {bad_len}")
        if {"post_kine_energy", "pre_kine_energy"}.issubset(step):
            increase = int((step["post_kine_energy"] > step["pre_kine_energy"] + 1e-9).sum())
            status = "WARN" if increase else "PASS"
            checks.add(status, "step kinetic energy monotonicity", f"post > pre entries: {increase}")
        if {"event", "track"}.issubset(track) and {"event", "track"}.issubset(step):
            track_pairs = unique_pairs(track["event"], track["track"])
            step_pairs = unique_pairs(step["event"], step["track"])
            orphan = step_pairs - track_pairs
            frac = len(orphan) / max(len(step_pairs), 1)
            if frac > 0.01:
                checks.add("FAIL", "step-track consistency", f"{len(orphan)} step track pairs absent from track tree")
            elif orphan:
                checks.add("WARN", "step-track consistency", f"{len(orphan)} step track pairs absent from track tree")
            else:
                checks.add("PASS", "step-track consistency", "all step (event, track) pairs are present in track tree")


def check_energy_closure(checks: CheckSet, event: dict[str, np.ndarray], step_tree: "RootTree | None") -> None:
    """Approximate event <-> step energy closure.

    The event (hit) tree stores raw summed sensitive-detector edep after a 2%
    Gaussian resolution smear and a 20 keV per-channel threshold
    (EventAction::EndOfEventAction / GausEnergy / IfThresholdTrigger).  The step
    tree stores the underlying raw per-step `de`, so summing step `de` over the
    sensitive crystal volumes per event must reproduce the recorded hit energy
    to within the smearing (sigma ~ 0.85% of E) plus the sub-threshold channel
    drop.  This exercises the whole SD -> hit accumulation path (task section
    6.1); a missing/double-counted hit, wrong unit, or channel-indexing bug
    breaks it.  A 5% relative tolerance with a 0.03 MeV absolute floor covers
    the smearing and the 20 keV threshold; the verdict is driven by the
    failed-event fraction so a handful of near-threshold events do not dominate.
    """
    if step_tree is None or not event or "event" not in event or "e" not in event:
        return
    step_totals = step_tree.sensitive_step_edep_by_event()
    if step_totals is None:
        checks.add("WARN", "event-step energy closure", "step tree lacks event/de/volume branches; closure not computed")
        return

    hit_totals: dict[int, float] = {}
    for e, energy in zip(event["event"].astype(int).tolist(), event["e"].astype(float).tolist()):
        hit_totals[e] = hit_totals.get(e, 0.0) + energy
    hit_events = {e: v for e, v in hit_totals.items() if v > 0.0}
    if not hit_events:
        checks.add("WARN", "event-step energy closure", "no positive-energy hit events to check")
        return

    rel_tol, abs_floor = 0.05, 0.03
    residuals = []
    n_fail = 0
    worst = (0.0, -1)
    for e, e_hit in hit_events.items():
        e_step = step_totals.get(e, 0.0)
        diff = e_step - e_hit
        residuals.append(diff / e_hit)
        if abs(diff) > max(rel_tol * e_hit, abs_floor):
            n_fail += 1
        if abs(diff) > abs(worst[0]):
            worst = (diff, e)
    residuals = np.asarray(residuals)
    n = len(hit_events)
    fail_frac = n_fail / n
    detail = (
        f"n={n} hit-events, mean_rel={float(np.mean(residuals)):+.4f} std_rel={float(np.std(residuals)):.4f} "
        f"fail_frac={fail_frac:.3e} (tol=max(5%,0.03MeV)); worst abs diff {worst[0]:+.4f} MeV at event {worst[1]}"
    )
    if fail_frac >= 1e-2:
        checks.add("FAIL", "event-step energy closure", detail, affected_entries=n_fail, related_branch="event.e vs step.de")
    elif fail_frac >= 1e-3:
        checks.add("WARN", "event-step energy closure", detail, affected_entries=n_fail)
    else:
        checks.add("PASS", "event-step energy closure", detail)


def check_event_id_alignment(
    checks: CheckSet,
    reaction: dict[str, np.ndarray],
    event: dict[str, np.ndarray],
    track: dict[str, np.ndarray],
    step: dict[str, np.ndarray],
) -> None:
    sets = {}
    for label, data in [("reaction", reaction), ("event", event), ("track", track), ("step", step)]:
        if data and "event" in data:
            values = data["event"].astype(int)
            neg = int((values < 0).sum())
            if neg:
                checks.add("FAIL", f"{label} event ids nonnegative", f"{neg} negative event ids")
            else:
                checks.add("PASS", f"{label} event ids nonnegative", f"{len(set(values.tolist()))} unique event ids")
            sets[label] = set(values.tolist())
    if "reaction" in sets and "event" in sets:
        event_without_reaction = sets["event"] - sets["reaction"]
        reaction_without_hit = sets["reaction"] - sets["event"]
        if event_without_reaction:
            frac = len(event_without_reaction) / max(len(sets["event"]), 1)
            status = "FAIL" if frac > 0.05 else "WARN"
            checks.add(status, "event ids mapped to reaction", f"{len(event_without_reaction)} hit events absent from reaction tree")
        else:
            checks.add("PASS", "event ids mapped to reaction", "all hit events have reaction rows")
        # Reaction events without a detector hit are expected physics: small
        # solid angle plus the 20 keV per-channel threshold means most reaction
        # events legitimately leave no above-threshold hit.  Only flag the real
        # failure mode -- a detector that produces no hits at all.
        n_reaction = max(len(sets["reaction"]), 1)
        hit_frac = 1.0 - len(reaction_without_hit) / n_reaction
        if hit_frac <= 0.0:
            checks.add("WARN", "reaction events with detector hits", "no reaction event produced an above-threshold detector hit")
        else:
            checks.add(
                "PASS",
                "reaction events with detector hits",
                f"{len(reaction_without_hit)}/{len(sets['reaction'])} reaction events have no above-threshold hit "
                f"(expected: small solid angle + 20 keV threshold); hit fraction={hit_frac:.3f}",
            )
    if "track" in sets and "step" in sets:
        missing = sets["step"] - sets["track"]
        checks.add("FAIL" if missing else "PASS", "step event ids mapped to track", f"orphan step event ids: {len(missing)}")


def check_reaction_content(checks: CheckSet, reaction: dict[str, np.ndarray], scenario: str) -> None:
    if not reaction or "reaction_channel" not in reaction:
        checks.add("FAIL", "reaction channel branch", "reaction_channel branch unavailable")
        return
    counts = Counter(reaction["reaction_channel"].astype(int).tolist())
    detail = ", ".join(f"{k}:{counts[k]}" for k in sorted(counts))
    if scenario == "alpha_step":
        status = "PASS" if counts.get(0, 0) > 0 else "FAIL"
        checks.add(status, "alpha reaction channel", detail)
    elif scenario == "gamma_step":
        status = "PASS" if counts.get(3, 0) > 0 else "FAIL"
        checks.add(status, "gamma reaction channel", detail)
    else:
        checks.add("PASS", "reaction channel inventory", detail)


def check_particles(checks: CheckSet, scenario: str, track_tree: RootTree | None, step_tree: RootTree | None) -> dict[str, object]:
    summary: dict[str, object] = {}
    for label, tree in [("track", track_tree), ("step", step_tree)]:
        if tree is None or "particle" not in tree.branches:
            continue
        particles = tree.strings("particle", limit=250000)
        counts = Counter(particles.tolist())
        summary[f"{label}_particles"] = dict(counts.most_common(20))
        if scenario == "alpha_step" and label == "track":
            n_alpha = int(sum(v for k, v in counts.items() if "alpha" in k))
            checks.add("PASS" if n_alpha else "FAIL", "alpha track presence", f"alpha track rows: {n_alpha}")
        if scenario == "gamma_step" and label == "track":
            n_gamma = int(counts.get("gamma", 0))
            checks.add("PASS" if n_gamma else "FAIL", "gamma track presence", f"gamma track rows: {n_gamma}")
    if step_tree is not None:
        for branch in ["volume", "process"]:
            if branch in step_tree.branches:
                values = step_tree.strings(branch, limit=250000)
                key = "step_processes" if branch == "process" else f"step_{branch}s"
                summary[key] = dict(Counter(values.tolist()).most_common(50))
    return summary


def write_csv(path: Path, tests: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["status", "name", "detail", "affected_entries", "related_branch"])
        writer.writeheader()
        for item in tests:
            writer.writerow({
                "status": item.get("status", ""),
                "name": item.get("name", ""),
                "detail": item.get("detail", ""),
                "affected_entries": item.get("affected_entries", ""),
                "related_branch": item.get("related_branch", ""),
            })


def write_markdown(path: Path, summary: dict[str, object]) -> None:
    lines = [
        f"# Detector response summary: {summary['scenario']}",
        "",
        f"- Overall status: `{summary['overall_status']}`",
        f"- Git commit: `{summary['git_commit']}`",
        "",
        "## Files",
    ]
    for label, info in summary["files"].items():
        lines.append(f"- {label}: `{info.get('path', '')}` entries=`{info.get('entries', 'n/a')}`")
    lines.extend(["", "## Checks", ""])
    for item in summary["tests"]:
        lines.append(f"- [{item['status']}] {item['name']}: {item.get('detail', '')}")
    lines.extend(["", "## Detector Summary", ""])
    det_summary = summary.get("detector_summary", {})
    if det_summary:
        lines.append("| detector_type | name | hits | events | total_edep | copy_min | copy_max |")
        lines.append("| --- | --- | ---: | ---: | ---: | ---: | ---: |")
        for key, row in sorted(det_summary.items(), key=lambda kv: int(kv[0])):
            lines.append(
                f"| {key} | {row['name']} | {row['hits']} | {row['events']} | "
                f"{row['total_edep']:.8g} | {row['copy_min']} | {row['copy_max']} |"
            )
    else:
        lines.append("No event detector summary available.")
    path.write_text("\n".join(lines) + "\n")


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "log", "--oneline", "-1"], text=True).strip()
    except Exception:
        return "unknown"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", required=True, choices=["alpha_step", "gamma_step", "event_stats"])
    parser.add_argument("--reaction", required=True)
    parser.add_argument("--event", required=True)
    parser.add_argument("--track", default=None)
    parser.add_argument("--step", default=None)
    parser.add_argument("--expect-track", action="store_true")
    parser.add_argument("--expect-step", action="store_true")
    parser.add_argument("--json", required=True)
    parser.add_argument("--summary-md", required=True)
    parser.add_argument("--checks-csv", required=True)
    parser.add_argument("--max-entries", type=int, default=None, help="debug limit for numeric ROOT branch reads")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    checks = CheckSet()

    reaction_tree = load_tree(args.reaction, checks, "reaction", required=True)
    event_tree = load_tree(args.event, checks, "event", required=True)
    track_tree = load_tree(args.track, checks, "track", required=args.expect_track)
    step_tree = load_tree(args.step, checks, "step", required=args.expect_step)

    require_branches(checks, "reaction", reaction_tree, ["event", "reaction_channel"])
    require_branches(checks, "event", event_tree, ["event", "detector_type", "copy_no", "ring_id", "module_id", "segment_id", "e", "time", "x", "y", "z", "pdg", "track_id", "parent_id"])
    if track_tree is not None:
        require_branches(checks, "track", track_tree, ["event", "track", "e", "x", "y", "z", "ts", "length", "volume", "particle"])
    if step_tree is not None:
        require_branches(checks, "step", step_tree, ["event", "track", "de", "pre_x", "pre_y", "pre_z", "pre_kine_energy", "post_x", "post_y", "post_z", "post_kine_energy", "length", "volume", "particle", "process"])

    reaction = reaction_tree.arrays(["event", "reaction_channel", "n_prompt_gammas", "gamma1_energy", "gamma2_energy", "e_alpha1", "e_alpha2", "e_alpha3"], args.max_entries) if reaction_tree else {}
    event = event_tree.arrays(["event", "detector_type", "array_id", "ring_id", "module_id", "segment_id", "copy_no", "ring", "sector", "e", "time", "x", "y", "z", "pdg", "track_id", "parent_id"], args.max_entries) if event_tree else {}
    track = track_tree.arrays(["event", "track", "e", "x", "y", "z", "ts", "length"], args.max_entries) if track_tree else {}
    step = step_tree.arrays(["event", "track", "de", "pre_x", "pre_y", "pre_z", "pre_total_energy", "pre_kine_energy", "post_x", "post_y", "post_z", "post_total_energy", "post_kine_energy", "length"], args.max_entries) if step_tree else {}

    numeric_sanity(checks, "reaction", reaction, reaction.keys())
    numeric_sanity(checks, "event", event, event.keys())
    numeric_sanity(checks, "track", track, track.keys())
    numeric_sanity(checks, "step", step, step.keys())

    check_reaction_content(checks, reaction, args.scenario)
    check_event_id_alignment(checks, reaction, event, track, step)
    check_ids(checks, event)
    check_event_hits(checks, event, args.scenario)
    check_track_step(checks, track, step)
    string_summary = check_particles(checks, args.scenario, track_tree, step_tree)

    check_energy_closure(checks, event, step_tree)

    files = {}
    for label, tree in [("reaction", reaction_tree), ("event", event_tree), ("track", track_tree), ("step", step_tree)]:
        files[label] = {
            "path": str(tree.path) if tree else "",
            "entries": tree.entries if tree else 0,
            "branches": tree.branches if tree else [],
        }

    summary = {
        "scenario": args.scenario,
        "overall_status": checks.overall,
        "git_commit": git_commit(),
        "files": files,
        "tests": checks.tests,
        "detector_summary": detector_summary(event),
        "string_summary": string_summary,
    }

    json_path = Path(args.json)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_markdown(Path(args.summary_md), summary)
    write_csv(Path(args.checks_csv), checks.tests)
    print(f"{args.scenario}: {checks.overall}")
    return 0 if checks.overall != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
