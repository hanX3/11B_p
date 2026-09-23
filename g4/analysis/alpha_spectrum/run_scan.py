#!/usr/bin/env python3
"""Run an energy scan with 1,000,000 incident protons per point.

Terminal: python3 analysis/alpha_spectrum/run_scan.py
Separate totals: python3 analysis/alpha_spectrum/run_scan.py --energies 50 200 10 --separate
The default 150:10:300 scan is merged; --separate keeps individual points only.
Existing validated points are reused. Worker ROOT files are removed by batch.sh.
Detailed execution logs and generated macros are kept in /tmp, not data/.
"""

import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
DATA = PROJECT / "data"
EVENTS = 1_000_000


def open_root(path):
    return uproot.open(path, handler=uproot.source.file.MemmapSource)


def validate_point(path, energy):
    fields = ["event", "projectile_kinetic_lab", "e_alpha1", "e_alpha2", "e_alpha3",
              "resonance_id", "h11b_reaction_channel", "scale_factor_162",
              "scale_factor_675", "enable_direct_decay", "cross_section_bias_factor",
              "configured_162_alpha0_branching_fraction", "sigma_gamma_total_barn",
              "h11b675_decay_model_used", "h11b675_strict_weight",
              "h11b675_strict_weight_max", "h11b675_strict_sampling_attempts"]
    with open_root(path) as f:
        a = f["reaction"].arrays(fields, library="np")
        seeds = np.unique(f["RunInfo"]["random_seed"].array(library="np"))
    n = len(a["event"])
    ids = a["event"]
    if not 0 < n <= EVENTS or len(np.unique(ids)) != n or ids.min() < 0 or ids.max() >= EVENTS:
        raise ValueError(f"{path.name}: invalid reaction event IDs/count: {n}")
    if len(seeds) != 1:
        raise ValueError("A single point must have one master random seed")
    for field, value in [("scale_factor_162", 1), ("scale_factor_675", 1),
                         ("enable_direct_decay", 0), ("sigma_gamma_total_barn", 0),
                         ("cross_section_bias_factor", 1e11),
                         ("configured_162_alpha0_branching_fraction", .05)]:
        if not np.allclose(a[field], value, rtol=0, atol=1e-12):
            raise ValueError(f"Unexpected configuration: {field}")
    proton = a["projectile_kinetic_lab"] * 1000
    if not np.all(np.isfinite(proton)) or proton.max() > energy + 1e-6 or proton.min() <= 0:
        raise ValueError(f"Unexpected reaction energies for {energy} keV")
    alpha = np.column_stack([a[f"e_alpha{i}"] for i in (1, 2, 3)])
    q = alpha.sum(axis=1) - proton / 1000
    if not np.all(np.isfinite(alpha)) or np.any(alpha <= 0) or np.ptp(q) > 1e-7:
        raise ValueError("Alpha energy conservation failed")
    channel = a["h11b_reaction_channel"]
    if not np.all(np.isin(channel, [1, 2, 4])):
        raise ValueError("Unexpected reaction channel")
    strict = a["resonance_id"] == 675
    if not np.all(a["h11b675_decay_model_used"][strict] == 2):
        raise ValueError("675-keV strict model was not used")
    if (np.any(a["h11b675_strict_sampling_attempts"][strict] >= 10000)
            or np.any(a["h11b675_strict_weight"][strict] > a["h11b675_strict_weight_max"][strict])):
        raise ValueError("675-keV rejection sampling exceeded its configured limits")
    print(f"Validated {energy} keV: {n:,} reactions, "
          f"gs={np.count_nonzero(channel == 1):,}, "
          f"162-exc={np.count_nonzero(channel == 2):,}, "
          f"675-exc={np.count_nonzero(channel == 4):,}; "
          f"mean proton={proton.mean():.6f} keV", flush=True)
    return n, int(seeds[0])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--energies", type=int, nargs=3, default=(150, 300, 10),
                        metavar=("START", "STOP", "STEP"), help="Inclusive LAB energy range in keV")
    parser.add_argument("--separate", action="store_true",
                        help="Save separate black total spectra and an overview; do not merge points")
    args = parser.parse_args(argv)
    first, last, step = args.energies
    if first <= 0 or last < first or step <= 0 or (last - first) % step:
        parser.error("Use positive energies/step and a stop reachable from start in whole steps")
    energies = tuple(range(first, last + 1, step))
    combined = DATA / f"spectrum_{first}_{last}_step{step}_lab_sum.root"
    logs = Path(tempfile.mkdtemp(prefix=f"11bp_scan_{first}_{last}_"))
    print(f"Logs and generated macros: {logs}", flush=True)
    template = (PROJECT / "macros/spectrum_scan_lab.mac").read_text()
    paths = []
    started = time.monotonic()
    for energy in energies:
        path = DATA / f"spectrum_{energy}keV_lab_162_675.root"
        paths.append(path)
        if path.exists():
            validate_point(path, energy)
            with open_root(path) as f:
                info = f["ScanInfo"].arrays(library="np")
                if int(info["energy_keV"][0]) != energy or int(info["incident_events"][0]) != EVENTS:
                    raise ValueError(f"Incorrect scan metadata in {path}")
            continue
        macro = logs / f"spectrum_{energy}.mac"
        macro.write_text(re.sub(r"^/gun/energy .*", f"/gun/energy {energy} keV", template, flags=re.M))
        before = set(DATA.glob("*.root"))
        log = logs / f"{energy}keV.log"
        print(f"Running {energy} keV, {EVENTS:,} incident events, 8 threads; log={log}", flush=True)
        with log.open("w") as stream:
            subprocess.run([str(PROJECT / "batch.sh"), str(macro), "8"], cwd=PROJECT,
                           env={**os.environ, "CLEAN_THREADS": "1"}, stdout=stream,
                           stderr=subprocess.STDOUT, check=True)
        text = log.read_text(errors="replace")
        if re.search(r"COMMAND NOT FOUND|parameter .*out of range|Batch is interrupted|FatalException", text, re.I):
            raise RuntimeError(f"Geant4 reported a macro/runtime error; see {log}")
        if not re.search(rf"The run was\s+{EVENTS}\s+events", text):
            raise RuntimeError(f"The requested incident event count was not completed; see {log}")
        produced = set(DATA.glob("*.root")) - before
        merged = [p for p in produced if p.name.endswith("_merged.root")]
        if len(merged) != 1 or len(produced) != 1:
            raise RuntimeError(f"Expected one merged ROOT and no worker files, found {produced}")
        result = merged[0]
        n, seed = validate_point(result, energy)
        with uproot.update(result) as f:
            info = f.mktree("ScanInfo", {"energy_keV": "int32", "incident_events": "int64",
                                          "reaction_entries": "int64", "random_seed": "uint64"})
            info.extend({"energy_keV": np.array([energy], dtype=np.int32),
                         "incident_events": np.array([EVENTS], dtype=np.int64),
                         "reaction_entries": np.array([n], dtype=np.int64),
                         "random_seed": np.array([seed], dtype=np.uint64)})
        result.rename(path)
        with (PROJECT / "data.log").open("a") as stream:
            stream.write(f"Renamed scan output: {result.name} -> {path.name}\n")
        print(f"Completed {energy} keV; elapsed {(time.monotonic()-started)/60:.1f} min", flush=True)
    if not args.separate and not combined.exists():
        # Sequential hadd keeps the reaction entry ranges aligned with ScanInfo.
        partial = combined.with_suffix(".partial.root")
        with (logs / "hadd_sum.log").open("w") as stream:
            subprocess.run(["hadd", "-f", str(partial), *map(str, paths)],
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
        with open_root(partial) as f:
            info = f["ScanInfo"].arrays(library="np")
            if (f["reaction"].num_entries != int(info["reaction_entries"].sum())
                    or not np.array_equal(info["energy_keV"], energies)
                    or not np.all(info["incident_events"] == EVENTS)):
                raise ValueError("Combined ROOT failed the scan completeness check")
        partial.rename(combined)
    if not args.separate:
        print(f"Combined ROOT: {combined}", flush=True)
    print(f"Total elapsed: {(time.monotonic()-started)/60:.1f} min", flush=True)
    if __package__:
        from .plot_lab_total import main as plot_main
    else:
        from plot_lab_total import main as plot_main
    if args.separate:
        plot_main([*map(str, paths), "--separate", "--overview"])
    else:
        plot_main([str(combined), "--out", str(HERE / f"alpha_lab_{first}_{last}_sum")])


if __name__ == "__main__":
    main()
