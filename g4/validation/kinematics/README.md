# Reaction Kinematics Validation

This directory contains reaction-generator validation tools for
`11B(p,3alpha)` and `11B(p,gamma)12C` kinematics.  The checks are intentionally
limited to the reaction event generator and ROOT reaction tree diagnostics; they
do not validate detector geometry, energy deposition, tracking cuts, or
visualization.

All kinematic calculations use MeV and MeV/c with `c = 1`.

## Channel Map

The channel IDs come from `include/H11BReaction.hh` and are recorded in
`include/DataStructure.hh`:

```text
reaction_channel = 0  -> 165-keV sequential alpha channel
reaction_channel = 1  -> 675-keV sequential alpha channel
reaction_channel = 2  -> direct 3-alpha phase-space channel
reaction_channel = 3  -> gamma-capture channel
```

Run `scripts/inspect_reaction_branches.py` on each ROOT file first when working
with new output.  It prints all available branches and the observed channel
counts.

## Mass Constants Used

The Python validators derive `m_p`, `m_11B`, and `m_alpha` from the ROOT
kinematic branches that the C++ generator itself filled.  This avoids false
residuals from using a Python mass table that differs from Geant4's ion table.

For gamma capture, the effective `12C` ground-state mass is reconstructed from
the recoil-aware primary-gamma two-body formula:

```text
E_gamma = (M_i^2 - M_f^2) / (2*M_i)
M_f = m_12C_ground + E_final_state
```

The current C++ implementation explicitly generates the recoil `12C`.  Cascade
branches generate the secondary gamma in the recoiling intermediate `12C` rest
frame.

## Manual Use

Inspect a reaction ROOT file:

```bash
python3 validation/kinematics/scripts/inspect_reaction_branches.py data/reaction_merged_*.root
```

Run the core 3-alpha checks:

```bash
python3 validation/kinematics/scripts/validate_four_momentum_closure.py \
  --root data/reaction_merged_xxx.root \
  --channel all \
  --outdir validation/kinematics/plots \
  --report validation/kinematics/reports/four_momentum.txt
```

Run channel-specific checks:

```bash
python3 validation/kinematics/scripts/validate_qvalue_energy_budget.py --root data/reaction_merged_xxx.root --channel all
python3 validation/kinematics/scripts/validate_sequential_decay_kinematics.py --root data/reaction_merged_xxx.root --channel all
python3 validation/kinematics/scripts/validate_direct_3alpha_phase_space.py --root data/reaction_merged_xxx.root --channel 2
python3 validation/kinematics/scripts/validate_gamma_capture_kinematics.py --root data/reaction_merged_xxx.root --channel 3
```

## One-Command Run

After building `build/HB`, run:

```bash
bash validation/kinematics/run_kinematics_validation.sh
```

Optional arguments:

```bash
bash validation/kinematics/run_kinematics_validation.sh --events 100000 --jobs 4
```

The script runs the validation macros, links the generated `reaction_merged_*`
ROOT files into `validation/kinematics/output/`, writes plots to
`validation/kinematics/plots/`, and writes
`validation/kinematics/reports/kinematics_validation_summary.md`.

## Pass Criteria

Default numeric tolerance is `1e-6 MeV` or `1e-6 MeV/c` for energy and momentum
closure.  The validators classify residuals as:

```text
PASS: failed-event fraction < 1e-4 and no systematic mean-level bias
WARN: failed-event fraction in [1e-4, 1e-2), or missing/non-applicable checks
FAIL: failed-event fraction >= 1e-2, or systematic residual bias
```

If ROOT output is written with lower precision in the future, use
`--tolerance 1e-4` and record that tolerance in the report.

## Output Layout

```text
validation/kinematics/
  macros/    validation-only Geant4 macros
  scripts/   Python validation scripts
  output/    symlinks to generated ROOT files
  plots/     residual and Dalitz plots
  reports/   text and Markdown validation reports
```
