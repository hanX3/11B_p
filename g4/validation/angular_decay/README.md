# Angular Distribution and Decay Model Validation

This directory validates reaction-level angular distributions and the 675-keV
strict three-alpha decay model.  It does not validate detector acceptance,
energy deposition, geometry, or cross-section normalization.

The validation uses reaction ROOT branches such as:

```text
cos_theta_primary_cm
cos_chi_secondary_8be
cos_theta_gamma_cm
primary_angular_mode
secondary_angular_model
h11b675_strict_weight
h11b675_strict_sampling_attempts
```

## Run

Use the one-command driver after `build/HB` exists:

```bash
bash validation/angular_decay/run_angular_decay_validation.sh
```

For a quick smoke test:

```bash
bash validation/angular_decay/run_angular_decay_validation.sh --events 1000 --strict-events 1000 --jobs 4
```

Outputs:

```text
validation/angular_decay/angular_decay_validation_report.md
validation/angular_decay/angular_decay_validation_summary.json
validation/angular_decay/plots/
validation/angular_decay/output/
```

## Important Macro Routing Note

The existing `macros/validation_angle_common.mac` enables direct 3-alpha and
sets both sequential fractions to zero.  Direct 3-alpha events do not exercise
the sequential primary/secondary angular samplers.  The run script therefore
uses the existing on/off macro names but generates temporary macros that force
sequential alpha routing for alpha Legendre validation:

```text
/h11b/enableDirectDecay false
/h11b/165SequentialDecayFraction 1.0
/h11b/675SequentialDecayFraction 1.0
```

For 675 secondary Legendre validation it also forces:

```text
/h11b/675AlphaDecayModel legacyLegendreA2A4
```

because the strict `symmetrizedCoherentL1L3` generator bypasses the legacy
secondary A2/A4 sampler by design.

## Checks

Legendre on/off tests compare:

```text
W(cos theta) = 1 + a1 P1(cos theta) + a2 P2(cos theta)
W(cos chi)   = 1 + A2 P2(cos chi) + A4 P4(cos chi)
```

against the ROOT reaction-level observables using histogram chi-square,
one-sample KS statistic, moments, ROOT coefficient branches, and on/off
two-sample KS effect size.

675 strict tests check:

```text
L1-only vs L3-only
coherent vs incoherent
phase response
permutation symmetrization
rejection-sampling diagnostics
```

The strict model stores kinetic energies in `e_3alpha_cm_alpha*`; Dalitz plots
use normalized kinetic-energy coordinates:

```text
x = sqrt(3) * (E2 - E3) / (E1 + E2 + E3)
y = (2E1 - E2 - E3) / (E1 + E2 + E3)
```

