# H + 11B Simulation Notes

## Default Detector Array Geometry

The default detector geometry is `CurrentChamberPortAligned`.  The target is at
`(0, 0, 170 mm)`, not at the chamber center, so detector positions and rotations
are generated relative to the target position.

Current default layout:

```text
Si:
  8 rectangular barrel modules inside ChamberBody vacuum
  barrel radius = 70 mm
  barrel length = 120 mm
  thickness = 0.5 mm
  1 backward annular module at z = TargetZPos - 120 mm

LaBr3:
  4 near-target side-port modules at z = +158.5 mm
  (+x, -x, +y, -y), radius = 225.6 mm

HPGe:
  4 backward side-port modules at z = -158.5 mm
  (+x, -x, +y, -y), radius = 214.75 mm
  detector axes are aligned through the corresponding backward side-port
  aperture to avoid intersecting the chamber shell
```

Si is placed inside the chamber vacuum because alpha particles cannot pass
through the stainless-steel chamber wall.  LaBr3 and HPGe are outside the
chamber and aligned to existing side ports.  The event ROOT tree keeps the old
`ring/sector` branches and adds module-level fields such as `detector_type`,
`ring_id`, `module_id`, `segment_id`, `copy_no`, `time`, `pdg`, `track_id`, and
`parent_id`.

## 675-keV Angular Distribution

For the 675-keV resonance, the primary alpha direction is sampled isotropically
in the 12C center-of-mass frame.  The relevant angular structure is the
internal angular correlation in the 8Be*(2+) -> alpha + alpha decay.

The implemented physics default is the Stave 2011 pure L=3 internal angular
correlation,

```text
W(chi) = 1 + (2/7) P2(cos chi) - (9/7) P4(cos chi)
```

where chi is the secondary-alpha angle relative to the 8Be* recoil direction.
The runtime command

```text
/h11b/675DecayModel isotropicSequential
/h11b/675DecayModel stave2011L3
```

selects the 675-keV alpha1 internal decay model.  In `stave2011L3`, the
8Be*(2+) excitation-energy sampling uses pure L=3 alpha+8Be penetrability for
the first breakup, consistent with the pure L=3 internal angular correlation.
In `isotropicSequential`, the internal 8Be*(2+) decay is isotropic and the
first-breakup penetrability keeps the previous L=1/L=3 mixture for backward
compatibility comparisons.

`Kuhlwein2022L1L3` is reserved but not implemented.  It is not exposed as a
valid `/h11b/675DecayModel` candidate for production runs; if requested
manually, the current model is left unchanged.

## Background 3-Alpha Channel

`Background3Alpha` represents the residual 11B(p,3alpha) yield not assigned to
the explicit 165-keV or 675-keV resonance components.  It is not elastic
scattering.

When `/h11b/backgroundMode alpha1Sequential` is used, the event is generated as
an alpha + 8Be*(2+) sequential decay and inherits the current 675-keV alpha1
internal decay model.  The reaction label records the actual internal model,
for example:

```text
background3alpha_alpha1_stave2011L3
background3alpha_alpha1_isotropicSequential
```

ROOT output also records the model through `h11b675_decay_model_used` and
`background_sequential_model`.

## Gamma Capture Channels

Gamma capture is implemented as an independent exit channel competing with the
3-alpha channels:

```text
p + 11B -> 12C* -> 12C + gamma
```

It is not part of `Background3Alpha`, and it is not generated after a 3-alpha
decay.  The 3-alpha evaluated decomposition remains:

```text
sigma_3alpha_eval = sigma_165_used + sigma_675_used + sigma_background_3alpha
```

The physical total cross section and Monte-Carlo sampling total are tracked
separately:

```text
sigma_total_physical_all = sigma_3alpha_eval
                         + sigma_gamma_165_0_physical
                         + sigma_gamma_165_1_physical
                         + sigma_gamma_675_physical

sigma_total_sampling_all = sigma_3alpha_sampling_total
                         + sigma_gamma_165_0_sampling
                         + sigma_gamma_165_1_sampling
                         + sigma_gamma_675_sampling
```

These are identical unless gamma or cross-section biasing is enabled.  Event
sampling uses `sigma_total_sampling_all`; physical normalization should use
`sigma_total_physical_all` together with the per-event `event_weight`.

The 165-keV gamma0/gamma1 physical cross sections are calculated with the same
Breit-Wigner form as the alpha resonances, using gamma partial widths
`Gamma_gamma0 = 0.66 eV` and `Gamma_gamma1 = 18.0 eV`.  The 165-keV gamma0
direction follows the Craig 1956 form:

```text
W(theta) = 1 - 0.19 cos(theta) + 0.21 cos(theta)^2
```

where theta is measured relative to the incident proton beam direction.  The
165-keV gamma1 branch is generated as an isotropic two-gamma cascade through
the 4.439-MeV 12C state.

The 675-keV gamma channel is a phenomenological sensitivity model, not a
partial-width calculation:

```text
sigma_gamma_675 = 1e-5 * sigma_675_model
```

The generated 675-keV gamma line is selected with the default
`RelativeLineTable` model.  The relative primary-line intensities are taken
from Zijderhand et al. 1990 Table 3 at theta = 55 deg and used as approximate
branching fractions under an isotropic-emission assumption:

```text
16.58 -> 0.000 MeV   : 15.7
16.58 -> 4.439 MeV   : 100.0
16.58 -> 12.710 MeV  : 6.8
16.58 -> 15.110 MeV  : 0.16
```

The 7.654-MeV Hoyle-state upper-limit line is excluded by default.  Only the
4.439-MeV final-state branch generates a secondary 4.439-MeV cascade gamma;
the other 675-keV branches generate only the primary gamma.  No angular
distribution correction is applied to the 675-keV gamma model.

The reaction tree records the selected gamma line through fields such as
`gamma_final_state_energy_MeV`, `gamma_primary_energy_MeV`,
`gamma_relative_intensity_used`, `gamma_branch_fraction_used`,
`gamma_branch_from_relative_table`, and `gamma_cascade_generated`.

## Beam Energy and Proton Step Limit

The default proton beam energy is `165 keV`.  Batch macros can override it
before `/run/beamOn` with the standard particle-gun command, for example:

```text
/gun/energy 675 keV
/run/beamOn 100000
```

To increase reaction statistics for mechanism and detector-response checks,
set a p + 11B cross-section sampling bias factor before `/run/beamOn`:

```text
/h11b/crossSectionBiasFactor 1000
/gun/energy 675 keV
/run/beamOn 100000
```

This multiplies the sampling cross section used by Geant4 to trigger the
reaction process.  It does not change the physical cross sections recorded in
the reaction tree; biased events carry `event_weight = physical / sampling`.
The reaction tree also records `cross_section_bias_factor` and the biased
sampling cross-section fields.

The target logical volume applies the `StepMax4Proton` maximum proton step
through `G4UserLimits`, and `ProtonStepLimiterPhysics` registers
`G4StepLimiter` for protons.  This keeps the proton energy loss sampled finely
inside the target so narrow p + 11B resonance regions are not skipped by long
tracking steps.
