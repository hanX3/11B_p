# cs_channel validation

Validates the p+11B cross-section model and reaction-channel selection logic
described in `gpt/codex_cs_channel_validation.md` against the actual code in
`include/H11BCrossSection.hh`, `src/H11BCrossSection.cc`,
`include/H11BReaction.hh`, `src/H11BReaction.cc`, and the ROOT branches written
by `src/RootIO.cc`. It does not touch detector geometry, sensitive detectors,
materials, or energy-deposition logic.

## Layout

```text
validation/cs_channel/
├── README.md                       (this file)
├── _common.py                      shared ROOT-loading / PASS-WARN-FAIL helpers
├── validate_cs_components.py       section 3 / 5.1: cross-section algebra
├── validate_channel_sampling.py    section 4 / 5.2: sampled channel counts vs probability
├── validate_gamma_branching.py     section 5.3: gamma resonance/branch selection
├── validate_event_weight.py        section 5.4: event_weight vs bias factors
├── run_cs_channel_validation.sh    runs the validation macros + all scripts above
├── plots/                          (currently unused; reserved for future diagnostic plots)
└── reports/                        per-macro PASS/WARN/FAIL text output + summary.txt
```

Run the whole suite from the project root:

```bash
bash validation/cs_channel/run_cs_channel_validation.sh
```

or run one script against an existing merged ROOT file:

```bash
python3 validation/cs_channel/validate_cs_components.py data/reaction_merged_xxx.root
```

## `channel_probability_*` normalization (see task doc section 9)

`H11BReaction.cc` fills two different families of "probability" branches with
**different denominators**; they must not be summed or compared directly:

```text
channel_probability_165 / _675 / _directdecay / _directdecay_165 / _directdecay_675
    = sigma_{...} / sigma_total                     (physical 3-alpha total only)

channel_probability_gamma / probability_gamma_165_0 / _165_1 / _675_*
    = sigma_{...} / sigma_total_sampling_all         (sampling total, includes bias factors)
```

Because the two groups are normalized against different totals, none of these
diagnostic branches can be added together to reconstruct a single 0-1
probability across all four reaction channels. For that, use the derived
quantities computed in `validate_channel_sampling.py`:

```text
p_sampling_165 / p_sampling_675 / p_sampling_directdecay / p_sampling_gamma
    = sigma_{...}_sampling_b / sigma_total_sampling_all_b
```

These four are the ones actually used by `SelectReactionChannel()` and are the
ones this validation suite checks against observed event counts.

We did not rename or add ROOT branches for this — the task doc says to do so
only if validation finds the existing branches are actually wrong, and they
are not; they are just easy to misread.

## Validation methodology note: per-event ratios, not mean/mean

Beam protons lose energy traversing the target before reacting, and the
165/675 keV cross sections are narrow resonances, so `sigma_*_sampling_b`
varies substantially event-to-event even at fixed beam energy. The expected
sampling probability for a channel is the **mean of the per-event ratio**
`sigma_channel_sampling_b[i] / sigma_total_sampling_all_b[i]`, not
`mean(sigma_channel_sampling_b) / mean(sigma_total_sampling_all_b)` — the two
differ by tens of statistical sigma in `validation_165gamma_only.mac`
(the naive mean/mean estimate gives 96.5% expected gamma fraction; the
per-event estimate and the observed count both agree at ~95.2%).
`validate_channel_sampling.py` uses the per-event form throughout, including
for the direct-decay 165-vs-675 split, which additionally must be conditioned
on the direct-decay subset (`reaction_channel == 2`), not averaged over all
events.

## Findings from this validation pass

Two real issues were found in the simulation code (not in this validation
suite). Both were confirmed by the FAIL entries below, then fixed in
`src/H11BCrossSection.cc` and `src/H11BReaction.cc`; re-running the full
suite after the fix turned both `[FAIL]`s into `[PASS]` with no new
failures (236 PASS / 0 FAIL / 10 WARN, up from 225 PASS / 7 FAIL / 12 WARN).
See `reports/summary.txt` and `reports/gamma165_cs_components.txt` /
`reports/gamma675_gamma_branching.txt` for the raw evidence from the run
that caught them.

### 1. `gamma_bias_factor` ROOT branch is stale for non-gamma-channel events (fixed)

`reaction_data.gamma_bias_factor` is only assigned inside
`H11BReaction::GenerateGammaCapture()` (`src/H11BReaction.cc:1681`). For events
where `reaction_channel` is 0, 1, or 2 (not gamma capture), the branch keeps
the `H11BReactionData` constructor default of `1.0`, even when
`/h11b/gammaBiasFactor` is set to something else (e.g. 10000 in
`validation_165gamma_only.mac`). The actual `sigma_gamma_*_sampling_b`
branches are unaffected (they are computed from the real global bias factor
inside `H11BCrossSection::CalculateComponents()` for every event, regardless
of channel), and `event_weight` is likewise unaffected since it is computed
from `selected_components` directly, not from this diagnostic branch. So this
is a diagnostics-only inconsistency: reconstructing
`sigma_gamma_*_physical_b * cross_section_bias_factor * gamma_bias_factor` for
a 3-alpha-channel event gives the wrong answer because `gamma_bias_factor` for
that event's row is stale.

Confirmed with `data/reaction_merged_20260702_142129.root`
(`validation_165gamma_only.mac`): `gamma_bias_factor` reads `1.0` for all 4158
`reaction_channel == 2` rows and `10000.0` for all 83259
`reaction_channel == 3` rows, though the macro sets a single global
`/h11b/gammaBiasFactor 10000` for the whole run.

**Fix:** `FillCrossSectionDiagnostics()` in `src/H11BReaction.cc` (called for
every event, from every reaction-channel code path) now sets
`data.gamma_bias_factor = H11BConfig::GetGammaBiasFactor()` unconditionally,
so the branch reflects the run's configured bias factor for every event
regardless of channel. The redundant per-event assignment inside
`GenerateGammaCapture()` was removed since it's now set uniformly upstream.

### 2. 675 keV gamma-capture cross section is coupled to the 675 sequential/direct split (fixed)

`H11BCrossSection::CalculateComponents()` computed:

```cpp
components.sigma_gamma_675_physical = Get675GammaTotalCrossSection(components.sigma_675);
```

where `components.sigma_675 = seq675 * components.sigma_675_total` — i.e. the
*sequential-only* portion of the 675 alpha channel, not `sigma_675_total` (or
`sigma_675_model`). `README.md`'s "Gamma Capture Channels" section documents
675 gamma capture as an independent exit channel with
`sigma_gamma_675 = 1e-5 * sigma_675_model`, decoupled from the alpha
sequential/direct split — the code did not match that description.

Effect observed: `validation_675gamma_only.mac` sets
`/h11b/675SequentialDecayFraction 0.0` (to isolate the direct 3-alpha phase
space) and `/h11b/enable675GammaCapture true` with `/h11b/gammaBiasFactor
10000`. Because `seq675 = 0`, `components.sigma_675 = 0`, so
`sigma_gamma_675_physical` was exactly zero for the entire 200000-event run and
**zero gamma-capture events were produced**, contradicting the macro's stated
purpose (task doc section 7.5: "reaction_channel = 3 存在, gamma_resonance =
675"). `validate_gamma_branching.py` reported this as `[FAIL] no
gamma-capture (reaction_channel=3) events in this file (macro_hint=675gamma_only)`
rather than silently treating it as an empty sample.

**Fix:** `H11BCrossSection.cc` now calls
`Get675GammaTotalCrossSection(components.sigma_675_model)`, matching the
documented `sigma_gamma_675 = 1e-5 * sigma_675_model` design and decoupling
675 gamma capture from `675SequentialDecayFraction`.

## Current summary status

As of the last run after the fixes above (`reports/summary.txt`), all
cross-section algebra, sequential/direct split, sampling algebra,
reaction-channel-count-vs-probability, direct-decay resonance split, 675
scale factor response, gamma resonance/branch selection, and event-weight
checks **PASS** for every macro — 236 PASS / 0 FAIL / 10 WARN (up from 225
PASS / 7 FAIL / 12 WARN before the fix). `validation_675gamma_only.mac` now
produces gamma-capture events and passes gamma branching / event-weight
checks. The remaining `[WARN]` entries are all the expected "this macro
doesn't exercise that channel" case (e.g. no direct-decay events in a
gamma-only run, no gamma events in a directdecay/scale675 run).
