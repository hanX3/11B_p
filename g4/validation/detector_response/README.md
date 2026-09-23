# Detector Response Validation

This directory validates the detector-response layer of the Geant4 application:
particle transport, energy deposition, sensitive-detector hits, and consistency
between the reaction, event, track, and step ROOT outputs.

It intentionally does not revalidate reaction kinematics, cross sections, or
angular distributions. Those checks live under the other `validation/`
subdirectories.

## Quick Run

```bash
bash validation/detector_response/run_detector_response_validation.sh --step-events 1000 --stats-events 20000 --jobs 4
```

For a faster smoke test:

```bash
bash validation/detector_response/run_detector_response_validation.sh --step-events 200 --stats-events 2000 --jobs 4 --skip-build
```

Outputs are written to:

- `validation/detector_response/results/`
- `validation/detector_response/plots/`
- `validation/detector_response/logs/`
- `validation/detector_response/detector_response_validation_report.md`

## Scenarios

- `alpha_step`: 165-keV alpha-rich sequential decay with reaction/event/track/step output enabled.
- `gamma_step`: 165-keV gamma-capture-rich run with reaction/event/track/step output enabled.
- `event_stats`: 675-keV mixed response run with reaction/event output only.

The optional no-reaction transport macro described in the task note is not
included because the current runtime command set does not expose a clean switch
to disable only the custom H11B reaction process while preserving the rest of
the same run configuration.

The current ROOT output uses one `tr` tree per output file class:

- reaction: generated reaction quantities and channel labels.
- event: detector hit records after sensitive-detector aggregation, resolution, and threshold.
- track: one post-tracking record per Geant4 track.
- step: one record per Geant4 step when `/output/saveStep true`.

`event` energy is not a raw sum over all steps; it is the sensitive-detector hit
energy after the detector response path in `EventAction`, including threshold
and resolution. The validator therefore checks step sanity and hit sanity
separately, and reports step/event closure as a documented limitation unless
raw sensitive-detector step deposits are added to the ROOT output.
