# ProtoDUNE-HD Analysis Diagnostics

Modular LArSoft analyzers for inspecting ProtoDUNE-HD beam instrumentation,
reconstructed objects, simulation truth, associations, and staged selections.
The package is designed to keep raw observables, derived quantities, and physics
decisions independently auditable in data and Monte Carlo.

## Available diagnostics

- `PDHDBeamInstrumentation` records beam-event provenance, trigger information,
  beamline tracks, spectrometer momentum, TOF, Cherenkov values, and PID
  candidates.
- `PDHDBeamSelectionStages` records individual reconstructed beam-selection
  decisions without using MC truth to accept events.
- Additional inventory, association, calorimetry, blip, cosmic, and truth
  modules are present as documented templates for incremental implementation.

Every module has separate data and MC FHiCL files in `job/`. Implemented
modules use runnable `run_pdhd_*_{data,mc}.fcl` entry points. Scaffold modules
use unscheduled `run_pdhd_*_{data,mc}.fcl` templates marked `template_only`; these
document intended inputs but must not be treated as runnable analyzers.

See [DIAGNOSTICS_CATALOG.md](DIAGNOSTICS_CATALOG.md) for every module's inputs,
outputs, assumptions, limitations, and validation state. Historical scaffold
details are preserved in
[IMPLEMENTATION_STATUS_AND_DESIGN.md](IMPLEMENTATION_STATUS_AND_DESIGN.md).

## Repository layout

```text
PDHDAnalysisDiagnostics/
├── Alg/             Reusable extraction and selection algorithms
├── DataProducts/    Shared diagnostic data products
├── Modules/         art analyzer plugins
├── job/             Shared and mode-specific FHiCL configurations
├── CMakeLists.txt
└── DIAGNOSTICS_CATALOG.md
```

## Running the beam-instrumentation diagnostic

Build the containing `protoduneana` checkout in its configured LArSoft
environment before running either job. Data and MC use separate entry points
because their beam products and services differ.

Data:

```bash
lar -c run_pdhd_beam_instrumentation_data.fcl \
  -n 1 \
  -s INPUT_DATA.root \
  -T pdhd_beam_instrumentation_data.root
```

Monte Carlo:

```bash
lar -c run_pdhd_beam_instrumentation_mc.fcl \
  -n 1 \
  -s INPUT_MC.root \
  -T pdhd_beam_instrumentation_mc.root
```

Omit `-n 1` to use the FHiCL default and process all input events. The data job
can rebuild `beamevent` from retained timing information and IFBeam. Its current
`SkipLLT` path stores the exact Cherenkov values consumed by the official PID
utility while separately recording whether their detector provenance is valid.

## Output contract

Diagnostic trees retain event identifiers, selected product tags, availability,
cardinality, method names, and validity fields. Missing information is kept
distinct from physical zero. Data-like reconstruction defines nominal
selections; MC truth is reserved for labeling and performance measurements.

## Validation status

Beam-instrumentation behavior has been exercised on a limited PDHD data sample,
but the most recent source and schema changes require a new user build and
runtime validation. Template modules must not be treated as implemented or
validated. Consult the catalog before using a branch in a physics result.

## References

- [DUNE software repositories](https://github.com/orgs/DUNE/repositories)
- [LArSoft code documentation](https://code-doc.larsoft.org/docs/latest/html/)
