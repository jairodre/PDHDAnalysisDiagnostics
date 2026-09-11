# PDHD Beamline Implementation Notes

Reference: [ProtoDUNE Beam Lines wiki, revision 44709](https://wiki.dunescience.org/w/index.php?title=ProtoDUNE_Beam_Lines&oldid=44709).

These notes preserve detector and API facts for future diagnostic modules. They
are not a selection prescription and do not validate a particular production.

## H4-VLE and `ProtoDUNEBeamEvent`

ProtoDUNE-HD is on H4-VLE. The official `BeamEvent` producer builds a
`beam::ProtoDUNEBeamEvent` from detector timing/DAQ products and IFBeam data;
the original full-reconstructed input need not retain that output product.
Consumers should retain vectors and cardinalities rather than silently choosing
one TOF, momentum, track, or fiber combination.

| Observable | Official origin | Important interpretation |
|---|---|---|
| Beam trigger | Timing trigger and event-to-beam match | Store both stages separately. |
| TOF | Upstream/downstream XBTF combinations | Multiple channel combinations are physical ambiguities. |
| Momentum | Three XBPF spectrometer planes and spill magnet current | Multiple active-fiber triplets give multiple momenta. |
| Beamline track | Downstream perpendicular XBPF monitor pairs | Tracks are projected to the TPC face; retain every candidate. |
| CKov status, PDHD | Matching DAQ HSI/LLT word to event timestamp | Preserve HSI availability, match quality, and status provenance. |
| CKov pressure | IFBeam conditions record | Store channel, timestamp/unit, source, and validity independently. |

`CKov0` is the low-pressure counter and `CKov1` the high-pressure counter. A
numeric pressure or status value is not sufficient evidence of a valid detector
measurement: absent IFBeam data and unavailable HSI frames must remain explicit.

## Accessing and qualifying BeamEvent output

Data jobs usually place `std::vector<beam::ProtoDUNEBeamEvent>` under the
`beamevent` module label. Consume that exact product explicitly and require a
valid, non-empty handle before inspecting an entry:

```cpp
auto const beamHandle =
  event.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");
if (beamHandle->empty()) {
  throw art::Exception(art::errors::ProductNotFound)
    << "No ProtoDUNEBeamEvent produced for this event";
}
```

The conventional event-quality components are kept separate:

```cpp
bool const isBeamTrigger = (beamEvent.GetTimingTrigger() == 12);
bool const isMatched = beamEvent.CheckIsMatched();
bool const passesBeamEventQuality = isBeamTrigger && isMatched;
```

`GetTimingTrigger() == 12` identifies the PDHD beam-trigger classification;
`CheckIsMatched()` says whether the detector event was matched to the beamline
event. Store both booleans, not only their conjunction, so trigger and matching
inefficiencies can be diagnosed independently.

Many consumers use `(*beamHandle)[0]` because the producer normally creates one
matched beam-event object. Do not make that index an unrecorded assumption in a
diagnostic or selection: store `beamHandle->size()`, document the selection
policy, and either require exactly one entry or write one output record per
entry. A missing, empty, or non-unique collection is a distinct state, not a
valid beam candidate.

## Complete `ProtoDUNEBeamEvent` payload for analysis

The following is the nominal payload to retain before developing PID or
beam-to-TPC selections. It intentionally keeps all reconstructed combinations;
the scalar convenience getters are not an analysis choice and must not be used
to hide multiplicity.

| Quantity | Official API | Units / shape | Later analysis use and required companion fields |
|---|---|---|---|
| Trigger class | `GetTimingTrigger()` | integer scalar | Keep raw code, `is_beam_trigger`, `CheckIsMatched()`, and their conjunction. |
| TOF | `GetTOFs()` | ns, vector | Keep every value aligned with `GetTOFChans()` and `n_tof`; use channel code when studying PID/resolution. |
| Convenience TOF | `GetTOF()`, `GetTOFChan()` | ns, integer scalar | Optional diagnostic only: they select a producer-defined first combination. |
| Spectrometer momentum | `GetRecoBeamMomenta()` | GeV/c, vector | Keep all values and `n_momenta`; do not call an arbitrary index without recording its choice. |
| Indexed momentum | `GetNRecoBeamMomenta()`, `GetRecoBeamMomentum(i)` | count, GeV/c | Use only after `i < n_momenta`; save the chosen index and selection method if a later analysis needs one nominal value. |
| CKov firing | `GetCKov0Status()`, `GetCKov1Status()` | short scalar | Keep each status, its HSI/LLT availability and match provenance. For PDHD, status is DAQ-derived, not a pressure proxy. |
| CKov pressure | `GetCKov0Pressure()`, `GetCKov1Pressure()` | producer scalar | Keep separately from firing status, with independent IFBeam record/channel/timestamp/validity provenance. |
| Fiber occupancy | `GetActiveFibers(name)` | vector of active fiber IDs/statuses | Keep the full vector or an explicitly documented lossless representation, plus `n_active_fibers` per monitor. |
| Beamline tracks | `GetBeamTracks()` | vector of `recob::Track` | Keep `n_beam_tracks` and each trajectory's TPC-face position/direction; do not retain only track zero. |

### TOF

The upstream and downstream XBTF stations each have A/B channels. The possible
physical combinations are upstream/downstream `A/A`, `B/A`, `A/B`, and `B/B`.
`tof[i]` and `tof_channel[i]` describe the same combination and must always be
written together. A later PID study can choose a calibrated channel or a
quality rule, but the initial diagnostics must preserve all pairs.

```cpp
auto const &tofs = beamEvent.GetTOFs();
auto const &tofChannels = beamEvent.GetTOFChans();
if (tofs.size() != tofChannels.size()) {
  throw art::Exception(art::errors::DataCorruption)
    << "ProtoDUNEBeamEvent TOF values and channel codes are misaligned";
}
```

### Momentum spectrometry

Each momentum is formed from a possible active-fiber triplet in the three
spectrometer monitors and the spill magnet current. Multiple active fibers
therefore give real reconstruction ambiguity, not duplicate bookkeeping. Save
the complete `GetRecoBeamMomenta()` vector, its count, the active-fiber content
of the three spectrometer monitors, and the configured nominal beam momentum
from run conditions when available. Only a subsequently documented rule may
choose a nominal reconstruction; examples include one unambiguous triplet or a
value within a declared window around the run's nominal setting.

### Fiber monitors and projected beamline tracks

For PDHD, `GetActiveFibers(name)` accepts these monitor names:

```text
XBPF022697  XBPF022698  XBPF022701  XBPF022702
XBPF022707  XBPF022708  XBPF022716  XBPF022717
```

The last four monitors form two perpendicular XY pairs downstream of the
spectrometer. BeamEvent combines their possible hit positions and projects each
line to the TPC face. Thus `GetBeamTracks()` is a vector even for a valid beam
event. For every track retained in an output tree, write start and end
coordinates, start and end directions, and the trajectory-point count. The
official projection has its terminal trajectory position at `z = 0 cm`; it is
the beamline-to-TPC matching seed, not the reconstructed TPC track.

### CKov and beamline PID inputs

Store the two CKov firing statuses, pressures, and their validity/provenance as
independent observables. At 1--2 GeV/c, the low-pressure counter separates
electrons from non-electrons and TOF separates pion/muon from proton; at 3
GeV/c, the two counter statuses are the useful separation; at 6--7 GeV/c, the
high-pressure counter separates kaon from proton. These are H4-VLE beamline
capabilities, not a ready-made PDHD analysis selection. They require the
run-period pressure settings, validated firing status, and validated TOF/momentum
calibrations before using them as PID cuts.

### Minimum analysis-ready branch groups

For a flat analysis tree, use separate branch groups rather than one overloaded
``beam PID'' value:

```text
beam_event:     product_valid, event_count, entry_index, timing_trigger,
                is_beam_trigger, is_matched, passes_beam_event_quality
tof:            values_ns[], channel_codes[], count
momentum:       values_GeV_c[], count, selected_index, selected_method
ckov:           ckov0_status, ckov1_status, status_valid/source,
                pressure_beamevent, pressure_ifbeam, pressure_valid/source
fibers:         monitor_name[], active_fiber_ids[] or documented full encoding,
                counts[]
beam_tracks:    count, per-track position, direction, and trajectory metadata
run_context:    run/subrun/event, nominal beam momentum, conditions source
```

Use invalid booleans, empty vectors, and NaN only according to documented
branch contracts; never turn a missing product, pressure, TOF, or track into a
physical zero. This layout supports later efficiency/purity studies and the
data-like beam selection without losing the raw BeamEvent evidence.

### Scope: beamline reconstruction versus TPC reconstruction and MC truth

`ProtoDUNEBeamEvent` describes the reconstructed external beamline and its
match to the detector event. It does **not** select a Pandora beam slice, choose
a TPC track, classify a particle species, or provide MC truth. Those are later,
separate stages. The recommended order is:

```text
BeamEvent availability and quality
  -> preserve all beamline TOF/momentum/track/PID inputs
  -> choose a documented beamline candidate only if the analysis needs one
  -> match that beamline candidate to reconstructed TPC objects
  -> on MC only, associate the chosen reconstructed object to truth
  -> measure efficiency, purity, and mis-tag rates
```

The data path must stop before the truth-association stage. For MC, retain the
same reconstructed BeamEvent-style and TPC-observable branches where available,
then add truth origin, PDG, and hit-based association branches in a separate
truth group. Never substitute generated momentum or PDG for an unavailable
data-like beamline measurement.

## IFBeam reference workflow

Beamline devices are recorded at CERN and transferred to the FNAL IFBeam
database. The browser requires FNAL authentication. For an independent manual
check, use event name `z,pdune`, the exact full device-field name, and a narrow
time interval around the run. CSV output is appropriate for an auditable local
snapshot; do not infer a device field from a shortened name.

An IFBeam bundle is a selected set of device fields. A bundle named
`*_ANALYSIS` is intended for analysis use, but that label alone is not a
validation: inspect the bundle and verify the required fields exist for the
run period. For 2024 PDHD, the active `BeamEvent.fcl` configuration is the
source of truth for the bundle and device names used by the running release.

Useful field convention: `GeneralTrigger:timestampCount` is a scalar count;
`GeneralTrigger:seconds[]` is an array. Arrays may contain padded zeros or a
device-specific record layout. Never use all entries blindly; use the matching
count and the official BeamEvent decoder.

## PDHD prerequisites and processing chain

The documented BeamEvent workflow is:

```text
PDHD timing decoder (RD timestamp)
  + trigger decoder (beam-trigger classification)
  + CTB decoder (HSI/LLT CKov firing bits)
  + IFBeam service (spill instrumentation and conditions)
      -> pdhd_beamevent
      -> vector<beam::ProtoDUNEBeamEvent>
```

The 2024 keepup files may omit products needed for one part of this chain. This
does not invalidate the remaining BeamEvent fields; it means every product's
availability must be reported separately. In particular, a job using
`SkipLLT: true` can still reconstruct beam profile, momentum, TOF, and IFBeam
conditions, but cannot claim the event-level CKov firing status from HSI words.

For PDHD, BeamEvent identifies a beam trigger from trigger-decoder information,
then matches the detector event to the beamline general trigger. The documented
PDHD convention assumes a zero spill offset; PDSP has additional timing-offset
and S11 fallback machinery that must not be copied into PDHD without a
detector-specific validation.

## Instrument-specific details

### TOF and general triggers

The H4-VLE upstream and downstream XBTF counters provide TOF. Multiple valid
upstream/downstream matches can occur, so `GetTOFs()` and `GetTOFChans()` must
remain aligned vectors. The channel codes represent the A/B combination in the
upstream and downstream planes. The BeamEvent configuration uses timing windows
for downstream-to-general-trigger and upstream-to-downstream matching; treat
these as release/configuration parameters, not analysis cuts.

### XBPF profile monitors

PDHD uses these eight H4-VLE monitor names:

```text
XBPF022697  XBPF022698  XBPF022701  XBPF022702
XBPF022707  XBPF022708  XBPF022716  XBPF022717
```

Each monitor has 192 one-millimetre fibers. Three spectrometer planes measure
the bending trajectory; BeamEvent combines their lateral positions with the
spill magnet current to form every possible momentum. The two downstream
horizontal/vertical monitor pairs form possible beamline tracks projected to
the TPC face at `z = 0 cm`. Adjacent fired fibers may be merged by the official
decoder. Preserve the resulting multiplicity rather than choosing a first
momentum or track.

Raw PDHD/PDVD XBPF IFBeam records use a 10-word layout. The timestamp uses the
seconds word and an 8-ns-tick word; the remaining six words encode 192 fibers.
This is useful for debugging an official decoder, not a reason to duplicate its
bit-level parsing in an analysis module.

### Cherenkov counters

For PDHD/PDVD, event-level CKov firing status comes from an HSI frame matched
within five Detector Timing System ticks of the event timestamp, then from the
appropriate LLT trigger bit. It is not interchangeable with pressure.

Pressure and spill-level conditions are IFBeam quantities. Their validity needs
the exact configured channel, a returned conditions record, its timestamp/unit,
and successful matching to the event/spill. Keep pressure validity independent
of HSI availability and preserve both sources in output.

## Reusable implementation rules

- Use `GetTOFs()` with `GetTOFChans()`, `GetRecoBeamMomenta()`,
  `GetActiveFibers(name)`, and `GetBeamTracks()` as vector-valued observables.
- Store `GetTimingTrigger()` and `CheckIsMatched()` independently before making
  an analysis-level good-trigger decision.
- For PDHD CKov status, do not replace unavailable HSI/LLT information with an
  IFBeam pressure lookup. They answer different questions.
- Keep the direct `ProtoDUNEBeamEvent` fields and any independently validated
  IFBeam conditions result in separate branches with explicit provenance.
- Verify the run-period IFBeam bundle and exact device fields in the IFBeam Data
  Browser before enabling a pressure/status provenance flag.
- Prefer the active official `BeamEvent` producer and `ProtoDUNEBeamEvent` API
  over local raw-IFBeam decoding. The wiki describes dunesw `v10_10_04d00`;
  compare its configuration/API with the active `v10_17_01d00` release before
  applying any numeric setting or low-level format detail.

## Pending validation for the current 2024 data path

1. Confirm the configured `pdhd_beamevent` bundle and CKov device fields return
   pressure records for a representative beam run.
2. Inspect HSI-frame availability and the LLT timestamp match before accepting
   `CKov0/1` firing status.
3. Record the confirmed IFBeam channel names and timestamp unit in data FHiCL;
   only then set the IFBeam provenance-valid flag.
4. If LLT status is needed, process an input retaining `ctbrawdecoder:daqLLT`
   or a compatible raw/keepup production; do not turn off `SkipLLT` on the
   current input that lacks that product.
