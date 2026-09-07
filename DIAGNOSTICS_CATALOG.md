# PDHD Analysis Diagnostics Catalog

Each module has separate data and MC FHiCL under `job/`. All mode-specific files
begin with `run_pdhd_`; scaffold configurations are unscheduled and marked
`ImplementationStatus: "template_only"` until their analyzers are implemented.

The catalog is both an implementation queue and an analysis contract. “Template”
means the source contains design requirements but no registered art plugin.

## Quick catalog

| Module | Main responsibility | Mode | Status |
|---|---|---|---|
| `PDHDEventProductInventory` | Event identity, products, provenance, quality | Both | Template |
| `PDHDMCTruthInventory` | Generator records and primaries | MC | Template |
| `PDHDGeantParticleInventory` | Complete MCParticle facts | MC | Template |
| `PDHDGeantProcessInventory` | Ancestry and process edges | MC | Template |
| `PDHDEnergyDepositInventory` | IDE/energy-deposit truth | MC | Template |
| `PDHDBeamInstrumentation` | Trigger, beam tracks, momentum, TOF, PID | Both | Implemented; rebuild pending |
| `PDHDTimingInventory` | Trigger timing and reconstructed T0 | Both | Template |
| `PDHDHitInventory` | Compact hit and CNN facts | Both | Template |
| `PDHDPFParticleInventory` | Pandora hierarchy and metadata | Both | Template |
| `PDHDTrackInventory` | Track geometry and containment | Both | Template |
| `PDHDShowerInventory` | Shower geometry and reconstruction | Both | Template |
| `PDHDTrackFeatureInventory` | Track PID, momentum, CNN, summaries | Both | Template |
| `PDHDShowerFeatureInventory` | Shower energy and classifiers | Both | Template |
| `PDHDCalorimetryPointInventory` | Per-point track/shower calorimetry | Both | Template |
| `PDHDRecoAssociation` | Normalized reconstructed-object edges | Both | Template |
| `PDHDRecoTruthMatch` | Bidirectional hit-based matching | MC | Template |
| `PDHDBlipInventory` | Blip geometry, charge, and energy | Both | Template |
| `PDHDBlipAssociation` | Blip-to-object and blip-to-truth edges | Both/MC | Template |
| `PDHDBeamSelectionStages` | Data-like beam decisions and MC metrics | Both | Implemented, unbuilt |
| `PDHDCosmicSelectionStages` | Explicit cosmic and containment stages | Both | Template |

## Module contracts

### PDHDEventProductInventory

Records event identity, data/MC mode, timestamps, configured-product presence,
sizes, resolved provenance, trigger/timing/beam availability, and event-quality
status. It assumes expected InputTags are explicit and does not silently find a
replacement producer. Output is an event table plus file-level configuration
metadata. Its first implementation should follow art provenance APIs and the
ProtoDUNE empty-event utility. Limitation: availability proves neither semantic
correctness nor association compatibility. Validation: scaffold only.

### PDHDMCTruthInventory

Records every configured MCTruth object and generated particle, preserving the
generator label, origin, status, four-vectors, and record/particle indices. It
assumes generator products are explicitly configured and may be plural. Output
uses stable event and generator keys. ParticleInventoryService supplies later
MCTruth-to-Geant linkage. Limitation: generator primaries are not automatically
Geant particles or reconstructed beam objects. Validation: scaffold only.

### PDHDGeantParticleInventory

Records every MCParticle with ancestry identifiers, processes, four-vectors,
trajectory facts, active-volume entry/exit, and geometry-clipped length. It uses
ParticleInventoryService and reviewed geometry helpers. Output preserves raw and
clipped facts separately. Limitation: Geant track IDs are event-local and absent
parents or truncated trajectories are valid diagnostic states. Validation:
scaffold only.

### PDHDGeantProcessInventory

Normalizes direct parent-child edges and available process/trajectory transitions
so complete chains can be reconstructed later. It assumes deterministic, bounded
graph traversal and retains literal Geant process strings. Output includes edge
validity and generation depth. Limitation: high-level pion, kaon, capture, or
inelastic categories require separately documented derived definitions.
Validation: scaffold only.

### PDHDEnergyDepositInventory

Records or explicitly aggregates SimEnergyDeposit/IDE information with energy,
charge carriers, position, time, detector context, and G4 ancestry. It assumes a
production-available configured truth product and provides thinning controls.
Output labels deposited energy distinctly from kinetic and reconstructed energy.
Limitation: unrestricted deposit rows can be extremely large. Validation:
scaffold only.

### PDHDBeamInstrumentation

Records raw beam trigger quality, beam-track multiplicity, spectrometer momentum,
TOF, Cherenkov/PID information, and extrapolated position/direction with source
and validity. It uses ProtoDUNEBeamlineUtils only after checking PDHD data and MC
compatibility. Output is selection-independent. The implementation uses
`beamevent` for data and the simulated beam event under `generator` for MC,
requires a unique beam-event object before filling detailed fields, and leaves
MC PID disabled by default. Limitation: data and MC may expose different beam
products and calibrations. The preceding revision compiled and ran on PDHD data,
correctly reporting that
the tested full-reconstructed files lacked a stored `beamevent`. The data job now
schedules `pdhd_beamevent` through IFBeam; the MC job directly reads `generator`.
These revised sources and split jobs require a new build and runtime validation.
Startup logging records the selected data/MC tags and evaluation settings;
bounded event messages identify missing, non-unique, successfully extracted, or
failed beam records; `endJob` prints complete processing and trigger counts.
For the tested 2024 full-reconstructed data, decoded LLT frames are absent. The
data job uses the producer's supported `SkipLLT` mode. It records the exact
Cherenkov values passed by `ProtoDUNEBeamEvent` to `GetPID`, so the PID result
can be reproduced without adding an analyzer-side cut.
The tree records the exact selected product tag, source description, build
method, whether the product was made in the current job, separate Cherenkov
status/pressure sources and validity, PID method, and whether non-nominal input
handling was used. For this data path,
`non_nominal_input_handling_used=true` explicitly
identifies `SkipLLT`; it is a diagnostic alternative and is not equivalent to
an LLT-complete beam event. The `*_value_available` branches describe whether
values exist, while `*_provenance_valid` states whether their detector origin
has been validated. This distinction is required because `SkipLLT` can leave
zero-valued Cherenkov fields that `GetPID` still consumes.
`HasPerfectBeamMomentum` stores the result returned by the official
`ProtoDUNEBeamlineUtils::HasPerfectBeamMomentum` method; it is not the nominal
beam-momentum setting.
Trigger branches follow their official method names. For PDHD,
`IsGoodBeamlineTrigger` is the supported decision and combines
`GetTimingTrigger()==12` with `CheckIsMatched`. `GetBITrigger` is retained only
as a legacy diagnostic because the v10.17 PDHD producer explicitly sets it to
-1.
For MC, `mc_beam_pdg` stores the first generator primary from a non-cosmic
`MCTruth` record under `MCTruthTag`. The validity, product source, truth-record
count, and qualifying-primary count are separate branches, so missing or
ambiguous truth remains visible. Data stores the -999 sentinel with
`mc_beam_pdg_valid=false`.

### PDHDTimingInventory

Records event/trigger timing and all configured T0 candidates and associations.
It preserves native time units, source, confidence, and any derived drift-X as a
separate value. Output supports containment and blip-X auditing. Limitation:
production-dependent timing labels and conventions must be established before
implementation. Validation: scaffold only.

### PDHDHitInventory

Records compact recob::Hit observables and association keys, with independent PFP
metadata and hit-CNN outputs. It assumes configurable selection/thinning and never
saves waveforms. Optional MC contributions have explicit validity. Output enables
audits of calorimetry, blips, classifiers, and truth matching. Limitation: even
compact all-hit output can dominate file size. Validation: scaffold only.

### PDHDPFParticleInventory

Records the full Pandora hierarchy, primary/beam-slice metadata, TrackScore,
clear-cosmic and beam/cosmic scores, vertices, and associated-object counts. It
uses ProtoDUNEPFParticleUtils while preserving metadata provenance. Output contains
all PFPs without a physics cut. Limitation: missing or multiple track/shower
associations must remain explicit. Validation: scaffold only.

### PDHDTrackInventory

Records track geometry, trajectory validity, original endpoints/directions,
length, containment, and optional derived orientation. It assumes geometry-service
boundaries and configured fiducial margins in cm. Output is independent of truth,
calorimetry, and selection. Limitation: endpoint containment and full-trajectory
containment are different quantities. Validation: scaffold only.

### PDHDShowerInventory

Records shower start, direction, length/opening angle, best plane, native energy
and dE/dx vectors, errors, containment, and PFP link. It assumes the selected
PDHD reconstruction exposes these fields and does not infer shower status merely
from absent tracks. Calibrated features are separate. Limitation: reconstruction
products may omit or invalidate individual fields. Validation: scaffold only.

### PDHDTrackFeatureInventory

Records range/MCS momentum, ParticleID chi-square and NDF, PFP TrackScore, separate
weighted/unweighted CNN outputs, and calorimetry summaries with algorithm source
and validity. It assumes hypotheses and units are reviewed for the local release.
Output is intended for later cut development, not implicit acceptance. Limitation:
missing feature products must not trigger fallback overwrites. Validation:
scaffold only.

### PDHDShowerFeatureInventory

Records raw/calibrated shower energy, dE/dx, charge, topology, hit/space-point
counts, and separate classifier scores. It assumes HD-compatible calibration and
detector conditions. Output keeps SCE/no-SCE and calibrated/uncalibrated variants
independent. Limitation: SP calibration constants cannot be reused without HD
validation. Validation: scaffold only.

### PDHDCalorimetryPointInventory

Records one row per object, plane, and calorimetry point: XYZ, wire/tick, pitch,
dQ/dx, dE/dx, residual range, electric field, calibration variant, and validity.
It assumes vector sizes and units are checked by CalorimetryAlg. Output preserves
ordering and rejected-point diagnostics. Limitation: detailed rows increase output
volume but are necessary to diagnose bad-pitch energy tails. Validation:
scaffold only.

### PDHDRecoAssociation

Records normalized edges among PFPs, tracks, showers, clusters, hits, space points,
vertices, T0s, calorimetry, PID, and blips. It assumes producer-qualified stable
keys and validates every art association before indexing. Output retains all
multiplicities instead of choosing the first object. Limitation: association
existence does not establish physics correctness. Validation: scaffold only.

### PDHDRecoTruthMatch

Records all retained reco-to-truth and truth-to-reco matches, shared evidence,
purity, completeness, ranks, ambiguity, delta-ray policy, and method provenance.
It prefers configured standard BackTrackerMatchingData associations, then a
reviewed ProtoDUNETruthUtils/BackTracker path. Output includes staged truth-beam
reconstruction flags. Limitation: MC-only; all denominators require explicit
definitions. Validation: scaffold only.

### PDHDBlipInventory

Records blip position/uncertainty, detector context, size, charge, reconstructed
energy, hit counts, quality, T0/X source, and containment. It assumes the local
DUNE blip product API and calibration are verified. Output separates raw and
corrected coordinates. Limitation: reconstructed blip energy is not truth energy
or parent-particle kinetic energy. Validation: scaffold only.

### PDHDBlipAssociation

Records producer-supplied and analysis-derived candidate links from blips to PFPs,
tracks, showers, hits, and MC particles. It stores evidence, distances, closest
trajectory point, rank, method, ambiguity, and validity. Output never treats plot
proximity as a producer association. Limitation: truth edges exist only on MC and
geometric candidates require later performance validation. Validation: scaffold
only.

### PDHDBeamSelectionStages

Records every component of the nominal data-like beam selection: trigger quality,
beamline-track multiplicity, Pandora beam primary/type, official position and
direction residuals, candidate rank, final decision, and exact cut source/values.
The same reconstructed logic runs on data and MC; MC truth supplies efficiency,
purity, wrong-candidate, and staged reconstruction metrics afterward. Limitation:
official SP cuts require demonstrated PDHD applicability. The implementation
uses the measured or simulated beam-instrumentation track for the nominal match,
loads the official data window for the configured momentum, writes every stage,
and never reads truth. Validation: static inspection only; not compiled or run.

### PDHDCosmicSelectionStages

Records Pandora clear-cosmic/BDT diagnostics and independent start, end, both-end,
and full-trajectory containment decisions before any named final selection. It
assumes geometry boundaries and margins are explicit and permits diagnostic-only
operation with no rejection. MC output uses mutually exclusive truth categories.
Limitation: truth PDG/origin cannot enter a data selection. Validation: scaffold
only.

## Reference foundations

- DUNE `protoduneana` utilities and general PDSP/PDHD analyzer configuration.
- LArSoft ParticleInventoryService and BackTrackerService for simulation linkage.
- Standard `art::Assns` and `anab::BackTrackerMatchingData` products when present.

Exact APIs must be checked against the locally active `v10_17_01d00` software;
the online “latest” LArSoft documentation may describe a different release.

## Shared algorithms

| Algorithm | Implemented responsibility | Validation |
|---|---|---|
| `ProductInventoryAlg` | Uniform required/optional product state | Static only |
| `GeometryContainmentAlg` | Point/trajectory containment and segment clipping | Static only |
| `GeantHierarchyAlg` | Deterministic ancestry with missing/duplicate/cycle states | Static only |
| `RecoAssociationAlg` | Association cardinality and unique-object policy | Static only |
| `RecoTruthMatchAlg` | Ranked purity/completeness from supplied evidence | Static only |
| `CalorimetryAlg` | Bounds-checked `anab::Calorimetry` point extraction | Static only |
| `BeamInstrumentationAlg` | Common measured/simulated beam-event extraction | Static only |
| `BeamSelectionAlg` | Explicit data-like beam matching and conjunction | Static only |
| `CosmicSelectionAlg` | Configurable cosmic and containment conjunction | Static only |

The algorithms deliberately do not discover products or write TTrees. Modules
must obtain official art/LArSoft associations, record their source, and pass the
observed facts to these deterministic calculations. `BeamSelectionAlg` consumes
official data cut values supplied by FHiCL rather than calling the truth-assisted
MC branch of `ProtoDUNEBeamCuts::IsBeamlike`. `RecoTruthMatchAlg` does not replace
BackTracker or `ProtoDUNETruthUtils`; it standardizes metrics after those tools
provide hit, charge, or energy contributions.
