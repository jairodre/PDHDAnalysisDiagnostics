/**
 * @file PDHDBeamSelectionStages_module.cc
 * @brief Applies one staged, data-like beam-candidate definition to data and
 * MC.
 *
 * Data reads the reconstructed external beamline from DataBeamTag. MC reads
 * the generated beam primary from MCTruthTag and its Geant trajectory from
 * MCParticleTag; a stored MCBeamTag is retained separately for simulated
 * instrumentation diagnostics. The nominal decision requires the
 * corresponding unique beam reference, a Pandora beam-slice primary, and one
 * unambiguous reconstructed track or shower.
 * Beamline-to-TPC position and direction compatibility are stored as
 * observables, not cut by an inherited SP table. Tracks use a local TPC-entry
 * segment; showers use their reconstructed start and initial direction. MC
 * truth labels are diagnostic and do not select a particle species.
 * CandidateTree stores each selection component and match observable;
 * EventTree stores multiplicities plus one normalized all-track display
 * inventory. A candidate is selected only when exactly one passes the named
 * reconstruction stages.
 */
#include "TTree.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamInstrumentationAlg.h"
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamSelectionAlg.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace pdhd::diagnostics {
namespace {
// Candidate and reference-method codes are stable ROOT schema values.
constexpr int kCandidateNone = 0;
constexpr int kCandidateTrack = 1;
constexpr int kCandidateShower = 2;
constexpr int kCandidateAmbiguous = 3;
constexpr int kTPCReferenceUnavailable = 0;
// Code 1 is retained for already-produced ROOT files where Track::Vertex()
// was the position residual reference. New track rows use code 3, making the
// lower-Z entry definition explicit without reinterpreting old outputs.
constexpr int kTPCReferenceTrackVertexLocalDirection = 1;
constexpr int kTPCReferenceShowerStart = 2;
constexpr int kTPCReferenceTrackEntryLocalDirection = 3;
constexpr int kTruthCategoryUnavailable = 0;
constexpr int kTruthCategoryIncidentBeam = 1;
constexpr int kTruthCategoryCosmic = 2;
constexpr int kTruthCategoryOtherBeamPrimary = 3;
constexpr int kTruthCategoryPionInelastic = 4;
constexpr int kTruthCategoryDecay = 5;
constexpr int kTruthCategoryOther = 6;

int ClassifySelectedTrackTruth(bool const matchesIncidentBeam,
                               bool const originValid, int const origin,
                               int const beamOrigin,
                               int const cosmicOrigin,
                               std::string const &process) {
  if (matchesIncidentBeam) return kTruthCategoryIncidentBeam;
  if (originValid && origin == cosmicOrigin) return kTruthCategoryCosmic;
  if (process == "pi+Inelastic" || process == "pi-Inelastic") {
    return kTruthCategoryPionInelastic;
  }
  if (process == "Decay") return kTruthCategoryDecay;
  if (originValid && origin == beamOrigin && process == "primary") {
    return kTruthCategoryOtherBeamPrimary;
  }
  return originValid || !process.empty() ? kTruthCategoryOther
                                         : kTruthCategoryUnavailable;
}

template <class T> int IndexOf(std::vector<T> const &objects, T const *object) {
  if (!object)
    return -1;
  for (std::size_t i = 0; i < objects.size(); ++i)
    if (&objects[i] == object)
      return i;
  return -1;
}

// TTree keeps these module-owned addresses for the job lifetime.
template <class T>
void AddBranch(TTree &tree, char const *name, T &value) {
  tree.Branch(name, std::addressof(value));
}

struct TPCEntryDirection {
  bool valid = false;
  // NaN preserves the fact that no trajectory endpoint was available.
  Point3D entry{std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
  Direction3D direction;
  double sampledLengthCm = std::numeric_limits<double>::quiet_NaN();
  // True when the lower-Z entry corresponds to Track::End().  This lets the
  // nominal endpoint tangent be oriented into the TPC without guessing from
  // the fitted Vertex() ordering.
  bool entryAtTrackEnd = false;
};

// MC keeps generated-particle identity and transported Geant reference
// separate: the former labels the beam species, the latter defines entry
// position/direction for reco-to-truth residuals.
struct TruthBeamReference {
  bool truthProductAvailable = false;
  bool generatedPrimaryUnique = false;
  bool geantProductAvailable = false;
  bool geantMatchUnique = false;
  bool referenceValid = false;
  int generatedPrimaryCount = 0;
  int pdg = 0;
  // Generator and Geant TrackIDs are independent identifiers in this chain.
  int geantTrackId = -1;
  int trackId = -1;
  double initialMomentumGeV = std::numeric_limits<double>::quiet_NaN();
  Point3D start{std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
  Point3D entry{std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
  // The transported endpoint supports the MC reach-TPC/background study.
  Point3D end{std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN()};
  Direction3D direction;
  // Finite transported-Geant points in the configured TPC comparison range.
  // They are stored without interpolation or extrapolation.
  std::vector<double> trajectoryXcm;
  std::vector<double> trajectoryYcm;
  std::vector<double> trajectoryZcm;
  // Full finite path, including the simulated beamline upstream of the TPC.
  // It is an event-display diagnostic, distinct from the Z-window trajectory.
  bool fullTrajectoryValid = false;
  std::vector<double> fullTrajectoryXcm;
  std::vector<double> fullTrajectoryYcm;
  std::vector<double> fullTrajectoryZcm;
};



bool IsFinite(Point3D const &point) {
  return std::isfinite(point.x) && std::isfinite(point.y) &&
         std::isfinite(point.z);
}

bool IsValidDirection(Direction3D const &direction) {
  double const magnitude = std::sqrt(direction.x * direction.x +
                                     direction.y * direction.y +
                                     direction.z * direction.z);
  return std::isfinite(magnitude) && magnitude > 0.;
}

bool StoreGeantTrajectoryInZRange(simb::MCParticle const &particle,
                                  double const minimumZcm,
                                  double const maximumZcm,
                                  std::vector<double> &x,
                                  std::vector<double> &y,
                                  std::vector<double> &z) {
  x.clear();
  y.clear();
  z.clear();
  for (std::size_t index = 0; index < particle.NumberTrajectoryPoints();
       ++index) {
    auto const position = particle.Position(index);
    Point3D const point{position.X(), position.Y(), position.Z()};
    if (!IsFinite(point) || point.z < minimumZcm || point.z > maximumZcm) {
      continue;
    }
    x.push_back(point.x);
    y.push_back(point.y);
    z.push_back(point.z);
  }
  return !x.empty() && x.size() == y.size() && x.size() == z.size();
}

// Event-display provenance only: retain every finite transported-Geant point,
// including the simulated beamline before the TPC. This is kept in EventTree,
// not repeated for every primary PFP candidate row.
bool StoreFullGeantTrajectory(simb::MCParticle const &particle,
                              std::vector<double> &x,
                              std::vector<double> &y,
                              std::vector<double> &z) {
  x.clear();
  y.clear();
  z.clear();
  for (std::size_t index = 0; index < particle.NumberTrajectoryPoints();
       ++index) {
    auto const position = particle.Position(index);
    Point3D const point{position.X(), position.Y(), position.Z()};
    if (!IsFinite(point)) continue;
    x.push_back(point.x);
    y.push_back(point.y);
    z.push_back(point.z);
  }
  return !x.empty() && x.size() == y.size() && x.size() == z.size();
}

// The lower-Z endpoint is the TPC entry convention for this beam geometry.
// The first requested path length after that entry defines an optional local
// chord diagnostic; the nominal cosine uses the fitted endpoint tangent.
TPCEntryDirection FindTPCEntryDirection(recob::Track const &track,
                                        double const requestedLengthCm) {
  TPCEntryDirection result;
  std::size_t first = track.NumberTrajectoryPoints();
  std::size_t last = track.NumberTrajectoryPoints();
  for (std::size_t index = 0; index < track.NumberTrajectoryPoints(); ++index) {
    if (!track.HasValidPoint(index)) {
      continue;
    }
    if (first == track.NumberTrajectoryPoints()) {
      first = index;
    }
    last = index;
  }
  if (first == track.NumberTrajectoryPoints() || first == last) {
    return result;
  }

  auto const firstLocation = track.LocationAtPoint(first);
  auto const lastLocation = track.LocationAtPoint(last);
  bool const enterFromLast = lastLocation.Z() < firstLocation.Z();
  std::size_t const entryIndex = enterFromLast ? last : first;
  int const step = enterFromLast ? -1 : 1;
  auto const entryLocation = track.LocationAtPoint(entryIndex);
  result.entry = {entryLocation.X(), entryLocation.Y(), entryLocation.Z()};
  if (!IsFinite(result.entry)) {
    return result;
  }
  result.entryAtTrackEnd = enterFromLast;
  // Entry position remains useful even if no finite local chord can be made.
  if (!std::isfinite(requestedLengthCm) || requestedLengthCm <= 0.) {
    return result;
  }

  Point3D previous = result.entry;
  double accumulatedLengthCm = 0.;
  int const stop =
      enterFromLast ? static_cast<int>(first) : static_cast<int>(last);
  for (int index = static_cast<int>(entryIndex) + step;
       enterFromLast ? index >= stop : index <= stop; index += step) {
    std::size_t const trajectoryIndex = static_cast<std::size_t>(index);
    if (!track.HasValidPoint(trajectoryIndex)) {
      continue;
    }
    auto const location = track.LocationAtPoint(trajectoryIndex);
    Point3D const current{location.X(), location.Y(), location.Z()};
    if (!IsFinite(current)) {
      return TPCEntryDirection{};
    }
    double const stepLengthCm = std::sqrt(
        (current.x - previous.x) * (current.x - previous.x) +
        (current.y - previous.y) * (current.y - previous.y) +
        (current.z - previous.z) * (current.z - previous.z));
    if (!std::isfinite(stepLengthCm)) {
      return TPCEntryDirection{};
    }
    accumulatedLengthCm += stepLengthCm;
    if (accumulatedLengthCm >= requestedLengthCm) {
      result.direction = {current.x - result.entry.x, current.y - result.entry.y,
                          current.z - result.entry.z};
      double const chordLengthCm = std::sqrt(
          result.direction.x * result.direction.x +
          result.direction.y * result.direction.y +
          result.direction.z * result.direction.z);
      if (!std::isfinite(chordLengthCm) || chordLengthCm <= 0.) {
        return TPCEntryDirection{};
      }
      result.sampledLengthCm = accumulatedLengthCm;
      result.valid = true;
      return result;
    }
    previous = current;
  }
  // A short track has a valid entry position but no local-chord direction.
  return result;
}

// Match the SP BeamCuts convention: orient the reconstructed fitted endpoint
// tangent from the lower-Z end into the TPC.  This is the nominal cosine
// direction; the local chord remains an independent diagnostic below.
Direction3D OrientedTrackEntryDirection(recob::Track const &track,
                                        bool const entryAtTrackEnd) {
  auto const fittedDirection =
      entryAtTrackEnd ? track.EndDirection() : track.StartDirection();
  double const sign = entryAtTrackEnd ? -1. : 1.;
  return {sign * fittedDirection.X(), sign * fittedDirection.Y(),
          sign * fittedDirection.Z()};
}

// Preserve the full reconstructed trajectory in the same lower-Z-entry order
// used for the local direction. Invalid or non-finite points are not plotted.
std::vector<Point3D> ExtractTrackTrajectoryFromEntry(recob::Track const &track) {
  std::vector<Point3D> points;
  std::size_t first = track.NumberTrajectoryPoints();
  std::size_t last = track.NumberTrajectoryPoints();
  for (std::size_t index = 0; index < track.NumberTrajectoryPoints(); ++index) {
    if (!track.HasValidPoint(index)) continue;
    if (first == track.NumberTrajectoryPoints()) first = index;
    last = index;
  }
  if (first == track.NumberTrajectoryPoints()) return points;
  bool const enterFromLast = track.LocationAtPoint(last).Z() <
                             track.LocationAtPoint(first).Z();
  int const step = enterFromLast ? -1 : 1;
  int const stop = enterFromLast ? static_cast<int>(first)
                                 : static_cast<int>(last);
  for (int index = enterFromLast ? static_cast<int>(last)
                                 : static_cast<int>(first);
       enterFromLast ? index >= stop : index <= stop; index += step) {
    auto const location = track.LocationAtPoint(static_cast<std::size_t>(index));
    Point3D const point{location.X(), location.Y(), location.Z()};
    if (IsFinite(point)) points.push_back(point);
  }
  return points;
}

TruthBeamReference FindTruthBeamReference(
    art::Handle<std::vector<simb::MCTruth>> const &truthHandle,
    art::Handle<std::vector<simb::MCParticle>> const &particleHandle,
    int const requiredOrigin, double const energyToleranceGeV,
    double const minimumZcm, double const maximumZcm) {
  TruthBeamReference result;
  result.truthProductAvailable = bool(truthHandle);
  result.geantProductAvailable = bool(particleHandle);
  if (!truthHandle || !std::isfinite(energyToleranceGeV) ||
      energyToleranceGeV < 0. || !std::isfinite(minimumZcm) ||
      !std::isfinite(maximumZcm) || minimumZcm > maximumZcm) {
    return result;
  }

  // Keep a non-owning pointer only while the MCTruth handle is in scope.  The
  // selected generated primary supplies the beam identity (PDG, generator
  // TrackID, initial momentum, and energy used for the Geant lookup below).
  // Count every qualifying primary first: choosing the last one would hide an
  // ambiguous generator record behind a plausible-looking reference.
  simb::MCParticle const *generated = nullptr;
  for (auto const &truth : *truthHandle) {
    if (static_cast<int>(truth.Origin()) != requiredOrigin) continue;
    for (int index = 0; index < truth.NParticles(); ++index) {
      auto const &particle = truth.GetParticle(index);
      // "primary" identifies the generator-level beam particle within the
      // configured beam-origin MCTruth record; daughters are not beam seeds.
      if (particle.Process() != "primary") continue;
      ++result.generatedPrimaryCount;
      generated = &particle;
    }
  }
  result.generatedPrimaryUnique = result.generatedPrimaryCount == 1;
  if (!result.generatedPrimaryUnique) return result;

  result.pdg = generated->PdgCode();
  result.trackId = generated->TrackId();
  if (generated->NumberTrajectoryPoints() > 0) {
    result.initialMomentumGeV = generated->Momentum(0).P();
  }
  // Preserve the generator identity even when the Geant collection is absent;
  // only the transported entry reference then remains unavailable.
  if (!particleHandle) return result;

  simb::MCParticle const *geant = nullptr;
  unsigned int geantMatches = 0;
  for (auto const &particle : *particleHandle) {
    // ProtoDUNETruthUtils matches the generated beam primary to Geant by
    // PDG and initial energy. TrackID is not preserved between these products.
    if (particle.PdgCode() != result.pdg ||
        std::abs(particle.E() - generated->E()) > energyToleranceGeV) {
      continue;
    }
    geant = &particle;
    ++geantMatches;
  }
  result.geantMatchUnique = geantMatches == 1;
  if (!result.geantMatchUnique) return result;
  result.geantTrackId = geant->TrackId();
  result.fullTrajectoryValid = StoreFullGeantTrajectory(
      *geant, result.fullTrajectoryXcm, result.fullTrajectoryYcm,
      result.fullTrajectoryZcm);
  StoreGeantTrajectoryInZRange(*geant, std::max(0., minimumZcm), maximumZcm,
                               result.trajectoryXcm, result.trajectoryYcm,
                               result.trajectoryZcm);
  if (geant->NumberTrajectoryPoints() > 0) {
    auto const position = geant->Position(0);
    result.start = {position.X(), position.Y(), position.Z()};
    // Keep the final finite trajectory point separately from the entry point:
    // an interaction before Z=0 is physically meaningful, not missing data.
    for (std::size_t index = geant->NumberTrajectoryPoints(); index-- > 0;) {
      auto const endPosition = geant->Position(index);
      Point3D const end{endPosition.X(), endPosition.Y(), endPosition.Z()};
      if (IsFinite(end)) {
        result.end = end;
        break;
      }
    }
  }

  // The configured Z window declares the common PDHD comparison surface. It
  // is not an inferred active-volume test and therefore stays visible in FHiCL.
  for (std::size_t index = 0; index < geant->NumberTrajectoryPoints(); ++index) {
    auto const position = geant->Position(index);
    auto const momentum = geant->Momentum(index);
    if (position.Z() < minimumZcm || position.Z() > maximumZcm) continue;
    result.entry = {position.X(), position.Y(), position.Z()};
    result.direction = {momentum.Px(), momentum.Py(), momentum.Pz()};
    if (IsFinite(result.entry) && IsValidDirection(result.direction)) {
      result.referenceValid = true;
    }
    return result;
  }
  return result;
}
} // namespace
class PDHDBeamSelectionStages : public art::EDAnalyzer {
public:
  explicit PDHDBeamSelectionStages(fhicl::ParameterSet const &p)
      : EDAnalyzer(p), dataBeamTag_(p.get<art::InputTag>(
                           "DataBeamTag", art::InputTag{"beamevent"})),
        mcBeamTag_(p.get<art::InputTag>("MCBeamTag",
                                        art::InputTag{"generator"})),
        mcTruthTag_(
            p.get<art::InputTag>("MCTruthTag", art::InputTag{"generator"})),
        mcParticleTag_(p.get<art::InputTag>("MCParticleTag",
                                            art::InputTag{"largeant"})),
        pfpTag_(
            p.get<art::InputTag>("PFParticleTag", art::InputTag{"pandora"})),
        pfpMetadataTag_(p.get<art::InputTag>("PFParticleMetadataTag",
                                             pfpTag_)),
        pfpSliceTag_(p.get<art::InputTag>("PFParticleSliceTag", pfpTag_)),
        trackTag_(
            p.get<art::InputTag>("TrackTag", art::InputTag{"pandoraTrack"})),
        showerTag_(
            p.get<art::InputTag>("ShowerTag", art::InputTag{"pandoraShower"})),
        truthHitTag_(
            p.get<art::InputTag>("TruthHitTag", art::InputTag{"gaushit"})),
        beamInstrumentationNominalMomentumGeV_(
            p.get<double>("BeamInstrumentationNominalMomentumGeV")),
        evaluateDataPid_(p.get<bool>("EvaluateDataPID", false)),
        requireDataTrigger_(p.get<bool>("RequireGoodDataTrigger", true)),
        // The default is the named simb beam-MC convention rather than its
        // current integer representation (4). FHiCL may still override it
        // for a production that uses a different MCTruth origin.
        mcTruthBeamOrigin_(p.get<int>(
            "MCTruthBeamOrigin", static_cast<int>(simb::kSingleParticle))),
        mcCosmicOrigin_(p.get<int>("MCCosmicOrigin", 2)),
        mcTruthGeantEnergyToleranceGeV_(
            p.get<double>("MCTruthGeantEnergyToleranceGeV", 1.e-5)),
        mcTruthReferenceMinimumZcm_(p.get<double>("MCTruthReferenceMinimumZCm", 0.)),
        mcTruthReferenceMaximumZcm_(p.get<double>("MCTruthReferenceMaximumZCm", 600.)),
        tpcEntryDirectionLengthCm_(
            p.get<double>("TPCEntryDirectionLengthCm", 5.)),
        enableMCTrackTruthMatching_(
            p.get<bool>("EnableMCTrackTruthMatching", false)),
        printSummaryToStdout_(p.get<bool>("PrintSummaryToStdout", true)),
        maxEventMessages_(p.get<unsigned int>("MaxEventMessages", 3U)),
        beamline_(p.get<fhicl::ParameterSet>("BeamlineUtils")) {
    consumes<std::vector<beam::ProtoDUNEBeamEvent>>(dataBeamTag_);
    consumes<std::vector<beam::ProtoDUNEBeamEvent>>(mcBeamTag_);
    consumes<std::vector<simb::MCTruth>>(mcTruthTag_);
    consumes<std::vector<simb::MCParticle>>(mcParticleTag_);
    consumes<std::vector<recob::PFParticle>>(pfpTag_);
    consumes<std::vector<recob::Track>>(trackTag_);
    consumes<std::vector<recob::Hit>>(truthHitTag_);
    consumes<std::vector<recob::Shower>>(showerTag_);
  }
  void beginJob() override;
  void analyze(art::Event const &) override;
  void endJob() override;

private:
  // Event-local ART handles and associations used only while building primary
  // candidates.  This deliberately keeps ART pointers out of `out_`, whose
  // address is held by CandidateTree for the whole job.
  struct RecoPFPInputs {
    art::Handle<std::vector<recob::PFParticle>> pfps;
    art::Handle<std::vector<recob::Track>> tracks;
    art::Handle<std::vector<recob::Shower>> showers;
    std::unique_ptr<art::FindManyP<larpandoraobj::PFParticleMetadata>>
        metadata;
    std::unique_ptr<art::FindOneP<recob::Slice>> slices;
    std::unique_ptr<art::FindManyP<recob::Track>> trackAssociations;
    std::unique_ptr<art::FindManyP<recob::Shower>> showerAssociations;
  };

  // One transient Candidate is built for every primary PFParticle.  It is both
  // the source of one BeamSelectionCandidate ROOT row and the bookkeeping
  // record used after the loop to identify the unique selected candidate.
  // It owns values and collection indices only: pointers to ART products are
  // local to analyze() and are never stored in the ROOT-bound `out_` buffer.
  struct Candidate {
    // Stable collection indices and the reconstructed-object classification.
    int pfp = -1;
    int track = -1;
    int shower = -1;
    int recoCandidateType = kCandidateNone;
    int tpcReferenceMethod = kTPCReferenceUnavailable;
    // Association facts. `beamMetadata` is direct IsTestBeam metadata, while
    // `beamSlicePrimary` is the nominal requirement: this primary belongs to
    // the single Pandora beam-tagged slice found for the event.
    bool primary = false;
    bool beamMetadata = false;
    bool beamSlicePrimary = false;
    bool metadataAssociationValid = false;
    bool sliceAssociationValid = false;
    bool trackAssociationValid = false;
    bool showerAssociationValid = false;
    bool hasRecoTrack = false;
    bool hasRecoShower = false;
    bool hasUnambiguousRecoObject = false;
    bool recoObjectLengthValid = false;
    // Reconstructed-object classification and diagnostic validity. Length and
    // beam--TPC matching are recorded studies, not nominal selection stages.
    bool positionMatchValid = false;
    bool directionMatchValid = false;
    bool matchValid = false;
    // Named selection gates repeated in each candidate row. Beam-event and
    // beam-track multiplicity are event-level facts; they are copied here so
    // candidate rows can be studied without joining EventTree.
    bool triggerPass = false;
    bool beamEventPass = false;
    bool beamTrackPass = false;
    // `passesAllSelectionGates` means this primary passes every named gate,
    // including the event-level unique-beam-event gate. `selected` additionally
    // requires that it is the sole passing primary in the event.
    bool passesAllSelectionGates = false;
    bool selected = false;
    bool selectedTrack = false;
    bool selectedShower = false;
    // MC-only selected-track truth diagnostics use the official utility
    // matches. They retain the established ROOT schema and never enter
    // is_selected; unavailable custom-energy quantities remain NaN.
    bool recoTruthAssociationValid = false;
    bool recoTruthUtilityMatchValid = false;
    bool recoTruthUtilityOriginValid = false;
    bool recoTruthBackgroundCategoryValid = false;
    bool recoTruthUtilityMatchesBeamGeant = false;
    bool recoTruthHitMatchValid = false;
    bool recoTruthHitMatchesBeamGeant = false;
    bool recoTruthPurityValid = false;
    bool recoTruthCompletenessValid = false;
    int recoTruthHitTrackId = -1;
    int recoTruthHitPdg = 0;
    int recoTruthUtilityTrackId = -1;
    int recoTruthUtilityPdg = 0;
    int recoTruthUtilityOrigin = -1;
    int recoTruthBackgroundCategory = kTruthCategoryUnavailable;
    std::string recoTruthUtilityProcess;
    // Geant trajectory of the MCParticle matched to this selected reco track.
    // Points are retained only in the configured TPC comparison Z range.
    bool recoTruthUtilityTrajectoryValid = false;
    std::vector<double> recoTruthUtilityTrajectoryXcm;
    std::vector<double> recoTruthUtilityTrajectoryYcm;
    std::vector<double> recoTruthUtilityTrajectoryZcm;
    // Raw reconstruction-association multiplicities diagnose ambiguous rows.
    unsigned int trackAssociationCount = 0;
    unsigned int showerAssociationCount = 0;
    // Validity of the physical TPC reference used by match diagnostics. These
    // are not selection gates: a short track may still retain a valid lower-Z
    // trajectory-entry residual even when it has no local-chord measurement.
    bool recoStartValid = false;
    bool tpcEntryPositionValid = false;
    // Nominal track cosine uses the fitted endpoint tangent, oriented inward
    // from the lower-Z entry as in the SP BeamCuts implementation.
    bool tpcEntryEndpointDirectionValid = false;
    // The local chord is retained independently for direction-systematic
    // studies; it never replaces the nominal endpoint tangent.
    bool tpcEntryDirectionValid = false;
    bool localChordDirectionMatchValid = false;
    // Beam--TPC diagnostic observables. Position residuals require a reference
    // and reco start; cosine additionally requires both directions.
    double dx = std::numeric_limits<double>::quiet_NaN();
    double dy = std::numeric_limits<double>::quiet_NaN();
    double dz = std::numeric_limits<double>::quiet_NaN();
    double startz = std::numeric_limits<double>::quiet_NaN();
    double cos = std::numeric_limits<double>::quiet_NaN();
    double localChordDirectionCosine = std::numeric_limits<double>::quiet_NaN();
    // Track/Shower::Length is the native reconstructed-object length in cm.
    double recoObjectLengthCm = std::numeric_limits<double>::quiet_NaN();
    unsigned int recoTruthHitSharedCount = 0;
    unsigned int recoTruthHitSharedDeltaRayCount = 0;
    double recoTruthPurity = std::numeric_limits<double>::quiet_NaN();
    double recoTruthCompleteness = std::numeric_limits<double>::quiet_NaN();
    // Beam reference coordinates: data stores the fitted BI track; MC stores
    // the transported Geant reference and labels that source explicitly.
    double beamlineEndXcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndYcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndZcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndDirectionX = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndDirectionY = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndDirectionZ = std::numeric_limits<double>::quiet_NaN();
    double beamlineStartXcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineStartYcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineStartZcm = std::numeric_limits<double>::quiet_NaN();
    // Legacy reconstructed-object position: Track::Vertex() for tracks and
    // ShowerStart() for showers. It is preserved for direct comparison with
    // the physical lower-Z entry used by the nominal beam--TPC residual.
    double recoStartXcm = std::numeric_limits<double>::quiet_NaN();
    double recoStartYcm = std::numeric_limits<double>::quiet_NaN();
    double recoStartZcm = std::numeric_limits<double>::quiet_NaN();
    // The lower-Z trajectory point is the track beam--TPC position reference
    // and anchors its local entry direction. For showers it is ShowerStart().
    double tpcEntryXcm = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryYcm = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryZcm = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryEndpointDirectionX =
        std::numeric_limits<double>::quiet_NaN();
    double tpcEntryEndpointDirectionY =
        std::numeric_limits<double>::quiet_NaN();
    double tpcEntryEndpointDirectionZ =
        std::numeric_limits<double>::quiet_NaN();
    // Local lower-Z-entry chord direction, retained for direct comparison to
    // the nominal fitted endpoint tangent above.
    double tpcEntryDirectionX = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryDirectionY = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryDirectionZ = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryDirectionSampledLengthCm =
        std::numeric_limits<double>::quiet_NaN();
  };
  // Transient result shared by the all-track event inventory and the richer
  // selected-track annotation. MCParticle pointers remain event-local.
  struct TrackTruthMatch {
    bool evaluated = false;
    simb::MCParticle const *utilityParticle = nullptr;
    simb::MCParticle const *hitParticle = nullptr;
    bool originValid = false;
    int origin = -1;
    std::size_t sharedHitCount = 0;
    std::size_t sharedDeltaRayHitCount = 0;
  };
  // Initializes every event summary/output member before any product lookup.
  // This prevents an unavailable current-event quantity from inheriting a
  // plausible value from the preceding event.
  void ResetEventState(art::Event const &evt);
  // Data obtains the fitted instrumentation track. MC obtains its nominal
  // reference from generated identity matched to transported Geant; the
  // stored MC BeamEvent remains diagnostics-only.
  void AcquireBeamReference(art::Event const &evt,
                            BeamInstrumentationRecord &dataBeamReference);
  // Product and association availability are diagnostic states.  They do not
  // create fallback candidates when one of the inputs is missing.
  RecoPFPInputs AcquireRecoPFPInputs(art::Event const &evt);
  // IsTestBeam metadata identifies a single Pandora slice; candidate building
  // subsequently considers every primary associated with that slice.
  void FindPandoraBeamSlice(RecoPFPInputs const &inputs);
  // Builds one row per primary and evaluates only its named nominal gates.
  // Match residuals, length, PID, and truth labels remain diagnostics.
  std::vector<Candidate> BuildPrimaryCandidates(
      RecoPFPInputs const &inputs,
      BeamInstrumentationRecord const &dataBeamReference);
  // Applies event-wide uniqueness after all candidate-local gates are known.
  void SelectUniqueEventCandidate(std::vector<Candidate> &candidates);
  std::vector<TrackTruthMatch> BuildRecoTrackInventory(
      art::Event const &, RecoPFPInputs const &);
  // MC-only post-selection diagnostic. It labels only the selected track and
  // cannot alter the already-complete nominal decision.
  void AnnotateSelectedTrackTruth(
      art::Event const &evt, RecoPFPInputs const &inputs,
      std::vector<TrackTruthMatch> const &trackTruthMatches,
      std::vector<Candidate> &candidates);
  void LogSelectionResult();
  // CandidateTree binds module-owned `out_`, so copy each completed transient
  // candidate immediately before Fill rather than rebinding any addresses.
  void FillOutputTrees(std::vector<Candidate> const &candidates);
  art::InputTag dataBeamTag_, mcBeamTag_, mcTruthTag_, mcParticleTag_, pfpTag_,
      pfpMetadataTag_, pfpSliceTag_, trackTag_, showerTag_,
      truthHitTag_;
  double beamInstrumentationNominalMomentumGeV_;
  bool evaluateDataPid_;
  bool requireDataTrigger_;
  int mcTruthBeamOrigin_;
  int mcCosmicOrigin_;
  double mcTruthGeantEnergyToleranceGeV_;
  double mcTruthReferenceMinimumZcm_, mcTruthReferenceMaximumZcm_;
  double tpcEntryDirectionLengthCm_;
  bool enableMCTrackTruthMatching_;
  bool printSummaryToStdout_;
  unsigned int maxEventMessages_;
  protoana::ProtoDUNEBeamlineUtils beamline_;

  // Job counters expose product availability and selection outcomes.
  std::uint64_t eventsProcessed_ = 0, dataEvents_ = 0, mcEvents_ = 0;
  std::uint64_t beamProductAvailableEvents_ = 0, uniqueBeamEventEvents_ = 0;
  std::uint64_t uniqueBeamTrackEvents_ = 0, pfpAvailableEvents_ = 0;
  std::uint64_t selectedEvents_ = 0, ambiguousEvents_ = 0;
  unsigned int eventMessagesEmitted_ = 0;
  TTree *eventTree_ = nullptr, *candidateTree_ = nullptr;
  unsigned int run_ = 0, subrun_ = 0, event_ = 0;
  bool isData_ = false;
  bool beamProductAvailable_ = false, beamEventUnique_ = false,
       beamTrackUnique_ = false, triggerEvaluated_ = false,
       goodTrigger_ = false;
  bool triggerStageApplied_ = false;
  // These are diagnostic beam-instrumentation observables. They never enter
  // the staged candidate selection.
  bool beamInstrumentationPidEvaluated_ = false;
  std::vector<double> beamInstrumentationMomentaGeV_;
  std::vector<int> beamInstrumentationPidCandidates_;
  // MC-only diagnostics extracted from the stored generator BeamEvent.  They
  // are intentionally separate from the Geant truth reference used by the
  // nominal selection and beam--TPC residuals.
  bool simulatedInstrumentationProductAvailable_ = false;
  bool simulatedInstrumentationBeamEventUnique_ = false;
  bool simulatedInstrumentationBeamTrackUnique_ = false;
  // Preserve the stored simulated-BI start independently for future
  // beamline/Geant studies. It is diagnostic only and not drawn by default.
  bool simulatedInstrumentationStartValid_ = false;
  bool simulatedInstrumentationReferenceValid_ = false;
  Point3D simulatedInstrumentationStart_;
  Point3D simulatedInstrumentationEnd_;
  Direction3D simulatedInstrumentationEndDirection_;
  std::vector<double> simulatedInstrumentationMomentaGeV_;
  bool truthBeamProductAvailable_ = false, truthBeamPrimaryUnique_ = false,
       truthBeamGeantProductAvailable_ = false, truthBeamGeantMatchUnique_ = false,
       truthBeamReferenceValid_ = false;
  int truthBeamPrimaryCount_ = 0, truthBeamPdg_ = 0, truthBeamTrackId_ = -1,
      truthBeamGeantTrackId_ = -1;
  double truthBeamInitialMomentumGeV_ = std::numeric_limits<double>::quiet_NaN();
  Point3D truthBeamStart_{std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()};
  Point3D truthBeamEntry_{std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()};
  Point3D truthBeamEnd_{std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN()};
  Direction3D truthBeamDirection_{std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN()};
  std::vector<double> truthBeamTrajectoryXcm_, truthBeamTrajectoryYcm_,
      truthBeamTrajectoryZcm_;
  // MC-only full transported trajectory for event-level candidate displays.
  bool truthBeamFullTrajectoryValid_ = false;
  std::vector<double> truthBeamFullTrajectoryXcm_, truthBeamFullTrajectoryYcm_,
      truthBeamFullTrajectoryZcm_;
  bool truthGeantParticleInventoryAvailable_ = false;
  std::vector<int> truthGeantTrackIds_, truthGeantMotherTrackIds_,
      truthGeantPdgs_;
  std::vector<unsigned long long> truthGeantTrajectoryOffsets_;
  std::vector<double> truthGeantTrajectoryXcm_, truthGeantTrajectoryYcm_,
      truthGeantTrajectoryZcm_;
  bool pfpMetadataAssociationAvailable_ = false;
  bool pfpSliceAssociationAvailable_ = false;
  bool pfpTrackAssociationAvailable_ = false;
  bool pfpShowerAssociationAvailable_ = false;
  bool pandoraBeamSliceFound_ = false, pandoraBeamSliceAmbiguous_ = false;
  int pandoraBeamSliceId_ = -1;
  unsigned int pandoraBeamSliceCount_ = 0;
  bool hasSelectedCandidate_ = false, selectionAmbiguous_ = false;
  int nPrimary_ = 0, nBeamPrimary_ = 0, nBeamSlicePrimary_ = 0,
      nPassing_ = 0, nPassingTrack_ = 0, nPassingShower_ = 0,
      selectedPfp_ = -1, selectedTrack_ = -1, selectedShower_ = -1,
      selectedCandidateType_ = kCandidateNone;
  // One normalized, event-level inventory owns every reco-track trajectory.
  // Offsets has Ntracks+1 entries; track i occupies [offsets[i], offsets[i+1]).
  std::vector<int> recoTrackIndices_, recoTrackInBeamSlice_;
  std::vector<int> recoTrackCosmicClassificationValid_, recoTrackIsCosmic_;
  std::vector<int> recoTrackTruthUtilityMatchValid_, recoTrackTruthUtilityPdg_;
  std::vector<int> recoTrackTruthHitMatchValid_, recoTrackTruthHitPdg_;
  std::vector<unsigned long long> recoTrackTrajectoryOffsets_;
  std::vector<double> recoTrackTrajectoryXcm_, recoTrackTrajectoryYcm_,
      recoTrackTrajectoryZcm_;
  Candidate out_;
  // Residuals are retained for HD cut development, not applied in this module.
  bool beamTpcMatchUsedForSelection_ = false;
  std::string referenceSource_ = "unavailable";
};
void PDHDBeamSelectionStages::beginJob() {
  mf::LogInfo("PDHDBeamSelectionStages")
      << "Starting staged beam selection. DataBeamTag='"
      << dataBeamTag_.encode() << "', MCTruthTag='" << mcTruthTag_.encode()
      << "', MCParticleTag='" << mcParticleTag_.encode()
      << "', PFParticleTag='" << pfpTag_.encode() << "', TrackTag='"
      << trackTag_.encode() << "', ShowerTag='" << showerTag_.encode()
      << "', PFParticleMetadataTag='" << pfpMetadataTag_.encode()
      << "', PFParticleSliceTag='" << pfpSliceTag_.encode()
      << "', BeamInstrumentationNominalMomentumGeV="
      << beamInstrumentationNominalMomentumGeV_
      << ", EvaluateDataPID=" << evaluateDataPid_
      << ", RequireGoodDataTrigger=" << requireDataTrigger_
      << ", MCTruthBeamOrigin=" << mcTruthBeamOrigin_
      << ", MCCosmicOrigin=" << mcCosmicOrigin_
      << ", MCTruthGeantEnergyToleranceGeV=" << mcTruthGeantEnergyToleranceGeV_
      << ", MCTruthReferenceZCm=[" << mcTruthReferenceMinimumZcm_ << ", "
      << mcTruthReferenceMaximumZcm_ << "]"
      << ", TPCEntryDirectionLengthCm=" << tpcEntryDirectionLengthCm_
      << ".";

  auto fs = art::ServiceHandle<art::TFileService>();
  eventTree_ =
      fs->make<TTree>("BeamSelectionEvent", "Beam selection event summary");
  // TFileService writes the complete trees at endJob. Disable ROOT's periodic
  // checkpoint keys so a successful job leaves one current tree cycle.
  eventTree_->SetAutoSave(0);
  AddBranch(*eventTree_, "run", run_);
  AddBranch(*eventTree_, "subrun", subrun_);
  AddBranch(*eventTree_, "event", event_);
  AddBranch(*eventTree_, "is_data", isData_);
  AddBranch(*eventTree_, "beam_product_available", beamProductAvailable_);
  AddBranch(*eventTree_, "beam_event_unique", beamEventUnique_);
  AddBranch(*eventTree_, "beam_track_unique", beamTrackUnique_);
  AddBranch(*eventTree_, "trigger_evaluated", triggerEvaluated_);
  AddBranch(*eventTree_, "good_trigger", goodTrigger_);
  AddBranch(*eventTree_, "trigger_stage_applied", triggerStageApplied_);
  // Data-only diagnostic values; vectors retain every beam-instrumentation
  // hypothesis instead of silently selecting one PID or momentum.
  AddBranch(*eventTree_, "beam_instrumentation_pid_evaluated",
            beamInstrumentationPidEvaluated_);
  AddBranch(*eventTree_, "beam_instrumentation_momenta_GeV",
            beamInstrumentationMomentaGeV_);
  AddBranch(*eventTree_, "beam_instrumentation_pid_candidates",
            beamInstrumentationPidCandidates_);
  // Stored MC BeamEvent diagnostics are kept distinct from the Geant truth
  // reference that defines the MC beam--TPC comparison.
  AddBranch(*eventTree_, "simulated_instrumentation_product_available",
            simulatedInstrumentationProductAvailable_);
  AddBranch(*eventTree_, "simulated_instrumentation_beam_event_unique",
            simulatedInstrumentationBeamEventUnique_);
  AddBranch(*eventTree_, "simulated_instrumentation_beam_track_unique",
            simulatedInstrumentationBeamTrackUnique_);
  AddBranch(*eventTree_, "simulated_instrumentation_reference_valid",
            simulatedInstrumentationReferenceValid_);
  AddBranch(*eventTree_, "simulated_instrumentation_start_valid",
            simulatedInstrumentationStartValid_);
  AddBranch(*eventTree_, "simulated_instrumentation_start_x_cm",
            simulatedInstrumentationStart_.x);
  AddBranch(*eventTree_, "simulated_instrumentation_start_y_cm",
            simulatedInstrumentationStart_.y);
  AddBranch(*eventTree_, "simulated_instrumentation_start_z_cm",
            simulatedInstrumentationStart_.z);
  AddBranch(*eventTree_, "simulated_instrumentation_end_x_cm",
            simulatedInstrumentationEnd_.x);
  AddBranch(*eventTree_, "simulated_instrumentation_end_y_cm",
            simulatedInstrumentationEnd_.y);
  AddBranch(*eventTree_, "simulated_instrumentation_end_z_cm",
            simulatedInstrumentationEnd_.z);
  AddBranch(*eventTree_, "simulated_instrumentation_end_direction_x",
            simulatedInstrumentationEndDirection_.x);
  AddBranch(*eventTree_, "simulated_instrumentation_end_direction_y",
            simulatedInstrumentationEndDirection_.y);
  AddBranch(*eventTree_, "simulated_instrumentation_end_direction_z",
            simulatedInstrumentationEndDirection_.z);
  AddBranch(*eventTree_, "simulated_instrumentation_momenta_GeV",
            simulatedInstrumentationMomentaGeV_);
  // MC truth fields are unavailable on data; their validity fields preserve
  // that distinction rather than encoding a data value as a truth result.
  AddBranch(*eventTree_, "truth_beam_product_available",
            truthBeamProductAvailable_);
  AddBranch(*eventTree_, "truth_beam_primary_unique",
            truthBeamPrimaryUnique_);
  AddBranch(*eventTree_, "truth_beam_geant_product_available",
            truthBeamGeantProductAvailable_);
  AddBranch(*eventTree_, "truth_beam_geant_match_unique",
            truthBeamGeantMatchUnique_);
  AddBranch(*eventTree_, "truth_beam_reference_valid",
            truthBeamReferenceValid_);
  AddBranch(*eventTree_, "truth_beam_primary_count", truthBeamPrimaryCount_);
  AddBranch(*eventTree_, "truth_beam_pdg", truthBeamPdg_);
  // The original branch is the generated primary TrackID; keep the Geant ID
  // separately so their intentionally different provenance is visible.
  AddBranch(*eventTree_, "truth_beam_track_id", truthBeamTrackId_);
  AddBranch(*eventTree_, "truth_beam_geant_track_id", truthBeamGeantTrackId_);
  AddBranch(*eventTree_, "truth_beam_initial_momentum_GeV",
            truthBeamInitialMomentumGeV_);
  AddBranch(*eventTree_, "truth_beam_start_x_cm", truthBeamStart_.x);
  AddBranch(*eventTree_, "truth_beam_start_y_cm", truthBeamStart_.y);
  AddBranch(*eventTree_, "truth_beam_start_z_cm", truthBeamStart_.z);
  AddBranch(*eventTree_, "truth_beam_entry_x_cm", truthBeamEntry_.x);
  AddBranch(*eventTree_, "truth_beam_entry_y_cm", truthBeamEntry_.y);
  AddBranch(*eventTree_, "truth_beam_entry_z_cm", truthBeamEntry_.z);
  AddBranch(*eventTree_, "truth_beam_end_x_cm", truthBeamEnd_.x);
  AddBranch(*eventTree_, "truth_beam_end_y_cm", truthBeamEnd_.y);
  AddBranch(*eventTree_, "truth_beam_end_z_cm", truthBeamEnd_.z);
  AddBranch(*eventTree_, "truth_beam_direction_x", truthBeamDirection_.x);
  AddBranch(*eventTree_, "truth_beam_direction_y", truthBeamDirection_.y);
  AddBranch(*eventTree_, "truth_beam_direction_z", truthBeamDirection_.z);
  AddBranch(*eventTree_, "truth_beam_full_trajectory_valid",
            truthBeamFullTrajectoryValid_);
  AddBranch(*eventTree_, "truth_beam_full_trajectory_x_cm",
            truthBeamFullTrajectoryXcm_);
  AddBranch(*eventTree_, "truth_beam_full_trajectory_y_cm",
            truthBeamFullTrajectoryYcm_);
  AddBranch(*eventTree_, "truth_beam_full_trajectory_z_cm",
            truthBeamFullTrajectoryZcm_);
  // Event-local Geant inventory for direct-daughter topology studies. This is
  // separate from the incident-beam reference above; offsets has N+1 entries.
  AddBranch(*eventTree_, "truth_geant_particle_inventory_available",
            truthGeantParticleInventoryAvailable_);
  AddBranch(*eventTree_, "truth_geant_track_ids", truthGeantTrackIds_);
  AddBranch(*eventTree_, "truth_geant_mother_track_ids",
            truthGeantMotherTrackIds_);
  AddBranch(*eventTree_, "truth_geant_pdgs", truthGeantPdgs_);
  AddBranch(*eventTree_, "truth_geant_trajectory_offsets",
            truthGeantTrajectoryOffsets_);
  AddBranch(*eventTree_, "truth_geant_trajectory_x_cm",
            truthGeantTrajectoryXcm_);
  AddBranch(*eventTree_, "truth_geant_trajectory_y_cm",
            truthGeantTrajectoryYcm_);
  AddBranch(*eventTree_, "truth_geant_trajectory_z_cm",
            truthGeantTrajectoryZcm_);
  AddBranch(*eventTree_, "pfp_metadata_association_available",
            pfpMetadataAssociationAvailable_);
  AddBranch(*eventTree_, "pfp_slice_association_available",
            pfpSliceAssociationAvailable_);
  AddBranch(*eventTree_, "pfp_track_association_available",
            pfpTrackAssociationAvailable_);
  AddBranch(*eventTree_, "pfp_shower_association_available",
            pfpShowerAssociationAvailable_);
  AddBranch(*eventTree_, "pandora_beam_slice_found", pandoraBeamSliceFound_);
  AddBranch(*eventTree_, "pandora_beam_slice_ambiguous",
            pandoraBeamSliceAmbiguous_);
  AddBranch(*eventTree_, "pandora_beam_slice_id", pandoraBeamSliceId_);
  AddBranch(*eventTree_, "pandora_beam_slice_count", pandoraBeamSliceCount_);
  AddBranch(*eventTree_, "n_primary_pfp", nPrimary_);
  // This counts direct IsTestBeam metadata; the next branch is the seed.
  AddBranch(*eventTree_, "n_beam_primary_pfp", nBeamPrimary_);
  AddBranch(*eventTree_, "n_beam_slice_primary_pfp", nBeamSlicePrimary_);
  AddBranch(*eventTree_, "n_passing", nPassing_);
  AddBranch(*eventTree_, "n_passing_track", nPassingTrack_);
  AddBranch(*eventTree_, "n_passing_shower", nPassingShower_);
  AddBranch(*eventTree_, "has_selected_candidate", hasSelectedCandidate_);
  AddBranch(*eventTree_, "selection_ambiguous", selectionAmbiguous_);
  AddBranch(*eventTree_, "selected_pfp_index", selectedPfp_);
  AddBranch(*eventTree_, "selected_track_index", selectedTrack_);
  AddBranch(*eventTree_, "selected_shower_index", selectedShower_);
  AddBranch(*eventTree_, "selected_candidate_type", selectedCandidateType_);
  AddBranch(*eventTree_, "beam_reference_source", referenceSource_);
  // This is a measurement definition, not a numerical selection threshold.
  AddBranch(*eventTree_, "tpc_entry_direction_length_cm",
            tpcEntryDirectionLengthCm_);
  AddBranch(*eventTree_, "beam_tpc_match_used_for_selection",
            beamTpcMatchUsedForSelection_);
  // Display-only all-track inventory. Trajectories are flattened so the
  // event tree owns one copy without requiring nested-vector dictionaries.
  AddBranch(*eventTree_, "reco_track_indices", recoTrackIndices_);
  AddBranch(*eventTree_, "reco_track_in_pandora_beam_slice",
            recoTrackInBeamSlice_);
  AddBranch(*eventTree_, "reco_track_cosmic_classification_valid",
            recoTrackCosmicClassificationValid_);
  AddBranch(*eventTree_, "reco_track_is_cosmic", recoTrackIsCosmic_);
  AddBranch(*eventTree_, "reco_track_truth_utility_match_valid",
            recoTrackTruthUtilityMatchValid_);
  AddBranch(*eventTree_, "reco_track_truth_utility_pdg",
            recoTrackTruthUtilityPdg_);
  AddBranch(*eventTree_, "reco_track_truth_hit_match_valid",
            recoTrackTruthHitMatchValid_);
  AddBranch(*eventTree_, "reco_track_truth_hit_pdg",
            recoTrackTruthHitPdg_);
  AddBranch(*eventTree_, "reco_track_trajectory_offsets",
            recoTrackTrajectoryOffsets_);
  AddBranch(*eventTree_, "reco_track_trajectory_x_cm",
            recoTrackTrajectoryXcm_);
  AddBranch(*eventTree_, "reco_track_trajectory_y_cm",
            recoTrackTrajectoryYcm_);
  AddBranch(*eventTree_, "reco_track_trajectory_z_cm",
            recoTrackTrajectoryZcm_);
  candidateTree_ = fs->make<TTree>("BeamSelectionCandidate",
                                   "One row per primary PFParticle");
  candidateTree_->SetAutoSave(0);
  AddBranch(*candidateTree_, "run", run_);
  AddBranch(*candidateTree_, "subrun", subrun_);
  AddBranch(*candidateTree_, "event", event_);
  AddBranch(*candidateTree_, "is_data", isData_);
  AddBranch(*candidateTree_, "beam_instrumentation_pid_evaluated",
            beamInstrumentationPidEvaluated_);
  AddBranch(*candidateTree_, "beam_instrumentation_momenta_GeV",
            beamInstrumentationMomentaGeV_);
  AddBranch(*candidateTree_, "beam_instrumentation_pid_candidates",
            beamInstrumentationPidCandidates_);
  AddBranch(*candidateTree_, "simulated_instrumentation_product_available",
            simulatedInstrumentationProductAvailable_);
  AddBranch(*candidateTree_, "simulated_instrumentation_beam_event_unique",
            simulatedInstrumentationBeamEventUnique_);
  AddBranch(*candidateTree_, "simulated_instrumentation_beam_track_unique",
            simulatedInstrumentationBeamTrackUnique_);
  AddBranch(*candidateTree_, "simulated_instrumentation_reference_valid",
            simulatedInstrumentationReferenceValid_);
  AddBranch(*candidateTree_, "simulated_instrumentation_start_valid",
            simulatedInstrumentationStartValid_);
  AddBranch(*candidateTree_, "simulated_instrumentation_start_x_cm",
            simulatedInstrumentationStart_.x);
  AddBranch(*candidateTree_, "simulated_instrumentation_start_y_cm",
            simulatedInstrumentationStart_.y);
  AddBranch(*candidateTree_, "simulated_instrumentation_start_z_cm",
            simulatedInstrumentationStart_.z);
  AddBranch(*candidateTree_, "simulated_instrumentation_end_x_cm",
            simulatedInstrumentationEnd_.x);
  AddBranch(*candidateTree_, "simulated_instrumentation_end_y_cm",
            simulatedInstrumentationEnd_.y);
  AddBranch(*candidateTree_, "simulated_instrumentation_end_z_cm",
            simulatedInstrumentationEnd_.z);
  AddBranch(*candidateTree_, "simulated_instrumentation_end_direction_x",
            simulatedInstrumentationEndDirection_.x);
  AddBranch(*candidateTree_, "simulated_instrumentation_end_direction_y",
            simulatedInstrumentationEndDirection_.y);
  AddBranch(*candidateTree_, "simulated_instrumentation_end_direction_z",
            simulatedInstrumentationEndDirection_.z);
  AddBranch(*candidateTree_, "simulated_instrumentation_momenta_GeV",
            simulatedInstrumentationMomentaGeV_);
  // Event-level generated beam identity is repeated per candidate so a single
  // candidate tree supports data/MC overlay studies without a fragile join.
  AddBranch(*candidateTree_, "truth_beam_product_available",
            truthBeamProductAvailable_);
  AddBranch(*candidateTree_, "truth_beam_primary_unique",
            truthBeamPrimaryUnique_);
  AddBranch(*candidateTree_, "truth_beam_geant_product_available",
            truthBeamGeantProductAvailable_);
  AddBranch(*candidateTree_, "truth_beam_geant_match_unique",
            truthBeamGeantMatchUnique_);
  AddBranch(*candidateTree_, "truth_beam_reference_valid",
            truthBeamReferenceValid_);
  AddBranch(*candidateTree_, "truth_beam_pdg", truthBeamPdg_);
  AddBranch(*candidateTree_, "truth_beam_track_id", truthBeamTrackId_);
  AddBranch(*candidateTree_, "truth_beam_geant_track_id", truthBeamGeantTrackId_);
  AddBranch(*candidateTree_, "truth_beam_initial_momentum_GeV",
            truthBeamInitialMomentumGeV_);
  AddBranch(*candidateTree_, "truth_beam_start_x_cm", truthBeamStart_.x);
  AddBranch(*candidateTree_, "truth_beam_start_y_cm", truthBeamStart_.y);
  AddBranch(*candidateTree_, "truth_beam_start_z_cm", truthBeamStart_.z);
  AddBranch(*candidateTree_, "truth_beam_entry_x_cm", truthBeamEntry_.x);
  AddBranch(*candidateTree_, "truth_beam_entry_y_cm", truthBeamEntry_.y);
  AddBranch(*candidateTree_, "truth_beam_entry_z_cm", truthBeamEntry_.z);
  AddBranch(*candidateTree_, "truth_beam_end_x_cm", truthBeamEnd_.x);
  AddBranch(*candidateTree_, "truth_beam_end_y_cm", truthBeamEnd_.y);
  AddBranch(*candidateTree_, "truth_beam_end_z_cm", truthBeamEnd_.z);
  AddBranch(*candidateTree_, "truth_beam_direction_x", truthBeamDirection_.x);
  AddBranch(*candidateTree_, "truth_beam_direction_y", truthBeamDirection_.y);
  AddBranch(*candidateTree_, "truth_beam_direction_z", truthBeamDirection_.z);
  AddBranch(*candidateTree_, "truth_beam_trajectory_x_cm",
            truthBeamTrajectoryXcm_);
  AddBranch(*candidateTree_, "truth_beam_trajectory_y_cm",
            truthBeamTrajectoryYcm_);
  AddBranch(*candidateTree_, "truth_beam_trajectory_z_cm",
            truthBeamTrajectoryZcm_);
  // MC-only selected-track truth diagnostics. Branch names are retained while
  // their matching source is the official ProtoDUNETruthUtils API.
  AddBranch(*candidateTree_, "reco_truth_association_valid",
            out_.recoTruthAssociationValid);
  AddBranch(*candidateTree_, "reco_truth_utility_match_valid",
            out_.recoTruthUtilityMatchValid);
  AddBranch(*candidateTree_, "reco_truth_utility_origin_valid",
            out_.recoTruthUtilityOriginValid);
  AddBranch(*candidateTree_, "reco_truth_utility_origin",
            out_.recoTruthUtilityOrigin);
  AddBranch(*candidateTree_, "reco_truth_utility_process",
            out_.recoTruthUtilityProcess);
  AddBranch(*candidateTree_, "reco_truth_background_category_valid",
            out_.recoTruthBackgroundCategoryValid);
  AddBranch(*candidateTree_, "reco_truth_background_category",
            out_.recoTruthBackgroundCategory);
  AddBranch(*candidateTree_, "reco_truth_utility_track_id",
            out_.recoTruthUtilityTrackId);
  AddBranch(*candidateTree_, "reco_truth_utility_pdg",
            out_.recoTruthUtilityPdg);
  AddBranch(*candidateTree_, "reco_truth_utility_trajectory_valid",
            out_.recoTruthUtilityTrajectoryValid);
  AddBranch(*candidateTree_, "reco_truth_utility_trajectory_x_cm",
            out_.recoTruthUtilityTrajectoryXcm);
  AddBranch(*candidateTree_, "reco_truth_utility_trajectory_y_cm",
            out_.recoTruthUtilityTrajectoryYcm);
  AddBranch(*candidateTree_, "reco_truth_utility_trajectory_z_cm",
            out_.recoTruthUtilityTrajectoryZcm);
  AddBranch(*candidateTree_, "reco_truth_utility_matches_beam_geant",
            out_.recoTruthUtilityMatchesBeamGeant);
  AddBranch(*candidateTree_, "reco_truth_hit_match_valid",
            out_.recoTruthHitMatchValid);
  AddBranch(*candidateTree_, "reco_truth_hit_track_id",
            out_.recoTruthHitTrackId);
  AddBranch(*candidateTree_, "reco_truth_hit_pdg", out_.recoTruthHitPdg);
  AddBranch(*candidateTree_, "reco_truth_hit_shared_count",
            out_.recoTruthHitSharedCount);
  AddBranch(*candidateTree_, "reco_truth_hit_shared_delta_ray_count",
            out_.recoTruthHitSharedDeltaRayCount);
  AddBranch(*candidateTree_, "reco_truth_hit_matches_beam_geant",
            out_.recoTruthHitMatchesBeamGeant);
  AddBranch(*candidateTree_, "reco_truth_purity_valid",
            out_.recoTruthPurityValid);
  AddBranch(*candidateTree_, "reco_truth_purity", out_.recoTruthPurity);
  AddBranch(*candidateTree_, "reco_truth_completeness_valid",
            out_.recoTruthCompletenessValid);
  AddBranch(*candidateTree_, "reco_truth_completeness",
            out_.recoTruthCompleteness);
  AddBranch(*candidateTree_, "pfp_index", out_.pfp);
  AddBranch(*candidateTree_, "track_index", out_.track);
  AddBranch(*candidateTree_, "shower_index", out_.shower);
  // 0=none, 1=track, 2=shower, 3=ambiguous associations.
  AddBranch(*candidateTree_, "reco_candidate_type", out_.recoCandidateType);
  // 0=unavailable, 1=legacy Track::Vertex/local direction,
  // 2=ShowerStart/Direction, 3=lower-Z track entry/fitted endpoint tangent.
  AddBranch(*candidateTree_, "tpc_reference_method", out_.tpcReferenceMethod);
  AddBranch(*candidateTree_, "is_primary", out_.primary);
  // Direct IsTestBeam metadata remains diagnostic; the slice is nominal.
  AddBranch(*candidateTree_, "is_pandora_beam_primary", out_.beamMetadata);
  AddBranch(*candidateTree_, "is_pandora_beam_slice_primary",
            out_.beamSlicePrimary);
  AddBranch(*candidateTree_, "metadata_association_valid",
            out_.metadataAssociationValid);
  AddBranch(*candidateTree_, "slice_association_valid",
            out_.sliceAssociationValid);
  AddBranch(*candidateTree_, "track_association_valid",
            out_.trackAssociationValid);
  AddBranch(*candidateTree_, "shower_association_valid",
            out_.showerAssociationValid);
  AddBranch(*candidateTree_, "track_association_count",
            out_.trackAssociationCount);
  AddBranch(*candidateTree_, "shower_association_count",
            out_.showerAssociationCount);
  AddBranch(*candidateTree_, "has_reco_track", out_.hasRecoTrack);
  AddBranch(*candidateTree_, "has_reco_shower", out_.hasRecoShower);
  AddBranch(*candidateTree_, "has_unambiguous_reco_object",
            out_.hasUnambiguousRecoObject);
  AddBranch(*candidateTree_, "reco_object_length_valid",
            out_.recoObjectLengthValid);
  AddBranch(*candidateTree_, "reco_object_length_cm", out_.recoObjectLengthCm);
  // This is the named reconstruction stage used by the nominal conjunction.
  AddBranch(*candidateTree_, "passes_unambiguous_reco_object",
            out_.hasUnambiguousRecoObject);
  AddBranch(*candidateTree_, "passes_pandora_beam_slice_primary",
            out_.beamSlicePrimary);
  AddBranch(*candidateTree_, "passes_trigger_stage", out_.triggerPass);
  AddBranch(*candidateTree_, "trigger_stage_applied", triggerStageApplied_);
  AddBranch(*candidateTree_, "passes_beam_event_unique", out_.beamEventPass);
  AddBranch(*candidateTree_, "passes_beam_track_unique", out_.beamTrackPass);
  AddBranch(*candidateTree_, "position_match_valid", out_.positionMatchValid);
  AddBranch(*candidateTree_, "direction_match_valid", out_.directionMatchValid);
  // `match_valid` is the stricter conjunction retained for compatibility.
  AddBranch(*candidateTree_, "match_valid", out_.matchValid);
  // These explicit coordinates make the beamline-end to lower-Z TPC-entry
  // residuals independently reproducible without inferring either point.
  AddBranch(*candidateTree_, "beamline_end_x_cm", out_.beamlineEndXcm);
  AddBranch(*candidateTree_, "beamline_end_y_cm", out_.beamlineEndYcm);
  AddBranch(*candidateTree_, "beamline_end_z_cm", out_.beamlineEndZcm);
  AddBranch(*candidateTree_, "beamline_end_direction_x",
            out_.beamlineEndDirectionX);
  AddBranch(*candidateTree_, "beamline_end_direction_y",
            out_.beamlineEndDirectionY);
  AddBranch(*candidateTree_, "beamline_end_direction_z",
            out_.beamlineEndDirectionZ);
  // The fitted beamline start/end define the measured incoming trajectory.
  AddBranch(*candidateTree_, "beamline_start_x_cm", out_.beamlineStartXcm);
  AddBranch(*candidateTree_, "beamline_start_y_cm", out_.beamlineStartYcm);
  AddBranch(*candidateTree_, "beamline_start_z_cm", out_.beamlineStartZcm);
  AddBranch(*candidateTree_, "reco_start_valid", out_.recoStartValid);
  AddBranch(*candidateTree_, "reco_start_x_cm", out_.recoStartXcm);
  AddBranch(*candidateTree_, "reco_start_y_cm", out_.recoStartYcm);
  AddBranch(*candidateTree_, "reco_start_z_cm", out_.recoStartZcm);
  // Track lower-Z entry is the nominal beam--TPC position reference; for
  // showers it coincides with ShowerStart(). `reco_start_*` remains the
  // independent Track::Vertex() diagnostic for tracks.
  AddBranch(*candidateTree_, "tpc_entry_position_valid",
            out_.tpcEntryPositionValid);
  AddBranch(*candidateTree_, "tpc_entry_x_cm", out_.tpcEntryXcm);
  AddBranch(*candidateTree_, "tpc_entry_y_cm", out_.tpcEntryYcm);
  AddBranch(*candidateTree_, "tpc_entry_z_cm", out_.tpcEntryZcm);
  // The nominal cosine uses the beamline direction and the oriented fitted
  // Track endpoint tangent (or Shower::Direction).  The local track chord is
  // retained in distinct branches for robustness comparisons.
  AddBranch(*candidateTree_, "tpc_entry_endpoint_direction_valid",
            out_.tpcEntryEndpointDirectionValid);
  AddBranch(*candidateTree_, "tpc_entry_endpoint_direction_x",
            out_.tpcEntryEndpointDirectionX);
  AddBranch(*candidateTree_, "tpc_entry_endpoint_direction_y",
            out_.tpcEntryEndpointDirectionY);
  AddBranch(*candidateTree_, "tpc_entry_endpoint_direction_z",
            out_.tpcEntryEndpointDirectionZ);
  AddBranch(*candidateTree_, "tpc_entry_direction_valid",
            out_.tpcEntryDirectionValid);
  AddBranch(*candidateTree_, "tpc_entry_direction_sampled_length_cm",
            out_.tpcEntryDirectionSampledLengthCm);
  AddBranch(*candidateTree_, "tpc_entry_direction_x",
            out_.tpcEntryDirectionX);
  AddBranch(*candidateTree_, "tpc_entry_direction_y",
            out_.tpcEntryDirectionY);
  AddBranch(*candidateTree_, "tpc_entry_direction_z",
            out_.tpcEntryDirectionZ);
  AddBranch(*candidateTree_, "local_chord_direction_match_valid",
            out_.localChordDirectionMatchValid);
  AddBranch(*candidateTree_, "delta_x_cm", out_.dx);
  AddBranch(*candidateTree_, "delta_y_cm", out_.dy);
  // Difference of the two measured/reconstructed endpoints, not absolute Z.
  AddBranch(*candidateTree_, "delta_z_cm", out_.dz);
  AddBranch(*candidateTree_, "entrance_z_cm", out_.startz);
  AddBranch(*candidateTree_, "direction_cosine", out_.cos);
  AddBranch(*candidateTree_, "local_chord_direction_cosine",
            out_.localChordDirectionCosine);
  AddBranch(*candidateTree_, "is_selected", out_.selected);
  AddBranch(*candidateTree_, "is_selected_track", out_.selectedTrack);
  AddBranch(*candidateTree_, "is_selected_shower", out_.selectedShower);
  AddBranch(*candidateTree_, "beam_reference_source", referenceSource_);

}
void PDHDBeamSelectionStages::ResetEventState(art::Event const &evt) {
  ++eventsProcessed_;
  run_ = evt.run();
  subrun_ = evt.subRun();
  event_ = evt.event();
  isData_ = evt.isRealData();
  if (isData_) {
    ++dataEvents_;
  } else {
    ++mcEvents_;
  }
  beamProductAvailable_ = beamEventUnique_ = beamTrackUnique_ = false;
  pfpMetadataAssociationAvailable_ = false;
  pfpSliceAssociationAvailable_ = false;
  pfpTrackAssociationAvailable_ = false;
  pfpShowerAssociationAvailable_ = false;
  pandoraBeamSliceFound_ = pandoraBeamSliceAmbiguous_ = false;
  pandoraBeamSliceId_ = -1;
  pandoraBeamSliceCount_ = 0;
  triggerEvaluated_ = goodTrigger_ = hasSelectedCandidate_ =
      selectionAmbiguous_ = false;
  triggerStageApplied_ = isData_ && requireDataTrigger_;
  beamInstrumentationPidEvaluated_ = false;
  beamInstrumentationMomentaGeV_.clear();
  beamInstrumentationPidCandidates_.clear();
  simulatedInstrumentationProductAvailable_ = false;
  simulatedInstrumentationBeamEventUnique_ = false;
  simulatedInstrumentationBeamTrackUnique_ = false;
  simulatedInstrumentationStartValid_ = false;
  simulatedInstrumentationReferenceValid_ = false;
  simulatedInstrumentationStart_ = {
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  simulatedInstrumentationEnd_ = {
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  simulatedInstrumentationEndDirection_ = {
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  simulatedInstrumentationMomentaGeV_.clear();
  nPrimary_ = nBeamPrimary_ = nBeamSlicePrimary_ = nPassing_ = 0;
  nPassingTrack_ = nPassingShower_ = 0;
  selectedPfp_ = selectedTrack_ = selectedShower_ = -1;
  selectedCandidateType_ = kCandidateNone;
  referenceSource_ = "unavailable";
  truthBeamProductAvailable_ = truthBeamPrimaryUnique_ = false;
  truthBeamGeantProductAvailable_ = truthBeamGeantMatchUnique_ = false;
  truthBeamReferenceValid_ = false;
  truthBeamPrimaryCount_ = 0;
  truthBeamPdg_ = 0;
  truthBeamTrackId_ = -1;
  truthBeamGeantTrackId_ = -1;
  truthBeamInitialMomentumGeV_ = std::numeric_limits<double>::quiet_NaN();
  truthBeamStart_ = {std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN()};
  truthBeamEntry_ = truthBeamStart_;
  truthBeamEnd_ = truthBeamStart_;
  truthBeamDirection_ = {std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN()};
  truthBeamTrajectoryXcm_.clear();
  truthBeamTrajectoryYcm_.clear();
  truthBeamTrajectoryZcm_.clear();
  truthBeamFullTrajectoryValid_ = false;
  truthBeamFullTrajectoryXcm_.clear();
  truthBeamFullTrajectoryYcm_.clear();
  truthBeamFullTrajectoryZcm_.clear();
  truthGeantParticleInventoryAvailable_ = false;
  truthGeantTrackIds_.clear();
  truthGeantMotherTrackIds_.clear();
  truthGeantPdgs_.clear();
  truthGeantTrajectoryOffsets_.clear();
  truthGeantTrajectoryXcm_.clear();
  truthGeantTrajectoryYcm_.clear();
  truthGeantTrajectoryZcm_.clear();
  recoTrackIndices_.clear();
  recoTrackInBeamSlice_.clear();
  recoTrackCosmicClassificationValid_.clear();
  recoTrackIsCosmic_.clear();
  recoTrackTruthUtilityMatchValid_.clear();
  recoTrackTruthUtilityPdg_.clear();
  recoTrackTruthHitMatchValid_.clear();
  recoTrackTruthHitPdg_.clear();
  recoTrackTrajectoryOffsets_.clear();
  recoTrackTrajectoryXcm_.clear();
  recoTrackTrajectoryYcm_.clear();
  recoTrackTrajectoryZcm_.clear();
}

void PDHDBeamSelectionStages::AcquireBeamReference(
    art::Event const &evt, BeamInstrumentationRecord &bi) {
  // Phase 1: obtain the event-level beam reference. Data uses the fitted BI
  // track; MC uses generated-primary plus transported-Geant truth. The stored
  // simulated BeamEvent is retained separately for diagnostics only.
  if (isData_) {
    auto const beamHandle =
        evt.getHandle<std::vector<beam::ProtoDUNEBeamEvent>>(dataBeamTag_);
    beamProductAvailable_ = bool(beamHandle);
    beamEventUnique_ = beamHandle && beamHandle->size() == 1;
    beamProductAvailableEvents_ += beamProductAvailable_;
    uniqueBeamEventEvents_ += beamEventUnique_;
    if (beamEventUnique_) {
      bi = BeamInstrumentationAlg::Extract(
          beamHandle->front(), beamline_, BeamReferenceSource::DataInstrumentation,
          beamInstrumentationNominalMomentumGeV_, requireDataTrigger_,
          evaluateDataPid_, {});
      triggerEvaluated_ = bi.triggerEvaluated;
      goodTrigger_ = bi.goodTrigger;
      beamInstrumentationPidEvaluated_ = evaluateDataPid_;
      beamInstrumentationMomentaGeV_ = bi.momentaGeV;
      beamInstrumentationPidCandidates_ = bi.pidCandidates;
      beamTrackUnique_ = bi.tracks.size() == 1;
      uniqueBeamTrackEvents_ += beamTrackUnique_;
      referenceSource_ = "data_instrumentation_fitted_end";
    }
  } else {
    // The generated BeamEvent is retained as a diagnostic-only simulated
    // instrumentation reference.  Do not use it to replace the Geant truth
    // reference or alter any nominal MC selection gate.
    auto const simulatedBeamHandle =
        evt.getHandle<std::vector<beam::ProtoDUNEBeamEvent>>(mcBeamTag_);
    simulatedInstrumentationProductAvailable_ = bool(simulatedBeamHandle);
    simulatedInstrumentationBeamEventUnique_ =
        simulatedBeamHandle && simulatedBeamHandle->size() == 1;
    if (simulatedInstrumentationBeamEventUnique_) {
      auto const simulatedInstrumentation = BeamInstrumentationAlg::Extract(
          simulatedBeamHandle->front(), beamline_,
          BeamReferenceSource::SimulatedInstrumentation,
          beamInstrumentationNominalMomentumGeV_, false, false, {});
      simulatedInstrumentationMomentaGeV_ =
          simulatedInstrumentation.momentaGeV;
      simulatedInstrumentationBeamTrackUnique_ =
          simulatedInstrumentation.tracks.size() == 1;
      if (simulatedInstrumentationBeamTrackUnique_) {
        simulatedInstrumentationStart_ =
            simulatedInstrumentation.tracks.front().start;
        simulatedInstrumentationEnd_ =
            simulatedInstrumentation.tracks.front().end;
        simulatedInstrumentationEndDirection_ =
            simulatedInstrumentation.tracks.front().endDirection;
        simulatedInstrumentationReferenceValid_ =
            IsFinite(simulatedInstrumentationEnd_);
        simulatedInstrumentationStartValid_ =
            IsFinite(simulatedInstrumentationStart_);
      }
    }
    auto const truthHandle = evt.getHandle<std::vector<simb::MCTruth>>(mcTruthTag_);
    auto const particleHandle =
        evt.getHandle<std::vector<simb::MCParticle>>(mcParticleTag_);
    truthGeantParticleInventoryAvailable_ = bool(particleHandle);
    truthGeantTrajectoryOffsets_.push_back(0);
    if (particleHandle) {
      for (auto const &particle : *particleHandle) {
        truthGeantTrackIds_.push_back(particle.TrackId());
        truthGeantMotherTrackIds_.push_back(particle.Mother());
        truthGeantPdgs_.push_back(particle.PdgCode());
        for (std::size_t index = 0;
             index < particle.NumberTrajectoryPoints(); ++index) {
          auto const position = particle.Position(index);
          Point3D const point{position.X(), position.Y(), position.Z()};
          if (!IsFinite(point)) continue;
          truthGeantTrajectoryXcm_.push_back(point.x);
          truthGeantTrajectoryYcm_.push_back(point.y);
          truthGeantTrajectoryZcm_.push_back(point.z);
        }
        truthGeantTrajectoryOffsets_.push_back(
            static_cast<unsigned long long>(truthGeantTrajectoryXcm_.size()));
      }
    }
    auto const truth = FindTruthBeamReference(
        truthHandle, particleHandle, mcTruthBeamOrigin_,
        mcTruthGeantEnergyToleranceGeV_,
        mcTruthReferenceMinimumZcm_, mcTruthReferenceMaximumZcm_);
    truthBeamProductAvailable_ = truth.truthProductAvailable;
    truthBeamPrimaryUnique_ = truth.generatedPrimaryUnique;
    truthBeamGeantProductAvailable_ = truth.geantProductAvailable;
    truthBeamGeantMatchUnique_ = truth.geantMatchUnique;
    truthBeamReferenceValid_ = truth.referenceValid;
    truthBeamPrimaryCount_ = truth.generatedPrimaryCount;
    truthBeamPdg_ = truth.pdg;
    truthBeamTrackId_ = truth.trackId;
    truthBeamGeantTrackId_ = truth.geantTrackId;
    truthBeamInitialMomentumGeV_ = truth.initialMomentumGeV;
    truthBeamStart_ = truth.start;
    truthBeamEntry_ = truth.entry;
    truthBeamEnd_ = truth.end;
    truthBeamDirection_ = truth.direction;
    truthBeamTrajectoryXcm_ = truth.trajectoryXcm;
    truthBeamTrajectoryYcm_ = truth.trajectoryYcm;
    truthBeamTrajectoryZcm_ = truth.trajectoryZcm;
    truthBeamFullTrajectoryValid_ = truth.fullTrajectoryValid;
    truthBeamFullTrajectoryXcm_ = truth.fullTrajectoryXcm;
    truthBeamFullTrajectoryYcm_ = truth.fullTrajectoryYcm;
    truthBeamFullTrajectoryZcm_ = truth.fullTrajectoryZcm;
    // The truth reference identifies a generated beam event; it never applies
    // a truth-PDG species cut to the reconstructed candidate.
    beamProductAvailable_ = truthBeamProductAvailable_;
    beamEventUnique_ = truthBeamPrimaryUnique_;
    beamTrackUnique_ = truthBeamReferenceValid_;
    // MC has no data trigger. The truth-primary requirement is exposed through
    // truth_beam_primary_unique and the generic beam-reference stages instead.
    triggerEvaluated_ = false;
    goodTrigger_ = false;
    beamProductAvailableEvents_ += beamProductAvailable_;
    uniqueBeamEventEvents_ += beamEventUnique_;
    uniqueBeamTrackEvents_ += beamTrackUnique_;
    if (truthBeamReferenceValid_) referenceSource_ = "mc_truth_geant_entry";
  }
}

PDHDBeamSelectionStages::RecoPFPInputs
PDHDBeamSelectionStages::AcquireRecoPFPInputs(art::Event const &evt) {
  // Phase 2: obtain reconstructed PFP collections and their association
  // helpers. Missing products remain explicit event/candidate status states.
  RecoPFPInputs inputs;
  // Short aliases make the ART retrieval below match the names used later in
  // the candidate builder. They refer to the same event-local storage.
  auto &ph = inputs.pfps;
  auto &th = inputs.tracks;
  auto &sh = inputs.showers;
  auto &metadata = inputs.metadata;
  auto &slices = inputs.slices;
  auto &tracks = inputs.trackAssociations;
  auto &showers = inputs.showerAssociations;
  ph = evt.getHandle<std::vector<recob::PFParticle>>(pfpTag_);
  th = evt.getHandle<std::vector<recob::Track>>(trackTag_);
  sh = evt.getHandle<std::vector<recob::Shower>>(showerTag_);
  pfpAvailableEvents_ += bool(ph);
  if ((!beamProductAvailable_ || !ph) &&
      eventMessagesEmitted_ < maxEventMessages_) {
    mf::LogWarning("PDHDBeamSelectionStages")
        << "Required input missing for run " << run_ << ", subrun " << subrun_
        << ", event " << event_ << ": beam reference product available="
        << beamProductAvailable_ << " from '"
        << (isData_ ? dataBeamTag_ : mcTruthTag_).encode()
        << "', PFParticle product available=" << bool(ph) << " from '"
        << pfpTag_.encode()
        << "'. The event summary is retained with no selected candidate.";
    ++eventMessagesEmitted_;
  }

  // Query each PFP association by collection index, not PFParticle::Self().
  // This follows Pandora's association products while remaining valid when
  // Self identifiers are not contiguous collection indices.
  if (ph) {
    try {
      metadata = std::make_unique<art::FindManyP<larpandoraobj::PFParticleMetadata>>(
          ph, evt, pfpMetadataTag_);
      pfpMetadataAssociationAvailable_ = true;
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_ < maxEventMessages_) {
        mf::LogWarning("PDHDBeamSelectionStages")
            << "No usable PFP-to-metadata association from tag '"
            << pfpMetadataTag_.encode() << "' for run " << run_ << ", subrun "
            << subrun_ << ", event " << event_ << ": " << error.what();
        ++eventMessagesEmitted_;
      }
    }
    try {
      slices = std::make_unique<art::FindOneP<recob::Slice>>(
          ph, evt, pfpSliceTag_);
      pfpSliceAssociationAvailable_ = true;
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_ < maxEventMessages_) {
        mf::LogWarning("PDHDBeamSelectionStages")
            << "No usable PFP-to-slice association from tag '"
            << pfpSliceTag_.encode() << "' for run " << run_ << ", subrun "
            << subrun_ << ", event " << event_ << ": " << error.what();
        ++eventMessagesEmitted_;
      }
    }
    try {
      if (th) {
        tracks = std::make_unique<art::FindManyP<recob::Track>>(
            ph, evt, trackTag_);
        pfpTrackAssociationAvailable_ = true;
      }
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_ < maxEventMessages_) {
        mf::LogWarning("PDHDBeamSelectionStages")
            << "No usable PFP-to-track association from tag '"
            << trackTag_.encode() << "' for run " << run_ << ", subrun "
            << subrun_ << ", event " << event_ << ": " << error.what();
        ++eventMessagesEmitted_;
      }
    }
    try {
      if (sh) {
        showers = std::make_unique<art::FindManyP<recob::Shower>>(
            ph, evt, showerTag_);
        pfpShowerAssociationAvailable_ = true;
      }
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_ < maxEventMessages_) {
        mf::LogWarning("PDHDBeamSelectionStages")
            << "No usable PFP-to-shower association from tag '"
            << showerTag_.encode() << "' for run " << run_ << ", subrun "
            << subrun_ << ", event " << event_ << ": " << error.what();
        ++eventMessagesEmitted_;
      }
    }
  }

  return inputs;
}

void PDHDBeamSelectionStages::FindPandoraBeamSlice(
    RecoPFPInputs const &inputs) {
  // Phase 3: identify the single Pandora beam slice. IsTestBeam metadata
  // identifies the beam slice. The nominal TPC candidate is every primary PFP
  // associated with that selected slice.
  auto const &ph = inputs.pfps;
  auto const &metadata = inputs.metadata;
  auto const &slices = inputs.slices;
  std::set<int> beamSliceIds;
  if (ph && metadata && slices) {
    for (std::size_t index = 0; index < ph->size(); ++index) {
      auto const &pfp = ph->at(index);
      if (!pfp.IsPrimary()) {
        continue;
      }
      try {
        auto const &metadataForPfp = metadata->at(index);
        auto const slice = slices->at(index);
        if (metadataForPfp.empty() || !slice) {
          continue;
        }
        auto const properties = metadataForPfp.front()->GetPropertiesMap();
        if (properties.find("IsTestBeam") != properties.end()) {
          beamSliceIds.insert(slice->ID());
        }
      } catch (std::exception const &error) {
        if (eventMessagesEmitted_ < maxEventMessages_) {
          mf::LogWarning("PDHDBeamSelectionStages")
              << "Could not query Pandora metadata/slice for PFP index "
              << index << " in run " << run_ << ", subrun " << subrun_
              << ", event " << event_ << ": " << error.what();
          ++eventMessagesEmitted_;
        }
      }
    }
  }
  pandoraBeamSliceCount_ = beamSliceIds.size();
  pandoraBeamSliceFound_ = pandoraBeamSliceCount_ == 1;
  pandoraBeamSliceAmbiguous_ = pandoraBeamSliceCount_ > 1;
  if (pandoraBeamSliceFound_) {
    pandoraBeamSliceId_ = *beamSliceIds.begin();
  }
}

std::vector<PDHDBeamSelectionStages::Candidate>
PDHDBeamSelectionStages::BuildPrimaryCandidates(
    RecoPFPInputs const &inputs,
    BeamInstrumentationRecord const &bi) {
  // Phase 4: build one diagnostic/selection record for each primary PFP. The
  // vector is intentionally retained until every primary has been counted: a
  // candidate cannot receive is_selected until event-wide uniqueness is known.
  auto const &ph = inputs.pfps;
  auto const &th = inputs.tracks;
  auto const &sh = inputs.showers;
  auto const &metadata = inputs.metadata;
  auto const &slices = inputs.slices;
  auto const &tracks = inputs.trackAssociations;
  auto const &showers = inputs.showerAssociations;
  std::vector<Candidate> candidates;
  if (ph) {
    for (std::size_t i = 0; i < ph->size(); ++i) {
      auto const &p = ph->at(i);
      if (!p.IsPrimary()) {
        continue;
      }

      // A Candidate starts from this primary PFP. It remains in `candidates`
      // even when an association is missing or ambiguous so cut-flow studies
      // can distinguish rejection from absent reconstruction information.
      Candidate c;
      c.pfp = i;
      c.primary = true;
      ++nPrimary_;
      // First establish the two related but distinct Pandora facts: direct
      // IsTestBeam metadata, and membership in the event's unique beam slice.
      if (metadata) {
        try {
          auto const &metadataForPfp = metadata->at(i);
          c.metadataAssociationValid = true;
          if (!metadataForPfp.empty()) {
            auto const properties = metadataForPfp.front()->GetPropertiesMap();
            c.beamMetadata = properties.find("IsTestBeam") != properties.end();
          }
        } catch (std::exception const &) {
          c.metadataAssociationValid = false;
        }
      }
      if (c.beamMetadata) {
        ++nBeamPrimary_;
      }
      if (slices) {
        try {
          auto const slice = slices->at(i);
          c.sliceAssociationValid = true;
          c.beamSlicePrimary =
              pandoraBeamSliceFound_ && slice &&
              static_cast<int>(slice->ID()) == pandoraBeamSliceId_;
        } catch (std::exception const &) {
          c.sliceAssociationValid = false;
        }
      }
      if (c.beamSlicePrimary) {
        ++nBeamSlicePrimary_;
      }
      // Resolve all associated object counts before choosing the single object
      // used for geometry. A first pointer is retained only to recover its
      // collection index; it never resolves an ambiguous association.
      recob::Track const *tr = nullptr;
      recob::Shower const *sw = nullptr;
      if (tracks) {
        try {
          auto const &associatedTracks = tracks->at(i);
          c.trackAssociationValid = true;
          c.trackAssociationCount = associatedTracks.size();
          if (!associatedTracks.empty()) {
            tr = associatedTracks.front().get();
          }
        } catch (std::exception const &) {
          c.trackAssociationValid = false;
        }
      }
      c.track = th ? IndexOf(*th, tr) : -1;
      if (showers) {
        try {
          auto const &associatedShowers = showers->at(i);
          c.showerAssociationValid = true;
          c.showerAssociationCount = associatedShowers.size();
          if (!associatedShowers.empty()) {
            sw = associatedShowers.front().get();
          }
        } catch (std::exception const &) {
          c.showerAssociationValid = false;
        }
      }
      c.shower = sh ? IndexOf(*sh, sw) : -1;
      c.hasRecoTrack = c.trackAssociationCount > 0;
      c.hasRecoShower = c.showerAssociationCount > 0;
      // Do not silently prioritize an object type or the first of many
      // associations: such PFPs remain explicit ambiguous diagnostic rows.
      if (c.trackAssociationCount == 1 && c.showerAssociationCount == 0) {
        c.recoCandidateType = kCandidateTrack;
        c.hasUnambiguousRecoObject = true;
        c.recoObjectLengthCm = tr->Length();
      } else if (c.trackAssociationCount == 0 &&
                 c.showerAssociationCount == 1) {
        c.recoCandidateType = kCandidateShower;
        c.hasUnambiguousRecoObject = true;
        c.recoObjectLengthCm = sw->Length();
      } else if (c.hasRecoTrack || c.hasRecoShower) {
        c.recoCandidateType = kCandidateAmbiguous;
      }
      // Length is useful for diagnostics, but a finite length is deliberately
      // not part of the unambiguous-track/shower selection definition.
      c.recoObjectLengthValid =
          c.hasUnambiguousRecoObject && std::isfinite(c.recoObjectLengthCm) &&
          c.recoObjectLengthCm >= 0.;
      // Fill the generic beam--TPC match input after the PFP's object type and
      // references are known. Its residuals are diagnostic observables only.
      BeamMatchInput mi;
      mi.beamReferenceValid = beamTrackUnique_;
      // Store the exact reference coordinates that define this row's residuals.
      // The source is data BI in data and transported Geant truth in MC.
      if (beamTrackUnique_) {
        if (isData_) {
          auto const &bt = bi.tracks.front();
          mi.beamPositionAtReference = bt.end;
          mi.beamDirection = bt.endDirection;
          c.beamlineEndXcm = bt.end.x;
          c.beamlineEndYcm = bt.end.y;
          c.beamlineEndZcm = bt.end.z;
          c.beamlineEndDirectionX = bt.endDirection.x;
          c.beamlineEndDirectionY = bt.endDirection.y;
          c.beamlineEndDirectionZ = bt.endDirection.z;
          c.beamlineStartXcm = bt.start.x;
          c.beamlineStartYcm = bt.start.y;
          c.beamlineStartZcm = bt.start.z;
          mi.referenceSource = BeamReferenceSource::DataInstrumentation;
        } else {
          // MC has no measured beamline fit. Store the transported truth entry
          // in the same reference-coordinate branches for overlay plots; the
          // explicit source branch prevents it being mistaken for data BI.
          mi.beamPositionAtReference = truthBeamEntry_;
          mi.beamDirection = truthBeamDirection_;
          c.beamlineEndXcm = truthBeamEntry_.x;
          c.beamlineEndYcm = truthBeamEntry_.y;
          c.beamlineEndZcm = truthBeamEntry_.z;
          c.beamlineEndDirectionX = truthBeamDirection_.x;
          c.beamlineEndDirectionY = truthBeamDirection_.y;
          c.beamlineEndDirectionZ = truthBeamDirection_.z;
          c.beamlineStartXcm = truthBeamStart_.x;
          c.beamlineStartYcm = truthBeamStart_.y;
          c.beamlineStartZcm = truthBeamStart_.z;
          // MC is compared to a transported Geant truth point, never to a
          // measured or reconstructed beam-instrumentation track.
          mi.referenceSource = BeamReferenceSource::TruthProjection;
        }
      }
      // Next derive the TPC reference from the unique reconstructed object.
      // No reference is invented for ambiguous or absent reco associations.
      if (c.recoCandidateType == kCandidateTrack) {
        auto const localDirection =
            FindTPCEntryDirection(*tr, tpcEntryDirectionLengthCm_);
        c.tpcReferenceMethod = kTPCReferenceTrackEntryLocalDirection;
        auto const vertex = tr->Vertex();
        c.recoStartXcm = vertex.X();
        c.recoStartYcm = vertex.Y();
        c.recoStartZcm = vertex.Z();
        c.recoStartValid =
            IsFinite({c.recoStartXcm, c.recoStartYcm, c.recoStartZcm});
        c.tpcEntryXcm = localDirection.entry.x;
        c.tpcEntryYcm = localDirection.entry.y;
        c.tpcEntryZcm = localDirection.entry.z;
        c.tpcEntryPositionValid = IsFinite(localDirection.entry);
        Direction3D endpointDirection{
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
        if (c.tpcEntryPositionValid) {
          endpointDirection =
              OrientedTrackEntryDirection(*tr, localDirection.entryAtTrackEnd);
          c.tpcEntryEndpointDirectionX = endpointDirection.x;
          c.tpcEntryEndpointDirectionY = endpointDirection.y;
          c.tpcEntryEndpointDirectionZ = endpointDirection.z;
          c.tpcEntryEndpointDirectionValid =
              IsValidDirection(endpointDirection);
        }
        c.tpcEntryDirectionValid = localDirection.valid;
        c.tpcEntryDirectionSampledLengthCm =
            localDirection.sampledLengthCm;
        // Use the same lower-Z physical entry for position residuals and the
        // fitted endpoint tangent. Track::Vertex() remains in reco_start_* as
        // a diagnostic of the reconstruction's vertex convention.
        mi.recoStartValid = c.tpcEntryPositionValid;
        mi.recoStart = {c.tpcEntryXcm, c.tpcEntryYcm, c.tpcEntryZcm};
        mi.recoDirectionValid = c.tpcEntryEndpointDirectionValid;
        if (c.tpcEntryEndpointDirectionValid) {
          mi.recoDirection = endpointDirection;
        }
        // A short track has no local-chord direction, but its lower-Z entry
        // residual and fitted endpoint tangent remain separately available.
        if (localDirection.valid) {
          c.tpcEntryDirectionX = localDirection.direction.x;
          c.tpcEntryDirectionY = localDirection.direction.y;
          c.tpcEntryDirectionZ = localDirection.direction.z;
        }
      } else if (c.recoCandidateType == kCandidateShower) {
        // ShowerStart/Direction are Pandora's initial shower quantities; no
        // arbitrary short-segment re-fit is substituted for a shower.
        auto const start = sw->ShowerStart();
        auto const direction = sw->Direction();
        Direction3D const initialDirection{direction.X(), direction.Y(),
                                           direction.Z()};
        c.tpcReferenceMethod = kTPCReferenceShowerStart;
        c.recoStartXcm = start.X();
        c.recoStartYcm = start.Y();
        c.recoStartZcm = start.Z();
        c.recoStartValid =
            IsFinite({c.recoStartXcm, c.recoStartYcm, c.recoStartZcm});
        c.tpcEntryXcm = start.X();
        c.tpcEntryYcm = start.Y();
        c.tpcEntryZcm = start.Z();
        c.tpcEntryPositionValid = c.recoStartValid;
        c.tpcEntryDirectionX = initialDirection.x;
        c.tpcEntryDirectionY = initialDirection.y;
        c.tpcEntryDirectionZ = initialDirection.z;
        c.tpcEntryDirectionValid = c.recoStartValid &&
                                   IsValidDirection(initialDirection);
        c.tpcEntryEndpointDirectionX = initialDirection.x;
        c.tpcEntryEndpointDirectionY = initialDirection.y;
        c.tpcEntryEndpointDirectionZ = initialDirection.z;
        c.tpcEntryEndpointDirectionValid = c.tpcEntryDirectionValid;
        mi.recoStartValid = c.tpcEntryPositionValid;
        mi.recoStart = {c.tpcEntryXcm, c.tpcEntryYcm, c.tpcEntryZcm};
        mi.recoDirectionValid = c.tpcEntryDirectionValid;
        if (c.tpcEntryDirectionValid) {
          mi.recoDirection = initialDirection;
        }
      }
      // Calculate all available match observables before evaluating selection.
      // Their validity and values are written even though they do not cut here.
      auto const m = BeamSelectionAlg::EvaluateInstrumentationMatch(mi);
      c.positionMatchValid = m.positionValid;
      c.directionMatchValid = m.directionValid;
      c.matchValid = m.valid;
      c.dx = m.deltaXcm;
      c.dy = m.deltaYcm;
      c.dz = m.deltaZcm;
      c.startz = m.entranceZcm;
      c.cos = m.directionCosine;
      if (c.recoCandidateType == kCandidateTrack &&
          c.tpcEntryDirectionValid) {
        BeamMatchInput localChordMatchInput = mi;
        localChordMatchInput.recoDirectionValid = true;
        localChordMatchInput.recoDirection = {
            c.tpcEntryDirectionX, c.tpcEntryDirectionY,
            c.tpcEntryDirectionZ};
        auto const localChordMatch =
            BeamSelectionAlg::EvaluateInstrumentationMatch(localChordMatchInput);
        c.localChordDirectionMatchValid = localChordMatch.directionValid;
        c.localChordDirectionCosine = localChordMatch.directionCosine;
      }
      c.triggerPass = isData_ ? (requireDataTrigger_
                                      ? (triggerEvaluated_ && goodTrigger_)
                                      : true)
                               : true;
      c.beamEventPass = beamEventUnique_;
      c.beamTrackPass = beamTrackUnique_;
      // Nominal gate equation for this primary:
      // trigger x unique beam track x Pandora beam-slice primary x exactly one
      // reco object. The separate unique-beam-event gate is event-level.
      // Position/direction residuals, object length, PID, and truth labels do
      // not enter is_selected.
      auto const nominalStageResult = BeamSelectionAlg::EvaluateCandidateStages(
          {c.triggerPass, c.beamTrackPass, c.beamSlicePrimary,
           c.hasUnambiguousRecoObject});
      bool const passesNominalCandidateStages = nominalStageResult.selected;
      bool const passesUniqueBeamEventGate = c.beamEventPass;
      c.passesAllSelectionGates =
          passesNominalCandidateStages && passesUniqueBeamEventGate;
      if (c.passesAllSelectionGates) {
        ++nPassing_;
        if (c.recoCandidateType == kCandidateTrack) {
          ++nPassingTrack_;
        } else if (c.recoCandidateType == kCandidateShower) {
          ++nPassingShower_;
        }
      }
      candidates.push_back(c);
    }
  }
  return candidates;
}

void PDHDBeamSelectionStages::SelectUniqueEventCandidate(
    std::vector<Candidate> &candidates) {
  // Phase 5: convert candidate-local pass flags into the event decision. A
  // single passing record is selected; zero records reject the event and two
  // or more records make it explicitly ambiguous.
  bool const hasUniquePassingCandidate = nPassing_ == 1;
  selectionAmbiguous_ = nPassing_ > 1;
  if (hasUniquePassingCandidate) {
    // Candidate::passesAllSelectionGates records the staged candidate outcome
    // and the separate unique-beam-event gate. is_selected is assigned only
    // after this event-level candidate-uniqueness requirement holds.
    auto const selectedCandidate = std::find_if(
        candidates.begin(), candidates.end(),
        [](Candidate const &candidate) {
          return candidate.passesAllSelectionGates;
        });
    if (selectedCandidate == candidates.end()) {
      throw cet::exception("PDHDBeamSelectionStages")
          << "Internal inconsistency: one passing candidate was counted but "
             "none was stored.";
    }
    selectedCandidate->selected = true;
    selectedCandidate->selectedTrack =
        selectedCandidate->recoCandidateType == kCandidateTrack;
    selectedCandidate->selectedShower =
        selectedCandidate->recoCandidateType == kCandidateShower;
    hasSelectedCandidate_ = true;
    selectedPfp_ = selectedCandidate->pfp;
    selectedTrack_ = selectedCandidate->track;
    selectedShower_ = selectedCandidate->shower;
    selectedCandidateType_ = selectedCandidate->recoCandidateType;
  }
  selectedEvents_ += hasSelectedCandidate_;
  ambiguousEvents_ += selectionAmbiguous_;
}

void PDHDBeamSelectionStages::AnnotateSelectedTrackTruth(
    art::Event const &evt, RecoPFPInputs const &inputs,
    std::vector<TrackTruthMatch> const &trackTruthMatches,
    std::vector<Candidate> &candidates) {
  // This post-selection diagnostic cannot modify is_selected or recover a
  // rejected event.
  if (isData_ || !enableMCTrackTruthMatching_ || !truthBeamGeantMatchUnique_) {
    return;
  }

  auto selected = std::find_if(
      candidates.begin(), candidates.end(), [](Candidate const &candidate) {
        return candidate.selectedTrack;
      });
  if (selected == candidates.end() || selected->track < 0) return;

  try {
    std::size_t const trackIndex = static_cast<std::size_t>(selected->track);
    if (!inputs.tracks || trackIndex >= inputs.tracks->size() ||
        trackIndex >= trackTruthMatches.size() ||
        !trackTruthMatches[trackIndex].evaluated) {
      return;
    }
    recob::Track const &track = inputs.tracks->at(trackIndex);
    TrackTruthMatch const &truthMatch = trackTruthMatches[trackIndex];

    auto const *utilityMatch = truthMatch.utilityParticle;
    if (utilityMatch) {
      selected->recoTruthUtilityMatchValid = true;
      selected->recoTruthUtilityTrackId = utilityMatch->TrackId();
      selected->recoTruthUtilityPdg = utilityMatch->PdgCode();
      selected->recoTruthUtilityMatchesBeamGeant =
          utilityMatch->TrackId() == truthBeamGeantTrackId_;
      selected->recoTruthUtilityProcess = utilityMatch->Process();
      selected->recoTruthUtilityOriginValid = truthMatch.originValid;
      selected->recoTruthUtilityOrigin = truthMatch.origin;
      selected->recoTruthBackgroundCategory = ClassifySelectedTrackTruth(
          selected->recoTruthUtilityMatchesBeamGeant,
          selected->recoTruthUtilityOriginValid,
          selected->recoTruthUtilityOrigin, mcTruthBeamOrigin_,
          mcCosmicOrigin_, selected->recoTruthUtilityProcess);
      selected->recoTruthBackgroundCategoryValid =
          selected->recoTruthBackgroundCategory != kTruthCategoryUnavailable;
      selected->recoTruthUtilityTrajectoryValid = StoreGeantTrajectoryInZRange(
          *utilityMatch, std::max(0., mcTruthReferenceMinimumZcm_),
          mcTruthReferenceMaximumZcm_,
          selected->recoTruthUtilityTrajectoryXcm,
          selected->recoTruthUtilityTrajectoryYcm,
          selected->recoTruthUtilityTrajectoryZcm);
    }

    // Standard shared-hit match preserves the established hit-match branches.
    if (truthMatch.hitParticle) {
      selected->recoTruthHitMatchValid = true;
      selected->recoTruthHitTrackId = truthMatch.hitParticle->TrackId();
      selected->recoTruthHitPdg = truthMatch.hitParticle->PdgCode();
      selected->recoTruthHitSharedCount = static_cast<unsigned int>(
          std::min(truthMatch.sharedHitCount,
                   static_cast<std::size_t>(
                       std::numeric_limits<unsigned int>::max())));
      selected->recoTruthHitSharedDeltaRayCount = static_cast<unsigned int>(
          std::min(truthMatch.sharedDeltaRayHitCount,
                   static_cast<std::size_t>(
                       std::numeric_limits<unsigned int>::max())));
      selected->recoTruthHitMatchesBeamGeant =
          truthMatch.hitParticle->TrackId() == truthBeamGeantTrackId_;
    }

    auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    protoana::ProtoDUNETruthUtils truthUtils;
    selected->recoTruthPurity =
        truthUtils.GetPurity(clockData, track, evt, trackTag_.encode());
    selected->recoTruthPurityValid = std::isfinite(selected->recoTruthPurity);
    selected->recoTruthCompleteness = truthUtils.GetCompleteness(
        clockData, track, evt, trackTag_.encode(), truthHitTag_.encode());
    selected->recoTruthCompletenessValid =
        std::isfinite(selected->recoTruthCompleteness);
    selected->recoTruthAssociationValid = true;
  } catch (std::exception const &error) {
    // Preserve the already-written selection decision and make the failed
    // optional diagnostic visible through its default invalid branches.
    if (eventMessagesEmitted_ < maxEventMessages_) {
      mf::LogWarning("PDHDBeamSelectionStages")
          << "Could not label selected MC track for run " << run_
          << ", subrun " << subrun_ << ", event " << event_ << ": "
          << error.what();
      ++eventMessagesEmitted_;
    }
  }
}

std::vector<PDHDBeamSelectionStages::TrackTruthMatch>
PDHDBeamSelectionStages::BuildRecoTrackInventory(
    art::Event const &evt, RecoPFPInputs const &inputs) {
  std::vector<TrackTruthMatch> truthMatches;
  recoTrackTrajectoryOffsets_.push_back(0ULL);
  if (!inputs.tracks) return truthMatches;
  truthMatches.resize(inputs.tracks->size());

  std::set<int> beamSliceTrackIndices;
  if (pandoraBeamSliceFound_ && inputs.pfps && inputs.slices &&
      inputs.trackAssociations) {
    for (std::size_t pfp = 0; pfp < inputs.pfps->size(); ++pfp) {
      try {
        auto const slice = inputs.slices->at(pfp);
        if (!slice || static_cast<int>(slice->ID()) != pandoraBeamSliceId_)
          continue;
        for (auto const &track : inputs.trackAssociations->at(pfp)) {
          int const index = IndexOf(*inputs.tracks, track.get());
          if (index >= 0) beamSliceTrackIndices.insert(index);
        }
      } catch (std::exception const &) {
      }
    }
  }
  std::unique_ptr<detinfo::DetectorClocksData> clockData;
  if (!isData_ && enableMCTrackTruthMatching_) {
    try {
      clockData = std::make_unique<detinfo::DetectorClocksData>(
          art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt));
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_ < maxEventMessages_) {
        mf::LogWarning("PDHDBeamSelectionStages")
            << "Could not obtain detector clocks for all-track truth labels "
            << "in run " << run_ << ", subrun " << subrun_ << ", event "
            << event_ << ": " << error.what();
        ++eventMessagesEmitted_;
      }
    }
  }
  protoana::ProtoDUNETruthUtils truthUtils;
  for (std::size_t index = 0; index < inputs.tracks->size(); ++index) {
    recoTrackIndices_.push_back(static_cast<int>(index));
    recoTrackInBeamSlice_.push_back(
        beamSliceTrackIndices.count(static_cast<int>(index)) ? 1 : 0);
    recoTrackCosmicClassificationValid_.push_back(0);
    recoTrackIsCosmic_.push_back(0);
    recoTrackTruthUtilityMatchValid_.push_back(0);
    recoTrackTruthUtilityPdg_.push_back(0);
    recoTrackTruthHitMatchValid_.push_back(0);
    recoTrackTruthHitPdg_.push_back(0);

    auto const trajectory =
        ExtractTrackTrajectoryFromEntry(inputs.tracks->at(index));
    for (auto const &point : trajectory) {
      recoTrackTrajectoryXcm_.push_back(point.x);
      recoTrackTrajectoryYcm_.push_back(point.y);
      recoTrackTrajectoryZcm_.push_back(point.z);
    }
    recoTrackTrajectoryOffsets_.push_back(
        static_cast<unsigned long long>(recoTrackTrajectoryXcm_.size()));

    if (clockData) {
      try {
        TrackTruthMatch &truthMatch = truthMatches[index];
        truthMatch.evaluated = true;
        auto const *match = truthUtils.GetMCParticleFromRecoTrack(
            *clockData, inputs.tracks->at(index), evt, trackTag_.encode());
        truthMatch.utilityParticle = match;
        if (match) {
          recoTrackTruthUtilityMatchValid_.back() = 1;
          recoTrackTruthUtilityPdg_.back() = match->PdgCode();
          art::ServiceHandle<cheat::ParticleInventoryService> inventory;
          auto const truth = inventory->TrackIdToMCTruth_P(match->TrackId());
          if (truth) {
            truthMatch.originValid = true;
            truthMatch.origin = static_cast<int>(truth->Origin());
            recoTrackCosmicClassificationValid_.back() = 1;
            recoTrackIsCosmic_.back() =
                truthMatch.origin == mcCosmicOrigin_ ? 1 : 0;
          }
        }
        auto const hitMatch = truthUtils.GetMCParticleByHits(
            *clockData, inputs.tracks->at(index), evt, trackTag_.encode(),
            truthHitTag_.encode());
        truthMatch.hitParticle = hitMatch.particle;
        truthMatch.sharedHitCount = hitMatch.nSharedHits;
        truthMatch.sharedDeltaRayHitCount = hitMatch.nSharedDeltaRayHits;
        if (hitMatch.particle) {
          recoTrackTruthHitMatchValid_.back() = 1;
          recoTrackTruthHitPdg_.back() = hitMatch.particle->PdgCode();
        }
      } catch (std::exception const &error) {
        truthMatches[index] = TrackTruthMatch{};
        recoTrackCosmicClassificationValid_.back() = 0;
        recoTrackIsCosmic_.back() = 0;
        recoTrackTruthUtilityMatchValid_.back() = 0;
        recoTrackTruthUtilityPdg_.back() = 0;
        recoTrackTruthHitMatchValid_.back() = 0;
        recoTrackTruthHitPdg_.back() = 0;
        if (eventMessagesEmitted_ < maxEventMessages_) {
          mf::LogWarning("PDHDBeamSelectionStages")
              << "Could not label reco track " << index << " in run "
              << run_ << ", subrun " << subrun_ << ", event " << event_
              << ": " << error.what();
          ++eventMessagesEmitted_;
        }
      }
    }
  }
  return truthMatches;
}

void PDHDBeamSelectionStages::LogSelectionResult() {
  if (eventMessagesEmitted_ < maxEventMessages_) {
    mf::LogInfo("PDHDBeamSelectionStages")
        << "Selection result for run " << run_ << ", subrun " << subrun_
        << ", event " << event_ << ": reference='" << referenceSource_
        << "', beam event unique=" << beamEventUnique_
        << ", beam track unique=" << beamTrackUnique_
        << ", Pandora beam slice found=" << pandoraBeamSliceFound_
        << " (id=" << pandoraBeamSliceId_
        << ", candidates=" << pandoraBeamSliceCount_ << ")"
        << ", primary PFPs=" << nPrimary_
        << ", direct IsTestBeam primaries=" << nBeamPrimary_
        << ", beam-slice primaries=" << nBeamSlicePrimary_
        << ", passing candidates=" << nPassing_ << " (track="
        << nPassingTrack_ << ", shower=" << nPassingShower_ << ")"
        << ", selected PFP index=" << selectedPfp_
        << ", selected type=" << selectedCandidateType_
        << ", ambiguous=" << selectionAmbiguous_ << ".";
    ++eventMessagesEmitted_;
  }
}

void PDHDBeamSelectionStages::FillOutputTrees(
    std::vector<Candidate> const &candidates) {
  // Phase 6: write every primary-PFP diagnostic row, including rejected and
  // ambiguous candidates, then write one event summary row with multiplicities.
  for (auto const &c : candidates) {
    out_ = c;
    candidateTree_->Fill();
  }
  eventTree_->Fill();
}

void PDHDBeamSelectionStages::analyze(art::Event const &evt) {
  ResetEventState(evt);

  BeamInstrumentationRecord dataBeamReference;
  AcquireBeamReference(evt, dataBeamReference);

  auto recoInputs = AcquireRecoPFPInputs(evt);
  FindPandoraBeamSlice(recoInputs);
  auto candidates = BuildPrimaryCandidates(recoInputs, dataBeamReference);

  SelectUniqueEventCandidate(candidates);
  auto const trackTruthMatches = BuildRecoTrackInventory(evt, recoInputs);
  AnnotateSelectedTrackTruth(evt, recoInputs, trackTruthMatches, candidates);
  LogSelectionResult();
  FillOutputTrees(candidates);
}

void PDHDBeamSelectionStages::endJob() {
  std::ostringstream summary;
  summary << "Beam-selection summary: processed=" << eventsProcessed_
          << " (data=" << dataEvents_ << ", MC=" << mcEvents_
          << "), beam product available=" << beamProductAvailableEvents_
          << ", exactly one beam event=" << uniqueBeamEventEvents_
          << ", exactly one beam track=" << uniqueBeamTrackEvents_
          << ", PFParticle product available=" << pfpAvailableEvents_
          << ", selected events=" << selectedEvents_
          << ", ambiguous selected sets=" << ambiguousEvents_ << ".";
  mf::LogInfo("PDHDBeamSelectionStages") << summary.str();
  // Match the proven calorimetry-module behavior for visible final counters.
  if (printSummaryToStdout_) {
    std::cout << "PDHDBeamSelectionStages: " << summary.str() << std::endl;
  }
}
} // namespace pdhd::diagnostics
DEFINE_ART_MODULE(pdhd::diagnostics::PDHDBeamSelectionStages)
