/**
 * @file PDHDBeamSelectionStages_module.cc
 * @brief Applies one staged, data-like beam-candidate definition to data and
 * MC.
 *
 * Data reads the reconstructed external beamline from DataBeamTag. MC reads
 * the generated beam primary from MCTruthTag and its Geant trajectory from
 * MCParticleTag. The nominal decision requires the corresponding unique beam
 * reference, a Pandora beam-slice primary, and one unambiguous reconstructed
 * track or shower.
 * Beamline-to-TPC position and direction compatibility are stored as
 * observables, not cut by an inherited SP table. Tracks use a local TPC-entry
 * segment; showers use their reconstructed start and initial direction. MC
 * truth labels are diagnostic and do not select a particle species.
 * CandidateTree stores each selection component and match
 * observable; EventTree stores multiplicities. A candidate is selected only
 * when exactly one passes the named reconstruction stages.
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
#include "cetlib_except/exception.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamInstrumentationAlg.h"
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamSelectionAlg.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
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
constexpr int kTPCReferenceTrackVertexLocalDirection = 1;
constexpr int kTPCReferenceShowerStart = 2;

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

// The lower-Z endpoint is the TPC entry convention for this beam geometry.
// Use only the first requested path length after that entry; this deliberately
// avoids both Track::StartDirection() and the full start-to-end chord.
TPCEntryDirection FindTPCEntryDirection(recob::Track const &track,
                                        double const requestedLengthCm) {
  TPCEntryDirection result;
  if (!std::isfinite(requestedLengthCm) || requestedLengthCm <= 0.) {
    return result;
  }

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
  // A short track does not silently use a different direction definition.
  return TPCEntryDirection{};
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

  simb::MCParticle const *generated = nullptr;
  for (auto const &truth : *truthHandle) {
    if (static_cast<int>(truth.Origin()) != requiredOrigin) continue;
    for (int index = 0; index < truth.NParticles(); ++index) {
      auto const &particle = truth.GetParticle(index);
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
        beamInstrumentationNominalMomentumGeV_(
            p.get<double>("BeamInstrumentationNominalMomentumGeV")),
        evaluateDataPid_(p.get<bool>("EvaluateDataPID", false)),
        requireDataTrigger_(p.get<bool>("RequireGoodDataTrigger", true)),
        mcTruthBeamOrigin_(p.get<int>("MCTruthBeamOrigin", 4)),
        mcTruthGeantEnergyToleranceGeV_(
            p.get<double>("MCTruthGeantEnergyToleranceGeV", 1.e-5)),
        mcTruthReferenceMinimumZcm_(p.get<double>("MCTruthReferenceMinimumZCm", 0.)),
        mcTruthReferenceMaximumZcm_(p.get<double>("MCTruthReferenceMaximumZCm", 600.)),
        tpcEntryDirectionLengthCm_(
            p.get<double>("TPCEntryDirectionLengthCm", 5.)),
        printSummaryToStdout_(p.get<bool>("PrintSummaryToStdout", true)),
        maxEventMessages_(p.get<unsigned int>("MaxEventMessages", 3U)),
        beamline_(p.get<fhicl::ParameterSet>("BeamlineUtils")) {
    consumes<std::vector<beam::ProtoDUNEBeamEvent>>(dataBeamTag_);
    consumes<std::vector<simb::MCTruth>>(mcTruthTag_);
    consumes<std::vector<simb::MCParticle>>(mcParticleTag_);
    consumes<std::vector<recob::PFParticle>>(pfpTag_);
    consumes<std::vector<recob::Track>>(trackTag_);
    consumes<std::vector<recob::Shower>>(showerTag_);
  }
  void beginJob() override;
  void analyze(art::Event const &) override;
  void endJob() override;

private:
  struct Candidate {
    int pfp = -1;
    int track = -1;
    int shower = -1;
    int recoCandidateType = kCandidateNone;
    int tpcReferenceMethod = kTPCReferenceUnavailable;
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
    bool positionMatchValid = false;
    bool directionMatchValid = false;
    bool matchValid = false;
    bool triggerPass = false;
    bool beamEventPass = false;
    bool beamTrackPass = false;
    bool recoStartValid = false;
    bool tpcEntryDirectionValid = false;
    // `passesAllSelectionGates` is candidate-local; `selected` additionally
    // means it is the event's unique passing candidate.
    bool passesAllSelectionGates = false;
    bool selected = false;
    bool selectedTrack = false;
    bool selectedShower = false;
    unsigned int trackAssociationCount = 0;
    unsigned int showerAssociationCount = 0;
    // Position residuals remain unavailable only when a declared reference or
    // reconstructed start is absent; cosine additionally needs both directions.
    double dx = std::numeric_limits<double>::quiet_NaN();
    double dy = std::numeric_limits<double>::quiet_NaN();
    double dz = std::numeric_limits<double>::quiet_NaN();
    double startz = std::numeric_limits<double>::quiet_NaN();
    double cos = std::numeric_limits<double>::quiet_NaN();
    // Track/Shower::Length is the native reconstructed-object length in cm.
    double recoObjectLengthCm = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndXcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndYcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndZcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndDirectionX = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndDirectionY = std::numeric_limits<double>::quiet_NaN();
    double beamlineEndDirectionZ = std::numeric_limits<double>::quiet_NaN();
    double beamlineStartXcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineStartYcm = std::numeric_limits<double>::quiet_NaN();
    double beamlineStartZcm = std::numeric_limits<double>::quiet_NaN();
    // This is Track::Vertex() for tracks and ShowerStart() for showers.  It
    // is the position used for the SP-compatible residual definitions.
    double recoStartXcm = std::numeric_limits<double>::quiet_NaN();
    double recoStartYcm = std::numeric_limits<double>::quiet_NaN();
    double recoStartZcm = std::numeric_limits<double>::quiet_NaN();
    // The lower-Z trajectory point remains separate: it anchors the local
    // track direction and is not silently substituted as the residual start.
    double tpcEntryXcm = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryYcm = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryZcm = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryDirectionX = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryDirectionY = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryDirectionZ = std::numeric_limits<double>::quiet_NaN();
    double tpcEntryDirectionSampledLengthCm =
        std::numeric_limits<double>::quiet_NaN();
    std::vector<double> trackTrajectoryXcm;
    std::vector<double> trackTrajectoryYcm;
    std::vector<double> trackTrajectoryZcm;
  };
  art::InputTag dataBeamTag_, mcTruthTag_, mcParticleTag_, pfpTag_, pfpMetadataTag_,
      pfpSliceTag_, trackTag_, showerTag_;
  double beamInstrumentationNominalMomentumGeV_;
  bool evaluateDataPid_;
  bool requireDataTrigger_;
  int mcTruthBeamOrigin_;
  double mcTruthGeantEnergyToleranceGeV_;
  double mcTruthReferenceMinimumZcm_, mcTruthReferenceMaximumZcm_;
  double tpcEntryDirectionLengthCm_;
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
  AddBranch(*candidateTree_, "pfp_index", out_.pfp);
  AddBranch(*candidateTree_, "track_index", out_.track);
  AddBranch(*candidateTree_, "shower_index", out_.shower);
  // 0=none, 1=track, 2=shower, 3=ambiguous associations.
  AddBranch(*candidateTree_, "reco_candidate_type", out_.recoCandidateType);
  // 0=unavailable, 1=Track::Vertex with a local entry direction,
  // 2=ShowerStart/Direction.
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
  // These explicit coordinates make the beamline-end to reconstructed-start
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
  // Track local entry is retained for its first-segment direction diagnostic.
  // For showers it coincides with ShowerStart().
  AddBranch(*candidateTree_, "tpc_entry_x_cm", out_.tpcEntryXcm);
  AddBranch(*candidateTree_, "tpc_entry_y_cm", out_.tpcEntryYcm);
  AddBranch(*candidateTree_, "tpc_entry_z_cm", out_.tpcEntryZcm);
  // The cosine uses the beamline fitted endpoint direction and the declared
  // local track chord or reconstructed initial shower direction.
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
  // Full valid recob::Track trajectory, ordered from lower-Z TPC entry.
  AddBranch(*candidateTree_, "track_trajectory_x_cm", out_.trackTrajectoryXcm);
  AddBranch(*candidateTree_, "track_trajectory_y_cm", out_.trackTrajectoryYcm);
  AddBranch(*candidateTree_, "track_trajectory_z_cm", out_.trackTrajectoryZcm);
  AddBranch(*candidateTree_, "delta_x_cm", out_.dx);
  AddBranch(*candidateTree_, "delta_y_cm", out_.dy);
  // Difference of the two measured/reconstructed endpoints, not absolute Z.
  AddBranch(*candidateTree_, "delta_z_cm", out_.dz);
  AddBranch(*candidateTree_, "entrance_z_cm", out_.startz);
  AddBranch(*candidateTree_, "direction_cosine", out_.cos);
  AddBranch(*candidateTree_, "is_selected", out_.selected);
  AddBranch(*candidateTree_, "is_selected_track", out_.selectedTrack);
  AddBranch(*candidateTree_, "is_selected_shower", out_.selectedShower);
  AddBranch(*candidateTree_, "beam_reference_source", referenceSource_);
}
void PDHDBeamSelectionStages::analyze(art::Event const &evt) {
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
  BeamInstrumentationRecord bi;
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
    auto const truthHandle = evt.getHandle<std::vector<simb::MCTruth>>(mcTruthTag_);
    auto const particleHandle =
        evt.getHandle<std::vector<simb::MCParticle>>(mcParticleTag_);
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
  auto ph = evt.getHandle<std::vector<recob::PFParticle>>(pfpTag_);
  auto th = evt.getHandle<std::vector<recob::Track>>(trackTag_);
  auto sh = evt.getHandle<std::vector<recob::Shower>>(showerTag_);
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
  std::unique_ptr<art::FindManyP<larpandoraobj::PFParticleMetadata>> metadata;
  std::unique_ptr<art::FindOneP<recob::Slice>> slices;
  std::unique_ptr<art::FindManyP<recob::Track>> tracks;
  std::unique_ptr<art::FindManyP<recob::Shower>> showers;
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

  // Pandora's IsTestBeam metadata identifies the beam slice. The nominal TPC
  // candidate is every primary PFP associated with that selected slice.
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

  std::vector<Candidate> candidates;
  if (ph) {
    for (std::size_t i = 0; i < ph->size(); ++i) {
      auto const &p = ph->at(i);
      if (!p.IsPrimary()) {
        continue;
      }
      Candidate c;
      c.pfp = i;
      c.primary = true;
      ++nPrimary_;
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
      c.recoObjectLengthValid =
          c.hasUnambiguousRecoObject && std::isfinite(c.recoObjectLengthCm) &&
          c.recoObjectLengthCm >= 0.;
      BeamMatchInput mi;
      mi.beamReferenceValid = beamTrackUnique_;
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
      if (c.recoCandidateType == kCandidateTrack) {
        auto const trajectory = ExtractTrackTrajectoryFromEntry(*tr);
        for (auto const &point : trajectory) {
          c.trackTrajectoryXcm.push_back(point.x);
          c.trackTrajectoryYcm.push_back(point.y);
          c.trackTrajectoryZcm.push_back(point.z);
        }
        auto const localDirection =
            FindTPCEntryDirection(*tr, tpcEntryDirectionLengthCm_);
        c.tpcReferenceMethod = kTPCReferenceTrackVertexLocalDirection;
        auto const vertex = tr->Vertex();
        c.recoStartXcm = vertex.X();
        c.recoStartYcm = vertex.Y();
        c.recoStartZcm = vertex.Z();
        c.recoStartValid =
            IsFinite({c.recoStartXcm, c.recoStartYcm, c.recoStartZcm});
        c.tpcEntryXcm = localDirection.entry.x;
        c.tpcEntryYcm = localDirection.entry.y;
        c.tpcEntryZcm = localDirection.entry.z;
        c.tpcEntryDirectionValid = localDirection.valid;
        c.tpcEntryDirectionSampledLengthCm =
            localDirection.sampledLengthCm;
        mi.recoStartValid = c.recoStartValid;
        mi.recoStart = {c.recoStartXcm, c.recoStartYcm, c.recoStartZcm};
        mi.recoDirectionValid = localDirection.valid;
        // A short track has no substituted direction, but its Vertex residual
        // remains available when the external beamline reference is valid.
        if (localDirection.valid) {
          c.tpcEntryDirectionX = localDirection.direction.x;
          c.tpcEntryDirectionY = localDirection.direction.y;
          c.tpcEntryDirectionZ = localDirection.direction.z;
          // SP TrackStart* residuals use Track::Vertex(), while the direction
          // remains a short local chord from the physical lower-Z entry.
          mi.recoDirection = localDirection.direction;
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
        c.tpcEntryDirectionX = initialDirection.x;
        c.tpcEntryDirectionY = initialDirection.y;
        c.tpcEntryDirectionZ = initialDirection.z;
        c.tpcEntryDirectionValid = c.recoStartValid &&
                                   IsValidDirection(initialDirection);
        mi.recoStartValid = c.recoStartValid;
        mi.recoStart = {c.recoStartXcm, c.recoStartYcm, c.recoStartZcm};
        mi.recoDirectionValid = c.tpcEntryDirectionValid;
        if (c.tpcEntryDirectionValid) {
          mi.recoDirection = initialDirection;
        }
      }
      auto const m = BeamSelectionAlg::EvaluateInstrumentationMatch(mi);
      c.positionMatchValid = m.positionValid;
      c.directionMatchValid = m.directionValid;
      c.matchValid = m.valid;
      c.dx = m.deltaXcm;
      c.dy = m.deltaYcm;
      c.dz = m.deltaZcm;
      c.startz = m.entranceZcm;
      c.cos = m.directionCosine;
      c.triggerPass = isData_ ? (requireDataTrigger_
                                      ? (triggerEvaluated_ && goodTrigger_)
                                      : true)
                               : true;
      c.beamEventPass = beamEventUnique_;
      c.beamTrackPass = beamTrackUnique_;
      // These are the four per-candidate nominal stages.  The beam-event
      // uniqueness gate is event-level and remains intentionally outside the
      // staged candidate algorithm.
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
  for (auto const &c : candidates) {
    out_ = c;
    candidateTree_->Fill();
  }
  eventTree_->Fill();
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
