/**
 * @file PDHDBeamInstrumentation_module.cc
 * @brief Writes raw ProtoDUNE-HD beam-instrumentation observables for data and
 * MC.
 *
 * Data uses DataBeamTag (normally beamevent); MC uses MCBeamTag (normally
 * generator). One event-tree row records product cardinality, trigger,
 * momentum, TOF, Cherenkov, PID, monitor occupancy, and every beamline track.
 * Mode-specific FHiCL records whether the product was stored or rebuilt and
 * whether LLT/Cherenkov status is valid. Missing or non-unique products remain
 * explicit; unavailable floating-point quantities use NaN, not physical zero.
 */
#include "TTree.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamInstrumentationAlg.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace pdhd::diagnostics {
namespace {
// ROOT branches must point to module-owned objects that outlive every Fill().
template <typename T> void addBranch(TTree &tree, char const *name, T &value) {
  tree.Branch(name, std::addressof(value));
}
} // namespace

class PDHDBeamInstrumentation : public art::EDAnalyzer {
public:
  explicit PDHDBeamInstrumentation(fhicl::ParameterSet const &p)
      : EDAnalyzer(p), dataTag_(p.get<art::InputTag>(
                           "DataBeamTag", art::InputTag{"beamevent"})),
        mcTag_(p.get<art::InputTag>("MCBeamTag", art::InputTag{"generator"})),
        mcTruthTag_(
            p.get<art::InputTag>("MCTruthTag", art::InputTag{"generator"})),
        nominalMomentum_(p.get<double>("NominalMomentumGeV", 1.)),
        evaluateDataPid_(p.get<bool>("EvaluateDataPID", false)),
        evaluateMcTrigger_(p.get<bool>("EvaluateMCTrigger", false)),
        evaluateMcPid_(p.get<bool>("EvaluateMCPID", false)),
        beamEventSourceDescription_(
            p.get<std::string>("BeamEventSourceDescription")),
        beamEventBuildMethod_(p.get<std::string>("BeamEventBuildMethod")),
        beamEventProducedInCurrentJob_(
            p.get<bool>("BeamEventProducedInCurrentJob")),
        cherenkovStatusSource_(p.get<std::string>("CherenkovStatusSource")),
        cherenkovStatusProvenanceValid_(
            p.get<bool>("CherenkovStatusProvenanceValid")),
        cherenkovPressureSource_(
            p.get<std::string>("CherenkovPressureSource")),
        cherenkovPressureProvenanceValid_(
            p.get<bool>("CherenkovPressureProvenanceValid")),
        nonNominalInputHandlingUsed_(
            p.get<bool>("NonNominalInputHandlingUsed")),
        maxEventMessages_(p.get<unsigned int>("MaxEventMessages", 3U)),
        monitors_(p.get<std::vector<std::string>>(
            "FiberMonitorNames", {"XBPF022697", "XBPF022701", "XBPF022702"})),
        beamline_(p.get<fhicl::ParameterSet>("BeamlineUtils")) {
    consumes<std::vector<beam::ProtoDUNEBeamEvent>>(dataTag_);
    consumes<std::vector<beam::ProtoDUNEBeamEvent>>(mcTag_);
    consumes<std::vector<simb::MCTruth>>(mcTruthTag_);
  }
  void beginJob() override;
  void analyze(art::Event const &) override;
  void endJob() override;

private:
  // Configuration records both the selected product and how it was produced.
  art::InputTag dataTag_, mcTag_, mcTruthTag_;
  double nominalMomentum_;
  bool evaluateDataPid_, evaluateMcTrigger_, evaluateMcPid_;
  std::string beamEventSourceDescription_, beamEventBuildMethod_;
  bool beamEventProducedInCurrentJob_;
  std::string cherenkovStatusSource_;
  bool cherenkovStatusProvenanceValid_;
  std::string cherenkovPressureSource_;
  bool cherenkovPressureProvenanceValid_;
  bool nonNominalInputHandlingUsed_;
  unsigned int maxEventMessages_;
  std::vector<std::string> monitors_;
  protoana::ProtoDUNEBeamlineUtils beamline_;

  // Job counters summarize missing products and evaluated diagnostics.
  std::uint64_t eventsProcessed_ = 0;
  std::uint64_t dataEvents_ = 0;
  std::uint64_t mcEvents_ = 0;
  std::uint64_t productAvailableEvents_ = 0;
  std::uint64_t uniqueBeamEventEvents_ = 0;
  std::uint64_t triggerEvaluatedEvents_ = 0;
  std::uint64_t pidEvaluatedEvents_ = 0;
  std::uint64_t goodTriggerEvents_ = 0;
  std::uint64_t extractionErrorEvents_ = 0;
  unsigned int eventMessagesEmitted_ = 0;

  // Stable per-event storage used directly by the output tree branches.
  TTree *tree_ = nullptr;
  unsigned int run_ = 0, subrun_ = 0, event_ = 0;
  bool isData_ = false, available_ = false, unique_ = false;
  bool pidEvaluated_ = false;
  bool mcBeamPdgValid_ = false;
  bool cherenkovStatusValueAvailable_ = false;
  bool cherenkovPressureValueAvailable_ = false;
  unsigned int beamEventCount_ = 0;
  unsigned int mcTruthRecordCount_ = 0, mcBeamPrimaryCount_ = 0;
  int mcBeamPdg_ = -999;
  std::string source_, selectedProductTag_, pidMethod_, mcBeamPdgSource_;
  bool goodTrigger_ = false, triggerEvaluated_ = false;
  bool hasPerfectBeamMomentum_ = false, triggersMatched_ = false;
  int timingTrigger_ = -1, beamTrigger_ = -1;
  int ckov0_ = -1, ckov1_ = -1;
  double ckov0Pressure_ = std::numeric_limits<double>::quiet_NaN();
  double ckov1Pressure_ = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> momenta_, tof_, trackStartX_, trackStartY_, trackStartZ_,
      trackEndX_, trackEndY_, trackEndZ_;
  std::vector<double> trackStartDirX_, trackStartDirY_, trackStartDirZ_,
      trackEndDirX_, trackEndDirY_, trackEndDirZ_;
  std::vector<int> tofChannels_, pidCandidates_, activeFiberCounts_;
  std::vector<std::string> monitorNames_;
  void reset();
  void fillMCBeamPdg(art::Event const &);
};
void PDHDBeamInstrumentation::beginJob() {
  mf::LogInfo("PDHDBeamInstrumentation")
      << "Starting beam-instrumentation diagnostics. DataBeamTag='"
      << dataTag_.encode() << "', MCBeamTag='" << mcTag_.encode()
      << "', MCTruthTag='" << mcTruthTag_.encode()
      << "', NominalMomentumGeV=" << nominalMomentum_
      << ", EvaluateDataPID=" << evaluateDataPid_
      << ", EvaluateMCTrigger=" << evaluateMcTrigger_
      << ", EvaluateMCPID=" << evaluateMcPid_
      << ", BeamEventBuildMethod='" << beamEventBuildMethod_
      << "', BeamEventProducedInCurrentJob="
      << beamEventProducedInCurrentJob_ << ", CherenkovStatusSource='"
      << cherenkovStatusSource_ << "', CherenkovStatusProvenanceValid="
      << cherenkovStatusProvenanceValid_ << ", CherenkovPressureSource='"
      << cherenkovPressureSource_ << "', CherenkovPressureProvenanceValid="
      << cherenkovPressureProvenanceValid_
      << ", NonNominalInputHandlingUsed=" << nonNominalInputHandlingUsed_
      << ", MaxEventMessages=" << maxEventMessages_
      << ". Data jobs must schedule a producer for DataBeamTag when the input "
         "file does not retain that product.";

  tree_ = art::ServiceHandle<art::TFileService>()->make<TTree>(
      "BeamInstrumentation", "Raw beam instrumentation");

  // Event identity and product provenance remain available even on failures.
  addBranch(*tree_, "run", run_);
  addBranch(*tree_, "subrun", subrun_);
  addBranch(*tree_, "event", event_);
  addBranch(*tree_, "is_data", isData_);
  addBranch(*tree_, "product_available", available_);
  addBranch(*tree_, "beam_event_unique", unique_);
  addBranch(*tree_, "beam_event_count", beamEventCount_);
  addBranch(*tree_, "source", source_);
  addBranch(*tree_, "beam_event_product_tag", selectedProductTag_);
  addBranch(*tree_, "beam_event_source_description",
            beamEventSourceDescription_);
  addBranch(*tree_, "beam_event_build_method", beamEventBuildMethod_);
  addBranch(*tree_, "beam_event_produced_in_current_job",
            beamEventProducedInCurrentJob_);
  addBranch(*tree_, "non_nominal_input_handling_used",
            nonNominalInputHandlingUsed_);

  // MC truth identity is independent of simulated beam-instrumentation fields.
  // Data keeps the sentinel and mc_beam_pdg_valid=false.
  addBranch(*tree_, "mc_beam_pdg", mcBeamPdg_);
  addBranch(*tree_, "mc_beam_pdg_valid", mcBeamPdgValid_);
  addBranch(*tree_, "mc_beam_pdg_source", mcBeamPdgSource_);
  addBranch(*tree_, "mc_truth_record_count", mcTruthRecordCount_);
  addBranch(*tree_, "mc_beam_primary_count", mcBeamPrimaryCount_);

  // Selection decisions carry explicit evaluation and method information.
  addBranch(*tree_, "IsGoodBeamlineTrigger", goodTrigger_);
  addBranch(*tree_, "trigger_evaluated", triggerEvaluated_);
  addBranch(*tree_, "pid_evaluated", pidEvaluated_);
  addBranch(*tree_, "pid_method", pidMethod_);
  addBranch(*tree_, "ckov_status_value_available",
            cherenkovStatusValueAvailable_);
  addBranch(*tree_, "ckov_status_provenance_valid",
            cherenkovStatusProvenanceValid_);
  addBranch(*tree_, "ckov_status_source", cherenkovStatusSource_);
  addBranch(*tree_, "ckov_pressure_value_available",
            cherenkovPressureValueAvailable_);
  addBranch(*tree_, "ckov_pressure_provenance_valid",
            cherenkovPressureProvenanceValid_);
  addBranch(*tree_, "ckov_pressure_source", cherenkovPressureSource_);

  // These are uncut beam-event observables; units are encoded in their names.
  addBranch(*tree_, "HasPerfectBeamMomentum", hasPerfectBeamMomentum_);
  addBranch(*tree_, "GetTimingTrigger", timingTrigger_);
  addBranch(*tree_, "GetBITrigger", beamTrigger_);
  addBranch(*tree_, "CheckIsMatched", triggersMatched_);
  addBranch(*tree_, "momenta_GeV", momenta_);
  addBranch(*tree_, "tof_ns", tof_);
  addBranch(*tree_, "tof_channels", tofChannels_);
  addBranch(*tree_, "pid_candidates", pidCandidates_);
  addBranch(*tree_, "ckov0_status", ckov0_);
  addBranch(*tree_, "ckov1_status", ckov1_);
  addBranch(*tree_, "ckov0_pressure", ckov0Pressure_);
  addBranch(*tree_, "ckov1_pressure", ckov1Pressure_);
  addBranch(*tree_, "monitor_names", monitorNames_);
  addBranch(*tree_, "active_fiber_counts", activeFiberCounts_);
  addBranch(*tree_, "track_start_x_cm", trackStartX_);
  addBranch(*tree_, "track_start_y_cm", trackStartY_);
  addBranch(*tree_, "track_start_z_cm", trackStartZ_);
  addBranch(*tree_, "track_end_x_cm", trackEndX_);
  addBranch(*tree_, "track_end_y_cm", trackEndY_);
  addBranch(*tree_, "track_end_z_cm", trackEndZ_);
  addBranch(*tree_, "track_start_dir_x", trackStartDirX_);
  addBranch(*tree_, "track_start_dir_y", trackStartDirY_);
  addBranch(*tree_, "track_start_dir_z", trackStartDirZ_);
  addBranch(*tree_, "track_end_dir_x", trackEndDirX_);
  addBranch(*tree_, "track_end_dir_y", trackEndDirY_);
  addBranch(*tree_, "track_end_dir_z", trackEndDirZ_);
}
void PDHDBeamInstrumentation::reset() {
  // Sentinels distinguish unavailable information from valid zero response.
  available_ = unique_ = goodTrigger_ = triggerEvaluated_ = pidEvaluated_ =
      hasPerfectBeamMomentum_ = triggersMatched_ = false;
  cherenkovStatusValueAvailable_ = false;
  cherenkovPressureValueAvailable_ = false;
  mcBeamPdgValid_ = false;
  beamEventCount_ = 0;
  mcTruthRecordCount_ = 0;
  mcBeamPrimaryCount_ = 0;
  mcBeamPdg_ = -999;
  mcBeamPdgSource_ = "not_evaluated_data";
  source_ = "unavailable";
  selectedProductTag_ = "unavailable";
  pidMethod_ = "not_evaluated";
  timingTrigger_ = beamTrigger_ = ckov0_ = ckov1_ = -1;
  ckov0Pressure_ = ckov1Pressure_ =
      std::numeric_limits<double>::quiet_NaN();
  momenta_.clear();
  tof_.clear();
  tofChannels_.clear();
  pidCandidates_.clear();
  monitorNames_.clear();
  activeFiberCounts_.clear();
  trackStartX_.clear();
  trackStartY_.clear();
  trackStartZ_.clear();
  trackEndX_.clear();
  trackEndY_.clear();
  trackEndZ_.clear();
  trackStartDirX_.clear();
  trackStartDirY_.clear();
  trackStartDirZ_.clear();
  trackEndDirX_.clear();
  trackEndDirY_.clear();
  trackEndDirZ_.clear();
}

void PDHDBeamInstrumentation::fillMCBeamPdg(art::Event const &evt) {
  auto const truthHandle =
      evt.getHandle<std::vector<simb::MCTruth>>(mcTruthTag_);
  if (!truthHandle) {
    mcBeamPdgSource_ = "missing_" + mcTruthTag_.encode();
    return;
  }

  mcTruthRecordCount_ = truthHandle->size();
  mcBeamPdgSource_ = mcTruthTag_.encode() + ":noncosmic_primary";

  // This mirrors the generator-primary step used by ProtoDUNETruthUtils.
  // Cosmic-overlay records are excluded before identifying the beam primary.
  for (auto const &truth : *truthHandle) {
    if (truth.Origin() == simb::kCosmicRay) {
      continue;
    }
    for (int index = 0; index < truth.NParticles(); ++index) {
      auto const &particle = truth.GetParticle(index);
      if (particle.Process() != "primary") {
        continue;
      }
      ++mcBeamPrimaryCount_;
      // Preserve the established first generator-primary convention while
      // exposing the total count so ambiguous events remain identifiable.
      if (!mcBeamPdgValid_) {
        mcBeamPdg_ = particle.PdgCode();
        mcBeamPdgValid_ = true;
      }
    }
  }
}

void PDHDBeamInstrumentation::analyze(art::Event const &evt) {
  reset();
  ++eventsProcessed_;
  run_ = evt.run();
  subrun_ = evt.subRun();
  event_ = evt.event();
  isData_ = evt.isRealData();
  if (isData_) {
    ++dataEvents_;
  } else {
    ++mcEvents_;
    fillMCBeamPdg(evt);
  }

  // Data and MC consume distinct products but share the output schema.
  auto const tag = isData_ ? dataTag_ : mcTag_;
  selectedProductTag_ = tag.encode();
  auto const h = evt.getHandle<std::vector<beam::ProtoDUNEBeamEvent>>(tag);
  if (!h) {
    if (eventMessagesEmitted_ < maxEventMessages_) {
      mf::LogWarning("PDHDBeamInstrumentation")
          << "No vector<beam::ProtoDUNEBeamEvent> for run " << run_
          << ", subrun " << subrun_ << ", event " << event_ << " using "
          << (isData_ ? "data" : "MC") << " tag '" << tag.encode()
          << "'. The event row is retained with product_available=false.";
      ++eventMessagesEmitted_;
    }
    tree_->Fill();
    return;
  }

  available_ = true;
  ++productAvailableEvents_;
  beamEventCount_ = h->size();

  // Multiple beam events are ambiguous, so detailed fields are not selected.
  unique_ = h->size() == 1;
  if (!unique_) {
    if (eventMessagesEmitted_ < maxEventMessages_) {
      mf::LogWarning("PDHDBeamInstrumentation")
          << "Expected exactly one beam event for run " << run_ << ", subrun "
          << subrun_ << ", event " << event_ << " from tag '" << tag.encode()
          << "', but found " << beamEventCount_
          << ". Detailed fields remain unset and the event row is retained.";
      ++eventMessagesEmitted_;
    }
    tree_->Fill();
    return;
  }

  ++uniqueBeamEventEvents_;
  try {
    // Official utilities define trigger quality, momentum quality, and PID.
    auto const src = isData_ ? BeamReferenceSource::DataInstrumentation
                             : BeamReferenceSource::SimulatedInstrumentation;
    auto const evaluatePid = isData_ ? evaluateDataPid_ : evaluateMcPid_;
    pidMethod_ = evaluatePid ? "ProtoDUNEBeamlineUtils::GetPID"
                             : "not_evaluated";
    auto const r = BeamInstrumentationAlg::Extract(
        h->front(), beamline_, src, nominalMomentum_,
        isData_ || evaluateMcTrigger_, evaluatePid, monitors_);
    source_ = isData_ ? "data_instrumentation" : "simulated_instrumentation";
    goodTrigger_ = r.goodTrigger;
    triggerEvaluated_ = r.triggerEvaluated;
    pidEvaluated_ = evaluatePid;
    hasPerfectBeamMomentum_ = r.hasPerfectBeamMomentum;
    timingTrigger_ = r.timingTrigger;
    beamTrigger_ = r.beamTrigger;
    triggersMatched_ = r.triggersMatched;
    momenta_ = r.momentaGeV;
    tof_ = r.tofNs;
    tofChannels_ = r.tofChannels;
    pidCandidates_ = r.pidCandidates;
    // Store the exact Cherenkov values seen by GetPID. Separate provenance
    // flags identify whether they came from validated detector information.
    ckov0_ = r.ckov0Status;
    ckov1_ = r.ckov1Status;
    ckov0Pressure_ = r.ckov0Pressure;
    ckov1Pressure_ = r.ckov1Pressure;
    cherenkovStatusValueAvailable_ = ckov0_ >= 0 && ckov1_ >= 0;
    cherenkovPressureValueAvailable_ = std::isfinite(ckov0Pressure_) &&
                                       std::isfinite(ckov1Pressure_);
    monitorNames_ = r.monitorNames;
    activeFiberCounts_ = r.activeFiberCounts;
    triggerEvaluatedEvents_ += triggerEvaluated_;
    pidEvaluatedEvents_ += pidEvaluated_;
    goodTriggerEvents_ += goodTrigger_;

    if (eventMessagesEmitted_ < maxEventMessages_) {
      mf::LogInfo("PDHDBeamInstrumentation")
          << "Extracted " << source_ << " for run " << run_ << ", subrun "
          << subrun_ << ", event " << event_ << " from tag '" << tag.encode()
          << "': tracks=" << r.tracks.size() << ", momenta=" << momenta_.size()
          << ", TOFs=" << tof_.size()
          << ", PID candidates=" << pidCandidates_.size()
          << ", PID evaluated=" << pidEvaluated_
          << ", trigger evaluated=" << triggerEvaluated_
          << ", good trigger=" << goodTrigger_ << ".";
      ++eventMessagesEmitted_;
    }

    // Parallel vectors preserve every beamline track without choosing one.
    for (auto const &t : r.tracks) {
      trackStartX_.push_back(t.start.x);
      trackStartY_.push_back(t.start.y);
      trackStartZ_.push_back(t.start.z);
      trackEndX_.push_back(t.end.x);
      trackEndY_.push_back(t.end.y);
      trackEndZ_.push_back(t.end.z);
      trackStartDirX_.push_back(t.startDirection.x);
      trackStartDirY_.push_back(t.startDirection.y);
      trackStartDirZ_.push_back(t.startDirection.z);
      trackEndDirX_.push_back(t.endDirection.x);
      trackEndDirY_.push_back(t.endDirection.y);
      trackEndDirZ_.push_back(t.endDirection.z);
    }
  } catch (std::exception const &e) {
    ++extractionErrorEvents_;
    if (eventMessagesEmitted_ < maxEventMessages_) {
      mf::LogWarning("PDHDBeamInstrumentation")
          << "Beam extraction failed for run " << run_ << ", subrun " << subrun_
          << ", event " << event_ << " from tag '" << tag.encode()
          << "': " << e.what();
      ++eventMessagesEmitted_;
    }
  }
  tree_->Fill();
}

void PDHDBeamInstrumentation::endJob() {
  mf::LogInfo("PDHDBeamInstrumentation")
      << "Beam-instrumentation summary: processed=" << eventsProcessed_
      << " (data=" << dataEvents_ << ", MC=" << mcEvents_
      << "), product available=" << productAvailableEvents_
      << ", exactly one beam event=" << uniqueBeamEventEvents_
      << ", trigger evaluated=" << triggerEvaluatedEvents_
      << ", PID evaluated=" << pidEvaluatedEvents_
      << ", good trigger=" << goodTriggerEvents_
      << ", extraction errors=" << extractionErrorEvents_
      << ". Detailed event rows were written for every processed event.";
}
} // namespace pdhd::diagnostics
DEFINE_ART_MODULE(pdhd::diagnostics::PDHDBeamInstrumentation)
