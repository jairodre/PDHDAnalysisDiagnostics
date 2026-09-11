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
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
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
        ifbeamCherenkovSource_(p.get<std::string>("IFBeamCherenkovSource")),
        ifbeamTimestampUnit_(p.get<std::string>("IFBeamTimestampUnit")),
        ifbeamCherenkovProvenanceValid_(
            p.get<bool>("IFBeamCherenkovProvenanceValid")),
        ifbeamCkov0ChannelName_(p.get<std::string>("IFBeamCKov0ChannelName")),
        ifbeamCkov1ChannelName_(p.get<std::string>("IFBeamCKov1ChannelName")),
        nonNominalInputHandlingUsed_(
            p.get<bool>("NonNominalInputHandlingUsed")),
        printSummaryToStdout_(p.get<bool>("PrintSummaryToStdout", true)),
        maxEventMessages_(p.get<unsigned int>("MaxEventMessages", 3U)),
        printEventDetailsToStdout_(
            p.get<bool>("PrintEventDetailsToStdout", true)),
        maxEventDetailPrintouts_(
            p.get<unsigned int>("MaxEventDetailPrintouts", 3U)),
        monitors_(p.get<std::vector<std::string>>(
            "FiberMonitorNames",
            {"XBPF022697", "XBPF022698", "XBPF022701", "XBPF022702",
             "XBPF022707", "XBPF022708", "XBPF022716", "XBPF022717"})),
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
  // The official BeamEvent producer can obtain these values from IFBeam. The
  // configuration must explicitly certify the channel mapping before this
  // derived view is considered a found IFBeam record.
  std::string ifbeamCherenkovSource_, ifbeamTimestampUnit_;
  bool ifbeamCherenkovProvenanceValid_;
  std::string ifbeamCkov0ChannelName_, ifbeamCkov1ChannelName_;
  bool nonNominalInputHandlingUsed_;
  bool printSummaryToStdout_;
  unsigned int maxEventMessages_;
  // Detailed stdout reports are bounded by default. A FHiCL value of zero
  // explicitly requests every eligible event, matching legacy BeamlineReco.
  bool printEventDetailsToStdout_;
  unsigned int maxEventDetailPrintouts_;
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
  unsigned int eventDetailPrintoutsEmitted_ = 0;

  // Stable per-event storage used directly by the output tree branches.
  TTree *tree_ = nullptr;
  unsigned int run_ = 0, subrun_ = 0, event_ = 0;
  bool isData_ = false, available_ = false, unique_ = false;
  bool pidEvaluated_ = false;
  bool mcBeamPdgValid_ = false;
  // Recorded means the BeamEvent getter returned a non-sentinel field.  A
  // recorded field becomes analysis-usable only when its configured source is
  // validated; SkipLLT zeroes must never be presented as detector data.
  bool cherenkovStatusRecorded_ = false;
  bool cherenkovPressureRecorded_ = false;
  bool cherenkovStatusValueAvailable_ = false;
  bool cherenkovPressureValueAvailable_ = false;
  bool ifbeamCkov0RecordFound_ = false;
  bool ifbeamCkov1RecordFound_ = false;
  bool ifbeamCkovRecordFound_ = false;
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
  // Direct product values and the certified IFBeam-derived view are retained
  // independently. The latter remains NaN until its mapping is validated.
  double ckov0PressureIfbeam_ = std::numeric_limits<double>::quiet_NaN();
  double ckov1PressureIfbeam_ = std::numeric_limits<double>::quiet_NaN();
  double ckov0IfbeamTimestampRaw_ = std::numeric_limits<double>::quiet_NaN();
  double ckov1IfbeamTimestampRaw_ = std::numeric_limits<double>::quiet_NaN();
  double generalTriggerSeconds_ = std::numeric_limits<double>::quiet_NaN();
  double generalTriggerNanoseconds_ =
      std::numeric_limits<double>::quiet_NaN();
  double magnetCurrent_ = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> momenta_, tof_, trackStartX_, trackStartY_, trackStartZ_,
      trackEndX_, trackEndY_, trackEndZ_;
  std::vector<double> trackStartDirX_, trackStartDirY_, trackStartDirZ_,
      trackEndDirX_, trackEndDirY_, trackEndDirZ_;
  std::vector<int> tofChannels_, pidCandidates_, activeFiberCounts_;
  std::vector<std::string> monitorNames_;
  std::vector<int> monitorAvailable_;
  std::vector<double> fiberTimestampRaw_;
  // Flat vectors plus monitor-aligned offsets avoid nested STL ROOT branches.
  std::vector<short> activeFiberIdsFlat_;
  std::vector<unsigned int> activeFiberOffsets_;
  std::vector<int> glitchFiberIndicesFlat_;
  std::vector<unsigned int> glitchFiberOffsets_;
  void reset();
  void fillMCBeamPdg(art::Event const &);
  void printEventDetails(BeamInstrumentationRecord const &);
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
      << ", IFBeamCherenkovSource='" << ifbeamCherenkovSource_
      << "', IFBeamCherenkovProvenanceValid="
      << ifbeamCherenkovProvenanceValid_
      << ", NonNominalInputHandlingUsed=" << nonNominalInputHandlingUsed_
      << ", MaxEventMessages=" << maxEventMessages_
      << ", PrintEventDetailsToStdout=" << printEventDetailsToStdout_
      << ", MaxEventDetailPrintouts=" << maxEventDetailPrintouts_
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
  addBranch(*tree_, "ckov_status_recorded", cherenkovStatusRecorded_);
  addBranch(*tree_, "ckov_status_provenance_valid",
            cherenkovStatusProvenanceValid_);
  addBranch(*tree_, "ckov_status_source", cherenkovStatusSource_);
  addBranch(*tree_, "ckov_pressure_value_available",
            cherenkovPressureValueAvailable_);
  addBranch(*tree_, "ckov_pressure_recorded", cherenkovPressureRecorded_);
  addBranch(*tree_, "ckov_pressure_provenance_valid",
            cherenkovPressureProvenanceValid_);
  addBranch(*tree_, "ckov_pressure_source", cherenkovPressureSource_);
  addBranch(*tree_, "ckov_ifbeam_source", ifbeamCherenkovSource_);
  addBranch(*tree_, "ckov_ifbeam_timestamp_unit", ifbeamTimestampUnit_);
  addBranch(*tree_, "ckov0_ifbeam_channel_name", ifbeamCkov0ChannelName_);
  addBranch(*tree_, "ckov1_ifbeam_channel_name", ifbeamCkov1ChannelName_);
  addBranch(*tree_, "ckov0_ifbeam_record_found", ifbeamCkov0RecordFound_);
  addBranch(*tree_, "ckov1_ifbeam_record_found", ifbeamCkov1RecordFound_);
  addBranch(*tree_, "ckov_ifbeam_record_found", ifbeamCkovRecordFound_);

  // These are uncut beam-event observables; units are encoded in their names.
  addBranch(*tree_, "HasPerfectBeamMomentum", hasPerfectBeamMomentum_);
  addBranch(*tree_, "GetTimingTrigger", timingTrigger_);
  addBranch(*tree_, "GetBITrigger", beamTrigger_);
  addBranch(*tree_, "CheckIsMatched", triggersMatched_);
  addBranch(*tree_, "general_trigger_seconds", generalTriggerSeconds_);
  addBranch(*tree_, "general_trigger_nanoseconds", generalTriggerNanoseconds_);
  addBranch(*tree_, "magnet_current", magnetCurrent_);
  addBranch(*tree_, "momenta_GeV", momenta_);
  addBranch(*tree_, "tof_ns", tof_);
  addBranch(*tree_, "tof_channels", tofChannels_);
  addBranch(*tree_, "pid_candidates", pidCandidates_);
  addBranch(*tree_, "ckov0_status", ckov0_);
  addBranch(*tree_, "ckov1_status", ckov1_);
  addBranch(*tree_, "ckov0_pressure", ckov0Pressure_);
  addBranch(*tree_, "ckov1_pressure", ckov1Pressure_);
  addBranch(*tree_, "ckov0_pressure_beamevent", ckov0Pressure_);
  addBranch(*tree_, "ckov1_pressure_beamevent", ckov1Pressure_);
  addBranch(*tree_, "ckov0_pressure_ifbeam", ckov0PressureIfbeam_);
  addBranch(*tree_, "ckov1_pressure_ifbeam", ckov1PressureIfbeam_);
  addBranch(*tree_, "ckov0_ifbeam_timestamp_raw",
            ckov0IfbeamTimestampRaw_);
  addBranch(*tree_, "ckov1_ifbeam_timestamp_raw",
            ckov1IfbeamTimestampRaw_);
  addBranch(*tree_, "monitor_names", monitorNames_);
  addBranch(*tree_, "monitor_available", monitorAvailable_);
  addBranch(*tree_, "active_fiber_counts", activeFiberCounts_);
  addBranch(*tree_, "fiber_timestamp_raw", fiberTimestampRaw_);
  addBranch(*tree_, "active_fiber_ids_flat", activeFiberIdsFlat_);
  addBranch(*tree_, "active_fiber_offsets", activeFiberOffsets_);
  addBranch(*tree_, "glitch_fiber_indices_flat", glitchFiberIndicesFlat_);
  addBranch(*tree_, "glitch_fiber_offsets", glitchFiberOffsets_);
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
  cherenkovStatusRecorded_ = false;
  cherenkovPressureRecorded_ = false;
  ifbeamCkov0RecordFound_ = false;
  ifbeamCkov1RecordFound_ = false;
  ifbeamCkovRecordFound_ = false;
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
  ckov0PressureIfbeam_ = ckov1PressureIfbeam_ =
      std::numeric_limits<double>::quiet_NaN();
  ckov0IfbeamTimestampRaw_ = ckov1IfbeamTimestampRaw_ =
      std::numeric_limits<double>::quiet_NaN();
  generalTriggerSeconds_ = generalTriggerNanoseconds_ = magnetCurrent_ =
      std::numeric_limits<double>::quiet_NaN();
  momenta_.clear();
  tof_.clear();
  tofChannels_.clear();
  pidCandidates_.clear();
  monitorNames_.clear();
  monitorAvailable_.clear();
  activeFiberCounts_.clear();
  fiberTimestampRaw_.clear();
  activeFiberIdsFlat_.clear();
  activeFiberOffsets_.clear();
  glitchFiberIndicesFlat_.clear();
  glitchFiberOffsets_.clear();
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

void PDHDBeamInstrumentation::printEventDetails(
    BeamInstrumentationRecord const &record) {
  if (!printEventDetailsToStdout_ ||
      (maxEventDetailPrintouts_ != 0 &&
       eventDetailPrintoutsEmitted_ >= maxEventDetailPrintouts_)) {
    return;
  }

  // This mirrors the legacy BeamlineReco inspection order, but keeps each
  // source and validity decision visible rather than implying a selection.
  std::ostringstream report;
  report << std::boolalpha << std::fixed << std::setprecision(3)
         << "\n================================================================\n"
         << " PDHDBeamInstrumentation | run " << run_ << "  subrun " << subrun_
         << "  event " << event_ << "\n"
         << "================================================================\n"
         << "[Input product]\n"
         << "  tag                : " << selectedProductTag_ << "\n"
         << "  entries            : " << beamEventCount_ << "\n"
         << "  available / unique : " << static_cast<bool>(available_) << " / "
         << static_cast<bool>(unique_) << "\n"
         << "\n[Beam-event quality]\n"
         << "  timing trigger     : " << timingTrigger_ << "\n"
         << "  BI trigger         : " << beamTrigger_ << "\n"
         << "  timing beam trigger: " << (timingTrigger_ == 12) << "\n"
         << "  beamline matched   : " << static_cast<bool>(triggersMatched_) << "\n"
         << "  good beam trigger  : " << static_cast<bool>(goodTrigger_) << "\n"
         << "  perfect momentum   : " << static_cast<bool>(hasPerfectBeamMomentum_)
         << "\n"
         << "\n[Spill and beamline]\n"
         << "  trigger time       : " << generalTriggerSeconds_ << " s + "
         << generalTriggerNanoseconds_ << " ns\n"
         << "  magnet current     : " << magnetCurrent_ << "\n"
         << "  reconstructed p    : ";
  if (momenta_.empty()) {
    report << "unavailable";
  } else {
    for (std::size_t index = 0; index < momenta_.size(); ++index) {
      report << (index == 0 ? "" : ", ") << momenta_[index] << " GeV/c";
    }
  }
  report << "\n  TOF measurements   : ";
  if (tof_.empty()) {
    report << "unavailable";
  }
  for (std::size_t index = 0; index < tof_.size(); ++index) {
    report << (index == 0 ? "" : "; ") << tof_[index] << " ns (channel ";
    if (index < tofChannels_.size()) {
      report << tofChannels_[index];
    } else {
      report << "missing";
    }
    report << ')';
  }
  if (tof_.size() != tofChannels_.size()) {
    report << " [invalid: TOF/channel counts differ]";
  }
  report << "\n\n[Cherenkov]\n"
         << "  BeamEvent status   : C0=" << ckov0_ << ", C1=" << ckov1_
         << "  (recorded=" << static_cast<bool>(cherenkovStatusRecorded_)
         << ", analysis usable="
         << static_cast<bool>(cherenkovStatusValueAvailable_)
         << ", provenance valid="
         << static_cast<bool>(cherenkovStatusProvenanceValid_) << ")\n"
         << "  BeamEvent pressure : C0=" << ckov0Pressure_ << ", C1="
         << ckov1Pressure_ << "  (recorded="
         << static_cast<bool>(cherenkovPressureRecorded_)
         << ", analysis usable="
         << static_cast<bool>(cherenkovPressureValueAvailable_)
         << ", provenance valid="
         << static_cast<bool>(cherenkovPressureProvenanceValid_) << ")\n"
         << "  BeamEvent source   : " << cherenkovStatusSource_ << "\n"
         << "  IFBeam records     : C0=" << static_cast<bool>(ifbeamCkov0RecordFound_)
         << ", C1=" << static_cast<bool>(ifbeamCkov1RecordFound_) << "\n"
         << "  IFBeam pressure    : C0=" << ckov0PressureIfbeam_ << ", C1="
         << ckov1PressureIfbeam_ << "\n"
         << "  IFBeam source      : " << ifbeamCherenkovSource_ << "\n"
         << "\n[PID candidates]\n"
         << "  method / evaluated : " << pidMethod_ << " / "
         << static_cast<bool>(pidEvaluated_) << "\n"
         << "  PDG candidates     : ";
  if (pidCandidates_.empty()) {
    report << "none";
  } else {
    for (std::size_t index = 0; index < pidCandidates_.size(); ++index) {
      report << (index == 0 ? "" : ", ") << pidCandidates_[index];
    }
  }
  report << "\n\n[Profile monitors: active fibers]\n";
  for (std::size_t index = 0; index < monitorNames_.size(); ++index) {
    report << "  " << std::left << std::setw(11) << monitorNames_[index]
           << " available=";
    if (index < monitorAvailable_.size()) {
      report << static_cast<bool>(monitorAvailable_[index]);
    } else {
      report << "missing_availability";
    }
    report << "  timestamp_raw=";
    if (index < fiberTimestampRaw_.size()) {
      report << fiberTimestampRaw_[index];
    } else {
      report << "missing_timestamp";
    }
    report << "  count=";
    if (index < activeFiberCounts_.size()) {
      report << activeFiberCounts_[index];
    } else {
      report << "missing_count";
    }
    report << "  IDs=[";
    if (index < record.activeFiberIds.size()) {
      auto const &fibers = record.activeFiberIds[index];
      for (std::size_t fiberIndex = 0; fiberIndex < fibers.size(); ++fiberIndex) {
        report << (fiberIndex == 0 ? "" : ", ") << fibers[fiberIndex];
      }
    }
    report << "]  glitches=[";
    if (index < record.glitchFiberIndices.size()) {
      auto const &glitches = record.glitchFiberIndices[index];
      for (std::size_t glitchIndex = 0; glitchIndex < glitches.size();
           ++glitchIndex) {
        report << (glitchIndex == 0 ? "" : ", ") << glitches[glitchIndex];
      }
    }
    report << "]\n";
  }
  report << "\n[Beamline tracks] count=" << record.tracks.size() << "\n";
  for (std::size_t index = 0; index < record.tracks.size(); ++index) {
    auto const &track = record.tracks[index];
    report << "  [" << index << "] start [cm] = (" << track.start.x << ", "
           << track.start.y << ", " << track.start.z << ")  dir = ("
           << track.startDirection.x << ',' << track.startDirection.y << ','
           << track.startDirection.z << ")\n"
           << "      end   [cm] = (" << track.end.x << ", " << track.end.y
           << ", " << track.end.z << ")  dir = ("
           << track.endDirection.x << ',' << track.endDirection.y << ','
           << track.endDirection.z << ")\n";
  }
  if (!isData_) {
    report << "\n[MC truth beam]\n"
           << "  PDG / valid         : " << mcBeamPdg_ << " / "
           << static_cast<bool>(mcBeamPdgValid_) << "\n"
           << "  source / primaries  : " << mcBeamPdgSource_ << " / "
           << mcBeamPrimaryCount_ << "\n";
  }
  report << "================================================================\n";
  std::cout << report.str() << std::endl;
  ++eventDetailPrintoutsEmitted_;
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
    generalTriggerSeconds_ = r.generalTriggerSeconds;
    generalTriggerNanoseconds_ = r.generalTriggerNanoseconds;
    magnetCurrent_ = r.magnetCurrent;
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
    cherenkovStatusRecorded_ = ckov0_ >= 0 && ckov1_ >= 0;
    cherenkovPressureRecorded_ = std::isfinite(ckov0Pressure_) &&
                                 std::isfinite(ckov1Pressure_);
    // SkipLLT emits default-looking values when LLT information is absent.
    // Keep those raw fields inspectable, but do not make them usable inputs to
    // PID or analysis until the channel source has been validated.
    cherenkovStatusValueAvailable_ =
        cherenkovStatusRecorded_ && cherenkovStatusProvenanceValid_;
    cherenkovPressureValueAvailable_ =
        cherenkovPressureRecorded_ && cherenkovPressureProvenanceValid_;
    // pdhd_beamevent is the only IFBeam access in this job. Do not label its
    // values as IFBeam records until FHiCL certifies the PDHD channel mapping
    // and the producer supplied detector timestamps for both Cherenkovs.
    ckov0IfbeamTimestampRaw_ = h->front().GetCKov0Time();
    ckov1IfbeamTimestampRaw_ = h->front().GetCKov1Time();
    ifbeamCkov0RecordFound_ =
        isData_ && beamEventProducedInCurrentJob_ &&
        ifbeamCherenkovProvenanceValid_ && std::isfinite(ckov0Pressure_) &&
        std::isfinite(ckov0IfbeamTimestampRaw_) &&
        ckov0IfbeamTimestampRaw_ > 0.;
    ifbeamCkov1RecordFound_ =
        isData_ && beamEventProducedInCurrentJob_ &&
        ifbeamCherenkovProvenanceValid_ && std::isfinite(ckov1Pressure_) &&
        std::isfinite(ckov1IfbeamTimestampRaw_) &&
        ckov1IfbeamTimestampRaw_ > 0.;
    ifbeamCkovRecordFound_ =
        ifbeamCkov0RecordFound_ && ifbeamCkov1RecordFound_;
    if (ifbeamCkov0RecordFound_) {
      ckov0PressureIfbeam_ = ckov0Pressure_;
    }
    if (ifbeamCkov1RecordFound_) {
      ckov1PressureIfbeam_ = ckov1Pressure_;
    }
    monitorNames_ = r.monitorNames;
    monitorAvailable_ = r.monitorAvailable;
    activeFiberCounts_ = r.activeFiberCounts;
    fiberTimestampRaw_ = r.fiberTimestampRaw;
    // IDs for monitor i occupy [offset[i], offset[i + 1]) in the flat branch.
    activeFiberOffsets_.push_back(0U);
    glitchFiberOffsets_.push_back(0U);
    for (std::size_t index = 0; index < r.monitorNames.size(); ++index) {
      if (index < r.activeFiberIds.size()) {
        activeFiberIdsFlat_.insert(activeFiberIdsFlat_.end(),
                                   r.activeFiberIds[index].begin(),
                                   r.activeFiberIds[index].end());
      }
      activeFiberOffsets_.push_back(
          static_cast<unsigned int>(activeFiberIdsFlat_.size()));
      if (index < r.glitchFiberIndices.size()) {
        glitchFiberIndicesFlat_.insert(glitchFiberIndicesFlat_.end(),
                                       r.glitchFiberIndices[index].begin(),
                                       r.glitchFiberIndices[index].end());
      }
      glitchFiberOffsets_.push_back(
          static_cast<unsigned int>(glitchFiberIndicesFlat_.size()));
    }
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

    // Write the complete human-readable inspection record before flattening
    // tracks into tree vectors; no candidate or track is selected here.
    printEventDetails(r);

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
  std::ostringstream summary;
  summary << "Beam-instrumentation summary: processed=" << eventsProcessed_
          << " (data=" << dataEvents_ << ", MC=" << mcEvents_
          << "), product available=" << productAvailableEvents_
          << ", exactly one beam event=" << uniqueBeamEventEvents_
          << ", trigger evaluated=" << triggerEvaluatedEvents_
          << ", PID evaluated=" << pidEvaluatedEvents_
          << ", good trigger=" << goodTriggerEvents_
          << ", extraction errors=" << extractionErrorEvents_
          << ". Detailed event rows were written for every processed event.";
  mf::LogInfo("PDHDBeamInstrumentation") << summary.str();
  // Preserve the summary when the configured message destination filters it.
  if (printSummaryToStdout_) {
    std::cout << "PDHDBeamInstrumentation: " << summary.str() << std::endl;
  }
}
} // namespace pdhd::diagnostics
DEFINE_ART_MODULE(pdhd::diagnostics::PDHDBeamInstrumentation)
