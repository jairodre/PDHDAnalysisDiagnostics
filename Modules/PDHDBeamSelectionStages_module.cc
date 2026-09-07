/**
 * @file PDHDBeamSelectionStages_module.cc
 * @brief Applies one staged, data-like beam-candidate definition to data and
 * MC.
 *
 * Data reads beam instrumentation from DataBeamTag and MC from MCBeamTag. Both
 * use the simulated/measured beam track, Pandora beam-primary metadata, and the
 * official data position/direction windows for NominalMomentum. No MC truth is
 * read. CandidateTree stores every primary PFP and every component decision;
 * EventTree stores multiplicities and the deterministic best passing candidate.
 */
#include "TTree.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamInstrumentationAlg.h"
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamSelectionAlg.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace pdhd::diagnostics {
namespace {
BeamMatchCuts LoadOfficialDataCuts(fhicl::ParameterSet const &all,
                                   std::string const &momentum) {
  for (auto const &p : all.get<std::vector<fhicl::ParameterSet>>("DataCuts"))
    if (p.get<std::string>("Momentum") == momentum) {
      auto const x = p.get<std::vector<double>>("TrackStartXCut");
      auto const y = p.get<std::vector<double>>("TrackStartYCut");
      auto const z = p.get<std::vector<double>>("TrackStartZCut");
      if (x.size() != 2 || y.size() != 2 || z.size() != 2)
        throw cet::exception("PDHDBeamSelectionStages")
            << "Beam cut windows require two values";
      return {{x[0], x[1]},
              {y[0], y[1]},
              {z[0], z[1]},
              p.get<double>("TrackDirCut"),
              "official_data"};
    }
  throw cet::exception("PDHDBeamSelectionStages")
      << "No official data cuts for momentum " << momentum;
}
template <class T> int IndexOf(std::vector<T> const &objects, T const *object) {
  if (!object)
    return -1;
  for (std::size_t i = 0; i < objects.size(); ++i)
    if (&objects[i] == object)
      return i;
  return -1;
}
} // namespace
class PDHDBeamSelectionStages : public art::EDAnalyzer {
public:
  explicit PDHDBeamSelectionStages(fhicl::ParameterSet const &p)
      : EDAnalyzer(p), dataBeamTag_(p.get<art::InputTag>(
                           "DataBeamTag", art::InputTag{"beamevent"})),
        mcBeamTag_(
            p.get<art::InputTag>("MCBeamTag", art::InputTag{"generator"})),
        pfpTag_(
            p.get<art::InputTag>("PFParticleTag", art::InputTag{"pandora"})),
        trackTag_(
            p.get<art::InputTag>("TrackTag", art::InputTag{"pandoraTrack"})),
        showerTag_(
            p.get<art::InputTag>("ShowerTag", art::InputTag{"pandoraShower"})),
        momentumLabel_(p.get<std::string>("NominalMomentum", "1")),
        acceptTracks_(p.get<bool>("AcceptTracks", true)),
        acceptShowers_(p.get<bool>("AcceptShowers", true)),
        requireDataTrigger_(p.get<bool>("RequireGoodDataTrigger", true)),
        requireMcTrigger_(p.get<bool>("RequireGoodMCTrigger", false)),
        beamline_(p.get<fhicl::ParameterSet>("BeamlineUtils")),
        cuts_(LoadOfficialDataCuts(p.get<fhicl::ParameterSet>("BeamCuts"),
                                   momentumLabel_)) {}
  void beginJob() override;
  void analyze(art::Event const &) override;

private:
  struct Candidate {
    int pfp = -1, track = -1, shower = -1, type = 0;
    bool primary = false, beam = false, typePass = false, matchValid = false;
    bool triggerPass = false, beamEventPass = false, beamTrackPass = false;
    bool x = false, y = false, z = false, dir = false, pass = false,
         selected = false;
    double dx = 0, dy = 0, startz = 0, cos = 0, score = 0;
  };
  art::InputTag dataBeamTag_, mcBeamTag_, pfpTag_, trackTag_, showerTag_;
  std::string momentumLabel_;
  bool acceptTracks_, acceptShowers_, requireDataTrigger_, requireMcTrigger_;
  protoana::ProtoDUNEBeamlineUtils beamline_;
  BeamMatchCuts cuts_;
  TTree *eventTree_ = nullptr, *candidateTree_ = nullptr;
  unsigned int run_ = 0, subrun_ = 0, event_ = 0;
  bool isData_ = false;
  bool beamProductAvailable_ = false, beamEventUnique_ = false,
       beamTrackUnique_ = false, triggerEvaluated_ = false,
       goodTrigger_ = false;
  bool hasSelectedCandidate_ = false, selectionAmbiguous_ = false;
  int nPrimary_ = 0, nBeamPrimary_ = 0, nPassing_ = 0, selectedPfp_ = -1,
      selectedTrack_ = -1, selectedShower_ = -1;
  Candidate out_;
  std::string cutSource_ = "official_data", referenceSource_ = "unavailable";
};
void PDHDBeamSelectionStages::beginJob() {
  auto fs = art::ServiceHandle<art::TFileService>();
  eventTree_ =
      fs->make<TTree>("BeamSelectionEvent", "Beam selection event summary");
#define E(n, v) eventTree_->Branch(n, &v)
  E("run", run_);
  E("subrun", subrun_);
  E("event", event_);
  E("is_data", isData_);
  E("beam_product_available", beamProductAvailable_);
  E("beam_event_unique", beamEventUnique_);
  E("beam_track_unique", beamTrackUnique_);
  E("trigger_evaluated", triggerEvaluated_);
  E("good_trigger", goodTrigger_);
  E("n_primary_pfp", nPrimary_);
  E("n_beam_primary_pfp", nBeamPrimary_);
  E("n_passing", nPassing_);
  E("has_selected_candidate", hasSelectedCandidate_);
  E("selection_ambiguous", selectionAmbiguous_);
  E("selected_pfp_index", selectedPfp_);
  E("selected_track_index", selectedTrack_);
  E("selected_shower_index", selectedShower_);
  E("nominal_momentum", momentumLabel_);
  E("beam_match_cut_source", cutSource_);
  E("beam_reference_source", referenceSource_);
  E("cut_delta_x_min_cm", cuts_.deltaXcm.minimum);
  E("cut_delta_x_max_cm", cuts_.deltaXcm.maximum);
  E("cut_delta_y_min_cm", cuts_.deltaYcm.minimum);
  E("cut_delta_y_max_cm", cuts_.deltaYcm.maximum);
  E("cut_entrance_z_min_cm", cuts_.entranceZcm.minimum);
  E("cut_entrance_z_max_cm", cuts_.entranceZcm.maximum);
  E("cut_min_direction_cosine", cuts_.minimumDirectionCosine);
#undef E
  candidateTree_ = fs->make<TTree>("BeamSelectionCandidate",
                                   "One row per primary PFParticle");
#define C(n, v) candidateTree_->Branch(n, &v)
  C("run", run_);
  C("subrun", subrun_);
  C("event", event_);
  C("is_data", isData_);
  C("pfp_index", out_.pfp);
  C("track_index", out_.track);
  C("shower_index", out_.shower);
  C("object_type", out_.type);
  C("is_primary", out_.primary);
  C("is_pandora_beam_primary", out_.beam);
  C("passes_object_type", out_.typePass);
  C("passes_trigger_stage", out_.triggerPass);
  C("passes_beam_event_unique", out_.beamEventPass);
  C("passes_beam_track_unique", out_.beamTrackPass);
  C("match_valid", out_.matchValid);
  C("delta_x_cm", out_.dx);
  C("delta_y_cm", out_.dy);
  C("entrance_z_cm", out_.startz);
  C("direction_cosine", out_.cos);
  C("match_score", out_.score);
  C("passes_delta_x", out_.x);
  C("passes_delta_y", out_.y);
  C("passes_entrance_z", out_.z);
  C("passes_direction", out_.dir);
  C("passes_all", out_.pass);
  C("is_selected", out_.selected);
#undef C
}
void PDHDBeamSelectionStages::analyze(art::Event const &evt) {
  run_ = evt.run();
  subrun_ = evt.subRun();
  event_ = evt.event();
  isData_ = evt.isRealData();
  beamProductAvailable_ = beamEventUnique_ = beamTrackUnique_ = false;
  triggerEvaluated_ = goodTrigger_ = hasSelectedCandidate_ =
      selectionAmbiguous_ = false;
  nPrimary_ = nBeamPrimary_ = nPassing_ = 0;
  selectedPfp_ = selectedTrack_ = selectedShower_ = -1;
  referenceSource_ = "unavailable";
  BeamInstrumentationRecord bi;
  auto bh = evt.getHandle<std::vector<beam::ProtoDUNEBeamEvent>>(
      isData_ ? dataBeamTag_ : mcBeamTag_);
  beamProductAvailable_ = bool(bh);
  beamEventUnique_ = bh && bh->size() == 1;
  if (beamEventUnique_) {
    bool const eval = isData_ ? requireDataTrigger_ : requireMcTrigger_;
    bi = BeamInstrumentationAlg::Extract(
        bh->front(), beamline_,
        isData_ ? BeamReferenceSource::DataInstrumentation
                : BeamReferenceSource::SimulatedInstrumentation,
        std::stod(momentumLabel_), eval, false, {});
    triggerEvaluated_ = bi.triggerEvaluated;
    goodTrigger_ = bi.goodTrigger;
    beamTrackUnique_ = bi.tracks.size() == 1;
  }
  if (beamEventUnique_)
    referenceSource_ =
        isData_ ? "data_instrumentation" : "simulated_instrumentation";
  auto ph = evt.getHandle<std::vector<recob::PFParticle>>(pfpTag_);
  auto th = evt.getHandle<std::vector<recob::Track>>(trackTag_);
  auto sh = evt.getHandle<std::vector<recob::Shower>>(showerTag_);
  std::vector<Candidate> candidates;
  protoana::ProtoDUNEPFParticleUtils pu;
  if (ph)
    for (std::size_t i = 0; i < ph->size(); ++i) {
      auto const &p = ph->at(i);
      if (!p.IsPrimary())
        continue;
      Candidate c;
      c.pfp = i;
      c.primary = true;
      ++nPrimary_;
      try {
        c.beam = pu.IsBeamParticle(p, evt, pfpTag_.encode());
      } catch (...) {
      }
      if (c.beam)
        ++nBeamPrimary_;
      recob::Track const *tr = nullptr;
      recob::Shower const *sw = nullptr;
      try {
        if (th)
          tr = pu.GetPFParticleTrack(p, evt, pfpTag_.encode(),
                                     trackTag_.encode());
      } catch (...) {
      }
      try {
        if (sh)
          sw = pu.GetPFParticleShower(p, evt, pfpTag_.encode(),
                                      showerTag_.encode());
      } catch (...) {
      }
      c.track = th ? IndexOf(*th, tr) : -1;
      c.shower = sh ? IndexOf(*sh, sw) : -1;
      BeamMatchInput mi;
      mi.beamReferenceValid = beamTrackUnique_;
      if (beamTrackUnique_) {
        auto const &bt = bi.tracks.front();
        mi.beamPositionAtReference = bt.end;
        mi.beamDirection = bt.endDirection;
        mi.referenceSource =
            isData_ ? BeamReferenceSource::DataInstrumentation
                    : BeamReferenceSource::SimulatedInstrumentation;
      }
      if (tr && acceptTracks_) {
        c.type = 1;
        c.typePass = true;
        bool rev = tr->End().Z() < tr->Start().Z();
        auto const &pos = rev ? tr->End() : tr->Start();
        auto const dir = rev ? tr->EndDirection() : tr->StartDirection();
        mi.recoObjectValid = true;
        mi.recoStart = {pos.X(), pos.Y(), pos.Z()};
        double s = rev ? -1. : 1.;
        mi.recoDirection = {s * dir.X(), s * dir.Y(), s * dir.Z()};
      } else if (sw && acceptShowers_) {
        c.type = 2;
        c.typePass = true;
        auto const pos = sw->ShowerStart(), dir = sw->Direction();
        double s = dir.Z() < 0. ? -1. : 1.;
        mi.recoObjectValid = true;
        mi.recoStart = {pos.X(), pos.Y(), pos.Z()};
        mi.recoDirection = {s * dir.X(), s * dir.Y(), s * dir.Z()};
      }
      auto const m = BeamSelectionAlg::Match(mi, cuts_);
      c.matchValid = m.valid;
      c.dx = m.deltaXcm;
      c.dy = m.deltaYcm;
      c.startz = m.entranceZcm;
      c.cos = m.directionCosine;
      c.score = m.matchScore;
      c.x = m.passesDeltaX;
      c.y = m.passesDeltaY;
      c.z = m.passesEntranceZ;
      c.dir = m.passesDirection;
      c.triggerPass = (isData_ ? requireDataTrigger_ : requireMcTrigger_)
                          ? (triggerEvaluated_ && goodTrigger_)
                          : true;
      c.beamEventPass = beamEventUnique_;
      c.beamTrackPass = beamTrackUnique_;
      c.pass = BeamSelectionAlg::Select(
                   {c.triggerPass, c.beamTrackPass, c.beam, c.typePass, m})
                   .selected &&
               c.beamEventPass;
      if (c.pass)
        ++nPassing_;
      candidates.push_back(c);
    }
  auto best = std::min_element(candidates.begin(), candidates.end(),
                               [](auto const &a, auto const &b) {
                                 if (a.pass != b.pass)
                                   return a.pass > b.pass;
                                 if (a.score != b.score)
                                   return a.score < b.score;
                                 return a.pfp < b.pfp;
                               });
  selectionAmbiguous_ = nPassing_ > 1;
  if (best != candidates.end() && best->pass) {
    best->selected = true;
    hasSelectedCandidate_ = true;
    selectedPfp_ = best->pfp;
    selectedTrack_ = best->track;
    selectedShower_ = best->shower;
  }
  for (auto const &c : candidates) {
    out_ = c;
    candidateTree_->Fill();
  }
  eventTree_->Fill();
}
} // namespace pdhd::diagnostics
DEFINE_ART_MODULE(pdhd::diagnostics::PDHDBeamSelectionStages)
