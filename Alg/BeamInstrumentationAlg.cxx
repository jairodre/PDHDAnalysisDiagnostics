/** @file BeamInstrumentationAlg.cxx @brief Implements common beam-event
 * extraction. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamInstrumentationAlg.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
namespace pdhd::diagnostics {
BeamInstrumentationRecord BeamInstrumentationAlg::Extract(
    beam::ProtoDUNEBeamEvent const &e, protoana::ProtoDUNEBeamlineUtils &utils,
    BeamReferenceSource source, double momentum, bool evaluateTrigger,
    bool evaluatePid, std::vector<std::string> const &monitors) {
  BeamInstrumentationRecord r;
  r.valid = true;
  r.source = source;
  r.timingTrigger = e.GetTimingTrigger();
  // Retained for diagnostics; the v10.17 PDHD producer sets this legacy field
  // to -1 and uses timing-trigger plus beam-spill matching instead.
  r.beamTrigger = e.GetBITrigger();
  r.triggersMatched = e.CheckIsMatched();
  r.triggerEvaluated = evaluateTrigger;
  // Official PDHD decision: GetTimingTrigger()==12 && CheckIsMatched().
  if (evaluateTrigger) {
    r.goodTrigger = utils.IsGoodBeamlineTrigger(e);
  }

  // This official quality flag is not the configured nominal momentum value.
  r.hasPerfectBeamMomentum = utils.HasPerfectBeamMomentum(e);
  r.momentaGeV = e.GetRecoBeamMomenta();
  r.tofNs = e.GetTOFs();
  r.tofChannels = e.GetTOFChans();
  r.ckov0Status = e.GetCKov0Status();
  r.ckov1Status = e.GetCKov1Status();
  r.ckov0Pressure = e.GetCKov0Pressure();
  r.ckov1Pressure = e.GetCKov1Pressure();

  // GetPID applies the official momentum-dependent TOF/XCET candidate table.
  if (evaluatePid) {
    r.pidCandidates = utils.GetPID(e, momentum);
  }

  // Retain all reconstructed beamline tracks; selection belongs downstream.
  for (auto const &t : e.GetBeamTracks()) {
    auto const &tr = t.Trajectory();
    auto const s = tr.Start(), en = tr.End();
    auto const sd = tr.StartDirection(), ed = tr.EndDirection();
    r.tracks.push_back({{s.X(), s.Y(), s.Z()},
                        {en.X(), en.Y(), en.Z()},
                        {sd.X(), sd.Y(), sd.Z()},
                        {ed.X(), ed.Y(), ed.Z()}});
  }
  // Monitor occupancy exposes beam-profile inputs without storing fiber hits.
  for (auto const &name : monitors) {
    r.monitorNames.push_back(name);
    r.activeFiberCounts.push_back(e.GetActiveFibers(name).size());
  }
  return r;
}
} // namespace pdhd::diagnostics
