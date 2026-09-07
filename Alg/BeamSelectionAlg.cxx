/**
 * @file BeamSelectionAlg.cxx
 * @brief Implements the reproducible conjunction used by nominal beam tagging.
 */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamSelectionAlg.h"
#include <algorithm>
#include <cmath>
namespace pdhd::diagnostics {
bool ClosedRange::Contains(double v) const noexcept {
  return std::isfinite(v) && std::isfinite(minimum) && std::isfinite(maximum) &&
         minimum <= maximum && v >= minimum && v <= maximum;
}
namespace {
bool Normalize(Direction3D const &in, Direction3D &out) {
  double const n = std::sqrt(in.x * in.x + in.y * in.y + in.z * in.z);
  if (!std::isfinite(n) || n <= 0.)
    return false;
  out = {in.x / n, in.y / n, in.z / n};
  return true;
}
} // namespace
BeamMatchResult BeamSelectionAlg::Match(BeamMatchInput const &in,
                                        BeamMatchCuts const &cuts) {
  BeamMatchResult r;
  r.cutSource = cuts.source;
  r.referenceSource = in.referenceSource;
  Direction3D beam, reco;
  if (!in.beamReferenceValid || !in.recoObjectValid ||
      !Normalize(in.beamDirection, beam) || !Normalize(in.recoDirection, reco))
    return r;
  r.valid = true;
  r.deltaXcm = in.recoStart.x - in.beamPositionAtReference.x;
  r.deltaYcm = in.recoStart.y - in.beamPositionAtReference.y;
  r.entranceZcm = in.recoStart.z;
  r.directionCosine =
      std::clamp(beam.x * reco.x + beam.y * reco.y + beam.z * reco.z, -1., 1.);
  r.passesDeltaX = cuts.deltaXcm.Contains(r.deltaXcm);
  r.passesDeltaY = cuts.deltaYcm.Contains(r.deltaYcm);
  r.passesEntranceZ = cuts.entranceZcm.Contains(r.entranceZcm);
  r.passesDirection = std::isfinite(cuts.minimumDirectionCosine) &&
                      r.directionCosine >= cuts.minimumDirectionCosine;
  r.passesAll = r.passesDeltaX && r.passesDeltaY && r.passesEntranceZ &&
                r.passesDirection;
  auto normalized = [](double value, ClosedRange const &range) {
    double const halfWidth = .5 * (range.maximum - range.minimum);
    return halfWidth > 0.
               ? std::abs(value - .5 * (range.minimum + range.maximum)) /
                     halfWidth
               : INFINITY;
  };
  r.matchScore = normalized(r.deltaXcm, cuts.deltaXcm) +
                 normalized(r.deltaYcm, cuts.deltaYcm) +
                 normalized(r.entranceZcm, cuts.entranceZcm) +
                 (1. - r.directionCosine);
  return r;
}
BeamSelectionResult BeamSelectionAlg::Select(BeamSelectionInput const &in) {
  BeamSelectionResult r;
  r.passesGoodBeamTrigger = in.goodBeamTrigger;
  r.passesBeamlineMultiplicity = in.exactlyOneBeamlineTrack;
  r.passesPandoraBeamPrimary = in.pandoraBeamPrimary;
  r.passesRecoObjectType = in.acceptedRecoObjectType;
  r.passesInstrumentationMatch = in.match.valid && in.match.passesAll;
  r.selected = r.passesGoodBeamTrigger && r.passesBeamlineMultiplicity &&
               r.passesPandoraBeamPrimary && r.passesRecoObjectType &&
               r.passesInstrumentationMatch;
  return r;
}
} // namespace pdhd::diagnostics
