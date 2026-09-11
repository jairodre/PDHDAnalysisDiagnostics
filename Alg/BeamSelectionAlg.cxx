/**
 * @file BeamSelectionAlg.cxx
 * @brief Implements the reproducible conjunction used by nominal beam tagging.
 */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamSelectionAlg.h"
#include <algorithm>
#include <cmath>
namespace pdhd::diagnostics {
namespace {
bool Normalize(Direction3D const &in, Direction3D &out) {
  double const n = std::sqrt(in.x * in.x + in.y * in.y + in.z * in.z);
  if (!std::isfinite(n) || n <= 0.)
    return false;
  out = {in.x / n, in.y / n, in.z / n};
  return true;
}
} // namespace
BeamMatchResult BeamSelectionAlg::EvaluateInstrumentationMatch(
    BeamMatchInput const &in) {
  BeamMatchResult r;
  r.referenceSource = in.referenceSource;
  Direction3D beam, reco;
  if (!in.beamReferenceValid || !in.recoObjectValid ||
      !Normalize(in.beamDirection, beam) || !Normalize(in.recoDirection, reco))
    return r;
  r.valid = true;
  r.deltaXcm = in.recoStart.x - in.beamPositionAtReference.x;
  r.deltaYcm = in.recoStart.y - in.beamPositionAtReference.y;
  r.entranceZcm = in.recoStart.z;
  // Round-off can move a normalized dot product just outside its physical range.
  double const dotProduct = beam.x * reco.x + beam.y * reco.y + beam.z * reco.z;
  r.directionCosine = std::max(-1., std::min(dotProduct, 1.));
  return r;
}
BeamSelectionResult BeamSelectionAlg::EvaluateCandidateStages(
    BeamSelectionInput const &in) {
  BeamSelectionResult r;
  r.passesGoodBeamTrigger = in.goodBeamTrigger;
  r.passesBeamlineMultiplicity = in.exactlyOneBeamlineTrack;
  r.passesPandoraBeamSlicePrimary = in.pandoraBeamSlicePrimary;
  r.passesUnambiguousRecoObject = in.hasUnambiguousRecoObject;
  r.selected = r.passesGoodBeamTrigger && r.passesBeamlineMultiplicity &&
               r.passesPandoraBeamSlicePrimary &&
               r.passesUnambiguousRecoObject;
  return r;
}
} // namespace pdhd::diagnostics
