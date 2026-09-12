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
  double const magnitude = std::sqrt(in.x * in.x + in.y * in.y + in.z * in.z);
  if (!std::isfinite(magnitude) || magnitude <= 0.)
    return false;
  out = {in.x / magnitude, in.y / magnitude, in.z / magnitude};
  return true;
}
} // namespace
BeamMatchResult
BeamSelectionAlg::EvaluateInstrumentationMatch(BeamMatchInput const &in) {
  BeamMatchResult result;
  result.referenceSource = in.referenceSource;
  Direction3D normalizedBeamDirection;
  Direction3D normalizedRecoDirection;
  result.positionValid = in.beamReferenceValid && in.recoStartValid;
  if (!result.positionValid) {
    // Keep direction-dependent observables unavailable without a reference.
    return result;
  }
  result.deltaXcm = in.recoStart.x - in.beamPositionAtReference.x;
  result.deltaYcm = in.recoStart.y - in.beamPositionAtReference.y;
  result.deltaZcm = in.recoStart.z - in.beamPositionAtReference.z;
  result.entranceZcm = in.recoStart.z;
  result.directionValid =
      in.recoDirectionValid &&
      Normalize(in.beamDirection, normalizedBeamDirection) &&
      Normalize(in.recoDirection, normalizedRecoDirection);
  if (!result.directionValid) {
    // Position residuals remain valid even when either direction cannot
    // normalize.
    return result;
  }
  // Round-off can move a normalized dot product just outside its physical
  // range.
  double const dotProduct =
      normalizedBeamDirection.x * normalizedRecoDirection.x +
      normalizedBeamDirection.y * normalizedRecoDirection.y +
      normalizedBeamDirection.z * normalizedRecoDirection.z;
  result.directionCosine = std::max(-1., std::min(dotProduct, 1.));
  result.valid = true;
  return result;
}
BeamSelectionResult
BeamSelectionAlg::EvaluateCandidateStages(BeamSelectionInput const &in) {
  BeamSelectionResult result;
  result.passesGoodBeamTrigger = in.goodBeamTrigger;
  result.passesBeamlineMultiplicity = in.exactlyOneBeamlineTrack;
  result.passesPandoraBeamSlicePrimary = in.pandoraBeamSlicePrimary;
  result.passesUnambiguousRecoObject = in.hasUnambiguousRecoObject;
  result.selected = result.passesGoodBeamTrigger &&
                    result.passesBeamlineMultiplicity &&
                    result.passesPandoraBeamSlicePrimary &&
                    result.passesUnambiguousRecoObject;
  return result;
}
} // namespace pdhd::diagnostics
