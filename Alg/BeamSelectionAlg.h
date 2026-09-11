/**
 * @file BeamSelectionAlg.h
 * @brief Evaluates explicit, data-like beam matching and selection components.
 *
 * Both data and MC pass a reconstructed beamline reference; truth origin never
 * enters this calculation. Beamline-to-TPC residuals are recorded as
 * observables; this algorithm deliberately contains no numerical cut table.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_BEAMSELECTIONALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_BEAMSELECTIONALG_H
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/GeometryContainmentAlg.h"
#include <limits>
namespace pdhd::diagnostics {
struct Direction3D {
  double x = 0., y = 0., z = 0.;
};
enum class BeamReferenceSource {
  Unavailable,
  DataInstrumentation,
  SimulatedInstrumentation,
  TruthProjection
};
struct BeamMatchInput {
  bool beamReferenceValid = false;
  Point3D beamPositionAtReference;
  Direction3D beamDirection;
  bool recoObjectValid = false;
  Point3D recoStart;
  Direction3D recoDirection;
  BeamReferenceSource referenceSource = BeamReferenceSource::Unavailable;
};
struct BeamMatchResult {
  bool valid = false;
  // NaN is an unavailable observable, never a physically meaningful zero.
  double deltaXcm = std::numeric_limits<double>::quiet_NaN();
  double deltaYcm = std::numeric_limits<double>::quiet_NaN();
  double entranceZcm = std::numeric_limits<double>::quiet_NaN();
  double directionCosine = std::numeric_limits<double>::quiet_NaN();
  BeamReferenceSource referenceSource = BeamReferenceSource::Unavailable;
};
struct BeamSelectionInput {
  bool goodBeamTrigger = false, exactlyOneBeamlineTrack = false;
  // The nominal TPC seed is a primary PFP in Pandora's beam-tagged slice.
  bool pandoraBeamSlicePrimary = false, hasUnambiguousRecoObject = false;
};
struct BeamSelectionResult {
  bool passesGoodBeamTrigger = false, passesBeamlineMultiplicity = false;
  bool passesPandoraBeamSlicePrimary = false;
  bool passesUnambiguousRecoObject = false;
  bool selected = false;
};
class BeamSelectionAlg {
public:
  // These evaluate named stages only; neither function ranks candidates.
  static BeamMatchResult EvaluateInstrumentationMatch(BeamMatchInput const &);
  static BeamSelectionResult EvaluateCandidateStages(BeamSelectionInput const &);
};
} // namespace pdhd::diagnostics
#endif
