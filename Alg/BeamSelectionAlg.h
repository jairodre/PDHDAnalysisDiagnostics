/**
 * @file BeamSelectionAlg.h
 * @brief Evaluates explicit, data-like beam matching and selection components.
 *
 * Official data-cut values are supplied by FHiCL. Both data and MC pass a
 * reconstructed beamline reference; truth origin never enters this calculation.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_BEAMSELECTIONALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_BEAMSELECTIONALG_H
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/GeometryContainmentAlg.h"
#include <string>
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
struct ClosedRange {
  double minimum = 0., maximum = 0.;
  bool Contains(double) const noexcept;
};
struct BeamMatchCuts {
  ClosedRange deltaXcm, deltaYcm, entranceZcm;
  double minimumDirectionCosine = 0.;
  std::string source = "official_data";
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
  double deltaXcm = 0., deltaYcm = 0., entranceZcm = 0., directionCosine = 0.;
  double matchScore = 0.;
  bool passesDeltaX = false, passesDeltaY = false, passesEntranceZ = false;
  bool passesDirection = false, passesAll = false;
  std::string cutSource;
  BeamReferenceSource referenceSource = BeamReferenceSource::Unavailable;
};
struct BeamSelectionInput {
  bool goodBeamTrigger = false, exactlyOneBeamlineTrack = false;
  bool pandoraBeamPrimary = false, acceptedRecoObjectType = false;
  BeamMatchResult match;
};
struct BeamSelectionResult {
  bool passesGoodBeamTrigger = false, passesBeamlineMultiplicity = false;
  bool passesPandoraBeamPrimary = false, passesRecoObjectType = false;
  bool passesInstrumentationMatch = false, selected = false;
};
class BeamSelectionAlg {
public:
  static BeamMatchResult Match(BeamMatchInput const &, BeamMatchCuts const &);
  static BeamSelectionResult Select(BeamSelectionInput const &);
};
} // namespace pdhd::diagnostics
#endif
