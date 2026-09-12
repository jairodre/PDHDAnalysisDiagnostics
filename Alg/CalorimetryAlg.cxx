/**
 * @file CalorimetryAlg.cxx
 * @brief Implements bounds-checked extraction of every calorimetry vector.
 *
 * The largest input-vector size defines row count. Missing elements remain
 * invalid rather than receiving physics-looking sentinel values.
 */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include <algorithm>
#include <cmath>
#include <utility>
namespace pdhd::diagnostics {
CalorimetryResult CalorimetryAlg::Extract(anab::Calorimetry const &calorimetry,
                                          CalorimetryVariant variant,
                                          bool calibrated, std::string producer,
                                          std::string dQdxUnits) {
  CalorimetryResult result;
  auto const &dEdxValues = calorimetry.dEdx();
  auto const &dQdxValues = calorimetry.dQdx();
  auto const &residualRanges = calorimetry.ResidualRange();
  auto const &trackPitches = calorimetry.TrkPitchVec();
  auto const &positions = calorimetry.XYZ();
  auto const &trajectoryPointIndices = calorimetry.TpIndices();
  auto const &electricFields = calorimetry.Efield();
  auto const &phiValues = calorimetry.Phi();
  std::size_t const pointCount = std::max(
      {dEdxValues.size(), dQdxValues.size(), residualRanges.size(),
       trackPitches.size(), positions.size(), trajectoryPointIndices.size(),
       electricFields.size(), phiValues.size()});

  result.points.reserve(pointCount);
  for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
    CalorimetryPoint point;
    point.pointIndex = pointIndex;
    if (pointIndex < dEdxValues.size() &&
        std::isfinite(dEdxValues[pointIndex])) {
      point.dEdxValid = true;
      point.dEdxMeVPerCm = dEdxValues[pointIndex];
    }
    if (pointIndex < dQdxValues.size() &&
        std::isfinite(dQdxValues[pointIndex])) {
      point.dQdxValid = true;
      point.dQdx = dQdxValues[pointIndex];
    }
    if (pointIndex < residualRanges.size() &&
        std::isfinite(residualRanges[pointIndex])) {
      point.residualRangeValid = true;
      point.residualRangeCm = residualRanges[pointIndex];
    }
    if (pointIndex < trackPitches.size() &&
        std::isfinite(trackPitches[pointIndex]) &&
        trackPitches[pointIndex] > 0.) {
      point.pitchValid = true;
      point.pitchCm = trackPitches[pointIndex];
      result.summary.sumPitchCm += trackPitches[pointIndex];
      result.summary.maximumPitchCm = std::max(
          result.summary.maximumPitchCm, double(trackPitches[pointIndex]));
    } else {
      ++result.summary.invalidPitchCount;
    }
    if (pointIndex < positions.size() &&
        std::isfinite(positions[pointIndex].X()) &&
        std::isfinite(positions[pointIndex].Y()) &&
        std::isfinite(positions[pointIndex].Z())) {
      point.positionValid = true;
      point.xCm = positions[pointIndex].X();
      point.yCm = positions[pointIndex].Y();
      point.zCm = positions[pointIndex].Z();
    }
    if (pointIndex < trajectoryPointIndices.size()) {
      point.trajectoryPointIndexValid = true;
      point.trajectoryPointIndex = trajectoryPointIndices[pointIndex];
    }
    if (pointIndex < electricFields.size() &&
        std::isfinite(electricFields[pointIndex])) {
      point.electricFieldValid = true;
      point.electricField = electricFields[pointIndex];
    }
    if (pointIndex < phiValues.size() && std::isfinite(phiValues[pointIndex])) {
      point.phiValid = true;
      point.phiDegrees = phiValues[pointIndex];
    }
    if (point.dEdxValid && point.dQdxValid && point.residualRangeValid &&
        point.pitchValid && point.positionValid) {
      ++result.summary.completePointCount;
    }
    result.points.push_back(point);
  }
  auto const planeId = calorimetry.PlaneID();
  result.summary.cryostat = planeId.Cryostat;
  result.summary.tpc = planeId.TPC;
  result.summary.plane = planeId.Plane;
  result.summary.kineticEnergyMeV = calorimetry.KineticEnergy();
  result.summary.rangeCm = calorimetry.Range();
  result.summary.pointCount = pointCount;
  result.summary.variant = variant;
  result.summary.calibrated = calibrated;
  result.summary.producer = std::move(producer);
  result.summary.dQdxUnits = std::move(dQdxUnits);
  result.summary.valid = planeId.isValid && pointCount > 0;
  return result;
}
} // namespace pdhd::diagnostics
