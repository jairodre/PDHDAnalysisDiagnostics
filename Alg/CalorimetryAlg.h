/**
 * @file CalorimetryAlg.h
 * @brief Extracts validated point and summary records from anab::Calorimetry.
 *
 * dE/dx is MeV/cm, pitch and residual range are cm, and coordinates are cm.
 * dQ/dx units depend on the producer and are carried as provenance metadata.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_CALORIMETRYALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_CALORIMETRYALG_H
#include <cstddef>
#include <string>
#include <vector>
namespace anab {
class Calorimetry;
}
namespace pdhd::diagnostics {
enum class CalorimetryVariant { Unknown, SCE, NoSCE };
struct CalorimetryPoint {
  std::size_t pointIndex = 0, trajectoryPointIndex = 0;
  bool trajectoryPointIndexValid = false, positionValid = false;
  bool dEdxValid = false, dQdxValid = false, residualRangeValid = false;
  bool pitchValid = false, electricFieldValid = false, phiValid = false;
  double xCm = 0., yCm = 0., zCm = 0., dEdxMeVPerCm = 0., dQdx = 0.;
  double residualRangeCm = 0., pitchCm = 0., electricField = 0.,
         phiDegrees = 0.;
};
struct CalorimetrySummary {
  bool valid = false;
  unsigned int cryostat = 0, tpc = 0, plane = 0;
  double kineticEnergyMeV = 0., rangeCm = 0., sumPitchCm = 0.,
         maximumPitchCm = 0.;
  std::size_t pointCount = 0, completePointCount = 0, invalidPitchCount = 0;
  CalorimetryVariant variant = CalorimetryVariant::Unknown;
  bool calibrated = false;
  std::string producer, dQdxUnits;
};
struct CalorimetryResult {
  CalorimetrySummary summary;
  std::vector<CalorimetryPoint> points;
};
class CalorimetryAlg {
public:
  static CalorimetryResult Extract(anab::Calorimetry const &,
                                   CalorimetryVariant, bool calibrated,
                                   std::string producer, std::string dQdxUnits);
};
} // namespace pdhd::diagnostics
#endif
