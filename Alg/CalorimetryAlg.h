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
  std::size_t pointIndex = 0;
  // Valid only when the corresponding Calorimetry vector had this element.
  std::size_t trajectoryPointIndex = 0;
  bool trajectoryPointIndexValid = false;
  bool positionValid = false;
  bool dEdxValid = false;
  bool dQdxValid = false;
  bool residualRangeValid = false;
  bool pitchValid = false;
  bool electricFieldValid = false;
  bool phiValid = false;
  // Coordinates, dE/dx, residual range, and pitch use cm-based units.
  double xCm = 0.;
  double yCm = 0.;
  double zCm = 0.;
  double dEdxMeVPerCm = 0.;
  // Units are carried by CalorimetrySummary::dQdxUnits.
  double dQdx = 0.;
  double residualRangeCm = 0.;
  double pitchCm = 0.;
  // Producer convention; this layer does not assign units to E field or phi.
  double electricField = 0.;
  double phiDegrees = 0.;
};
struct CalorimetrySummary {
  bool valid = false;
  unsigned int cryostat = 0;
  unsigned int tpc = 0;
  unsigned int plane = 0;
  double kineticEnergyMeV = 0.;
  double rangeCm = 0.;
  double sumPitchCm = 0.;
  double maximumPitchCm = 0.;
  std::size_t pointCount = 0;
  std::size_t completePointCount = 0;
  std::size_t invalidPitchCount = 0;
  CalorimetryVariant variant = CalorimetryVariant::Unknown;
  bool calibrated = false;
  std::string producer;
  std::string dQdxUnits;
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
