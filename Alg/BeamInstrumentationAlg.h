/**
 * @file BeamInstrumentationAlg.h
 * @brief Extracts one ProtoDUNEBeamEvent into a common data/MC diagnostic
 * record.
 *
 * The caller selects the data or simulated product. PID is normally evaluated
 * only on data because simulated Cherenkov/TOF response may be incomplete.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_BEAMINSTRUMENTATIONALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_BEAMINSTRUMENTATIONALG_H
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/BeamSelectionAlg.h"
#include <limits>
#include <string>
#include <vector>
namespace beam {
class ProtoDUNEBeamEvent;
}
namespace protoana {
class ProtoDUNEBeamlineUtils;
}
namespace pdhd::diagnostics {
struct BeamTrackRecord {
  // Beamline-extrapolated endpoints are in cm; directions are unitless.
  Point3D start, end;
  Direction3D startDirection, endDirection;
};
struct BeamInstrumentationRecord {
  // Trigger/PID booleans are meaningful only when their evaluated flag is set.
  bool valid = false, goodTrigger = false, triggerEvaluated = false,
       hasPerfectBeamMomentum = false;
  // PDHD BeamEvent sets the legacy BI trigger to -1; use the timing/match
  // fields through IsGoodBeamlineTrigger for the supported trigger decision.
  int timingTrigger = -1, beamTrigger = -1;
  bool triggersMatched = false;
  // BeamEvent preserves the matched general-trigger timestamp as seconds plus
  // nanoseconds. Fiber timestamps retain the producer's native convention.
  double generalTriggerSeconds = std::numeric_limits<double>::quiet_NaN();
  double generalTriggerNanoseconds = std::numeric_limits<double>::quiet_NaN();
  double magnetCurrent = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> momentaGeV, tofNs;
  std::vector<int> tofChannels, pidCandidates;
  // Raw Cherenkov values are validated by the calling mode configuration.
  // Keep unavailable fields visibly distinct from physical off (status 0) or
  // on (status 1) responses. Extract() overwrites these with the direct
  // ProtoDUNEBeamEvent getters when a beam event is available.
  int ckov0Status = -1, ckov1Status = -1;
  double ckov0Pressure = std::numeric_limits<double>::quiet_NaN();
  double ckov1Pressure = std::numeric_limits<double>::quiet_NaN();
  std::vector<BeamTrackRecord> tracks;
  std::vector<std::string> monitorNames;
  std::vector<int> monitorAvailable;
  std::vector<int> activeFiberCounts;
  std::vector<double> fiberTimestampRaw;
  std::vector<std::vector<short>> activeFiberIds;
  std::vector<std::vector<int>> glitchFiberIndices;
  BeamReferenceSource source = BeamReferenceSource::Unavailable;
};
class BeamInstrumentationAlg {
public:
  static BeamInstrumentationRecord
  Extract(beam::ProtoDUNEBeamEvent const &, protoana::ProtoDUNEBeamlineUtils &,
          BeamReferenceSource, double nominalMomentumGeV, bool evaluateTrigger,
          bool evaluatePid, std::vector<std::string> const &monitors);
};
} // namespace pdhd::diagnostics
#endif
