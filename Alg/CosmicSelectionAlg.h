/**
 * @file CosmicSelectionAlg.h
 * @brief Combines explicit reconstructed cosmic diagnostics without hidden
 * cuts.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_COSMICSELECTIONALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_COSMICSELECTIONALG_H
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/GeometryContainmentAlg.h"
namespace pdhd::diagnostics {
struct CosmicSelectionConfig {
  bool requireClearCosmic = false;
  bool requireBdt = false;
  bool requireStartContained = false;
  bool requireEndContained = false;
  bool requireBothEndpointsContained = false;
  bool requireFullTrajectoryContained = false;
  double minimumCosmicBdtScore = 0.;
};
struct CosmicSelectionInput {
  // Invalid classifier inputs fail only when their corresponding cut is used.
  bool clearCosmicValid = false;
  bool isClearCosmic = false;
  bool bdtValid = false;
  double cosmicBdtScore = 0.;
  ContainmentResult containment;
};
struct CosmicSelectionResult {
  // Disabled stages intentionally pass, so every result records its full gate.
  bool passesClearCosmic = true;
  bool passesBdt = true;
  bool passesStartContained = true;
  bool passesEndContained = true;
  bool passesBothEndpointsContained = true;
  bool passesFullTrajectoryContained = true;
  bool selected = false;
};
class CosmicSelectionAlg {
public:
  static CosmicSelectionResult Select(CosmicSelectionInput const &,
                                      CosmicSelectionConfig const &);
};
} // namespace pdhd::diagnostics
#endif
