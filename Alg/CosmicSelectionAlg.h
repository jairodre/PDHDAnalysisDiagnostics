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
  bool requireClearCosmic = false, requireBdt = false;
  bool requireStartContained = false, requireEndContained = false;
  bool requireBothEndpointsContained = false,
       requireFullTrajectoryContained = false;
  double minimumCosmicBdtScore = 0.;
};
struct CosmicSelectionInput {
  bool clearCosmicValid = false, isClearCosmic = false;
  bool bdtValid = false;
  double cosmicBdtScore = 0.;
  ContainmentResult containment;
};
struct CosmicSelectionResult {
  bool passesClearCosmic = true, passesBdt = true, passesStartContained = true;
  bool passesEndContained = true, passesBothEndpointsContained = true;
  bool passesFullTrajectoryContained = true, selected = false;
};
class CosmicSelectionAlg {
public:
  static CosmicSelectionResult Select(CosmicSelectionInput const &,
                                      CosmicSelectionConfig const &);
};
} // namespace pdhd::diagnostics
#endif
