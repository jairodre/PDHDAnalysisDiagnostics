/** @file CosmicSelectionAlg.cxx @brief Implements independently configurable
 * cosmic stages. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/CosmicSelectionAlg.h"
#include <cmath>
namespace pdhd::diagnostics {
CosmicSelectionResult
CosmicSelectionAlg::Select(CosmicSelectionInput const &in,
                           CosmicSelectionConfig const &cfg) {
  CosmicSelectionResult r;
  r.passesClearCosmic =
      !cfg.requireClearCosmic || (in.clearCosmicValid && in.isClearCosmic);
  r.passesBdt =
      !cfg.requireBdt || (in.bdtValid && std::isfinite(in.cosmicBdtScore) &&
                          in.cosmicBdtScore >= cfg.minimumCosmicBdtScore);
  r.passesStartContained =
      !cfg.requireStartContained ||
      (in.containment.volumeValid && in.containment.startContained);
  r.passesEndContained =
      !cfg.requireEndContained ||
      (in.containment.volumeValid && in.containment.endContained);
  r.passesBothEndpointsContained =
      !cfg.requireBothEndpointsContained ||
      (in.containment.volumeValid && in.containment.bothEndpointsContained);
  r.passesFullTrajectoryContained =
      !cfg.requireFullTrajectoryContained ||
      (in.containment.volumeValid && in.containment.allSampledPointsContained);
  r.selected = r.passesClearCosmic && r.passesBdt && r.passesStartContained &&
               r.passesEndContained && r.passesBothEndpointsContained &&
               r.passesFullTrajectoryContained;
  return r;
}
} // namespace pdhd::diagnostics
