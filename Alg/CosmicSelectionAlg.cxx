/** @file CosmicSelectionAlg.cxx @brief Implements independently configurable
 * cosmic stages. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/CosmicSelectionAlg.h"
#include <cmath>
namespace pdhd::diagnostics {
CosmicSelectionResult
CosmicSelectionAlg::Select(CosmicSelectionInput const &input,
                           CosmicSelectionConfig const &config) {
  CosmicSelectionResult result;

  // A disabled stage passes; an enabled stage requires valid input.
  result.passesClearCosmic = !config.requireClearCosmic ||
                             (input.clearCosmicValid && input.isClearCosmic);
  result.passesBdt = !config.requireBdt ||
                     (input.bdtValid && std::isfinite(input.cosmicBdtScore) &&
                      input.cosmicBdtScore >= config.minimumCosmicBdtScore);

  bool const containmentValid = input.containment.volumeValid;
  result.passesStartContained =
      !config.requireStartContained ||
      (containmentValid && input.containment.startContained);
  result.passesEndContained =
      !config.requireEndContained ||
      (containmentValid && input.containment.endContained);
  result.passesBothEndpointsContained =
      !config.requireBothEndpointsContained ||
      (containmentValid && input.containment.bothEndpointsContained);
  result.passesFullTrajectoryContained =
      !config.requireFullTrajectoryContained ||
      (containmentValid && input.containment.allSampledPointsContained);

  result.selected = result.passesClearCosmic && result.passesBdt &&
                    result.passesStartContained && result.passesEndContained &&
                    result.passesBothEndpointsContained &&
                    result.passesFullTrajectoryContained;
  return result;
}
} // namespace pdhd::diagnostics
