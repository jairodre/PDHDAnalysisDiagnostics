/** @file ProductInventoryAlg.cxx @brief Implements uniform product status
 * classification. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/ProductInventoryAlg.h"
#include <utility>
namespace pdhd::diagnostics {
ProductStatus ProductInventoryAlg::Evaluate(ProductObservation observation,
                                            ProductRequirement requirement) {
  ProductStatus status;
  static_cast<ProductObservation &>(status) = std::move(observation);
  status.requirement = requirement;

  bool const configurationComplete = !status.logicalName.empty() &&
                                     !status.configuredInputTag.empty() &&
                                     !status.expectedType.empty();
  if (!configurationComplete) {
    status.state = ProductState::Invalid;
  } else if (status.handleValid) {
    status.state = ProductState::Present;
    status.usable = true;
  } else {
    status.state = requirement == ProductRequirement::Required
                       ? ProductState::MissingRequired
                       : ProductState::MissingOptional;
  }
  return status;
}
} // namespace pdhd::diagnostics
