/** @file ProductInventoryAlg.cxx @brief Implements uniform product status
 * classification. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/ProductInventoryAlg.h"
#include <utility>
namespace pdhd::diagnostics {
ProductStatus ProductInventoryAlg::Evaluate(ProductObservation observation,
                                            ProductRequirement requirement) {
  ProductStatus r;
  static_cast<ProductObservation &>(r) = std::move(observation);
  r.requirement = requirement;
  if (r.logicalName.empty() || r.configuredInputTag.empty() ||
      r.expectedType.empty()) {
    r.state = ProductState::Invalid;
  } else if (r.handleValid) {
    r.state = ProductState::Present;
    r.usable = true;
  } else {
    r.state = requirement == ProductRequirement::Required
                  ? ProductState::MissingRequired
                  : ProductState::MissingOptional;
  }
  return r;
}
} // namespace pdhd::diagnostics
