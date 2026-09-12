/**
 * @file ProductInventoryAlg.h
 * @brief Standardizes product-availability and provenance status records.
 *
 * The analyzer performs typed art lookups and passes observed facts here. This
 * keeps all modules consistent without expensive event-wide product discovery.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_PRODUCTINVENTORYALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_PRODUCTINVENTORYALG_H
#include <cstddef>
#include <string>
namespace pdhd::diagnostics {
enum class ProductRequirement { Optional, Required };
enum class ProductState { Present, MissingOptional, MissingRequired, Invalid };
struct ProductObservation {
  std::string logicalName;
  std::string configuredInputTag;
  std::string expectedType;
  bool handleValid = false;
  std::size_t collectionSize = 0;
  std::string resolvedModule;
  std::string resolvedInstance;
  std::string resolvedProcess;
  // Populated by the caller when a typed lookup fails.
  std::string failureReason;
};
struct ProductStatus : ProductObservation {
  ProductRequirement requirement = ProductRequirement::Optional;
  ProductState state = ProductState::Invalid;
  bool usable = false;
};
class ProductInventoryAlg {
public:
  static ProductStatus Evaluate(ProductObservation, ProductRequirement);
};
} // namespace pdhd::diagnostics
#endif
