/** @file RecoAssociationAlg.cxx @brief Implements reco-association cardinality
 * summaries. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/RecoAssociationAlg.h"
#include <utility>
namespace pdhd::diagnostics {
AssociationSummary RecoAssociationAlg::Summarize(std::size_t count, bool valid,
                                                 std::string source) {
  AssociationSummary summary;
  summary.count = count;
  summary.associationProductValid = valid;
  summary.source = std::move(source);

  if (!valid) {
    summary.cardinality = AssociationCardinality::Invalid;
  } else if (count == 0) {
    summary.cardinality = AssociationCardinality::None;
  } else if (count == 1) {
    summary.cardinality = AssociationCardinality::Unique;
  } else {
    summary.cardinality = AssociationCardinality::Ambiguous;
  }
  return summary;
}
} // namespace pdhd::diagnostics
