/** @file RecoAssociationAlg.cxx @brief Implements reco-association cardinality
 * summaries. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/RecoAssociationAlg.h"
#include <utility>
namespace pdhd::diagnostics {
AssociationSummary RecoAssociationAlg::Summarize(std::size_t count, bool valid,
                                                 std::string source) {
  AssociationSummary r;
  r.count = count;
  r.associationProductValid = valid;
  r.source = std::move(source);
  if (!valid)
    r.cardinality = AssociationCardinality::Invalid;
  else if (count == 0)
    r.cardinality = AssociationCardinality::None;
  else if (count == 1)
    r.cardinality = AssociationCardinality::Unique;
  else
    r.cardinality = AssociationCardinality::Ambiguous;
  return r;
}
} // namespace pdhd::diagnostics
