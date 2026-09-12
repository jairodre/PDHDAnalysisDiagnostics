/**
 * @file RecoAssociationAlg.h
 * @brief Common cardinality and ambiguity handling for reco associations.
 *
 * Official utilities or art association helpers obtain objects; this layer
 * preserves multiplicity and prevents implicit selection of element zero.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_RECOASSOCIATIONALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_RECOASSOCIATIONALG_H
#include <cstddef>
#include <string>
#include <vector>
namespace pdhd::diagnostics {
enum class AssociationCardinality { None, Unique, Ambiguous, Invalid };
struct AssociationSummary {
  AssociationCardinality cardinality = AssociationCardinality::Invalid;
  std::size_t count = 0;
  bool associationProductValid = false;
  std::string source;
};
struct RecoAssociationEdge {
  std::size_t sourceIndex = 0;
  std::size_t targetIndex = 0;
  std::string sourceCollection;
  std::string targetCollection;
  std::string associationType;
};
class RecoAssociationAlg {
public:
  static AssociationSummary Summarize(std::size_t, bool,
                                      std::string source = {});
  template <typename T>
  static T const *UniqueOrNull(std::vector<T const *> const &objects,
                               bool valid) {
    return valid && objects.size() == 1 ? objects.front() : nullptr;
  }
};
} // namespace pdhd::diagnostics
#endif
