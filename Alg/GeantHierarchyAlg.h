/**
 * @file GeantHierarchyAlg.h
 * @brief Builds deterministic Geant particle ancestry from event-local facts.
 *
 * ParticleInventoryService supplies MCParticle/MCTruth facts to the analyzer;
 * this pure graph layer detects missing parents, duplicate IDs, cycles, and
 * truncated chains without owning service or event pointers.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_GEANTHIERARCHYALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_GEANTHIERARCHYALG_H
#include <cstddef>
#include <map>
#include <string>
#include <vector>
namespace pdhd::diagnostics {
struct GeantParticleFact {
  int trackId = 0, motherTrackId = 0, pdg = 0;
  std::string process, endProcess;
};
enum class AncestryState {
  Complete,
  MissingParticle,
  MissingMother,
  DuplicateTrackId,
  Cycle,
  DepthLimit
};
struct AncestryResult {
  AncestryState state = AncestryState::MissingParticle;
  std::vector<int>
      trackIds; // requested particle first, oldest resolved ancestor last
};
struct GeantHierarchy {
  std::map<int, GeantParticleFact> particles;
  std::map<int, std::vector<int>> daughters;
  std::vector<int> duplicateTrackIds;
};
class GeantHierarchyAlg {
public:
  static GeantHierarchy Build(std::vector<GeantParticleFact> const &);
  static AncestryResult Ancestors(GeantHierarchy const &, int trackId,
                                  std::size_t maximumDepth = 1000);
  static bool HasAncestor(GeantHierarchy const &, int trackId,
                          int ancestorTrackId, std::size_t maximumDepth = 1000);
};
} // namespace pdhd::diagnostics
#endif
