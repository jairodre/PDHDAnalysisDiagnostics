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
  // Minimal event-local particle fact; `motherTrackId == 0` marks a primary.
  int trackId = 0;
  // Zero denotes an event-local primary particle.
  int motherTrackId = 0;
  int pdg = 0;
  std::string process;
  std::string endProcess;
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
  // Partial `trackIds` are retained on failure to expose where traversal ended.
  AncestryState state = AncestryState::MissingParticle;
  std::vector<int>
      trackIds; // requested particle first, oldest resolved ancestor last
};
struct GeantHierarchy {
  // Lookup graph built from facts. Duplicate IDs are retained as invalidity
  // evidence rather than silently selecting one particle.
  std::map<int, GeantParticleFact> particles;
  std::map<int, std::vector<int>> daughters;
  std::vector<int> duplicateTrackIds;
};
class GeantHierarchyAlg {
public:
  // Builds stable ID-indexed parent/daughter maps; duplicate input IDs are not
  // overwritten for ancestry purposes.
  static GeantHierarchy Build(std::vector<GeantParticleFact> const &);
  // Traverses child to primary, returning a terminal state for missing,
  // duplicate, cyclic, and depth-limited histories.
  static AncestryResult Ancestors(GeantHierarchy const &, int trackId,
                                  std::size_t maximumDepth = 1000);
  // Tests resolved ancestors only: the requested track is not its own ancestor.
  static bool HasAncestor(GeantHierarchy const &, int trackId,
                          int ancestorTrackId, std::size_t maximumDepth = 1000);
};
} // namespace pdhd::diagnostics
#endif
