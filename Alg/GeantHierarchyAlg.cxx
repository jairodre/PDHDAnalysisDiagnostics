/** @file GeantHierarchyAlg.cxx @brief Implements guarded event-local Geant
 * ancestry traversal. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/GeantHierarchyAlg.h"
#include <algorithm>
#include <set>
namespace pdhd::diagnostics {
GeantHierarchy
GeantHierarchyAlg::Build(std::vector<GeantParticleFact> const &input) {
  GeantHierarchy hierarchy;
  for (GeantParticleFact const &particle : input) {
    bool const inserted =
        hierarchy.particles.emplace(particle.trackId, particle).second;
    if (!inserted) {
      hierarchy.duplicateTrackIds.push_back(particle.trackId);
    }
  }
  std::sort(hierarchy.duplicateTrackIds.begin(),
            hierarchy.duplicateTrackIds.end());
  hierarchy.duplicateTrackIds.erase(
      std::unique(hierarchy.duplicateTrackIds.begin(),
                  hierarchy.duplicateTrackIds.end()),
      hierarchy.duplicateTrackIds.end());

  for (auto const &[trackId, particle] : hierarchy.particles) {
    if (particle.motherTrackId != 0 && particle.motherTrackId != trackId) {
      hierarchy.daughters[particle.motherTrackId].push_back(trackId);
    }
  }
  for (auto &motherAndDaughters : hierarchy.daughters) {
    std::sort(motherAndDaughters.second.begin(),
              motherAndDaughters.second.end());
  }
  return hierarchy;
}
AncestryResult GeantHierarchyAlg::Ancestors(GeantHierarchy const &h,
                                            int trackId,
                                            std::size_t maximumDepth) {
  AncestryResult result;
  if (std::binary_search(h.duplicateTrackIds.begin(), h.duplicateTrackIds.end(),
                         trackId)) {
    result.state = AncestryState::DuplicateTrackId;
    return result;
  }

  std::set<int> visited;
  int currentTrackId = trackId;
  for (std::size_t depth = 0; depth <= maximumDepth; ++depth) {
    auto const particle = h.particles.find(currentTrackId);
    if (particle == h.particles.end()) {
      result.state = result.trackIds.empty() ? AncestryState::MissingParticle
                                             : AncestryState::MissingMother;
      return result;
    }
    if (!visited.insert(currentTrackId).second) {
      result.state = AncestryState::Cycle;
      return result;
    }
    result.trackIds.push_back(currentTrackId);

    int const motherTrackId = particle->second.motherTrackId;
    if (motherTrackId == 0) {
      result.state = AncestryState::Complete;
      return result;
    }
    if (std::binary_search(h.duplicateTrackIds.begin(),
                           h.duplicateTrackIds.end(), motherTrackId)) {
      result.state = AncestryState::DuplicateTrackId;
      return result;
    }
    currentTrackId = motherTrackId;
  }
  result.state = AncestryState::DepthLimit;
  return result;
}
bool GeantHierarchyAlg::HasAncestor(GeantHierarchy const &h, int trackId,
                                    int ancestorTrackId,
                                    std::size_t maximumDepth) {
  auto const result = Ancestors(h, trackId, maximumDepth);
  auto const firstAncestor = result.trackIds.begin() +
                             std::min<std::size_t>(1, result.trackIds.size());
  return std::find(firstAncestor, result.trackIds.end(), ancestorTrackId) !=
         result.trackIds.end();
}
} // namespace pdhd::diagnostics
