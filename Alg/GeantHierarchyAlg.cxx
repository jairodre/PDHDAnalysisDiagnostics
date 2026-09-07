/** @file GeantHierarchyAlg.cxx @brief Implements guarded event-local Geant
 * ancestry traversal. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/GeantHierarchyAlg.h"
#include <algorithm>
#include <set>
namespace pdhd::diagnostics {
GeantHierarchy
GeantHierarchyAlg::Build(std::vector<GeantParticleFact> const &input) {
  GeantHierarchy h;
  for (auto const &p : input) {
    auto const inserted = h.particles.emplace(p.trackId, p).second;
    if (!inserted)
      h.duplicateTrackIds.push_back(p.trackId);
  }
  std::sort(h.duplicateTrackIds.begin(), h.duplicateTrackIds.end());
  h.duplicateTrackIds.erase(
      std::unique(h.duplicateTrackIds.begin(), h.duplicateTrackIds.end()),
      h.duplicateTrackIds.end());
  for (auto const &[id, p] : h.particles) {
    if (p.motherTrackId != 0 && p.motherTrackId != id)
      h.daughters[p.motherTrackId].push_back(id);
  }
  for (auto &[mother, daughters] : h.daughters)
    std::sort(daughters.begin(), daughters.end());
  return h;
}
AncestryResult GeantHierarchyAlg::Ancestors(GeantHierarchy const &h,
                                            int trackId,
                                            std::size_t maximumDepth) {
  AncestryResult r;
  if (std::binary_search(h.duplicateTrackIds.begin(), h.duplicateTrackIds.end(),
                         trackId)) {
    r.state = AncestryState::DuplicateTrackId;
    return r;
  }
  std::set<int> visited;
  int current = trackId;
  for (std::size_t depth = 0; depth <= maximumDepth; ++depth) {
    auto const it = h.particles.find(current);
    if (it == h.particles.end()) {
      r.state = r.trackIds.empty() ? AncestryState::MissingParticle
                                   : AncestryState::MissingMother;
      return r;
    }
    if (!visited.insert(current).second) {
      r.state = AncestryState::Cycle;
      return r;
    }
    r.trackIds.push_back(current);
    int const mother = it->second.motherTrackId;
    if (mother == 0) {
      r.state = AncestryState::Complete;
      return r;
    }
    if (std::binary_search(h.duplicateTrackIds.begin(),
                           h.duplicateTrackIds.end(), mother)) {
      r.state = AncestryState::DuplicateTrackId;
      return r;
    }
    current = mother;
  }
  r.state = AncestryState::DepthLimit;
  return r;
}
bool GeantHierarchyAlg::HasAncestor(GeantHierarchy const &h, int trackId,
                                    int ancestorTrackId,
                                    std::size_t maximumDepth) {
  auto const result = Ancestors(h, trackId, maximumDepth);
  return std::find(result.trackIds.begin() +
                       std::min<std::size_t>(1, result.trackIds.size()),
                   result.trackIds.end(),
                   ancestorTrackId) != result.trackIds.end();
}
} // namespace pdhd::diagnostics
