/** @file RecoTruthMatchAlg.cxx @brief Implements deterministic truth-match
 * metrics. */
#include "protoduneana/PDHDAnalysisDiagnostics/Alg/RecoTruthMatchAlg.h"
#include <algorithm>
#include <cmath>
#include <utility>
namespace pdhd::diagnostics {
TruthMatchResult
RecoTruthMatchAlg::Rank(std::vector<TruthContribution> const &input,
                        MatchEvidence evidence, std::string source,
                        double tieTolerance) {
  TruthMatchResult r;
  r.evidence = evidence;
  r.source = std::move(source);
  if (!std::isfinite(tieTolerance) || tieTolerance < 0.)
    return r;
  for (auto const &c : input) {
    if (!std::isfinite(c.sharedEvidence) || c.sharedEvidence < 0.)
      continue;
    TruthMatch m;
    m.truthTrackId = c.truthTrackId;
    m.sharedEvidence = c.sharedEvidence;
    m.sharedHits = c.sharedHits;
    m.sharedDeltaRayHits = c.sharedDeltaRayHits;
    if (std::isfinite(c.recoTotalEvidence) && c.recoTotalEvidence > 0.) {
      m.purity = c.sharedEvidence / c.recoTotalEvidence;
      m.purityValid = true;
    }
    if (std::isfinite(c.truthTotalEvidence) && c.truthTotalEvidence > 0.) {
      m.completeness = c.sharedEvidence / c.truthTotalEvidence;
      m.completenessValid = true;
    }
    r.matches.push_back(m);
  }
  std::stable_sort(r.matches.begin(), r.matches.end(),
                   [](auto const &a, auto const &b) {
                     if (a.sharedEvidence != b.sharedEvidence)
                       return a.sharedEvidence > b.sharedEvidence;
                     return a.truthTrackId < b.truthTrackId;
                   });
  for (std::size_t i = 0; i < r.matches.size(); ++i)
    r.matches[i].rank = i;
  r.valid = !r.matches.empty();
  r.ambiguousBest = r.matches.size() > 1 &&
                    std::abs(r.matches[0].sharedEvidence -
                             r.matches[1].sharedEvidence) <= tieTolerance;
  return r;
}
} // namespace pdhd::diagnostics
