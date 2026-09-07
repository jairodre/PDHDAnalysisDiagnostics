/**
 * @file RecoTruthMatchAlg.h
 * @brief Computes auditable ranks, purity, and completeness from truth
 * contributions.
 *
 * BackTrackerMatchingData or ProtoDUNETruthUtils supplies contribution
 * evidence; this class performs deterministic, source-independent metric
 * calculation.
 */
#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_RECOTRUTHMATCHALG_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_RECOTRUTHMATCHALG_H
#include <string>
#include <vector>
namespace pdhd::diagnostics {
enum class MatchEvidence { HitCount, Charge, Energy };
struct TruthContribution {
  int truthTrackId = 0;
  double sharedEvidence = 0.;
  double recoTotalEvidence = 0.;
  double truthTotalEvidence = 0.;
  std::size_t sharedHits = 0;
  std::size_t sharedDeltaRayHits = 0;
};
struct TruthMatch {
  int truthTrackId = 0;
  unsigned int rank = 0;
  double sharedEvidence = 0.;
  double purity = 0.;
  double completeness = 0.;
  std::size_t sharedHits = 0, sharedDeltaRayHits = 0;
  bool purityValid = false, completenessValid = false;
};
struct TruthMatchResult {
  std::vector<TruthMatch> matches;
  bool valid = false, ambiguousBest = false;
  MatchEvidence evidence = MatchEvidence::HitCount;
  std::string source;
};
class RecoTruthMatchAlg {
public:
  static TruthMatchResult Rank(std::vector<TruthContribution> const &,
                               MatchEvidence, std::string source,
                               double tieTolerance = 1.e-12);
};
} // namespace pdhd::diagnostics
#endif
