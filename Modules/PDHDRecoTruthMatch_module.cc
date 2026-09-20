/**
 * @file PDHDRecoTruthMatch_module.cc
 * @brief Writes one event-level tree relating every reconstructed track to MC truth.
 *
 * Inputs are the configured Pandora PFParticles, slices, tracks, and hits.  In
 * MC, ProtoDUNETruthUtils supplies energy- and hit-ranked contributors, while
 * ParticleInventoryService resolves the configured Geant hierarchy.  The tree
 * retains all reconstructed tracks, marks membership in Pandora's unique
 * IsTestBeam slice, and stores flattened reconstructed/Geant trajectories.
 * Data events keep the reconstructed inventory and explicitly mark every
 * truth-only quantity unavailable; no truth service is called for data.
 */

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

template <class T>
void AddBranch(TTree &tree, char const *name, T &value) {
  // ROOT keeps a pointer to value.  Before each tree_->Fill(), the module
  // updates that same scalar or vector for the current event.
  tree.Branch(name, std::addressof(value));
}

unsigned int ToUnsigned(std::size_t value) {
  // ROOT branch types use unsigned int for stored counts.  Saturate rather
  // than silently narrowing a larger C++ container size.
  return static_cast<unsigned int>(std::min(
      value, static_cast<std::size_t>(std::numeric_limits<unsigned int>::max())));
}

} // namespace

namespace pdhd::diagnostics {

class PDHDRecoTruthMatch : public art::EDAnalyzer {
public:
  explicit PDHDRecoTruthMatch(fhicl::ParameterSet const &);

private:
  void beginJob() override;
  void analyze(art::Event const &) override;

  void ResetEventState(art::Event const &);
  // Restores every truth-only field for the most recently appended reco track.
  // This keeps scalar facts consistent with rolled-back flattened ranges.
  void ResetLatestTrackTruth();
  void AppendRecoTrackTrajectory(recob::Track const &);
  int AppendGeantHierarchy(int mcTrackId, int parentIndex,
                           unsigned int depth,
                           std::unordered_set<int> &activePath,
                           bool &truncated,
                           cheat::ParticleInventoryService &inventory);

  // FHiCL-configured labels for the event products and associations read in
  // analyze().
  art::InputTag pfpTag_;
  art::InputTag pfpMetadataTag_;
  art::InputTag pfpSliceTag_;
  art::InputTag trackTag_;
  art::InputTag hitTag_;
  art::InputTag mcParticleTag_;
  bool enableMCTruthMatching_;
  unsigned int maximumHierarchyDepth_;
  unsigned int maxEventMessages_;
  unsigned int eventMessagesEmitted_ = 0;

  // ROOT output object and event-level branches.
  TTree *tree_ = nullptr;

  unsigned int run_ = 0;
  unsigned int subrun_ = 0;
  unsigned int event_ = 0;
  bool isData_ = false;
  bool pfpProductAvailable_ = false;
  bool trackProductAvailable_ = false;
  bool pfpMetadataAssociationAvailable_ = false;
  bool pfpSliceAssociationAvailable_ = false;
  bool pfpTrackAssociationAvailable_ = false;
  bool pandoraBeamSliceFound_ = false;
  bool pandoraBeamSliceAmbiguous_ = false;
  int pandoraBeamSliceId_ = -1;
  unsigned int pandoraBeamSliceCount_ = 0;
  bool mcParticleProductAvailable_ = false;
  bool truthHitProductAvailable_ = false;
  bool geantHierarchyDuplicateFree_ = false;
  unsigned int geantHierarchyDuplicateTrackIdCount_ = 0;
  bool mcTruthMatchingEvaluated_ = false;
  std::string energyMatchMethod_;
  std::string hitMatchMethod_;

  // One element per reco track.  All vectors in this group share the same
  // reco-track index, so e.g. reco_track_id[i] and reco_track_in_selected_beam_slice[i]
  // describe the same recob::Track.
  std::vector<int> recoTrackCollectionIndex_;
  std::vector<int> recoTrackId_;
  std::vector<unsigned char> recoTrackInSelectedBeamSlice_;
  std::vector<unsigned char> recoTrackPfpAssociationValid_;
  std::vector<unsigned int> recoTrackPfpAssociationCount_;
  std::vector<unsigned int> recoTrackSelectedBeamSliceMembershipCount_;
  std::vector<unsigned long long> recoTrackTrajectoryOffsets_;
  std::vector<double> recoTrackTrajectoryXcm_;
  std::vector<double> recoTrackTrajectoryYcm_;
  std::vector<double> recoTrackTrajectoryZcm_;

  // Leading MC match and scalar truth properties for each reco track.
  std::vector<unsigned char> recoTrackEnergyMatchValid_;
  std::vector<unsigned char> recoTrackEnergyMatchAmbiguous_;
  std::vector<int> recoTrackEnergyBestMcTrackId_;
  std::vector<int> recoTrackEnergyBestPdg_;
  std::vector<unsigned char> recoTrackEnergyBestTrajectoryValid_;
  std::vector<std::string> recoTrackEnergyBestProcess_;
  std::vector<std::string> recoTrackEnergyBestEndProcess_;
  std::vector<double> recoTrackEnergyBestStartXcm_;
  std::vector<double> recoTrackEnergyBestStartYcm_;
  std::vector<double> recoTrackEnergyBestStartZcm_;
  std::vector<double> recoTrackEnergyBestEndXcm_;
  std::vector<double> recoTrackEnergyBestEndYcm_;
  std::vector<double> recoTrackEnergyBestEndZcm_;
  std::vector<double> recoTrackEnergyBestInitialPxGeV_;
  std::vector<double> recoTrackEnergyBestInitialPyGeV_;
  std::vector<double> recoTrackEnergyBestInitialPzGeV_;
  std::vector<int> recoTrackEnergyBestMotherMcTrackId_;
  std::vector<unsigned char> recoTrackEnergyBestMotherValid_;
  std::vector<int> recoTrackEnergyBestMotherPdg_;
  std::vector<unsigned int> recoTrackEnergyBestDaughterCount_;
  std::vector<unsigned char> recoTrackHitMatchValid_;
  std::vector<unsigned char> recoTrackHitMatchAmbiguous_;
  std::vector<int> recoTrackHitBestMcTrackId_;
  std::vector<int> recoTrackHitBestPdg_;
  std::vector<unsigned int> recoTrackHitBestSharedHits_;
  std::vector<unsigned int> recoTrackHitBestSharedDeltaRayHits_;
  std::vector<unsigned char> recoTrackPurityValid_;
  std::vector<double> recoTrackPurity_;
  std::vector<unsigned char> recoTrackCompletenessValid_;
  std::vector<double> recoTrackCompleteness_;
  std::vector<unsigned char> recoTrackTruthOriginValid_;
  std::vector<int> recoTrackTruthOrigin_;

  // Flattened all-contributor lists.  Offsets map a reco-track index to its
  // range: [offsets[i], offsets[i + 1]).
  std::vector<unsigned long long> recoTrackEnergyMatchOffsets_;
  std::vector<int> energyMatchMcTrackId_;
  std::vector<int> energyMatchPdg_;
  std::vector<double> energyMatchEnergyFraction_;
  std::vector<unsigned int> energyMatchRank_;

  std::vector<unsigned long long> recoTrackHitMatchOffsets_;
  std::vector<int> hitMatchMcTrackId_;
  std::vector<int> hitMatchPdg_;
  std::vector<unsigned int> hitMatchSharedHits_;
  std::vector<unsigned int> hitMatchSharedDeltaRayHits_;
  std::vector<unsigned int> hitMatchRank_;

  // Flattened Geant family trees rooted at the leading energy-matched
  // MCParticle for every reco track.
  std::vector<unsigned char> recoTrackHierarchyValid_;
  std::vector<unsigned char> recoTrackHierarchyTruncated_;
  std::vector<int> recoTrackHierarchyRootMcTrackId_;
  std::vector<unsigned long long> recoTrackHierarchyOffsets_;
  std::vector<int> hierarchyRecoTrackIndex_;
  std::vector<int> hierarchyMcTrackId_;
  std::vector<int> hierarchyParentIndex_;
  std::vector<int> hierarchyMotherMcTrackId_;
  std::vector<int> hierarchyPdg_;
  std::vector<unsigned int> hierarchyDepth_;
  std::vector<std::string> hierarchyProcess_;
  std::vector<std::string> hierarchyEndProcess_;
  std::vector<unsigned long long> hierarchyDaughterBegin_;
  std::vector<unsigned int> hierarchyDaughterCount_;
  std::vector<int> hierarchyDaughterIndices_;
  std::vector<unsigned long long> hierarchyTrajectoryOffsets_;
  std::vector<double> hierarchyTrajectoryXcm_;
  std::vector<double> hierarchyTrajectoryYcm_;
  std::vector<double> hierarchyTrajectoryZcm_;
  std::vector<double> hierarchyTrajectoryPxGeV_;
  std::vector<double> hierarchyTrajectoryPyGeV_;
  std::vector<double> hierarchyTrajectoryPzGeV_;
};

PDHDRecoTruthMatch::PDHDRecoTruthMatch(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      pfpTag_(p.get<art::InputTag>("PFParticleTag", art::InputTag{"pandora"})),
      pfpMetadataTag_(p.get<art::InputTag>("PFParticleMetadataTag", pfpTag_)),
      pfpSliceTag_(p.get<art::InputTag>("PFParticleSliceTag", pfpTag_)),
      trackTag_(p.get<art::InputTag>("TrackTag", art::InputTag{"pandoraTrack"})),
      hitTag_(p.get<art::InputTag>("HitTag", art::InputTag{"gaushit"})),
      mcParticleTag_(p.get<art::InputTag>("MCParticleTag", art::InputTag{"largeant"})),
      enableMCTruthMatching_(p.get<bool>("EnableMCTruthMatching", false)),
      maximumHierarchyDepth_(p.get<unsigned int>("MaximumHierarchyDepth", 1000U)),
      maxEventMessages_(p.get<unsigned int>("MaxEventMessages", 3U)) {
  // Declare the products this EDAnalyzer may read.  The InputTags above say
  // which producer label to use; consumes<T>() declares the C++ product type.
  consumes<std::vector<recob::PFParticle>>(pfpTag_);
  consumes<std::vector<recob::Track>>(trackTag_);
  consumes<std::vector<recob::Hit>>(hitTag_);
  consumes<std::vector<simb::MCParticle>>(mcParticleTag_);
}

void PDHDRecoTruthMatch::beginJob() {
  auto fileService = art::ServiceHandle<art::TFileService>();
  tree_ = fileService->make<TTree>("RecoTruthMatch",
                                   "reco tracks, truth matches, and Geant hierarchies");
  // TFileService writes the completed tree at job end.  Disable ROOT's
  // intermediate autosaves so the output has one RecoTruthMatch key.
  tree_->SetAutoSave(0);

  // Event provenance and product/association availability.
  AddBranch(*tree_, "run", run_);
  AddBranch(*tree_, "subrun", subrun_);
  AddBranch(*tree_, "event", event_);
  AddBranch(*tree_, "is_data", isData_);
  AddBranch(*tree_, "pfp_product_available", pfpProductAvailable_);
  AddBranch(*tree_, "track_product_available", trackProductAvailable_);
  AddBranch(*tree_, "pfp_metadata_association_available", pfpMetadataAssociationAvailable_);
  AddBranch(*tree_, "pfp_slice_association_available", pfpSliceAssociationAvailable_);
  AddBranch(*tree_, "pfp_track_association_available", pfpTrackAssociationAvailable_);
  AddBranch(*tree_, "pandora_beam_slice_found", pandoraBeamSliceFound_);
  AddBranch(*tree_, "pandora_beam_slice_ambiguous", pandoraBeamSliceAmbiguous_);
  AddBranch(*tree_, "pandora_beam_slice_id", pandoraBeamSliceId_);
  AddBranch(*tree_, "pandora_beam_slice_count", pandoraBeamSliceCount_);
  AddBranch(*tree_, "mcparticle_product_available", mcParticleProductAvailable_);
  AddBranch(*tree_, "truth_hit_product_available", truthHitProductAvailable_);
  AddBranch(*tree_, "geant_hierarchy_duplicate_free", geantHierarchyDuplicateFree_);
  AddBranch(*tree_, "geant_hierarchy_duplicate_track_id_count",
            geantHierarchyDuplicateTrackIdCount_);
  AddBranch(*tree_, "mc_truth_matching_evaluated", mcTruthMatchingEvaluated_);
  AddBranch(*tree_, "energy_match_method", energyMatchMethod_);
  AddBranch(*tree_, "hit_match_method", hitMatchMethod_);

  // Reconstructed-track inventory and Pandora beam-slice membership.
  AddBranch(*tree_, "reco_track_collection_index", recoTrackCollectionIndex_);
  AddBranch(*tree_, "reco_track_id", recoTrackId_);
  AddBranch(*tree_, "reco_track_in_selected_beam_slice", recoTrackInSelectedBeamSlice_);
  AddBranch(*tree_, "reco_track_pfp_association_valid", recoTrackPfpAssociationValid_);
  AddBranch(*tree_, "reco_track_pfp_association_count", recoTrackPfpAssociationCount_);
  AddBranch(*tree_, "reco_track_selected_beam_slice_membership_count",
            recoTrackSelectedBeamSliceMembershipCount_);
  AddBranch(*tree_, "reco_track_trajectory_offsets", recoTrackTrajectoryOffsets_);
  AddBranch(*tree_, "reco_track_trajectory_x_cm", recoTrackTrajectoryXcm_);
  AddBranch(*tree_, "reco_track_trajectory_y_cm", recoTrackTrajectoryYcm_);
  AddBranch(*tree_, "reco_track_trajectory_z_cm", recoTrackTrajectoryZcm_);

  // Leading MC matches and derived truth quantities, one entry per reco track.
  AddBranch(*tree_, "reco_track_energy_match_valid", recoTrackEnergyMatchValid_);
  AddBranch(*tree_, "reco_track_energy_match_ambiguous", recoTrackEnergyMatchAmbiguous_);
  AddBranch(*tree_, "reco_track_energy_best_mc_track_id", recoTrackEnergyBestMcTrackId_);
  AddBranch(*tree_, "reco_track_energy_best_pdg", recoTrackEnergyBestPdg_);
  AddBranch(*tree_, "reco_track_energy_best_trajectory_valid",
            recoTrackEnergyBestTrajectoryValid_);
  AddBranch(*tree_, "reco_track_energy_best_process", recoTrackEnergyBestProcess_);
  AddBranch(*tree_, "reco_track_energy_best_end_process", recoTrackEnergyBestEndProcess_);
  AddBranch(*tree_, "reco_track_energy_best_start_x_cm", recoTrackEnergyBestStartXcm_);
  AddBranch(*tree_, "reco_track_energy_best_start_y_cm", recoTrackEnergyBestStartYcm_);
  AddBranch(*tree_, "reco_track_energy_best_start_z_cm", recoTrackEnergyBestStartZcm_);
  AddBranch(*tree_, "reco_track_energy_best_end_x_cm", recoTrackEnergyBestEndXcm_);
  AddBranch(*tree_, "reco_track_energy_best_end_y_cm", recoTrackEnergyBestEndYcm_);
  AddBranch(*tree_, "reco_track_energy_best_end_z_cm", recoTrackEnergyBestEndZcm_);
  AddBranch(*tree_, "reco_track_energy_best_initial_px_GeV",
            recoTrackEnergyBestInitialPxGeV_);
  AddBranch(*tree_, "reco_track_energy_best_initial_py_GeV",
            recoTrackEnergyBestInitialPyGeV_);
  AddBranch(*tree_, "reco_track_energy_best_initial_pz_GeV",
            recoTrackEnergyBestInitialPzGeV_);
  AddBranch(*tree_, "reco_track_energy_best_mother_mc_track_id",
            recoTrackEnergyBestMotherMcTrackId_);
  AddBranch(*tree_, "reco_track_energy_best_mother_valid",
            recoTrackEnergyBestMotherValid_);
  AddBranch(*tree_, "reco_track_energy_best_mother_pdg",
            recoTrackEnergyBestMotherPdg_);
  AddBranch(*tree_, "reco_track_energy_best_daughter_count",
            recoTrackEnergyBestDaughterCount_);
  AddBranch(*tree_, "reco_track_hit_match_valid", recoTrackHitMatchValid_);
  AddBranch(*tree_, "reco_track_hit_match_ambiguous", recoTrackHitMatchAmbiguous_);
  AddBranch(*tree_, "reco_track_hit_best_mc_track_id", recoTrackHitBestMcTrackId_);
  AddBranch(*tree_, "reco_track_hit_best_pdg", recoTrackHitBestPdg_);
  AddBranch(*tree_, "reco_track_hit_best_shared_hits", recoTrackHitBestSharedHits_);
  AddBranch(*tree_, "reco_track_hit_best_shared_delta_ray_hits",
            recoTrackHitBestSharedDeltaRayHits_);
  AddBranch(*tree_, "reco_track_purity_valid", recoTrackPurityValid_);
  AddBranch(*tree_, "reco_track_purity", recoTrackPurity_);
  AddBranch(*tree_, "reco_track_completeness_valid", recoTrackCompletenessValid_);
  AddBranch(*tree_, "reco_track_completeness", recoTrackCompleteness_);
  AddBranch(*tree_, "reco_track_truth_origin_valid", recoTrackTruthOriginValid_);
  AddBranch(*tree_, "reco_track_truth_origin", recoTrackTruthOrigin_);

  // Complete energy- and hit-ranked MC contributor lists, flattened across
  // reco tracks and addressed by their offset branches.
  AddBranch(*tree_, "reco_track_energy_match_offsets", recoTrackEnergyMatchOffsets_);
  AddBranch(*tree_, "energy_match_mc_track_id", energyMatchMcTrackId_);
  AddBranch(*tree_, "energy_match_pdg", energyMatchPdg_);
  AddBranch(*tree_, "energy_match_energy_fraction", energyMatchEnergyFraction_);
  AddBranch(*tree_, "energy_match_rank", energyMatchRank_);
  AddBranch(*tree_, "reco_track_hit_match_offsets", recoTrackHitMatchOffsets_);
  AddBranch(*tree_, "hit_match_mc_track_id", hitMatchMcTrackId_);
  AddBranch(*tree_, "hit_match_pdg", hitMatchPdg_);
  AddBranch(*tree_, "hit_match_shared_hits", hitMatchSharedHits_);
  AddBranch(*tree_, "hit_match_shared_delta_ray_hits", hitMatchSharedDeltaRayHits_);
  AddBranch(*tree_, "hit_match_rank", hitMatchRank_);

  // Recursive Geant MCParticle hierarchy and its flattened trajectories.
  AddBranch(*tree_, "reco_track_hierarchy_valid", recoTrackHierarchyValid_);
  AddBranch(*tree_, "reco_track_hierarchy_truncated", recoTrackHierarchyTruncated_);
  AddBranch(*tree_, "reco_track_hierarchy_root_mc_track_id", recoTrackHierarchyRootMcTrackId_);
  AddBranch(*tree_, "reco_track_hierarchy_offsets", recoTrackHierarchyOffsets_);
  AddBranch(*tree_, "hierarchy_reco_track_index", hierarchyRecoTrackIndex_);
  AddBranch(*tree_, "hierarchy_mc_track_id", hierarchyMcTrackId_);
  AddBranch(*tree_, "hierarchy_parent_index", hierarchyParentIndex_);
  AddBranch(*tree_, "hierarchy_mother_mc_track_id", hierarchyMotherMcTrackId_);
  AddBranch(*tree_, "hierarchy_pdg", hierarchyPdg_);
  AddBranch(*tree_, "hierarchy_depth", hierarchyDepth_);
  AddBranch(*tree_, "hierarchy_process", hierarchyProcess_);
  AddBranch(*tree_, "hierarchy_end_process", hierarchyEndProcess_);
  AddBranch(*tree_, "hierarchy_daughter_begin", hierarchyDaughterBegin_);
  AddBranch(*tree_, "hierarchy_daughter_count", hierarchyDaughterCount_);
  AddBranch(*tree_, "hierarchy_daughter_indices", hierarchyDaughterIndices_);
  AddBranch(*tree_, "hierarchy_trajectory_offsets", hierarchyTrajectoryOffsets_);
  AddBranch(*tree_, "hierarchy_trajectory_x_cm", hierarchyTrajectoryXcm_);
  AddBranch(*tree_, "hierarchy_trajectory_y_cm", hierarchyTrajectoryYcm_);
  AddBranch(*tree_, "hierarchy_trajectory_z_cm", hierarchyTrajectoryZcm_);
  AddBranch(*tree_, "hierarchy_trajectory_px_GeV", hierarchyTrajectoryPxGeV_);
  AddBranch(*tree_, "hierarchy_trajectory_py_GeV", hierarchyTrajectoryPyGeV_);
  AddBranch(*tree_, "hierarchy_trajectory_pz_GeV", hierarchyTrajectoryPzGeV_);
}

void PDHDRecoTruthMatch::ResetEventState(art::Event const &evt) {
  // Event identity and availability flags.
  run_ = evt.run();
  subrun_ = evt.subRun();
  event_ = evt.id().event();
  isData_ = evt.isRealData();
  eventMessagesEmitted_ = 0;
  pfpProductAvailable_ = false;
  trackProductAvailable_ = false;
  pfpMetadataAssociationAvailable_ = false;
  pfpSliceAssociationAvailable_ = false;
  pfpTrackAssociationAvailable_ = false;
  pandoraBeamSliceFound_ = false;
  pandoraBeamSliceAmbiguous_ = false;
  pandoraBeamSliceId_ = -1;
  pandoraBeamSliceCount_ = 0;
  mcParticleProductAvailable_ = false;
  truthHitProductAvailable_ = false;
  geantHierarchyDuplicateFree_ = false;
  geantHierarchyDuplicateTrackIdCount_ = 0;
  mcTruthMatchingEvaluated_ = false;
  energyMatchMethod_ = "ProtoDUNETruthUtils::GetMCParticleListFromRecoTrack";
  hitMatchMethod_ = "ProtoDUNETruthUtils::GetMCParticleListByHits";

  // Reco-track vectors and their flattened trajectory.
  recoTrackCollectionIndex_.clear();
  recoTrackId_.clear();
  recoTrackInSelectedBeamSlice_.clear();
  recoTrackPfpAssociationValid_.clear();
  recoTrackPfpAssociationCount_.clear();
  recoTrackSelectedBeamSliceMembershipCount_.clear();
  recoTrackTrajectoryOffsets_.assign(1, 0ULL);
  recoTrackTrajectoryXcm_.clear();
  recoTrackTrajectoryYcm_.clear();
  recoTrackTrajectoryZcm_.clear();

  // Per-reco-track truth defaults.  A false validity flag and NaN numerical
  // value distinguish unavailable MC information from a physical zero.
  recoTrackEnergyMatchValid_.clear();
  recoTrackEnergyMatchAmbiguous_.clear();
  recoTrackEnergyBestMcTrackId_.clear();
  recoTrackEnergyBestPdg_.clear();
  recoTrackEnergyBestTrajectoryValid_.clear();
  recoTrackEnergyBestProcess_.clear();
  recoTrackEnergyBestEndProcess_.clear();
  recoTrackEnergyBestStartXcm_.clear();
  recoTrackEnergyBestStartYcm_.clear();
  recoTrackEnergyBestStartZcm_.clear();
  recoTrackEnergyBestEndXcm_.clear();
  recoTrackEnergyBestEndYcm_.clear();
  recoTrackEnergyBestEndZcm_.clear();
  recoTrackEnergyBestInitialPxGeV_.clear();
  recoTrackEnergyBestInitialPyGeV_.clear();
  recoTrackEnergyBestInitialPzGeV_.clear();
  recoTrackEnergyBestMotherMcTrackId_.clear();
  recoTrackEnergyBestMotherValid_.clear();
  recoTrackEnergyBestMotherPdg_.clear();
  recoTrackEnergyBestDaughterCount_.clear();
  recoTrackHitMatchValid_.clear();
  recoTrackHitMatchAmbiguous_.clear();
  recoTrackHitBestMcTrackId_.clear();
  recoTrackHitBestPdg_.clear();
  recoTrackHitBestSharedHits_.clear();
  recoTrackHitBestSharedDeltaRayHits_.clear();
  recoTrackPurityValid_.clear();
  recoTrackPurity_.clear();
  recoTrackCompletenessValid_.clear();
  recoTrackCompleteness_.clear();
  recoTrackTruthOriginValid_.clear();
  recoTrackTruthOrigin_.clear();

  // All-contributor lists and their per-track offsets.
  recoTrackEnergyMatchOffsets_.assign(1, 0ULL);
  energyMatchMcTrackId_.clear();
  energyMatchPdg_.clear();
  energyMatchEnergyFraction_.clear();
  energyMatchRank_.clear();
  recoTrackHitMatchOffsets_.assign(1, 0ULL);
  hitMatchMcTrackId_.clear();
  hitMatchPdg_.clear();
  hitMatchSharedHits_.clear();
  hitMatchSharedDeltaRayHits_.clear();
  hitMatchRank_.clear();

  // Flattened Geant hierarchy and its per-particle trajectory.
  recoTrackHierarchyValid_.clear();
  recoTrackHierarchyTruncated_.clear();
  recoTrackHierarchyRootMcTrackId_.clear();
  recoTrackHierarchyOffsets_.assign(1, 0ULL);
  hierarchyRecoTrackIndex_.clear();
  hierarchyMcTrackId_.clear();
  hierarchyParentIndex_.clear();
  hierarchyMotherMcTrackId_.clear();
  hierarchyPdg_.clear();
  hierarchyDepth_.clear();
  hierarchyProcess_.clear();
  hierarchyEndProcess_.clear();
  hierarchyDaughterBegin_.clear();
  hierarchyDaughterCount_.clear();
  hierarchyDaughterIndices_.clear();
  hierarchyTrajectoryOffsets_.assign(1, 0ULL);
  hierarchyTrajectoryXcm_.clear();
  hierarchyTrajectoryYcm_.clear();
  hierarchyTrajectoryZcm_.clear();
  hierarchyTrajectoryPxGeV_.clear();
  hierarchyTrajectoryPyGeV_.clear();
  hierarchyTrajectoryPzGeV_.clear();
}

void PDHDRecoTruthMatch::ResetLatestTrackTruth() {
  recoTrackEnergyMatchValid_.back() = 0U;
  recoTrackEnergyMatchAmbiguous_.back() = 0U;
  recoTrackEnergyBestMcTrackId_.back() = 0;
  recoTrackEnergyBestPdg_.back() = 0;
  recoTrackEnergyBestTrajectoryValid_.back() = 0U;
  recoTrackEnergyBestProcess_.back().clear();
  recoTrackEnergyBestEndProcess_.back().clear();
  double const nan = std::numeric_limits<double>::quiet_NaN();
  recoTrackEnergyBestStartXcm_.back() = nan;
  recoTrackEnergyBestStartYcm_.back() = nan;
  recoTrackEnergyBestStartZcm_.back() = nan;
  recoTrackEnergyBestEndXcm_.back() = nan;
  recoTrackEnergyBestEndYcm_.back() = nan;
  recoTrackEnergyBestEndZcm_.back() = nan;
  recoTrackEnergyBestInitialPxGeV_.back() = nan;
  recoTrackEnergyBestInitialPyGeV_.back() = nan;
  recoTrackEnergyBestInitialPzGeV_.back() = nan;
  recoTrackEnergyBestMotherMcTrackId_.back() = 0;
  recoTrackEnergyBestMotherValid_.back() = 0U;
  recoTrackEnergyBestMotherPdg_.back() = 0;
  recoTrackEnergyBestDaughterCount_.back() = 0U;
  recoTrackHitMatchValid_.back() = 0U;
  recoTrackHitMatchAmbiguous_.back() = 0U;
  recoTrackHitBestMcTrackId_.back() = 0;
  recoTrackHitBestPdg_.back() = 0;
  recoTrackHitBestSharedHits_.back() = 0U;
  recoTrackHitBestSharedDeltaRayHits_.back() = 0U;
  recoTrackPurityValid_.back() = 0U;
  recoTrackPurity_.back() = nan;
  recoTrackCompletenessValid_.back() = 0U;
  recoTrackCompleteness_.back() = nan;
  recoTrackTruthOriginValid_.back() = 0U;
  recoTrackTruthOrigin_.back() = -1;
  recoTrackHierarchyValid_.back() = 0U;
  recoTrackHierarchyTruncated_.back() = 1U;
  recoTrackHierarchyRootMcTrackId_.back() = 0;
}

void PDHDRecoTruthMatch::AppendRecoTrackTrajectory(recob::Track const &track) {
  // recob::Track can contain invalid fitted points.  Save only valid positions
  // and close this track's flattened range with the next offset.
  for (std::size_t point = 0; point < track.NumberTrajectoryPoints(); ++point) {
    if (!track.HasValidPoint(point)) continue;
    auto const position = track.LocationAtPoint(point);
    recoTrackTrajectoryXcm_.push_back(position.X());
    recoTrackTrajectoryYcm_.push_back(position.Y());
    recoTrackTrajectoryZcm_.push_back(position.Z());
  }
  recoTrackTrajectoryOffsets_.push_back(
      static_cast<unsigned long long>(recoTrackTrajectoryXcm_.size()));
}

int PDHDRecoTruthMatch::AppendGeantHierarchy(
    int const mcTrackId, int const parentIndex, unsigned int const depth,
    std::unordered_set<int> &activePath, bool &truncated,
    cheat::ParticleInventoryService &inventory) {
  if (depth > maximumHierarchyDepth_ || activePath.count(mcTrackId) != 0U) {
    truncated = true;
    return -1;
  }

  // ParticleInventoryService resolves a Geant TrackID directly to the
  // corresponding MCParticle.  The MCParticle itself then provides its PDG,
  // processes, trajectory points, and direct daughter TrackIDs.
  auto const *particle = inventory.TrackIdToParticle_P(mcTrackId);
  if (!particle) {
    truncated = true;
    return -1;
  }

  activePath.insert(mcTrackId);
  // Save this particle before recursion.  Each daughter below records this
  // hierarchyIndex as its parent_index, which reconstructs the family tree in
  // the flat ROOT vectors.
  int const hierarchyIndex = static_cast<int>(hierarchyMcTrackId_.size());
  hierarchyRecoTrackIndex_.push_back(static_cast<int>(recoTrackId_.size() - 1U));
  hierarchyMcTrackId_.push_back(mcTrackId);
  hierarchyParentIndex_.push_back(parentIndex);
  hierarchyMotherMcTrackId_.push_back(particle->Mother());
  hierarchyPdg_.push_back(particle->PdgCode());
  hierarchyDepth_.push_back(depth);
  hierarchyProcess_.push_back(particle->Process());
  hierarchyEndProcess_.push_back(particle->EndProcess());
  hierarchyDaughterBegin_.push_back(0ULL);
  hierarchyDaughterCount_.push_back(0U);

  for (std::size_t point = 0; point < particle->NumberTrajectoryPoints(); ++point) {
    auto const position = particle->Position(point);
    auto const momentum = particle->Momentum(point);
    hierarchyTrajectoryXcm_.push_back(position.X());
    hierarchyTrajectoryYcm_.push_back(position.Y());
    hierarchyTrajectoryZcm_.push_back(position.Z());
    hierarchyTrajectoryPxGeV_.push_back(momentum.Px());
    hierarchyTrajectoryPyGeV_.push_back(momentum.Py());
    hierarchyTrajectoryPzGeV_.push_back(momentum.Pz());
  }
  hierarchyTrajectoryOffsets_.push_back(
      static_cast<unsigned long long>(hierarchyTrajectoryXcm_.size()));

  // MCParticle::Daughter(i) is a Geant TrackID.  Resolve and recurse through
  // every direct daughter, retaining only nodes that were successfully saved.
  std::vector<int> daughterIndices;
  for (int daughter = 0; daughter < particle->NumberDaughters(); ++daughter) {
    int const daughterTrackId = particle->Daughter(daughter);
    int const daughterIndex = AppendGeantHierarchy(
        daughterTrackId, hierarchyIndex, depth + 1U, activePath, truncated,
        inventory);
    if (daughterIndex >= 0) daughterIndices.push_back(daughterIndex);
  }
  // Store the direct-child index list as a second flattened array.
  hierarchyDaughterBegin_.at(static_cast<std::size_t>(hierarchyIndex)) =
      static_cast<unsigned long long>(hierarchyDaughterIndices_.size());
  hierarchyDaughterIndices_.insert(hierarchyDaughterIndices_.end(),
                                   daughterIndices.begin(), daughterIndices.end());
  hierarchyDaughterCount_.at(static_cast<std::size_t>(hierarchyIndex)) =
      ToUnsigned(daughterIndices.size());
  activePath.erase(mcTrackId);
  return hierarchyIndex;
}

void PDHDRecoTruthMatch::analyze(art::Event const &evt) {
  ResetEventState(evt);

  // 1. Read the reconstructed products configured in FHiCL.  A handle is the
  // event-local view of one ART product; it may be invalid when that product
  // was not produced in this input file.
  auto const trackHandle = evt.getHandle<std::vector<recob::Track>>(trackTag_);
  auto const pfpHandle = evt.getHandle<std::vector<recob::PFParticle>>(pfpTag_);
  auto const hitHandle = evt.getHandle<std::vector<recob::Hit>>(hitTag_);
  trackProductAvailable_ = static_cast<bool>(trackHandle);
  pfpProductAvailable_ = static_cast<bool>(pfpHandle);
  truthHitProductAvailable_ = static_cast<bool>(hitHandle);

  // 2. Construct the Pandora associations needed to answer two different
  // questions: which Slice contains a PFP, and which reco Tracks belong to a
  // PFP.  The associations are indexed by the PFParticle collection index.
  std::unique_ptr<art::FindManyP<larpandoraobj::PFParticleMetadata>> metadata;
  std::unique_ptr<art::FindOneP<recob::Slice>> slices;
  std::unique_ptr<art::FindManyP<recob::Track>> pfpTracks;
  if (pfpHandle) {
    try {
      metadata = std::make_unique<art::FindManyP<larpandoraobj::PFParticleMetadata>>(
          pfpHandle, evt, pfpMetadataTag_);
      pfpMetadataAssociationAvailable_ = true;
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_++ < maxEventMessages_)
        mf::LogWarning("PDHDRecoTruthMatch")
            << "PFP metadata association is unavailable for run " << run_
            << ", subrun " << subrun_ << ", event " << event_ << ": "
            << error.what();
    }
    try {
      slices = std::make_unique<art::FindOneP<recob::Slice>>(
          pfpHandle, evt, pfpSliceTag_);
      pfpSliceAssociationAvailable_ = true;
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_++ < maxEventMessages_)
        mf::LogWarning("PDHDRecoTruthMatch")
            << "PFP slice association is unavailable for run " << run_
            << ", subrun " << subrun_ << ", event " << event_ << ": "
            << error.what();
    }
    if (trackHandle) {
      try {
        pfpTracks = std::make_unique<art::FindManyP<recob::Track>>(
            pfpHandle, evt, trackTag_);
        pfpTrackAssociationAvailable_ = true;
      } catch (std::exception const &error) {
        if (eventMessagesEmitted_++ < maxEventMessages_)
          mf::LogWarning("PDHDRecoTruthMatch")
              << "PFP track association is unavailable for run " << run_
              << ", subrun " << subrun_ << ", event " << event_ << ": "
              << error.what();
      }
    }
  }

  // 3. Find the unique Pandora beam slice.  A primary PFP with the
  // IsTestBeam metadata property identifies the selected beam slice.
  std::set<int> beamSliceIds;
  if (pfpHandle && metadata && slices) {
    protoana::ProtoDUNEPFParticleUtils pfpUtils;
    for (std::size_t pfpIndex = 0; pfpIndex < pfpHandle->size(); ++pfpIndex) {
      auto const &pfp = pfpHandle->at(pfpIndex);
      if (!pfp.IsPrimary()) continue;
      try {
        auto const slice = slices->at(pfpIndex);
        // Use the version-matched ProtoDUNE utility for the IsTestBeam
        // convention. The direct associations above remain necessary to
        // retain their availability and the selected Slice ID in this tree.
        bool const isTestBeam =
            pfpUtils.IsBeamParticle(pfp, evt, pfpTag_.encode());
        if (slice && isTestBeam) {
          beamSliceIds.insert(static_cast<int>(slice->ID()));
        }
      } catch (std::exception const &error) {
        if (eventMessagesEmitted_++ < maxEventMessages_)
          mf::LogWarning("PDHDRecoTruthMatch")
              << "Could not read the Slice or metadata for primary PFParticle "
              << pfpIndex << " in run " << run_ << ", subrun " << subrun_
              << ", event " << event_ << ": " << error.what();
      }
    }
  }
  pandoraBeamSliceCount_ = ToUnsigned(beamSliceIds.size());
  pandoraBeamSliceFound_ = pandoraBeamSliceCount_ == 1U;
  pandoraBeamSliceAmbiguous_ = pandoraBeamSliceCount_ > 1U;
  if (pandoraBeamSliceFound_) pandoraBeamSliceId_ = *beamSliceIds.begin();

  // 4. Convert the PFP-to-Track association into per-track information.  The
  // output tree stores every reco track, then records whether it is associated
  // with the unique selected beam slice.
  std::vector<unsigned int> pfpAssociationCounts;
  std::vector<unsigned int> selectedBeamSliceMembershipCounts;
  if (trackHandle) {
    pfpAssociationCounts.assign(trackHandle->size(), 0U);
    selectedBeamSliceMembershipCounts.assign(trackHandle->size(), 0U);
  }
  if (pfpHandle && trackHandle && slices && pfpTracks) {
    // FindManyP returns art::Ptr<recob::Track>.  Map its raw Track pointer
    // back to this track collection's index so the per-track vectors align.
    std::unordered_map<recob::Track const *, std::size_t> trackIndices;
    for (std::size_t index = 0; index < trackHandle->size(); ++index)
      trackIndices.emplace(std::addressof(trackHandle->at(index)), index);

    for (std::size_t pfpIndex = 0; pfpIndex < pfpHandle->size(); ++pfpIndex) {
      try {
        auto const slice = slices->at(pfpIndex);
        bool const inSelectedBeamSlice =
            pandoraBeamSliceFound_ && slice &&
            static_cast<int>(slice->ID()) == pandoraBeamSliceId_;
        for (auto const &track : pfpTracks->at(pfpIndex)) {
          auto const trackIndex = trackIndices.find(track.get());
          if (trackIndex == trackIndices.end()) continue;
          ++pfpAssociationCounts.at(trackIndex->second);
          if (inSelectedBeamSlice)
            ++selectedBeamSliceMembershipCounts.at(trackIndex->second);
        }
      } catch (std::exception const &error) {
        if (eventMessagesEmitted_++ < maxEventMessages_)
          mf::LogWarning("PDHDRecoTruthMatch")
              << "Could not read the Slice or Track associations for PFParticle "
              << pfpIndex << " in run " << run_ << ", subrun " << subrun_
              << ", event " << event_ << ": " << error.what();
      }
    }
  }

  // 5. For MC, always record whether the configured Geant MCParticle product
  // exists. This remains a product-provenance fact even when optional truth
  // matching is disabled. ParticleInventoryService below resolves individual
  // TrackIDs; this product loop also records duplicate TrackID evidence.
  bool mcParticleProductAvailableForHierarchy = false;
  if (!isData_) {
    auto const mcParticleHandle = evt.getHandle<std::vector<simb::MCParticle>>(mcParticleTag_);
    mcParticleProductAvailable_ = static_cast<bool>(mcParticleHandle);
    if (mcParticleHandle) {
      std::set<int> seenTrackIds;
      std::set<int> duplicateTrackIds;
      for (auto const &particle : *mcParticleHandle) {
        if (!seenTrackIds.insert(particle.TrackId()).second)
          duplicateTrackIds.insert(particle.TrackId());
      }
      geantHierarchyDuplicateTrackIdCount_ =
          ToUnsigned(duplicateTrackIds.size());
      geantHierarchyDuplicateFree_ = duplicateTrackIds.empty();
      // Daughter traversal below uses MCParticle::Daughter directly.  A
      // duplicate elsewhere in the product is recorded above, but does not
      // prevent an otherwise resolvable matched particle from being saved.
      mcParticleProductAvailableForHierarchy = enableMCTruthMatching_;
    }
  }

  if (!trackHandle) {
    // Keep one event entry even when the configured track product is absent;
    // its availability branch tells the offline reader why track vectors are empty.
    tree_->Fill();
    return;
  }

  // 6. Prepare the services required by ProtoDUNETruthUtils.  Do not request
  // truth services for data: all truth branches retain their unavailable
  // defaults in that case.
  std::unique_ptr<detinfo::DetectorClocksData> clockData;
  cheat::ParticleInventoryService *inventory = nullptr;
  if (!isData_ && enableMCTruthMatching_) {
    try {
      clockData = std::make_unique<detinfo::DetectorClocksData>(
          art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt));
      art::ServiceHandle<cheat::ParticleInventoryService> inventoryService;
      inventory = inventoryService.operator->();
      mcTruthMatchingEvaluated_ = true;
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_++ < maxEventMessages_)
        mf::LogWarning("PDHDRecoTruthMatch")
            << "Detector clocks or ParticleInventoryService are unavailable for MC truth "
            << "matching in run " << run_ << ", subrun " << subrun_ << ", event "
            << event_ << ": " << error.what();
    }
  }

  // 7. Loop over the complete configured reco Track collection.  This module
  // does not discard non-beam-slice tracks; reco_track_in_selected_beam_slice
  // is the downstream selection flag.
  protoana::ProtoDUNETruthUtils truthUtils;
  for (std::size_t trackIndex = 0; trackIndex < trackHandle->size(); ++trackIndex) {
    auto const &track = trackHandle->at(trackIndex);
    auto const energyBegin = energyMatchMcTrackId_.size();
    auto const hitBegin = hitMatchMcTrackId_.size();
    auto const hierarchyBegin = hierarchyMcTrackId_.size();
    auto const daughterIndexBegin = hierarchyDaughterIndices_.size();
    // 7a. Save the reco-only information first.  These branches are equally
    // meaningful for data and MC.
    recoTrackCollectionIndex_.push_back(static_cast<int>(trackIndex));
    recoTrackId_.push_back(track.ID());
    recoTrackInSelectedBeamSlice_.push_back(
        selectedBeamSliceMembershipCounts.at(trackIndex) > 0U ? 1U : 0U);
    recoTrackPfpAssociationValid_.push_back(
        pfpTrackAssociationAvailable_ ? 1U : 0U);
    recoTrackPfpAssociationCount_.push_back(pfpAssociationCounts.at(trackIndex));
    recoTrackSelectedBeamSliceMembershipCount_.push_back(
        selectedBeamSliceMembershipCounts.at(trackIndex));
    AppendRecoTrackTrajectory(track);

    recoTrackEnergyMatchValid_.push_back(0U);
    recoTrackEnergyMatchAmbiguous_.push_back(0U);
    recoTrackEnergyBestMcTrackId_.push_back(0);
    recoTrackEnergyBestPdg_.push_back(0);
    recoTrackEnergyBestTrajectoryValid_.push_back(0U);
    recoTrackEnergyBestProcess_.emplace_back();
    recoTrackEnergyBestEndProcess_.emplace_back();
    recoTrackEnergyBestStartXcm_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestStartYcm_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestStartZcm_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestEndXcm_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestEndYcm_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestEndZcm_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestInitialPxGeV_.push_back(
        std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestInitialPyGeV_.push_back(
        std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestInitialPzGeV_.push_back(
        std::numeric_limits<double>::quiet_NaN());
    recoTrackEnergyBestMotherMcTrackId_.push_back(0);
    recoTrackEnergyBestMotherValid_.push_back(0U);
    recoTrackEnergyBestMotherPdg_.push_back(0);
    recoTrackEnergyBestDaughterCount_.push_back(0U);
    recoTrackHitMatchValid_.push_back(0U);
    recoTrackHitMatchAmbiguous_.push_back(0U);
    recoTrackHitBestMcTrackId_.push_back(0);
    recoTrackHitBestPdg_.push_back(0);
    recoTrackHitBestSharedHits_.push_back(0U);
    recoTrackHitBestSharedDeltaRayHits_.push_back(0U);
    recoTrackPurityValid_.push_back(0U);
    recoTrackPurity_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackCompletenessValid_.push_back(0U);
    recoTrackCompleteness_.push_back(std::numeric_limits<double>::quiet_NaN());
    recoTrackTruthOriginValid_.push_back(0U);
    recoTrackTruthOrigin_.push_back(-1);
    recoTrackHierarchyValid_.push_back(0U);
    recoTrackHierarchyTruncated_.push_back(0U);
    recoTrackHierarchyRootMcTrackId_.push_back(0);

    // Truth defaults above are appended for every reco track before any MC
    // call.  That preserves one-to-one vector alignment if matching fails.

    // 7b. In MC, use this reco Track and the configured Track label to obtain
    // its contributing MCParticles.  The explicit code below saves both the
    // full contributor lists and the leading contributor from each method.
    if (clockData && inventory) {
      try {
        // ProtoDUNETruthUtils returns MC contributors ordered by decreasing
        // deposited-energy fraction.  The first non-null entry is therefore
        // the leading energy match for this reconstructed track.
        auto const mcParticleMatches = truthUtils.GetMCParticleListFromRecoTrack(
            *clockData, track, evt, trackTag_.encode());
        simb::MCParticle const *leadingEnergyMatchedParticle = nullptr;
        double leadingEnergyFraction = 0.;
        unsigned int energyRank = 0U;
        for (auto const &match : mcParticleMatches) {
          // GetMCParticleListFromRecoTrack returns this pair:
          //   match.first  -> the contributing Geant MCParticle
          //   match.second -> its deposited-energy fraction in this reco track.
          simb::MCParticle const *particle = match.first;
          double const energyFraction = match.second;
          if (!particle || !std::isfinite(energyFraction) || energyFraction < 0.)
            continue;
          if (energyRank == 0U) {
            leadingEnergyMatchedParticle = particle;
            leadingEnergyFraction = energyFraction;
          } else if (energyRank == 1U &&
                     std::abs(energyFraction - leadingEnergyFraction) <= 1.e-12) {
            recoTrackEnergyMatchAmbiguous_.back() = 1U;
          }
          energyMatchMcTrackId_.push_back(particle->TrackId());
          energyMatchPdg_.push_back(particle->PdgCode());
          energyMatchEnergyFraction_.push_back(energyFraction);
          energyMatchRank_.push_back(energyRank);
          ++energyRank;
        }
        recoTrackEnergyMatchValid_.back() =
            leadingEnergyMatchedParticle ? 1U : 0U;

        // The hit utility independently returns contributors in descending
        // shared-hit count.  This ranking is intentionally separate from the
        // deposited-energy ranking above.
        auto const hitContributors = truthUtils.GetMCParticleListByHits(
            *clockData, track, evt, trackTag_.encode(), hitTag_.encode());
        simb::MCParticle const *bestHitParticle = nullptr;
        std::size_t leadingHitEvidence = 0U;
        std::size_t bestHitSharedHits = 0U;
        std::size_t bestHitSharedDeltaRayHits = 0U;
        unsigned int hitRank = 0U;
        for (auto const &match : hitContributors) {
          simb::MCParticle const *particle = match.particle;
          if (!particle) continue;
          std::size_t const hitEvidence =
              match.nSharedHits + match.nSharedDeltaRayHits;
          if (hitRank == 0U) {
            bestHitParticle = particle;
            leadingHitEvidence = hitEvidence;
            bestHitSharedHits = match.nSharedHits;
            bestHitSharedDeltaRayHits = match.nSharedDeltaRayHits;
          } else if (hitRank == 1U && hitEvidence == leadingHitEvidence) {
            recoTrackHitMatchAmbiguous_.back() = 1U;
          }
          hitMatchMcTrackId_.push_back(particle->TrackId());
          hitMatchPdg_.push_back(particle->PdgCode());
          hitMatchSharedHits_.push_back(ToUnsigned(match.nSharedHits));
          hitMatchSharedDeltaRayHits_.push_back(ToUnsigned(match.nSharedDeltaRayHits));
          hitMatchRank_.push_back(hitRank);
          ++hitRank;
        }

        // Save the leading hit-count match separately from the leading
        // deposited-energy match.  They can legitimately identify different
        // MCParticles for the same reconstructed track.
        if (bestHitParticle) {
          recoTrackHitMatchValid_.back() = 1U;
          recoTrackHitBestMcTrackId_.back() = bestHitParticle->TrackId();
          recoTrackHitBestPdg_.back() = bestHitParticle->PdgCode();
          recoTrackHitBestSharedHits_.back() = ToUnsigned(bestHitSharedHits);
          recoTrackHitBestSharedDeltaRayHits_.back() =
              ToUnsigned(bestHitSharedDeltaRayHits);
        }

        // 7c. Read the leading matched MCParticle directly.  These are the
        // standard MCParticle accessors a learner can reuse: TrackId, PDG,
        // process, endpoints, momentum, mother, daughters, and MCTruth origin.
        if (leadingEnergyMatchedParticle) {
          recoTrackEnergyBestMcTrackId_.back() =
              leadingEnergyMatchedParticle->TrackId();
          recoTrackEnergyBestPdg_.back() =
              leadingEnergyMatchedParticle->PdgCode();
          recoTrackEnergyBestProcess_.back() =
              leadingEnergyMatchedParticle->Process();
          recoTrackEnergyBestEndProcess_.back() =
              leadingEnergyMatchedParticle->EndProcess();
          if (leadingEnergyMatchedParticle->NumberTrajectoryPoints() > 0U) {
            recoTrackEnergyBestTrajectoryValid_.back() = 1U;
            recoTrackEnergyBestStartXcm_.back() = leadingEnergyMatchedParticle->Vx();
            recoTrackEnergyBestStartYcm_.back() = leadingEnergyMatchedParticle->Vy();
            recoTrackEnergyBestStartZcm_.back() = leadingEnergyMatchedParticle->Vz();
            recoTrackEnergyBestEndXcm_.back() = leadingEnergyMatchedParticle->EndX();
            recoTrackEnergyBestEndYcm_.back() = leadingEnergyMatchedParticle->EndY();
            recoTrackEnergyBestEndZcm_.back() = leadingEnergyMatchedParticle->EndZ();
            recoTrackEnergyBestInitialPxGeV_.back() = leadingEnergyMatchedParticle->Px();
            recoTrackEnergyBestInitialPyGeV_.back() = leadingEnergyMatchedParticle->Py();
            recoTrackEnergyBestInitialPzGeV_.back() = leadingEnergyMatchedParticle->Pz();
          }
          recoTrackEnergyBestMotherMcTrackId_.back() =
              leadingEnergyMatchedParticle->Mother();
          recoTrackEnergyBestDaughterCount_.back() =
              ToUnsigned(leadingEnergyMatchedParticle->NumberDaughters());
          if (leadingEnergyMatchedParticle->Mother() != 0) {
            auto const *mother = inventory->TrackIdToParticle_P(
                leadingEnergyMatchedParticle->Mother());
            if (mother) {
              recoTrackEnergyBestMotherValid_.back() = 1U;
              recoTrackEnergyBestMotherPdg_.back() = mother->PdgCode();
            }
          }
          auto const truth = inventory->TrackIdToMCTruth_P(
              leadingEnergyMatchedParticle->TrackId());
          if (truth) {
            recoTrackTruthOriginValid_.back() = 1U;
            recoTrackTruthOrigin_.back() = static_cast<int>(truth->Origin());
          }
          // Purity and completeness are separate utility calculations using
          // the same reco Track, detector clocks, and configured labels.
          auto const purity = truthUtils.GetPurity(*clockData, track, evt, trackTag_.encode());
          if (std::isfinite(purity)) {
            recoTrackPurityValid_.back() = 1U;
            recoTrackPurity_.back() = purity;
          }
          auto const completeness = truthUtils.GetCompleteness(
              *clockData, track, evt, trackTag_.encode(), hitTag_.encode());
          if (std::isfinite(completeness)) {
            recoTrackCompletenessValid_.back() = 1U;
            recoTrackCompleteness_.back() = completeness;
          }

          // 7d. Starting from this same MCParticle, save its full direct
          // daughter hierarchy for offline reaction-chain studies.
          recoTrackHierarchyRootMcTrackId_.back() =
              leadingEnergyMatchedParticle->TrackId();
          if (mcParticleProductAvailableForHierarchy) {
            std::unordered_set<int> activePath;
            bool truncated = false;
            int const rootIndex = AppendGeantHierarchy(
                leadingEnergyMatchedParticle->TrackId(), -1, 0U, activePath,
                truncated, *inventory);
            recoTrackHierarchyValid_.back() = rootIndex >= 0 ? 1U : 0U;
            recoTrackHierarchyTruncated_.back() = truncated ? 1U : 0U;
          }
        }

      } catch (std::exception const &error) {
        energyMatchMcTrackId_.resize(energyBegin);
        energyMatchPdg_.resize(energyBegin);
        energyMatchEnergyFraction_.resize(energyBegin);
        energyMatchRank_.resize(energyBegin);
        hitMatchMcTrackId_.resize(hitBegin);
        hitMatchPdg_.resize(hitBegin);
        hitMatchSharedHits_.resize(hitBegin);
        hitMatchSharedDeltaRayHits_.resize(hitBegin);
        hitMatchRank_.resize(hitBegin);
        hierarchyRecoTrackIndex_.resize(hierarchyBegin);
        hierarchyMcTrackId_.resize(hierarchyBegin);
        hierarchyParentIndex_.resize(hierarchyBegin);
        hierarchyMotherMcTrackId_.resize(hierarchyBegin);
        hierarchyPdg_.resize(hierarchyBegin);
        hierarchyDepth_.resize(hierarchyBegin);
        hierarchyProcess_.resize(hierarchyBegin);
        hierarchyEndProcess_.resize(hierarchyBegin);
        hierarchyDaughterBegin_.resize(hierarchyBegin);
        hierarchyDaughterCount_.resize(hierarchyBegin);
        hierarchyDaughterIndices_.resize(daughterIndexBegin);
        hierarchyTrajectoryOffsets_.resize(hierarchyBegin + 1U);
        hierarchyTrajectoryXcm_.resize(hierarchyTrajectoryOffsets_.back());
        hierarchyTrajectoryYcm_.resize(hierarchyTrajectoryOffsets_.back());
        hierarchyTrajectoryZcm_.resize(hierarchyTrajectoryOffsets_.back());
        hierarchyTrajectoryPxGeV_.resize(hierarchyTrajectoryOffsets_.back());
        hierarchyTrajectoryPyGeV_.resize(hierarchyTrajectoryOffsets_.back());
        hierarchyTrajectoryPzGeV_.resize(hierarchyTrajectoryOffsets_.back());
        ResetLatestTrackTruth();
        if (eventMessagesEmitted_++ < maxEventMessages_)
          mf::LogWarning("PDHDRecoTruthMatch")
              << "Truth matching failed for reco track " << track.ID()
              << " in run " << run_ << ", subrun " << subrun_ << ", event "
              << event_ << ": " << error.what();
      }
    }
    // Close this reco track's three flattened ranges after either successful
    // matching or the rollback in the exception handler.
    recoTrackEnergyMatchOffsets_.push_back(
        static_cast<unsigned long long>(energyMatchMcTrackId_.size()));
    recoTrackHitMatchOffsets_.push_back(
        static_cast<unsigned long long>(hitMatchMcTrackId_.size()));
    recoTrackHierarchyOffsets_.push_back(
        static_cast<unsigned long long>(hierarchyMcTrackId_.size()));
  }

  // 8. One tree entry represents this event and contains aligned vectors for
  // every reconstructed track collected above.
  tree_->Fill();
}

} // namespace pdhd::diagnostics

DEFINE_ART_MODULE(pdhd::diagnostics::PDHDRecoTruthMatch)
