/**
 * @file PDHDRecoTruthMatch_module.cc
 * @brief Writes one event-level tree relating every reconstructed track to MC truth.
 *
 * Inputs are the configured Pandora PFParticles, slices, tracks, and hits.  In
 * MC, ProtoDUNETruthUtils supplies energy- and hit-ranked contributors, while
 * ParticleInventoryService resolves MCTruth origin for the configured,
 * event-wide Geant MCParticle inventory.  The tree retains all reconstructed
 * tracks, marks membership in Pandora's unique IsTestBeam slice, and stores
 * flattened reconstructed/MCParticle XYZ trajectories.
 * Data events keep the reconstructed inventory and explicitly mark every
 * truth-only quantity unavailable; no truth service is called for data.
 *
 * ART products and services used (FHiCL key -> default label):
 *   std::vector<recob::PFParticle>     PFParticleTag -> "pandora"
 *   PFParticle<->PFParticleMetadata    PFParticleMetadataTag (default: PFParticleTag)
 *   PFParticle<->Slice                 PFParticleSliceTag    (default: PFParticleTag)
 *   PFParticle<->Track, Track<->Hit    TrackTag      -> "pandoraTrack"
 *   std::vector<recob::Hit>            HitTag        -> "gaushit"  (completeness only)
 *   std::vector<simb::MCParticle>      MCParticleTag -> "largeant" (MC only)
 *   cheat::BackTrackerService          hit -> sim::TrackIDE (used inside the utilities)
 *   cheat::ParticleInventoryService    TrackID -> MCParticle / MCTruth
 *   detinfo::DetectorClocksService     timing needed to backtrack hits
 *
 * Output: one TTree entry per art event.  Branches come in three groups:
 *   - event scalars (run/subrun/event, availability flags, beam-slice ID);
 *   - "reco_track_*" vectors, one element per recob::Track;
 *   - "mc_particle_*" vectors, one element per retained simb::MCParticle
 *     (low-level EM secondaries are dropped, see
 *     ExcludeElectromagneticSecondaryParticle, unless truth matching ran and
 *     they contributed to a reco track).
 * Variable-length lists (trajectory points, match contributors) are
 * "flattened": all objects' entries are concatenated into one vector, and an
 * "*_offsets" vector with N+1 entries gives object i the index range
 * [offsets[i], offsets[i+1]).
 *
 * Truth matching (EnableMCTruthMatching, default false) is optional because it
 * backtracks every hit of every track and is the slowest part of the job.
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
#include <cstdlib>
#include <exception>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// File-local helpers (anonymous namespace: not visible outside this file).
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

// Exclude configured EM secondary particles from the stored inventory.  Their
// metadata and trajectory coordinates are both omitted to reduce ROOT size.
// Process() is the Geant4 name of the process that created the particle
// (e.g. "compt" = Compton scattering, "eIoni"/"muIoni" = ionisation delta
// rays, "phot" = photoelectric effect, "conv" = pair production, "eBrem" =
// bremsstrahlung).  std::abs(pdg) == 11 covers both e- and e+.
// The filter only affects this output: truth matching still sees every
// particle.  An EM secondary listed as an energy or hit contributor of any
// reco track is kept anyway (see appendRetainedMCParticles in analyze()), so
// with EnableMCTruthMatching false every matching EM secondary is dropped.
// Only the contributor itself is kept, not its ancestors: a kept Compton
// electron can still have a dropped photon as its mother.
bool ExcludeElectromagneticSecondaryParticle(
    simb::MCParticle const &particle) {
  int const pdg = particle.PdgCode();
  std::string const &process = particle.Process();
  bool const electron = std::abs(pdg) == 11 &&
                        (process == "compt" || process == "eIoni" ||
                         process == "muIoni" || process == "hIoni" ||
                         process == "ionIoni" || process == "phot" ||
                         process == "conv");
  bool const photon =
      pdg == 22 && (process == "annihil" || process == "conv" ||
                    process == "eBrem");
  return electron || photon;
}

} // namespace

namespace pdhd::diagnostics {

// An EDAnalyzer only reads the event; it cannot add products to it (that is
// what an EDProducer does).  art calls beginJob() once, then analyze() once
// per event.
class PDHDRecoTruthMatch : public art::EDAnalyzer {
public:
  explicit PDHDRecoTruthMatch(fhicl::ParameterSet const &);

private:
  void beginJob() override;
  void analyze(art::Event const &) override;

  // Clears every branch vector and sets event scalars to their "not
  // available" defaults, so no value can leak from the previous event.
  void ResetEventState(art::Event const &);
  // Restores every truth-only field for the most recently appended reco track.
  // This keeps scalar facts consistent with rolled-back flattened ranges.
  void ResetLatestTrackTruth();
  // Appends one track's valid 3D points to the flattened reco trajectory.
  void AppendRecoTrackTrajectory(recob::Track const &);
  // Appends one MCParticle to the event-wide inventory.  inventory may be
  // nullptr, in which case the truth origin is left marked invalid.
  void AppendMCParticle(simb::MCParticle const &, cheat::ParticleInventoryService *);

  // FHiCL-configured labels for the event products and associations read in
  // analyze().  An art::InputTag is "module_label:instance:process"; usually
  // only the module label is given (e.g. "pandoraTrack").
  art::InputTag pfpTag_;
  art::InputTag pfpMetadataTag_;
  art::InputTag pfpSliceTag_;
  art::InputTag trackTag_;
  art::InputTag hitTag_;
  art::InputTag mcParticleTag_;
  bool enableMCTruthMatching_;
  // Caps the number of warnings printed per event so a broken input file
  // does not flood the job log.
  unsigned int maxEventMessages_;
  unsigned int eventMessagesEmitted_ = 0;

  // ROOT output object and event-level branches.  The TTree is owned by
  // TFileService; this raw pointer must not be deleted here.
  TTree *tree_ = nullptr;

  unsigned int run_ = 0;
  unsigned int subrun_ = 0;
  unsigned int event_ = 0;
  bool isData_ = false;
  // "*_available" flags record whether a product or association was found in
  // this event, so downstream code can tell "missing input" from "empty".
  bool pfpProductAvailable_ = false;
  bool trackProductAvailable_ = false;
  bool pfpMetadataAssociationAvailable_ = false;
  bool pfpSliceAssociationAvailable_ = false;
  bool pfpTrackAssociationAvailable_ = false;
  // Beam slice: found means exactly one slice had an IsTestBeam primary PFP;
  // ambiguous means more than one did (then no slice is selected and the ID
  // stays -1).
  bool pandoraBeamSliceFound_ = false;
  bool pandoraBeamSliceAmbiguous_ = false;
  int pandoraBeamSliceId_ = -1;
  unsigned int pandoraBeamSliceCount_ = 0;
  bool mcParticleProductAvailable_ = false;
  // Records only whether the HitTag product exists; the utilities read the
  // hits themselves.
  bool truthHitProductAvailable_ = false;
  // Geant TrackIDs should be unique within an event.  Duplicates would make
  // TrackID-based lookups (mother, match counts) ambiguous, so they are
  // counted and only the first retained occurrence is used for lookups.
  bool mcParticleDuplicateFree_ = false;
  unsigned int mcParticleDuplicateTrackIdCount_ = 0;
  // True only when truth matching was requested and its prerequisites
  // (MC, ParticleInventoryService, detector clocks) were all available.
  bool mcTruthMatchingEvaluated_ = false;
  // Names of the utility functions used, written into the file so the
  // output documents how the matches were made.
  std::string energyMatchMethod_;
  std::string hitMatchMethod_;

  // One element per reco track.  All vectors in this group share the same
  // reco-track index, so e.g. reco_track_id[i] and reco_track_in_selected_beam_slice[i]
  // describe the same recob::Track.  unsigned char is used for booleans
  // (0/1) because std::vector<bool> cannot be written as a ROOT branch.
  // Collection index = position in the std::vector<recob::Track> product;
  // track ID = recob::Track::ID().  They are usually equal for pandoraTrack
  // but are conceptually different, so both are stored.
  std::vector<int> recoTrackCollectionIndex_;
  std::vector<int> recoTrackId_;
  std::vector<unsigned char> recoTrackInSelectedBeamSlice_;
  std::vector<unsigned char> recoTrackPfpAssociationValid_;
  std::vector<unsigned int> recoTrackPfpAssociationCount_;
  std::vector<unsigned int> recoTrackSelectedBeamSliceMembershipCount_;
  // Flattened fitted track positions [cm], detector coordinates.
  std::vector<unsigned long long> recoTrackTrajectoryOffsets_;
  std::vector<double> recoTrackTrajectoryXcm_;
  std::vector<double> recoTrackTrajectoryYcm_;
  std::vector<double> recoTrackTrajectoryZcm_;

  // Leading MC match and scalar truth properties for each reco track.
  //   energy match: the MCParticle that deposited the largest fraction of the
  //                 energy seen by the track's hits;
  //   hit match:    the MCParticle sharing the most hits with the track.
  // "ambiguous" means the first and second candidates are tied.  Every value
  // has a "*_valid" flag; when it is 0 the value is a default (0, -1 or NaN),
  // not a measurement.  Units: positions in cm, momenta in GeV/c.
  // Truth origin is simb::Origin_t: 0 unknown, 1 beam neutrino, 2 cosmic,
  // 3 supernova neutrino, 4 single particle (particle gun / beam generator).
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
  // Purity and completeness refer to the leading *energy* match:
  //   purity       = shared hits / hits on the reco track;
  //   completeness = shared hits / all HitTag hits of that MCParticle.
  std::vector<unsigned char> recoTrackPurityValid_;
  std::vector<double> recoTrackPurity_;
  std::vector<unsigned char> recoTrackCompletenessValid_;
  std::vector<double> recoTrackCompleteness_;
  std::vector<unsigned char> recoTrackTruthOriginValid_;
  std::vector<int> recoTrackTruthOrigin_;

  // Flattened all-contributor lists.  Offsets map a reco-track index to its
  // range: [offsets[i], offsets[i + 1]).  Within a range, rank 0 is the
  // leading contributor.
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

  // One compact entry for every retained MCParticle. The configured EM filter
  // intentionally removes selected particles, except those that contributed
  // to a reco track, so mother TrackIDs can refer to a particle absent from
  // this size-reduced inventory.  EM secondaries kept for that reason can be
  // recognised by their PDG/process together with mc_particle_has_reco_track
  // == 1; they are not a complete sample of that population. Only trajectory
  // coordinates are retained; per-point momenta are intentionally not written.
  std::vector<int> mcParticleTrackId_;
  std::vector<int> mcParticleMotherTrackId_;
  std::vector<int> mcParticlePdg_;
  std::vector<std::string> mcParticleProcess_;
  std::vector<std::string> mcParticleEndProcess_;
  std::vector<unsigned int> mcParticleDaughterCount_;
  std::vector<unsigned int> mcParticleTrajectoryPointCount_;
  std::vector<unsigned char> mcParticleHasDrawableTrajectory_;
  std::vector<unsigned char> mcParticleTruthOriginValid_;
  std::vector<int> mcParticleTruthOrigin_;
  // Reverse direction of the matching: for each MCParticle, how many reco
  // tracks list it anywhere in their contributor lists (not only as the
  // leading match).  Filled only when truth matching ran.
  std::vector<unsigned char> mcParticleHasRecoTrack_;
  std::vector<unsigned int> mcParticleEnergyMatchedRecoTrackCount_;
  std::vector<unsigned int> mcParticleHitMatchedRecoTrackCount_;
  // Flattened Geant trajectory points [cm].  These are the points LArG4 kept,
  // which may be a thinned subset of the steps Geant4 actually took.
  std::vector<unsigned long long> mcParticleTrajectoryOffsets_;
  std::vector<double> mcParticleTrajectoryXcm_;
  std::vector<double> mcParticleTrajectoryYcm_;
  std::vector<double> mcParticleTrajectoryZcm_;
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
      maxEventMessages_(p.get<unsigned int>("MaxEventMessages", 3U)) {
  // p.get<T>(key, default) reads a FHiCL parameter and falls back to the
  // default when the key is absent from the job configuration.
  //
  // Declare the products this EDAnalyzer may read.  The InputTags above say
  // which producer label to use; consumes<T>() declares the C++ product type.
  // The associations read through FindManyP/FindOneP and the products read
  // inside the ProtoDUNE utilities are not declared here.
  consumes<std::vector<recob::PFParticle>>(pfpTag_);
  consumes<std::vector<recob::Track>>(trackTag_);
  consumes<std::vector<recob::Hit>>(hitTag_);
  consumes<std::vector<simb::MCParticle>>(mcParticleTag_);
}

void PDHDRecoTruthMatch::beginJob() {
  // TFileService owns the job's analysis ROOT file (lar -T <file>, or
  // services.TFileService.fileName).  make<T>() creates the object inside it,
  // in a directory named after this module's FHiCL label.
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
  AddBranch(*tree_, "mc_particle_duplicate_free", mcParticleDuplicateFree_);
  AddBranch(*tree_, "mc_particle_duplicate_track_id_count",
            mcParticleDuplicateTrackIdCount_);
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

  // Event-wide MCParticle inventory.  MCParticle::Mother() reconstructs the
  // genealogy; trajectory momenta are intentionally omitted to control size.
  // Filtered EM secondaries (those not matched to any reco track) have no
  // row, so a mother TrackID may not be found among mc_particle_track_id.
  AddBranch(*tree_, "mc_particle_track_id", mcParticleTrackId_);
  AddBranch(*tree_, "mc_particle_mother_track_id", mcParticleMotherTrackId_);
  AddBranch(*tree_, "mc_particle_pdg", mcParticlePdg_);
  AddBranch(*tree_, "mc_particle_process", mcParticleProcess_);
  AddBranch(*tree_, "mc_particle_end_process", mcParticleEndProcess_);
  AddBranch(*tree_, "mc_particle_daughter_count", mcParticleDaughterCount_);
  AddBranch(*tree_, "mc_particle_trajectory_point_count",
            mcParticleTrajectoryPointCount_);
  AddBranch(*tree_, "mc_particle_has_drawable_trajectory",
            mcParticleHasDrawableTrajectory_);
  AddBranch(*tree_, "mc_particle_truth_origin_valid", mcParticleTruthOriginValid_);
  AddBranch(*tree_, "mc_particle_truth_origin", mcParticleTruthOrigin_);
  AddBranch(*tree_, "mc_particle_has_reco_track", mcParticleHasRecoTrack_);
  AddBranch(*tree_, "mc_particle_energy_matched_reco_track_count",
            mcParticleEnergyMatchedRecoTrackCount_);
  AddBranch(*tree_, "mc_particle_hit_matched_reco_track_count",
            mcParticleHitMatchedRecoTrackCount_);
  AddBranch(*tree_, "mc_particle_trajectory_offsets", mcParticleTrajectoryOffsets_);
  AddBranch(*tree_, "mc_particle_trajectory_x_cm", mcParticleTrajectoryXcm_);
  AddBranch(*tree_, "mc_particle_trajectory_y_cm", mcParticleTrajectoryYcm_);
  AddBranch(*tree_, "mc_particle_trajectory_z_cm", mcParticleTrajectoryZcm_);
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
  mcParticleDuplicateFree_ = false;
  mcParticleDuplicateTrackIdCount_ = 0;
  mcTruthMatchingEvaluated_ = false;
  energyMatchMethod_ = "ProtoDUNETruthUtils::GetMCParticleListFromRecoTrack";
  hitMatchMethod_ = "ProtoDUNETruthUtils::GetMCParticleListByHits";

  // Reco-track vectors and their flattened trajectory.  Each offsets vector
  // starts as {0}; every appended object pushes its end index, so N objects
  // give N+1 offsets.
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

  // Event-wide MCParticle inventory and flattened XYZ trajectories.
  mcParticleTrackId_.clear();
  mcParticleMotherTrackId_.clear();
  mcParticlePdg_.clear();
  mcParticleProcess_.clear();
  mcParticleEndProcess_.clear();
  mcParticleDaughterCount_.clear();
  mcParticleTrajectoryPointCount_.clear();
  mcParticleHasDrawableTrajectory_.clear();
  mcParticleTruthOriginValid_.clear();
  mcParticleTruthOrigin_.clear();
  mcParticleHasRecoTrack_.clear();
  mcParticleEnergyMatchedRecoTrackCount_.clear();
  mcParticleHitMatchedRecoTrackCount_.clear();
  mcParticleTrajectoryOffsets_.assign(1, 0ULL);
  mcParticleTrajectoryXcm_.clear();
  mcParticleTrajectoryYcm_.clear();
  mcParticleTrajectoryZcm_.clear();
}

void PDHDRecoTruthMatch::ResetLatestTrackTruth() {
  // These values must stay identical to the defaults pushed for each track in
  // analyze() step 8a.
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
}

void PDHDRecoTruthMatch::AppendRecoTrackTrajectory(recob::Track const &track) {
  // recob::Track can contain invalid fitted points.  Save only valid positions
  // and close this track's flattened range with the next offset.
  for (std::size_t point = 0; point < track.NumberTrajectoryPoints(); ++point) {
    if (!track.HasValidPoint(point)) continue;
    // LocationAtPoint returns a 3D point (recob::tracking::Point_t) in cm.
    auto const position = track.LocationAtPoint(point);
    recoTrackTrajectoryXcm_.push_back(position.X());
    recoTrackTrajectoryYcm_.push_back(position.Y());
    recoTrackTrajectoryZcm_.push_back(position.Z());
  }
  recoTrackTrajectoryOffsets_.push_back(
      static_cast<unsigned long long>(recoTrackTrajectoryXcm_.size()));
}

void PDHDRecoTruthMatch::AppendMCParticle(
    simb::MCParticle const &particle,
    cheat::ParticleInventoryService *inventory) {
  // TrackId() is the Geant4 track number, unique within the event.  Mother()
  // is the parent's TrackId (0 for a generator primary).  NumberDaughters()
  // counts all Geant daughters, including any the EM filter drops here.
  mcParticleTrackId_.push_back(particle.TrackId());
  mcParticleMotherTrackId_.push_back(particle.Mother());
  mcParticlePdg_.push_back(particle.PdgCode());
  mcParticleProcess_.push_back(particle.Process());
  mcParticleEndProcess_.push_back(particle.EndProcess());
  mcParticleDaughterCount_.push_back(ToUnsigned(particle.NumberDaughters()));
  std::size_t const pointCount = particle.NumberTrajectoryPoints();
  mcParticleTrajectoryPointCount_.push_back(ToUnsigned(pointCount));
  mcParticleHasDrawableTrajectory_.push_back(pointCount >= 2U ? 1U : 0U);
  mcParticleTruthOriginValid_.push_back(0U);
  mcParticleTruthOrigin_.push_back(-1);
  mcParticleHasRecoTrack_.push_back(0U);
  mcParticleEnergyMatchedRecoTrackCount_.push_back(0U);
  mcParticleHitMatchedRecoTrackCount_.push_back(0U);

  // TrackIdToMCTruth_P returns the generator record (simb::MCTruth) this
  // particle descends from; its Origin() says beam, cosmic, single particle...
  if (inventory) {
    auto const truth = inventory->TrackIdToMCTruth_P(particle.TrackId());
    if (truth) {
      mcParticleTruthOriginValid_.back() = 1U;
      mcParticleTruthOrigin_.back() = static_cast<int>(truth->Origin());
    }
  }

  for (std::size_t point = 0; point < pointCount; ++point) {
    // Position() is a TLorentzVector (x, y, z in cm; t in ns); only the
    // spatial part is stored.
    auto const position = particle.Position(point);
    mcParticleTrajectoryXcm_.push_back(position.X());
    mcParticleTrajectoryYcm_.push_back(position.Y());
    mcParticleTrajectoryZcm_.push_back(position.Z());
  }
  mcParticleTrajectoryOffsets_.push_back(
      static_cast<unsigned long long>(mcParticleTrajectoryXcm_.size()));
}

void PDHDRecoTruthMatch::analyze(art::Event const &evt) {
  ResetEventState(evt);

  // 1. Read the reconstructed products configured in FHiCL.  A handle is the
  // event-local view of one ART product; it may be invalid when that product
  // was not produced in this input file.  getHandle() does not throw when the
  // product is missing (getValidHandle() would), so each handle is tested
  // before use.  The hit handle is read only to record availability; the
  // utilities read the hits they need themselves.
  auto const trackHandle = evt.getHandle<std::vector<recob::Track>>(trackTag_);
  auto const pfpHandle = evt.getHandle<std::vector<recob::PFParticle>>(pfpTag_);
  auto const hitHandle = evt.getHandle<std::vector<recob::Hit>>(hitTag_);
  trackProductAvailable_ = static_cast<bool>(trackHandle);
  pfpProductAvailable_ = static_cast<bool>(pfpHandle);
  truthHitProductAvailable_ = static_cast<bool>(hitHandle);

  // 2. Construct the Pandora associations needed to answer two different
  // questions: which Slice contains a PFP, and which reco Tracks belong to a
  // PFP.  The associations are indexed by the PFParticle collection index:
  // slices->at(i) is the Slice of pfpHandle->at(i).
  //   FindOneP  -> at most one associated object (an art::Ptr, may be null);
  //   FindManyP -> a std::vector<art::Ptr<T>>, possibly empty.
  // Their constructors throw when the association product is absent, which
  // is why each one is built inside try/catch.
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
  // IsTestBeam metadata property identifies the selected beam slice.  A
  // primary PFP has no parent PFP: it is the top of one Pandora hierarchy.
  std::set<int> beamSliceIds;
  if (pfpHandle && metadata && slices) {
    protoana::ProtoDUNEPFParticleUtils pfpUtils;
    for (std::size_t pfpIndex = 0; pfpIndex < pfpHandle->size(); ++pfpIndex) {
      auto const &pfp = pfpHandle->at(pfpIndex);
      if (!pfp.IsPrimary()) continue;
      try {
        auto const slice = slices->at(pfpIndex);
        // IsBeamParticle reads the PFP's metadata itself and returns true when
        // the "IsTestBeam" key is present.  The local `metadata` association
        // is therefore used only as an availability check; `slices` supplies
        // the Slice ID stored in this tree.
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
  // with the unique selected beam slice.  Counts rather than booleans are kept
  // so that a track associated with more than one PFP is visible downstream.
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

  // 5. Read the complete MCParticle product once. It is the source for the
  // size-reduced inventory, independent of reco/truth matching.  The MCParticle
  // product is written by Geant4 simulation (LArG4, label "largeant") and does
  // not exist in data, so it is not requested for data events.
  std::vector<simb::MCParticle> const *mcParticles = nullptr;
  // TrackID -> row in the stored compact inventory; filled after reco/truth
  // contributor lists establish which EM particles must be retained.
  std::unordered_map<int, std::size_t> mcParticleIndexByTrackId;
  if (!isData_) {
    auto const mcParticleHandle = evt.getHandle<std::vector<simb::MCParticle>>(mcParticleTag_);
    mcParticleProductAvailable_ = static_cast<bool>(mcParticleHandle);
    if (mcParticleHandle) {
      mcParticles = mcParticleHandle.product();
      std::set<int> seenTrackIds;
      std::set<int> duplicateTrackIds;
      for (std::size_t index = 0; index < mcParticles->size(); ++index) {
        auto const &particle = mcParticles->at(index);
        if (!seenTrackIds.insert(particle.TrackId()).second)
          duplicateTrackIds.insert(particle.TrackId());
      }
      mcParticleDuplicateTrackIdCount_ = ToUnsigned(duplicateTrackIds.size());
      mcParticleDuplicateFree_ = duplicateTrackIds.empty();
    }
  }

  // 6. Obtain the MCTruth-origin lookup now; compact-inventory construction is
  // deferred until reco/truth contributors have been collected below.
  cheat::ParticleInventoryService *inventory = nullptr;
  if (mcParticles) {
    try {
      // A ServiceHandle gives access to a job-wide service configured in
      // FHiCL.  The raw pointer stays valid for the whole job.
      art::ServiceHandle<cheat::ParticleInventoryService> inventoryService;
      inventory = inventoryService.operator->();
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_++ < maxEventMessages_)
        mf::LogWarning("PDHDRecoTruthMatch")
            << "ParticleInventoryService is unavailable for MCParticle origins in run "
            << run_ << ", subrun " << subrun_ << ", event "
            << event_ << ": " << error.what();
    }
  }

  // Builds the stored inventory in largeant order.  A particle is written
  // unless the EM filter selects it AND it is absent from the reco-track
  // contributor TrackIDs.  Called exactly once per event: with an empty set
  // when there are no tracks, otherwise after the track loop.  Contributor
  // TrackIDs come from ParticleInventoryService's particle list, so they only
  // match here when that list and MCParticleTag refer to the same product.
  auto appendRetainedMCParticles =
      [&](std::set<int> const &recoContributorTrackIds) {
        if (!mcParticles) return;
        for (auto const &particle : *mcParticles) {
          bool const retainedByRecoTrack =
              recoContributorTrackIds.find(particle.TrackId()) !=
              recoContributorTrackIds.end();
          if (ExcludeElectromagneticSecondaryParticle(particle) &&
              !retainedByRecoTrack) {
            continue;
          }
          std::size_t const storedIndex = mcParticleTrackId_.size();
          AppendMCParticle(particle, inventory);
          mcParticleIndexByTrackId.emplace(particle.TrackId(), storedIndex);
        }
      };

  if (!trackHandle) {
    // Without reconstructed tracks there are no truth contributors to retain.
    appendRetainedMCParticles({});
    tree_->Fill();
    return;
  }

  // 7. Detector clocks are required only by the reco/truth matching utilities:
  // the BackTracker converts a hit's TPC tick into the simulated time window
  // in which to look for energy deposits.
  std::unique_ptr<detinfo::DetectorClocksData> clockData;
  if (!isData_ && enableMCTruthMatching_ && inventory) {
    try {
      clockData = std::make_unique<detinfo::DetectorClocksData>(
          art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt));
      mcTruthMatchingEvaluated_ = true;
    } catch (std::exception const &error) {
      if (eventMessagesEmitted_++ < maxEventMessages_)
        mf::LogWarning("PDHDRecoTruthMatch")
            << "Detector clocks are unavailable for MC truth matching in run "
            << run_ << ", subrun " << subrun_ << ", event " << event_ << ": "
            << error.what();
    }
  }

  // 8. Loop over the complete configured reco Track collection.  This module
  // does not discard non-beam-slice tracks; reco_track_in_selected_beam_slice
  // is the downstream selection flag.
  protoana::ProtoDUNETruthUtils truthUtils;
  for (std::size_t trackIndex = 0; trackIndex < trackHandle->size(); ++trackIndex) {
    auto const &track = trackHandle->at(trackIndex);
    auto const energyBegin = energyMatchMcTrackId_.size();
    auto const hitBegin = hitMatchMcTrackId_.size();
    // 8a. Save the reco-only information first.  These branches are equally
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

    // Truth defaults above are appended for every reco track before any MC
    // call.  That preserves one-to-one vector alignment if matching fails.

    // 8b. In MC, use this reco Track and the configured Track label to obtain
    // its contributing MCParticles.  The explicit code below saves both the
    // full contributor lists and the leading contributor from each method.
    // clockData is non-null only for MC with EnableMCTruthMatching true.
    if (clockData && inventory) {
      try {
        // ProtoDUNETruthUtils returns MC contributors ordered by decreasing
        // deposited-energy fraction.  The first non-null entry is therefore
        // the leading energy match for this reconstructed track.
        // What the utility does internally:
        //   1. gets the track's hits from the Track<->Hit association
        //      created by the TrackTag producer;
        //   2. BackTrackerService::HitToTrackIDEs(hit) -> sim::TrackIDE
        //      (Geant TrackID + energy deposited in that hit);
        //   3. ParticleInventoryService::TrackIdToParticle_P(TrackID)
        //      -> simb::MCParticle;
        //   4. sums energy per MCParticle and divides by the total energy
        //      on the track's hits.
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
          // Rank 0 is the leading match; a rank-1 fraction equal to it within
          // floating-point precision marks the match as ambiguous.
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
        // The candidates are the same MCParticles found by the energy method;
        // for each one the utility counts the track's hits that backtrack to
        //   nSharedHits:         its TrackID;
        //   nSharedDeltaRayHits: minus its TrackID.  LArG4 labels energy from
        //                        EM secondaries it did not keep as separate
        //                        MCParticles with the parent's negative ID.
        // Ranking uses the sum of both.  The hit label is passed through, but
        // the track's hits themselves come from the Track<->Hit association.
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

        // 8c. Read the leading matched MCParticle directly.  These are the
        // standard MCParticle accessors a learner can reuse: TrackId, PDG,
        // process, endpoints, momentum, mother, daughters, and MCTruth origin.
        // Vx/Vy/Vz and Px/Py/Pz are taken at the first trajectory point and
        // EndX/EndY/EndZ at the last one, which is why the trajectory must
        // not be empty.  Momentum is in GeV/c.
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
          // Mother() == 0 means a generator primary, which has no mother.
          // The lookup goes through ParticleInventoryService, so it works
          // even when the mother was dropped from this module's output.
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
          // Both repeat the energy matching internally and use its leading
          // MCParticle (the same one as above):
          //   purity       = hits shared with it / hits on the track;
          //   completeness = hits shared with it / all its hits in HitTag.
          // Completeness backtracks every hit in the event, so it is the
          // most expensive call in this loop.
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

        }

      } catch (std::exception const &error) {
        // Roll back: drop this track's partially written contributors and
        // restore its truth defaults, so every vector stays aligned.
        energyMatchMcTrackId_.resize(energyBegin);
        energyMatchPdg_.resize(energyBegin);
        energyMatchEnergyFraction_.resize(energyBegin);
        energyMatchRank_.resize(energyBegin);
        hitMatchMcTrackId_.resize(hitBegin);
        hitMatchPdg_.resize(hitBegin);
        hitMatchSharedHits_.resize(hitBegin);
        hitMatchSharedDeltaRayHits_.resize(hitBegin);
        hitMatchRank_.resize(hitBegin);
        ResetLatestTrackTruth();
        if (eventMessagesEmitted_++ < maxEventMessages_)
          mf::LogWarning("PDHDRecoTruthMatch")
              << "Truth matching failed for reco track " << track.ID()
              << " in run " << run_ << ", subrun " << subrun_ << ", event "
              << event_ << ": " << error.what();
      }
    }
    // Close this reco track's two flattened ranges after either successful
    // matching or the rollback in the exception handler.
    recoTrackEnergyMatchOffsets_.push_back(
        static_cast<unsigned long long>(energyMatchMcTrackId_.size()));
    recoTrackHitMatchOffsets_.push_back(
        static_cast<unsigned long long>(hitMatchMcTrackId_.size()));
  }

  // Keep every Geant particle used by either reco-track truth method. The
  // contributor lists are complete only after all transactional rollbacks.
  std::set<int> recoContributorTrackIds(energyMatchMcTrackId_.begin(),
                                        energyMatchMcTrackId_.end());
  recoContributorTrackIds.insert(hitMatchMcTrackId_.begin(),
                                 hitMatchMcTrackId_.end());
  appendRetainedMCParticles(recoContributorTrackIds);

  // Derive event-wide MCParticle-to-reco flags only after all per-reco ranges
  // are complete. This makes the counts transactional if one track's truth
  // utility throws and its flattened contributors are rolled back.
  // For each reco track, each distinct contributor TrackID adds 1 to that
  // MCParticle's count. All listed contributors are retained above when their
  // TrackID resolves in the configured MCParticle product.
  auto countMatchedRecoTracks = [&](std::vector<unsigned long long> const &offsets,
                                    std::vector<int> const &matchTrackIds,
                                    std::vector<unsigned int> &counts) {
    for (std::size_t reco = 0; reco + 1U < offsets.size(); ++reco) {
      std::size_t const begin = static_cast<std::size_t>(offsets.at(reco));
      std::size_t const end = static_cast<std::size_t>(offsets.at(reco + 1U));
      if (begin > end || end > matchTrackIds.size()) continue;
      std::set<int> uniqueContributors;
      for (std::size_t match = begin; match < end; ++match) {
        int const trackId = matchTrackIds.at(match);
        if (!uniqueContributors.insert(trackId).second) continue;
        auto const found = mcParticleIndexByTrackId.find(trackId);
        if (found != mcParticleIndexByTrackId.end()) ++counts.at(found->second);
      }
    }
  };
  countMatchedRecoTracks(recoTrackEnergyMatchOffsets_, energyMatchMcTrackId_,
                         mcParticleEnergyMatchedRecoTrackCount_);
  countMatchedRecoTracks(recoTrackHitMatchOffsets_, hitMatchMcTrackId_,
                         mcParticleHitMatchedRecoTrackCount_);
  for (std::size_t particle = 0; particle < mcParticleHasRecoTrack_.size(); ++particle) {
    bool const hasRecoTrack =
        mcParticleEnergyMatchedRecoTrackCount_.at(particle) > 0U ||
        mcParticleHitMatchedRecoTrackCount_.at(particle) > 0U;
    mcParticleHasRecoTrack_.at(particle) = hasRecoTrack ? 1U : 0U;
  }

  // 9. One tree entry represents this event and contains aligned vectors for
  // every reconstructed track and every retained MCParticle.
  tree_->Fill();
}

} // namespace pdhd::diagnostics

// Registers the class with art so a FHiCL job can use it as
// module_type: PDHDRecoTruthMatch.
DEFINE_ART_MODULE(pdhd::diagnostics::PDHDRecoTruthMatch)
