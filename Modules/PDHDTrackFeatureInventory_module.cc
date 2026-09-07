/**
 * @file PDHDTrackFeatureInventory_module.cc
 * @brief Template for reconstructed track PID, momentum, CNN, and summary features.
 *
 * Inputs: tracks, ParticleID, calorimetry, PFP metadata and hit-level CNN output.
 * Outputs: range momenta under proton/muon/pion hypotheses, MCS momentum and fit
 * status, chi2 and NDF for proton/pion/muon/kaon when available, PFP TrackScore,
 * separate charge-weighted/unweighted CNN class summaries, and per-plane
 * calorimetry summaries. Every algorithm/source has its own validity and label;
 * no fallback score may overwrite another score. These are observables for later
 * cut development, not automatic beam or cosmic acceptance.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Confirm ParticleID score conventions and TrackMomentumCalculator support.
