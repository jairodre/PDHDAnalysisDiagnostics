/**
 * @file PDHDRecoTruthMatch_module.cc
 * @brief Template for persistent bidirectional reconstructed-to-truth diagnostics.
 *
 * Inputs: tracks/showers/PFPs and their hits, standard BackTrackerMatchingData
 * associations when available, otherwise the configured official backtracking
 * utilities. Outputs: all retained match candidates plus best-match summaries,
 * shared hits/charge/energy, purity, completeness, rank, ambiguity, delta-ray
 * policy, and association-method source. Also emit staged truth-beam reconstruction
 * flags from truth PFP match through track/shower and nominal beam selection.
 * Data events receive explicit unavailable status and no truth-service calls.
 *
 * Status: structure-only; MC-only content and no art plugin is registered yet.
 */

// TODO: Fix denominators and existing-association precedence before implementation.
