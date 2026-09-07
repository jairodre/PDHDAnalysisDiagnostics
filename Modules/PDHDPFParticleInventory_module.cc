/**
 * @file PDHDPFParticleInventory_module.cc
 * @brief Template for Pandora PFParticle hierarchy and metadata inventory.
 *
 * Inputs: configured PFParticle collection and official ProtoDUNEPFParticleUtils
 * associations. Outputs: collection index/self/parent/daughter keys, primary and
 * beam-slice status, hierarchy depth, PDG hypothesis, metadata, TrackScore,
 * clear-cosmic flag, beam/cosmic score, vertices, and associated object counts.
 * Keep PFP metadata scores distinct from hit-level CNN scores. Store every PFP,
 * including cosmics and ambiguous objects; selection decisions are separate.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Fix stable hierarchy-key and missing-metadata policies before coding.
