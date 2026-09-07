/**
 * @file PDHDCosmicSelectionStages_module.cc
 * @brief Template for independently inspectable cosmic-selection stages.
 *
 * Inputs: PFP/track inventories, Pandora clear-cosmic metadata/BDT scores, timing,
 * geometry, blip associations, and optional MC truth matches. Outputs: every raw
 * diagnostic and boolean stage, including independently named start-contained,
 * end-contained, both-contained, full-trajectory-contained, clear-cosmic, BDT,
 * and final selections. Analyses requesting endpoint containment must choose the
 * exact endpoint flag; the module must not silently strengthen it. MC truth PDG
 * categories should use mutually exclusive muon/proton/pion/kaon/other definitions.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Define named selection configurations without imposing cuts by default.
