/**
 * @file PDHDMCTruthInventory_module.cc
 * @brief Template for generator-level MCTruth and primary-particle inventory.
 *
 * Inputs: one or more explicitly configured std::vector<simb::MCTruth> products.
 * Outputs: one row per MCTruth record and one row per generated particle, keyed
 * by event, generator InputTag, MCTruth index, and particle index. Preserve
 * origin, generator status, PDG, four-position, four-momentum, mother/daughter
 * fields available at this stage, and beam-record diagnostics. Do not conflate
 * a generator primary with a Geant4 particle or a reconstructed beam candidate.
 * On real data, record module/product unavailability without calling MC services.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Implement after generator labels and beam-origin definitions are fixed.
