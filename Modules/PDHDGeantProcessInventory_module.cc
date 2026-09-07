/**
 * @file PDHDGeantProcessInventory_module.cc
 * @brief Template for Geant4 ancestry edges and interaction-process diagnostics.
 *
 * Inputs: the MCParticle graph resolved by ParticleInventoryService.
 * Outputs: direct parent-child edge rows and process/trajectory transition rows
 * with event keys, parent and child track IDs, generation depth, process names,
 * positions, momenta, and traversal-validity flags. The base output should keep
 * literal Geant4 information; pion/kaon/muon topology labels are derived views
 * and must document their definitions rather than replacing process strings.
 * Cycles, absent mothers, and particles outside the active volume remain visible.
 *
 * Status: structure-only; MC-only and not registered as a plugin.
 */

// TODO: Define whether trajectory transitions are available in current products.
