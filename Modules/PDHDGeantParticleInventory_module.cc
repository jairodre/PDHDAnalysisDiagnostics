/**
 * @file PDHDGeantParticleInventory_module.cc
 * @brief Template for the complete Geant4 MCParticle inventory.
 *
 * Inputs: ParticleInventoryService and its configured G4/MCTruth products.
 * Outputs: one row per MCParticle containing stable event/track keys, PDG,
 * status, mother ID, direct-daughter IDs/count, process/end process, start/end
 * four-vectors, trajectory size, active-volume entry/exit, path length, and
 * MCTruth/origin linkage validity. Preserve raw endpoints as well as geometry-
 * clipped quantities. Ancestry traversal belongs in GeantHierarchyAlg so loops,
 * missing parents, and deterministic ordering have one implementation.
 *
 * Status: structure-only; MC-only and not registered as a plugin.
 */

// TODO: Implement base facts before derived interaction classifications.
