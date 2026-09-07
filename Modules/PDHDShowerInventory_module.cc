/**
 * @file PDHDShowerInventory_module.cc
 * @brief Template for geometry-only reconstructed-shower inventory.
 *
 * Inputs: configured recob::Shower collections and PFP associations. Outputs:
 * stable shower/PFP keys, start, direction, length/open angle, best plane, plane
 * energy/dE/dx vectors as supplied by reconstruction, errors, containment, and
 * validity per optional quantity. Do not treat a missing track as proof that a
 * PFP is a shower; preserve track/shower multiplicities and ambiguity explicitly.
 * Calibrated shower features and truth matching remain separate extensions.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Verify which shower fields and associations exist in PDHD reconstruction.
