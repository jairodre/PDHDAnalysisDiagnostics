/**
 * @file PDHDBlipInventory_module.cc
 * @brief Template for reconstructed blip geometry, charge, and energy inventory.
 *
 * Inputs: configured DUNE blip products plus hit/cluster and timing associations.
 * Outputs: stable blip key, position and uncertainty, TPC/volume, size/extent,
 * charge/electron count, reconstructed energy with calibration definition, hit
 * counts, quality flags, T0/X-correction source, and detector/fiducial containment.
 * Keep raw and corrected coordinates separate. Reconstructed blip energy must
 * never be labeled as truth deposited or particle kinetic energy.
 *
 * Status: template-only contract; intentionally unbuilt, unscheduled, and
 * without an art plugin registration.
 */

// TODO: Confirm the local Blip data-product API and calibration units first.
