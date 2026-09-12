/**
 * @file PDHDEnergyDepositInventory_module.cc
 * @brief Template for MC energy-deposit/IDE truth records linked to G4 particles.
 *
 * Inputs: configured SimEnergyDeposit and/or SimChannel/IDE products plus the
 * particle inventory. Outputs should be one row per retained deposit or a clearly
 * documented aggregation, with G4 track ID, energy/electrons/photons, midpoint or
 * start/end position, time, TPC/plane/channel context when meaningful, active and
 * fiducial volume flags, and ancestry linkage. Keep deposited energy distinct
 * from kinetic energy and reconstructed calorimetric energy. Provide configurable
 * thinning/aggregation because unbounded deposit-level output can be very large.
 *
 * Status: template-only contract; MC-only, intentionally unbuilt, unscheduled,
 * and without an art plugin registration.
 */

// TODO: Choose the production-available truth product before fixing the schema.
