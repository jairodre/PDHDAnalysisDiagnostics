/**
 * @file PDHDShowerFeatureInventory_module.cc
 * @brief Template for calibrated shower and classifier feature output.
 *
 * Inputs: showers, PFPs, hits, shower-calorimetry products, detector conditions,
 * and optional CNN output. Outputs: raw and calibrated energy per plane, dE/dx,
 * charge, hit/space-point counts, conversion-gap/topology variables, and separate
 * CNN scores with source and validity. SCE/no-SCE and calibrated/uncalibrated
 * quantities remain independent. MC truth labels are attached only by the truth
 * matcher and never participate in reconstructed feature calculation.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Verify HD calibration support before defining calibrated energy branches.
