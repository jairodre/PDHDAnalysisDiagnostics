/**
 * @file PDHDRecoAssociation_module.cc
 * @brief Template for normalized edges among reconstructed object collections.
 *
 * Inputs: configured PFP, track, shower, cluster, hit, space-point, vertex, T0,
 * calorimetry, PID, and blip association products. Outputs: one edge per actual
 * association with event key, source/target collection identities and indices,
 * association type, metadata where present, and validity/ambiguity information.
 * This provides lossless joins between modular trees and prevents each analyzer
 * from independently choosing a first associated object. No selection is applied.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Define compact typed edge tables rather than opaque string-only edges.
