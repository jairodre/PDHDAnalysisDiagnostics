/**
 * @file PDHDBlipAssociation_module.cc
 * @brief Template for blip-to-reconstructed-object and blip-to-truth associations.
 *
 * Inputs: blips/hits, reconstructed tracks/showers/PFPs, timing, and optional MC
 * hit backtracking. Outputs: candidate edges with shared-hit evidence where
 * available, spatial/endpoint distances, closest trajectory point, longitudinal
 * and transverse coordinates, association rank/method, ambiguity, and validity.
 * Track association must not be inferred solely from a plotting distance. MC
 * truth linkage is validation-only and unavailable on data by construction.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Separate producer associations from analysis-level geometric candidates.
