/**
 * @file PDHDTrackInventory_module.cc
 * @brief Template for geometry-only reconstructed-track inventory.
 *
 * Inputs: configured recob::Track collections and PFP associations. Outputs:
 * stable track/PFP keys, trajectory-point counts and validity, raw start/end,
 * start/end directions, length, vertex/end covariance when available, containment
 * flags with explicit detector/fiducial margins, and orientation diagnostics.
 * Preserve original orientation; any beam-oriented or Z-oriented endpoints must
 * be additional derived columns. Calorimetry, PID, truth, and selection stay in
 * their dedicated modules so a geometry inventory works identically on data/MC.
 *
 * Status: template-only contract; intentionally unbuilt, unscheduled, and
 * without an art plugin registration.
 */

// TODO: Define geometry service boundaries and fiducial configuration in cm.
