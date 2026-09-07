/**
 * @file PDHDHitInventory_module.cc
 * @brief Template for compact reconstructed-hit and hit-classifier output.
 *
 * Inputs: recob::Hit plus optional wire, space-point, PFParticle/track/shower/blip,
 * CNN, and MC backtracking associations. Outputs: one row per selected hit with
 * channel, wire/plane/TPC/cryostat, peak/start/end tick, RMS, integral/amplitude,
 * multiplicity/local index, association keys, and separately named weighted and
 * unweighted CNN outputs. Optional MC columns store contribution provenance and
 * validity. Do not save raw waveforms here; provide configurable row limits or
 * object-based filtering to control output size without biasing default summaries.
 *
 * Status: structure-only; no art plugin is registered yet.
 */

// TODO: Define hit selection and CNN label configuration before enabling output.
