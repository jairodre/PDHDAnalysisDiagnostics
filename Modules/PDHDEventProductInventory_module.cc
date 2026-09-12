/**
 * @file PDHDEventProductInventory_module.cc
 * @brief Template for one-row-per-event product, provenance, and quality output.
 *
 * Inputs: configurable expected InputTags plus EventAuxiliary/provenance.
 * Outputs: EventTree and a once-per-file configuration/provenance tree.
 * Implement run/subrun/event, isRealData, timestamp, product presence and size,
 * resolved producer/process, beam-instrumentation availability, timing/trigger
 * availability, and empty-event status. Required and optional products need
 * separate policies. This module must never infer a physical zero from a
 * missing product and should be implemented before all object inventories.
 *
 * Status: template-only contract; intentionally unbuilt, unscheduled, and
 * without an art plugin registration.
 */

// TODO: Implement EDAnalyzer after the common event key, status enumeration,
// configuration snapshot, and HD data/MC product lists are reviewed.
