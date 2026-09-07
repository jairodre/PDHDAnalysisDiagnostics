/**
 * @file DiagnosticRecords.h
 * @brief Shared row and association records for PDHD diagnostic ntuples.
 *
 * Populate this header only with stable, reusable records shared by more than
 * one analyzer. Every record must document units, validity flags, sentinel
 * policy, and stable event/object keys. Keep reconstructed, truth, calibrated,
 * uncalibrated, SCE, and no-SCE quantities distinct. ROOT-only branch buffers
 * that are private to one analyzer should remain in that analyzer instead.
 * If records become persistent art products, add the required dictionaries and
 * schema-evolution policy before enabling their producers.
 */

#ifndef PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_DIAGNOSTICRECORDS_H
#define PROTODUNEANA_PDHDANALYSISDIAGNOSTICS_DIAGNOSTICRECORDS_H

namespace pdhd::diagnostics {

// TODO: Add reviewed records incrementally. Prefer explicit validity members
// over physics-looking sentinel values and include run/subrun/event plus the
// relevant collection index in every independently joined table.

} // namespace pdhd::diagnostics

#endif
