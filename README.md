# PDHD Analysis Diagnostics

This package is a scaffold for small, composable ProtoDUNE-HD diagnostic
analyzers. It separates event products, truth particles, reconstructed objects,
association edges, derived features, and selection decisions so each stage can
be tested without depending on a single monolithic ntuple module.

## Current status

The directory structure, implementation contracts, documentation catalog, and
base FHiCL tables exist. The reusable algorithms are implemented and registered
as `dune_PDHDAnalysisDiagnosticsAlg`. A preceding beam-instrumentation revision
compiled and ran on data; that test established that the selected reconstructed
files do not store `beamevent`. Separate data and MC jobs now handle this product
difference explicitly, but the revised sources still require a user build and
runtime validation. The two beam analyzer plugins are registered; all other
module templates intentionally perform no event processing. Add modules one at a
time only after their product labels, schemas, units, validity behavior, and
data/MC contract have been reviewed.

## Design rules

- Every table uses run, subrun, and event keys plus stable collection indices.
- Missing products are represented by availability/status fields, not physical
  zeros or silently substituted producers.
- Reconstructed observables define nominal data-like selections. MC truth is
  used to label and measure their performance, not to accept an event.
- Raw and derived quantities remain separate, including original/reoriented
  tracks, raw/T0-corrected coordinates, SCE/no-SCE calorimetry, and weighted/
  unweighted classifier scores.
- Object inventories preserve facts; association modules preserve graph edges;
  selection modules preserve all component decisions and final decisions.
- Official DUNE and LArSoft utilities are preferred after checking compatibility
  with the locally active release and PDHD production products.

See `DIAGNOSTICS_CATALOG.md` for the responsibility and implementation contract
of every planned analyzer.
