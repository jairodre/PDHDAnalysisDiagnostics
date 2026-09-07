# Shared calculations are built as one small library. These classes do not own
# ROOT trees or art modules; analyzers supply event products and persist results.
art_make_library(
  LIBRARY_NAME dune_PDHDAnalysisDiagnosticsAlg
  SOURCE
    ProductInventoryAlg.cxx
    GeometryContainmentAlg.cxx
    GeantHierarchyAlg.cxx
    RecoAssociationAlg.cxx
    RecoTruthMatchAlg.cxx
    CalorimetryAlg.cxx
    BeamInstrumentationAlg.cxx
    BeamSelectionAlg.cxx
    CosmicSelectionAlg.cxx
  LIBRARIES
    lardataobj::AnalysisBase
    dunecore::DuneObj
    ProtoDUNEUtilities
)

install_headers()
install_source()
