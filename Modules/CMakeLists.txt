# Only the two implemented beam analyzers are compiled and registered. The
# remaining source files are template-only contracts and are not discovered
# implicitly.
set(BEAM_MODULE_LIBRARIES
  dune_PDHDAnalysisDiagnosticsAlg
  ProtoDUNEUtilities
  dunecore::DuneObj
  nusimdata::SimulationBase
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Tree
)
cet_build_plugin(PDHDBeamInstrumentation art::Module
  SOURCE PDHDBeamInstrumentation_module.cc LIBRARIES ${BEAM_MODULE_LIBRARIES})
cet_build_plugin(PDHDBeamSelectionStages art::Module
  SOURCE PDHDBeamSelectionStages_module.cc LIBRARIES ${BEAM_MODULE_LIBRARIES})
install_source()
