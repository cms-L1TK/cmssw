# configuration for AnalyzerDTC

import FWCore.ParameterSet.Config as cms

TrackerDTCAnalyzer_params = cms.PSet (

  UseMCTruth              = cms.bool( True ),                                    # eneables analyze of TPs
  InputTagReconstructable = cms.InputTag( "StubAssociator", "Reconstructable" ), #
  InputTagSelection       = cms.InputTag( "StubAssociator", "UseForAlgEff"    ), #
  InputTagDTC             = cms.InputTag( "ProducerDTC",    "StubAccepted"    )  # DTC stubs

)
