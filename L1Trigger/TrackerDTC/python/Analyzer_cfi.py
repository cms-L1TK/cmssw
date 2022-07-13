import FWCore.ParameterSet.Config as cms

TrackerDTCAnalyzer_params = cms.PSet (

  InputTagAccepted        = cms.InputTag( "TrackerDTCProducer", "StubAccepted"     ), # dtc passed stubs selection
  InputTagLost            = cms.InputTag( "TrackerDTCProducer", "StubLost"         ), # dtc lost stubs selection
  InputTagReconstructable = cms.InputTag( "StubAssociator",     "Reconstructable"  ), #
  InputTagSelection       = cms.InputTag( "StubAssociator",     "UseForAlgEff"     ), #
  UseMCTruth              = cms.bool( True )                                          # eneables analyze of TPs                                                    # eneables analyze of TPs

)
