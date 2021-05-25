import FWCore.ParameterSet.Config as cms

StubAssociator_params = cms.PSet (
  UseTTStubAssMap         = cms.bool    ( False ),                                                   #
  InputTagTTStubDetSetVec = cms.InputTag( "TTStubsFromPhase2TrackerDigis",     "StubAccepted"    ), #
  InputTagTTClusterAssMap = cms.InputTag( "TTClusterAssociatorFromPixelDigis", "ClusterAccepted" ), #
  InputTagTTStubAssMap    = cms.InputTag( "TTStubAssociatorFromPixelDigis",    "StubAccepted"    ), #
  BranchReconstructable   = cms.string  ( "Reconstructable" ),                                      #
  BranchSelection         = cms.string  ( "UseForAlgEff"    )                                       #
)