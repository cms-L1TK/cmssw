# configuration for hybrid track reconstruction chain analyzer

import FWCore.ParameterSet.Config as cms

TrackTriggerAnalyzerMC_params = cms.PSet (

  StubAssociation         = cms.string  ( "StubAssociator"  ),
  BranchReconstructable   = cms.string  ( "Reconstructable" ),
  BranchSelection         = cms.string  ( "UseForAlgEff"    ),
  InputTagTTStubs         = cms.InputTag( "TTStubsFromPhase2TrackerDigis", "StubAccepted"        ),
  InputTagTTClusterAssMap = cms.InputTag( "TTClusterAssociatorFromPixelDigis", "ClusterAccepted" )

)
