# configuration for hybrid track reconstruction chain analyzer

import FWCore.ParameterSet.Config as cms

TrackTriggerAnalyzerMC_params = cms.PSet (

  StubAssociation = cms.string  ( "StubAssociator"  ),
  BranchFake      = cms.string  ( "UseForFake"      ),
  BranchEff       = cms.string  ( "UseForEff"       ),
  InputTagTTStubs = cms.InputTag( "TTStubsFromPhase2TrackerDigis", "StubAccepted" ),
  InputTagMC      = cms.InputTag( "Cleaner", "AtLeastOneCluster" )

)
