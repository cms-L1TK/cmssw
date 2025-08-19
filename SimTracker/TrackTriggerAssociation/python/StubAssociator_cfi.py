import FWCore.ParameterSet.Config as cms

StubAssociator_params = cms.PSet (

  InputTagTTStubDetSetVec = cms.InputTag( "TTStubsFromPhase2TrackerDigis",     "StubAccepted"      ), #
  InputTagTTClusterAssMap = cms.InputTag( "TTClusterAssociatorFromPixelDigis", "ClusterAccepted"   ), #
  #InputTagTTClusterAssMap = cms.InputTag( "Cleaner", "AtLeastOneCluster" ),
  BranchReconstructable   = cms.string  ( "Reconstructable" ),                                        # name of StubAssociation collection made with reconstractable TPs
  BranchSelection         = cms.string  ( "UseForAlgEff"    ),                                        # name of StubAssociation collection used for tracking efficiency 

  MinPt           = cms.double(  2.   ), # pt cut in GeV
  MaxEta0         = cms.double(  2.4  ), # max eta for TP with z0 = 0
  MaxZ0           = cms.double( 15.   ), # half lumi region size in cm
  MaxD0           = cms.double(  5.   ), # cut on impact parameter in cm
  MaxVertR        = cms.double(  1.   ), # cut on vertex pos r in cm
  MaxVertZ        = cms.double( 30.   ), # cut on vertex pos z in cm

)
