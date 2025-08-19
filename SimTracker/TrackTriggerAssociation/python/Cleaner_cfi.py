import FWCore.ParameterSet.Config as cms

Cleaner_params = cms.PSet (

  InputTag = cms.InputTag( "TTClusterAssociatorFromPixelDigis", "ClusterAccepted"   ), #
  Branch   = cms.string  ( "AtLeastOneCluster" )                                       # 

)
