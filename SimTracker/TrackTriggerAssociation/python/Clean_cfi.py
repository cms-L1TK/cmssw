import FWCore.ParameterSet.Config as cms

Clean_params = cms.PSet (
  InputTagTTClusterAssMap = cms.InputTag( "TTClusterAssociatorFromPixelDigis", "ClusterAccepted"   ), #
  InputTagTVs             = cms.InputTag( "CleanTV",                           "AtLeastOneCluster" ), #
  InputTagTPs             = cms.InputTag( "CleanTP",                           "AtLeastOneCluster" ), #
  Branch                  = cms.string  ( "AtLeastOneCluster" )                                       # 
)