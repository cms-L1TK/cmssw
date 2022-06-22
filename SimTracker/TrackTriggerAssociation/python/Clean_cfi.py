import FWCore.ParameterSet.Config as cms

Clean_params = cms.PSet (
  InputTagTTClusterAssMap = cms.InputTag( "TTClusterAssociatorFromPixelDigis", "ClusterAccepted" ), #
  InputTagTVs             = cms.InputTag( "CleanTV",                           "ClusterAccepted" ), #
  InputTagTPs             = cms.InputTag( "CleanTP",                           "ClusterAccepted" ), #
  Branch                  = cms.string  ( "ClusterAccepted" )                                       # 
)