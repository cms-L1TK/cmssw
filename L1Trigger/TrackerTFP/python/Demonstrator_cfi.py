# configuration of Demonstrator.
import FWCore.ParameterSet.Config as cms

TrackTriggerDemonstrator_params = cms.PSet (

  LabelIn  = cms.string( "TrackerTFPProducerHT"             ), #
  LabelOut = cms.string( "TrackerTFPProducerMHT"            ), #
  DirIPBB  = cms.string( "/heplnw039/tschuh/work/proj/mht/" ), # path to ipbb proj area
  RunTime  = cms.double( 4.5 )                                 # runtime in us

)