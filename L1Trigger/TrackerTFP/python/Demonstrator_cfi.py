# configuration of Demonstrator.
import FWCore.ParameterSet.Config as cms

TrackTriggerDemonstrator_params = cms.PSet (

  LabelIn  = cms.string( "TrackerTFPProducerZHT"             ), #
  LabelOut = cms.string( "TrackerTFPProducerKFin"            ), #
  DirIPBB  = cms.string( "/heplnw039/tschuh/work/proj/zhtkfin/" ), # path to ipbb proj area
  RunTime  = cms.double( 4.5 )                                 # runtime in us

)