# configuration of Demonstrator. This is used to compare FW with SW for the subset fo the chain between LabelIn & LabelOut. FW must be wrapped by EMP & compiled with IPBB.
import FWCore.ParameterSet.Config as cms

TrackTriggerDemonstrator_params = cms.PSet (

  LabelIn  = cms.string( "TrackFindingTrackletProducerKFin"  ), #
  LabelOut = cms.string( "TrackFindingTrackletProducerKF" ), #
  DirIPBB  = cms.string( "/heplnw039/tschuh/work/proj/kfinkf/" ), # path to ipbb proj area
  RunTime  = cms.double( 8.0 )                                  # runtime in us

)