# configuration of Demonstrator. This is used to compare FW with SW for the subset fo the chain between LabelIn & LabelOut. FW must be wrapped by EMP & compiled with IPBB.
import FWCore.ParameterSet.Config as cms

# these parameters a for ModelSim runs of FW
TrackTriggerDemonstrator_params = cms.PSet (

  LabelIn  = cms.string( "TrackFindingTrackletProducerTM"    ), #
  LabelOut = cms.string( "TrackFindingTrackletProducerDR"      ), #
  DirIPBB  = cms.string( "/heplnw039/tschuh/work/proj/tmdr/" ), # path to ipbb proj area
  RunTime  = cms.double( 4.5 ),                                   # runtime in us

  LinkMappingIn  = cms.vint32(  ),
  LinkMappingOut = cms.vint32(  )

)