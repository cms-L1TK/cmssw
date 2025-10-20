# configuration of Demonstrator. This is used to compare FW with SW for the subset fo the chain between LabelIn & LabelOut. FW must be wrapped by EMP & compiled with IPBB.
import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackFindingTracklet.Producer_cfi import TrackFindingTrackletProducer_params
from L1Trigger.TrackFindingTracklet.Analyzer_cfi import TrackFindingTrackletAnalyzer_params

# these parameters a for ModelSim runs of FW
TrackTriggerDemonstrator_params = cms.PSet (

  LabelIn  = TrackFindingTrackletProducer_params.InputLabelKF,  #
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelKF, #
  DirIPBB  = cms.string( "/home/hep/am2023/CMSSW_15_1_0_pre4/src/L1Trigger/TrackFindingTracklet/data/" ), # path to ipbb proj area
  RunTime  = cms.double( 7.5 ),                                   # runtime in us

  LinkMappingIn  = cms.vint32(  ),
  LinkMappingOut = cms.vint32(  )

)
