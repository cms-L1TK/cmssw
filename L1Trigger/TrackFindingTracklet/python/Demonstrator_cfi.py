# configuration of Demonstrator. This is used to compare FW with SW for the subset fo the chain between LabelIn & LabelOut. FW must be wrapped by EMP & compiled with IPBB.
import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackFindingTracklet.Producer_cfi import TrackFindingTrackletProducer_params
from L1Trigger.TrackFindingTracklet.Analyzer_cfi import TrackFindingTrackletAnalyzer_params

# these parameters a for ModelSim runs of FW
TrackTriggerDemonstrator_params = cms.PSet (

  LabelIn  = TrackFindingTrackletAnalyzer_params.OutputLabelKF,  #
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelTQ, #
  DirIPBB  = cms.string( "/home/hep/am2023/tschuh_new_pr_nov_2026/CMSSW_20_1_0_pre1/src/" ), # path to ipbb proj area
  RunTime  = cms.double( 6. ),                                   # runtime in us

  LinkMappingIn  = cms.vint32(  ),
  LinkMappingOut = cms.vint32(  )

)
