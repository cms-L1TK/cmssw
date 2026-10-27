# ESProducer providing the algorithm to run input data through modelsim and to compares results with expected output data
# and EDAnalyzer running the ESProduct produced by above ESProducer

import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackFindingTracklet.Demonstrator_cfi import TrackTriggerDemonstrator_params
from L1Trigger.TrackFindingTracklet.Producer_cfi import TrackFindingTrackletProducer_params
from L1Trigger.TrackFindingTracklet.Analyzer_cfi import TrackFindingTrackletAnalyzer_params

TrackFindingTrackletDemonstrator = cms.EDAnalyzer("trklet::Demonstrator", TrackTriggerDemonstrator_params, TrackFindingTrackletProducer_params)

TTDemoTM = TrackFindingTrackletDemonstrator.clone(
  LabelIn  = TrackFindingTrackletProducer_params.InputLabelTM,
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelTM,
  DirIPBB  = "/data/tschuh/work/proj/tbtm/"
)
TTDemoDR = TrackFindingTrackletDemonstrator.clone(
  LabelIn  = TrackFindingTrackletProducer_params.InputLabelDR,
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelDR,
  DirIPBB  = "/data/tschuh/work/proj/tmdr/"
)
TTDemoKF = TrackFindingTrackletDemonstrator.clone(
  LabelIn  = TrackFindingTrackletProducer_params.InputLabelKF,
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelKF,
  DirIPBB  = "/data/tschuh/work/proj/drkf/"
)
TTDemoTQ = TrackFindingTrackletDemonstrator.clone(
  LabelIn  = TrackFindingTrackletProducer_params.InputLabelTQ,
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelTQ,
  DirIPBB  = "/data/tschuh/work/proj/kftq/"
)
TTDemoTFP = TrackFindingTrackletDemonstrator.clone(
  LabelIn  = TrackFindingTrackletProducer_params.InputLabelTFP,
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelTFP,
  DirIPBB  = "/data/tschuh/work/proj/tqtfp/"
)
