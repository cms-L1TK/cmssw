import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackerDTC.Analyzer_cfi import TrackerDTCAnalyzer_params
from L1Trigger.TrackerDTC.DTC_cfi import TrackerDTC_params

TrackerDTCAnalyzer = cms.EDAnalyzer('trackerDTC::Analyzer', TrackerDTCAnalyzer_params, TrackerDTC_params)