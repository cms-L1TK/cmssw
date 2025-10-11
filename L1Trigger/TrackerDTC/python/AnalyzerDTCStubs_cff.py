import FWCore.ParameterSet.Config as cms

#=== Import default values for all parameters & define EDProducer.

from L1Trigger.TrackerDTC.AnalyzerDTCStubs_cfi import TrackerDTC_params
from L1Trigger.TrackTrigger.Setup_cff import TrackTriggerSetup

ProducerDTC = cms.EDAnalyzer('trackerDTC::AnalyzerDTCStubs', TrackerDTC_params)
