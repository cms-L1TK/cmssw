import FWCore.ParameterSet.Config as cms

#=== Import default values for all parameters & define EDProducer.

from L1Trigger.TrackerDTC.AnalyzerDTCStubs_cfi import TrackerDTC_params

ProducerDTC = cms.EDAnalyzer('trackerDTC::AnalyzerDTCStubs', TrackerDTC_params)
