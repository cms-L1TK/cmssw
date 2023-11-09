import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackerTFP.Analyzer_cfi import TrackerTFPAnalyzer_params
from L1Trigger.TrackerTFP.Producer_cfi import TrackerTFPProducer_params

TrackerTFPAnalyzerGP = cms.EDAnalyzer( 'trackerTFP::AnalyzerGP', TrackerTFPAnalyzer_params, TrackerTFPProducer_params )
TrackerTFPAnalyzerHT = cms.EDAnalyzer( 'trackerTFP::AnalyzerHT', TrackerTFPAnalyzer_params, TrackerTFPProducer_params )
TrackerTFPAnalyzerCTB = cms.EDAnalyzer( 'trackerTFP::AnalyzerCTB', TrackerTFPAnalyzer_params, TrackerTFPProducer_params )
TrackerTFPAnalyzerKF = cms.EDAnalyzer( 'trackerTFP::AnalyzerKF', TrackerTFPAnalyzer_params, TrackerTFPProducer_params )
TrackerTFPAnalyzerDR = cms.EDAnalyzer( 'trackerTFP::AnalyzerDR', TrackerTFPAnalyzer_params, TrackerTFPProducer_params )