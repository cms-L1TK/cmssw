# N.B., DUE TO THE CHANGE IN STUB WINDOW SIZES WITH CMSSW 14_2_0_PRE2, THIS JOB HAS BEEN NODIFIED TO
# RECREATE THE STUBS, WHICH IS NECESSARY WHEN RUNNING ON MONTE CARLO GENERATED WITH OLDER VERSIONS.

############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os
process = cms.Process("L1TrackNtuple")

############################################################
# edit options here
############################################################

# D88 was used for CMSSW_12_6 datasets, and D98 recommended for more recent ones.
#GEOMETRY = "D88"
GEOMETRY = "D98"

WRITE_DATA = False

############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.L1track = dict(limit = -1)
process.MessageLogger.Tracklet = dict(limit = -1)
process.MessageLogger.TrackTriggerHPH = dict(limit = -1)

print("using geometry " + GEOMETRY + " (tilted)")

process.load('Configuration.Geometry.GeometryExtendedRun4' + GEOMETRY + 'Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4' + GEOMETRY +'_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
# Change needed to run with D98 geometry in recent CMSSW versions.
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun4_realistic_v1', '')

############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

inputMC = ["/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/0b2b0b0b-f312-48a8-9d46-ccbadc69bbfd.root"]

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*inputMC))
process.TFileService = cms.Service("TFileService", fileName = cms.string('L1TrkNtuple.root'), closeFileFast = cms.untracked.bool(True))
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

############################################################
# L1 tracking: stubs / DTC emulation
############################################################

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')

# remake stubs?
#from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
#process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")

#from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
#TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

#process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
#process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)


# load code that associates stubs with mctruth
process.load( 'SimTracker.TrackTriggerAssociation.StubAssociator_cff' )
# DTC emulation
process.load('L1Trigger.TrackerDTC.DTC_cff')

# load code that analyzes DTCStubs
process.load('L1Trigger.TrackerDTC.Analyzer_cff')

# modify default cuts
#process.TrackTriggerSetup.FrontEnd.BendCut = 5.0
#process.TrackTriggerSetup.Hybrid.MinPt = 1.0

process.dtc = cms.Path(process.StubAssociator + process.ProducerDTC + process.AnalyzerDTC)

############################################################
# L1 tracking
############################################################

process.load ("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
process.load ('L1Trigger.TrackFindingTracklet.Producer_cff')
process.load ('L1Trigger.TrackFindingTracklet.Analyzer_cff')
L1TRK_NAME  = process.TrackFindingTrackletAnalyzer_params.OutputLabelTFP.value()
L1TRK_LABEL = process.TrackFindingTrackletProducer_params.BranchTTTracks.value()
L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"

process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag(cms.InputTag(L1TRK_NAME, L1TRK_LABEL))
process.HybridNewKF = cms.Sequence( process.L1THybridTracks + 
                                    process.ProducerTM      + 
                                    process.ProducerDR      + 
                                    process.ProducerKF      + 
                                    process.ProducerTQ      + 
                                    process.ProducerTFP)

process.TTTracksEmulation = cms.Path(process.HybridNewKF)

process.load('SimTracker.TrackTriggerAssociation.StubAssociator_cff')
process.TTTracksEmulationWithTruth = cms.Path( process.HybridNewKF                   +  
                                               process.TrackTriggerAssociatorTracks  + 
                                               process.StubAssociator                +  
                                               process.AnalyzerTracklet              + 
                                               process.AnalyzerTM                    + 
                                               process.AnalyzerDR                    + 
                                               process.AnalyzerKF                    + 
                                               process.AnalyzerTQ                    + 
                                               process.AnalyzerTFP)

from L1Trigger.TrackFindingTracklet.Customize_cff import *
fwConfig(process)
# Needed by L1TrackNtupleMaker
process.HitPatternHelperSetup.useNewKF = True
