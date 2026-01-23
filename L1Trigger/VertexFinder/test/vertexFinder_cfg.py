import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os
process = cms.Process("L1TrackNtuple")

GEOMETRY = "D110"

## ------------------ L1 Track Finder ------------------ ##
L1TRKALGO = 'HYBRID_DISPLACED_NEWKF_KILL'

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
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

inputMC = ["/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/0f0bcfd3-dafe-4dda-8d39-9765f6eae68e.root"]
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*inputMC))
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))
process.TFileService = cms.Service("TFileService", fileName = cms.string('VertexNTupler.root'), closeFileFast = cms.untracked.bool(True))

## ------------------ L1 TRG Algorithms ------------------ ##

process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackSelectionProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tVertexNTupler_cfi')
process.load('L1Trigger.VertexFinder.l1tInputDataProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tTPStubValueMapProducer_cfi')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.p = cms.Path(process.l1tGTTInputProducer + 
                     process.l1tTrackSelectionProducer + 
                     process.l1tVertexProducer)

