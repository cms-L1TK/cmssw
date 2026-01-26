import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

## ------------------ Generics ------------------ ##

process = cms.Process("L1TrackNtuple")

GEOMETRY = "D110"

inputMC = [
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/4f8965fd-fda6-414c-bc3e-92598ba7b251.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/8eb00500-1562-4ef0-bec0-d95ce1ce3c02.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/6d36a59f-b984-493d-a750-dfeb6b2515c0.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/3be594c4-9067-4ee7-b2c7-39f8894328e2.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/82f7377e-b0a1-4fee-8ce4-4b45cbe23182.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/abce35f4-2c95-4edd-aad3-1010eac9a426.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/0f0bcfd3-dafe-4dda-8d39-9765f6eae68e.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/e06e0c3d-b694-4a43-87f9-4c2a9bebaa17.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/33c9fb5a-93ea-4a29-abad-ec75900e7c55.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/74dcce3e-e636-4502-8dd9-64eecbb5e31d.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/612d219f-fc9c-4c05-b5f9-70649771f569.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/5ec71aac-41ad-4366-8229-f4b09aabd6a5.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/30f8edb9-b6a9-4596-a0ab-fccb07611910.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/a07eea4e-d869-41e9-8514-f1bce62fdbd2.root",
    "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/7e10bf8d-6b5d-4d60-bf44-3f7e6cc831a2.root",
]

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*inputMC))
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))
process.TFileService = cms.Service("TFileService", fileName = cms.string('VertexNTupler.root'), closeFileFast = cms.untracked.bool(True))

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
# Drop previously reconstructed L1 tracks + 
# their truth association to avoid risk of analysing 
# them instead of new tracks created by this job.
process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
process.source.inputCommands = cms.untracked.vstring()
process.source.inputCommands.append('keep *_*_*_*')
process.source.inputCommands.append('drop  *_*_*Level1TTTracks*_*')

## ------------------ L1 Track Finder ------------------ ##
L1TRKALGO = 'HYBRID_NEWKF'  # HYBRID, HYBRID_NEWKF, HYBRID_REDUCED

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.StubAssociator_cff')
process.load('L1Trigger.TrackTrigger.AnalyzerMC_cff' )
process.load('L1Trigger.TrackerDTC.DTC_cff')
process.load('L1Trigger.TrackerDTC.Analyzer_cff')

process.dtc = cms.Path(process.StubAssociator + 
                       process.AnalyzerMC + 
                       process.ProducerDTC + 
                       process.AnalyzerDTC)
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")

if (L1TRKALGO == 'HYBRID'):
    process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
    process.TTTracksEmulation = cms.Path(process.L1THybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1THybridTracksWithAssociators + process.AnalyzerTracklet)
    NHELIXPAR = 4
    L1TRK_NAME  = "l1tTTTracksFromTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"
elif (L1TRKALGO == 'HYBRID_NEWKF' or L1TRKALGO == 'HYBRID_REDUCED'):
    process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
    process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
    NHELIXPAR = 4
    L1TRK_NAME  = process.TrackFindingTrackletAnalyzer_params.OutputLabelTFP.value()
    L1TRK_LABEL = process.TrackFindingTrackletProducer_params.BranchTTTracks.value()
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"
    process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag(L1TRK_NAME, L1TRK_LABEL) )
    process.HybridNewKF = cms.Sequence(process.L1THybridTracks + process.ProducerTM + process.ProducerDR + process.ProducerKF + process.ProducerTQ + process.ProducerTFP)
    process.TTTracksEmulation = cms.Path(process.HybridNewKF)
    #process.TTTracksEmulationWithTruth = cms.Path(process.HybridNewKF +  process.TrackTriggerAssociatorTracks)
    # Optionally include code producing performance plots & end-of-job summary.
    process.TTTracksEmulationWithTruth = cms.Path(process.HybridNewKF +  process.TrackTriggerAssociatorTracks +  process.AnalyzerTracklet + process.AnalyzerTM + process.AnalyzerDR + process.AnalyzerKF + process.AnalyzerTFP )
    from L1Trigger.TrackFindingTracklet.Customize_cff import *
    if (L1TRKALGO == 'HYBRID_NEWKF'):
        fwConfig( process )
        # cheats to get good performance
        process.TrackTriggerSetup.KalmanFilter.UseTTStubResiduals = True
        process.TrackTriggerSetup.KalmanFilter.UseTTStubParameters = True
    if (L1TRKALGO == 'HYBRID_REDUCED'):
        reducedConfig( process )
    # Needed by L1TrackNtupleMaker
    process.HitPatternHelperSetup.useNewKF = True


## ------------------ L1 TRG Algorithms ------------------ ##

process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.l1tGTTInputProducer.l1TracksInputTag = cms.InputTag(L1TRK_NAME, L1TRK_LABEL)

# Making Sure Vertex Finder Uses New KF Tracks
tag = process.l1tGTTInputProducer.l1TracksInputTag
print(f"\033[1;32mSelected L1 Track Trigger Algorithm -> {L1TRKALGO}\033[0m")
print("\033[1;32m=== TTracks Constructed By Kalman Filter ===\033[0m")
print("Module :", tag.getModuleLabel())
print("Product:", tag.getProductInstanceLabel())
process.load('L1Trigger.L1TTrackMatch.l1tTrackSelectionProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackVertexAssociationProducer_cfi')
process.load('L1Trigger.VertexFinder.l1TrackZ0ResolutionProducer_cfi')

# Ignore 
# process.load('L1Trigger.VertexFinder.l1tVertexNTupler_cfi')
# process.load('L1Trigger.VertexFinder.l1tInputDataProducer_cfi')
# process.load('L1Trigger.VertexFinder.l1tTPStubValueMapProducer_cfi')

process.load('L1Trigger.VertexFinder.l1VertexAnalyzer_cfi')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.p = cms.Path(process.l1tGTTInputProducer + 
                     process.l1tTrackSelectionProducer + 
                     process.l1tTrackVertexAssociationProducer +
                     process.l1tVertexProducer + 
                     process.l1VertexAnalyzer + 
                     process.l1TrackZ0NtupleProducer)
