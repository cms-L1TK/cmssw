#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_1_0_pre4/bin/el9_amd64_gcc12/cmsRun
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.Utilities.FileUtils as FileUtils
import os

## ------------------ Generics ------------------ ##

def get_input_mc_line(dataset_database, line_number):
    with open(dataset_database, 'r') as file:
        lines = file.readlines()
        if line_number < 0 or line_number >= len(lines):
            raise IndexError("Line number out of range")
        return lines[line_number].strip()

# Set up VarParsing options
options = VarParsing.VarParsing('analysis')

# Add custom command-line arguments
options.register('cluster',
                 0, # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Cluster ID from HTCondor")

options.register('process',
                 0, # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Process ID from HTCondor")

# Parse command-line arguments
options.parseArguments()

# Print the ClusterID and ProcessID
print(f'~ Cluster ID: {options.cluster}')
print(f'~ Process ID: {options.process}')

DatasetDatabase = "/home/hep/am2023/ian_fixes_for_14/CMSSW_15_1_0_pre4/src/RelValTTbar_14TeV_CMSSW_15_1_0_pre5_GEN-SIM-DIGI-RAW_files.txt"

try:
    InputMC = [get_input_mc_line(DatasetDatabase, options.process)]
except Exception as e:
    print(f"Error: {e}")
    InputMC = []

for i in range(len(InputMC)):
    print(f"InputMC[{i}]: {InputMC[i]}")

process = cms.Process("L1TrackNtuple")

GEOMETRY = "D110"

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*InputMC))
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))
process.TFileService = cms.Service("TFileService", fileName = cms.string('VertexNTupler_' + str(options.process) + '.root'), closeFileFast = cms.untracked.bool(True))

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
L1TRKALGO = 'HYBRID_NEWKF'

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
#process.load('L1Trigger.VertexFinder.l1TrackZ0ResolutionAnalyzer_cfi')
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

