# this compares event by event the output of the C++ emulation with the ModelSim simulation of the firmware
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load( 'FWCore.MessageService.MessageLogger_cfi' )
process.load( 'Configuration.EventContent.EventContent_cff' )
process.load( 'Configuration.Geometry.GeometryExtendedRun4D110Reco_cff' ) 
process.load( 'Configuration.Geometry.GeometryExtendedRun4D110_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
process.load( 'L1Trigger.TrackFindingTracklet.Demonstrator_cff' )
from L1Trigger.TrackFindingTracklet.Customize_cff import *
#reducedConfig( process )
fwConfig( process )

# build schedule
process.TrackProcessorEmulation = cms.Sequence (  process.ProducerDTC
                                                + process.L1THybridTracks
                                                + process.ProducerTM
                                                + process.ProducerDR
                                                + process.ProducerKF
                                                + process.ProducerTQ
                                                + process.ProducerTFP
                                                )


Samples = ["/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/0f0bcfd3-dafe-4dda-8d39-9765f6eae68e.root"]

# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/4f8965fd-fda6-414c-bc3e-92598ba7b251.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/8eb00500-1562-4ef0-bec0-d95ce1ce3c02.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/6d36a59f-b984-493d-a750-dfeb6b2515c0.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/3be594c4-9067-4ee7-b2c7-39f8894328e2.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/82f7377e-b0a1-4fee-8ce4-4b45cbe23182.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/bbce3dfb-cbde-4a7d-a7d3-ea8eafedfcd7.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/e06e0c3d-b694-4a43-87f9-4c2a9bebaa17.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/33c9fb5a-93ea-4a29-abad-ec75900e7c55.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/74dcce3e-e636-4502-8dd9-64eecbb5e31d.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/612d219f-fc9c-4c05-b5f9-70649771f569.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/5ec71aac-41ad-4366-8229-f4b09aabd6a5.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/30f8edb9-b6a9-4596-a0ab-fccb07611910.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/a07eea4e-d869-41e9-8514-f1bce62fdbd2.root",
# "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/7e10bf8d-6b5d-4d60-bf44-3f7e6cc831a2.root"]

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring( Samples ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' ),
  skipEvents = cms.untracked.uint32( 0 ),
)
process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )
L1TRK_NAME  = process.TrackFindingTrackletAnalyzer_params.OutputLabelTFP.value()
L1TRK_LABEL = process.TrackFindingTrackletProducer_params.BranchTTTracks.value()
process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.l1tGTTInputProducer.l1TracksInputTag = cms.InputTag(L1TRK_NAME, L1TRK_LABEL)

process.load('L1Trigger.L1TTrackMatch.l1tTrackSelectionProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
process.load('L1Trigger.DemonstratorTools.l1tGTTFileWriter_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackVertexAssociationProducer_cfi')

# GTT File Writer
process.l1tGTTFileWriter.tracks = cms.untracked.InputTag(L1TRK_NAME, L1TRK_LABEL)
process.l1tGTTFileWriter.vertices = cms.untracked.InputTag("l1tVertexFinderEmulator", "L1VerticesEmulation")
process.l1tGTTFileWriter.selectedTracks = cms.untracked.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelectedEmulation")
process.l1tGTTFileWriter.vertexAssociatedTracks = cms.untracked.InputTag("l1tTrackVertexAssociationProducer", "Level1TTTracksSelectedAssociatedEmulation")
process.l1tGTTFileWriter.format = cms.untracked.string("EMPv2")

# Vertex Finder Emulator
process.l1tVertexFinderEmulator.VertexReconstruction.VxMinTrackPt = cms.double(0.0)
process.l1tVertexFinderEmulator.VertexReconstruction.Algorithm = cms.string("NNEmulation")

process.demo = cms.Path( process.TrackProcessorEmulation + 
                         process.l1tGTTInputProducer + 
                         process.l1tTrackSelectionProducer + 
                         process.l1tVertexFinderEmulator + 
                         process.l1tTrackVertexAssociationProducer +
                         process.l1tGTTFileWriter +
                         process.TrackerTFPDemonstrator)