# this compares event by event the output of the C++ emulation with the ModelSim simulation of the firmware
import FWCore.ParameterSet.Config as cms

process = cms.Process( "Demo" )
process.load( 'FWCore.MessageService.MessageLogger_cfi' )
process.load( 'Configuration.EventContent.EventContent_cff' )
process.load( 'Configuration.Geometry.GeometryExtendedRun4D110Reco_cff' ) 
process.load( 'Configuration.Geometry.GeometryExtendedRun4D110_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
# L1 tracking => hybrid emulation 
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
# load code that fits hybrid tracks
process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
#--- Load code that compares s/w with f/w
process.load( 'L1Trigger.TrackFindingTracklet.Demonstrator_cff' )
from L1Trigger.TrackFindingTracklet.Customize_cff import *
#reducedConfig( process )
fwConfig( process )

# build schedule
process.emu = cms.Sequence (  process.ProducerDTC
                            + process.L1THybridTracks
                            + process.ProducerTM
                            + process.ProducerDR
                            + process.ProducerKF
                            + process.ProducerTQ
                            + process.ProducerTFP
                           )


Samples = [
  "/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/0f0bcfd3-dafe-4dda-8d39-9765f6eae68e.root"
]

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
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

process.l1tGTTFileWriter.tracks = cms.untracked.InputTag(L1TRK_NAME, L1TRK_LABEL)
process.l1tGTTFileWriter.vertices = cms.untracked.InputTag("l1tVertexFinderEmulator", "L1VerticesEmulation")
process.l1tGTTFileWriter.selectedTracks = cms.untracked.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelectedEmulation")
process.l1tGTTFileWriter.vertexAssociatedTracks = cms.untracked.InputTag("l1tTrackVertexAssociationProducer", "Level1TTTracksSelectedAssociatedEmulation")
process.l1tGTTFileWriter.format = cms.untracked.string("EMPv2")
process.l1tVertexFinderEmulator.VertexReconstruction.VxMinTrackPt = cms.double(0.0)

process.demo = cms.Path(process.emu + 
                        #process.TrackerTFPDemonstrator + 
                        process.l1tGTTInputProducer + 
                        process.l1tTrackSelectionProducer + 
                        process.l1tVertexFinderEmulator + 
                        process.l1tTrackVertexAssociationProducer +
                        process.l1tGTTFileWriter)