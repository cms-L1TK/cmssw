import FWCore.ParameterSet.Config as cms

process = cms.Process( "Demo" )
process.load( 'FWCore.MessageService.MessageLogger_cfi' )
process.load( 'Configuration.EventContent.EventContent_cff' )
process.load( 'Configuration.Geometry.GeometryExtendedRun4D121Reco_cff' ) 
process.load( 'Configuration.Geometry.GeometryExtendedRun4D121_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# load code that associates stubs with mctruth
process.load( 'SimTracker.TrackTriggerAssociation.StubAssociator_cff' )
# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
# L1 tracking => hybrid emulation 
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
# Load code that analyzes hybrid emulation 
process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
# load code that fits hybrid tracks
process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
from L1Trigger.TrackFindingTracklet.Customize_cff import *
fwConfig( process )

# build schedule
process.mc       = cms.Sequence( process.StubAssociator  )
process.dtc      = cms.Sequence( process.ProducerDTC     )
process.tracklet = cms.Sequence( process.L1THybridTracks )
process.tm       = cms.Sequence( process.ProducerTM      )
process.dr       = cms.Sequence( process.ProducerDR      )
process.kf       = cms.Sequence( process.ProducerKF      )
process.tq       = cms.Sequence( process.ProducerTQ      )
process.data     = cms.Sequence( process.NtupleMakerTQ   )
process.tt       = cms.Path( process.mc + process.dtc + process.tracklet + process.tm + process.dr + process.kf + process.tq + process.data )
process.schedule = cms.Schedule( process.tt )

# create options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing( 'analysis' )


Samples = [
  #"/store/relval/CMSSW_20_0_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1/2590000/0033230b-a131-453a-95c0-fe14d5027d1f.root",
  "file:/data/store/CMSSW_20_0_0_pre1_RelValTTbar_14TeV_GEN-SIM-DIGI-RAW_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1_2590000.root"
]

options.register( 'inputMC', Samples, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed" )
options.register( 'Events',-1,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Number of Events to analyze" )
options.parseArguments()

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )
process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring( options.inputMC ),
  #skipEvents = cms.untracked.uint32( 1 ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)
process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )
process.MessageLogger.cerr.enableStatistics = False
process.MessageLogger.L1track = dict(limit = -1)
process.TFileService = cms.Service( "TFileService", fileName = cms.string( "tq.root" ) )
