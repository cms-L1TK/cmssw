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

# load code that cleans TP & TV collection and associates clean TPs with clusters
process.load( 'SimTracker.TrackTriggerAssociation.Cleaner_cff' )
# load code that created association maps of stubs with cleaned TPs
process.load( 'SimTracker.TrackTriggerAssociation.StubAssociator_cff' )
# load code that analyzes TPs
process.load( 'L1Trigger.TrackTrigger.AnalyzerMC_cff' )

#process.StubAssociatorC = process.StubAssociator.clone( InputTagTTClusterAssMap = cms.InputTag( "Cleaner", "AtLeastOneCluster" ) )

# build schedule
#process.path     = cms.Path( process.StubAssociator + process.Cleaner + process.StubAssociatorC + process.AnalyzerMC )
process.path     = cms.Path( process.AnalyzerMC )
process.schedule = cms.Schedule( process.path )

# create options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing( 'analysis' )
# specify input MC
Samples = [
  #"/store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2/130000/00c7f40e-b44e-4eea-a86b-def8f7d82b0e.root"
  "/store/trimmed.root"
]
options.register( 'inputMC', Samples, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed" )
# specify number of events to process.
options.register( 'Events',100,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Number of Events to analyze" )
options.parseArguments()

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )
process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring( options.inputMC ),
  #skipEvents = cms.untracked.uint32( 3537 ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)
