################################################################################################
# Run bit-accurate TMTT L1 tracking emulation. 
#
# To run execute do
# cmsRun L1Trigger/L1TTrackerTFP/test/test_cfg.py
# where the arguments take default values if you don't specify them. You can change defaults below.
#################################################################################################

import FWCore.ParameterSet.Config as cms

process = cms.Process( "CleanUp" )
process.load( 'FWCore.MessageService.MessageLogger_cfi' )
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'Configuration.EventContent.EventContent_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# load code that creates clean TVs
process.load( 'SimTracker.TrackTriggerAssociation.CleanTV_cff' )
# load code that creates clean TPs
process.load( 'SimTracker.TrackTriggerAssociation.CleanTP_cff' )
# load code that associates TTStubs with clean TPs
process.load( 'SimTracker.TrackTriggerAssociation.CleanAssoc_cff' )

# build schedule
process.clean = cms.Path( process.CleanTV + process.CleanTP + process.CleanAssoc )
process.schedule = cms.Schedule( process.clean )

# create options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing( 'analysis' )
# specify input MC
Samples = [
  '/store/relvalOrig/CMSSW_12_6_0_pre4/RelValSingleMuPt100/GEN-SIM-DIGI-RAW/125X_mcRun4_realistic_v2_2026D88noPU-v1/2580000/b46f31f3-ba67-421f-8ad8-c9ff53f487ed.root'
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
  #skipEvents = cms.untracked.uint32( 64+16+4+2 ),
  noEventSort = cms.untracked.bool( True ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)
#process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )
process.MessageLogger.cerr.enableStatistics = False

if True :
  process.writeDataset = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/heplnw039/store/relvalTrimmed/CMSSW_12_6_0_pre4/RelValSingleMuPt100/GEN-SIM-DIGI-RAW/125X_mcRun4_realistic_v2_2026D88noPU-v1/2580000/b46f31f3-ba67-421f-8ad8-c9ff53f487ed.root'),
    dataset = cms.untracked.PSet(
      filterName = cms.untracked.string(''),
      dataTier = cms.untracked.string('GEN-SIM')
    )
  )
  process.writeDataset.outputCommands.append('drop  *_*_*_*')
  process.writeDataset.outputCommands.append('keep  *_CleanTV_AtLeastOneCluster_*')
  process.writeDataset.outputCommands.append('keep  *_CleanTP_AtLeastOneCluster_*')
  process.writeDataset.outputCommands.append('keep  *_TTStubsFromPhase2TrackerDigis_ClusterAccepted_*')
  process.writeDataset.outputCommands.append('keep  *_TTStubsFromPhase2TrackerDigis_StubAccepted_*')
  process.writeDataset.outputCommands.append('keep  *_CleanAssoc_AtLeastOneCluster_*')

  process.pd = cms.EndPath(process.writeDataset)
  process.schedule.append(process.pd)