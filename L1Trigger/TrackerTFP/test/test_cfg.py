################################################################################################
# Run bit-accurate TMTT L1 tracking emulation. 
#
# To run execute do
# cmsRun L1Trigger/L1TTrackerTFP/test/test_cfg.py
# where the arguments take default values if you don't specify them. You can change defaults below.
#################################################################################################

import FWCore.ParameterSet.Config as cms

process = cms.Process( "Demo" )
process.load( 'FWCore.MessageService.MessageLogger_cfi' )
process.load( 'Configuration.Geometry.GeometryExtended2026D88Reco_cff' ) 
process.load( 'Configuration.Geometry.GeometryExtended2026D88_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# load code that associates stubs with mctruth
process.load( 'SimTracker.TrackTriggerAssociation.StubAssociator_cff' )
# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
# load code that analyzes DTCStubs
process.load( 'L1Trigger.TrackerDTC.Analyzer_cff' )
# cosutmize TT algorithm
from L1Trigger.TrackerDTC.Customize_cff import *
producerUseTMTT( process )
analyzerUseTMTT( process )
#--- Load code that produces tfp Stubs
process.load( 'L1Trigger.TrackerTFP.Producer_cff' )
from L1Trigger.TrackerTFP.Customize_cff import *
setupUseTMTT( process )
#--- Load code that analyzes tfp Stubs
process.load( 'L1Trigger.TrackerTFP.Analyzer_cff' )

# build schedule
process.mc = cms.Sequence( process.StubAssociator )
process.dtc = cms.Sequence( process.TrackerDTCProducer + process.TrackerDTCAnalyzer )
process.pp = cms.Sequence( process.TrackerTFPProducerPP )
process.gp = cms.Sequence( process.TrackerTFPProducerGP + process.TrackerTFPAnalyzerGP )
process.ht = cms.Sequence( process.TrackerTFPProducerHT + process.TrackerTFPAnalyzerHT )
process.ctb = cms.Sequence( process.TrackerTFPProducerCTB + process.TrackerTFPAnalyzerCTB )
process.kf = cms.Sequence( process.TrackerTFPProducerKF + process.TrackerTFPAnalyzerKF )
process.dr = cms.Sequence( process.TrackerTFPProducerDR + process.TrackerTFPAnalyzerDR )
process.tt = cms.Path( process.mc + process.dtc + process.pp + process.gp + process.ht + process.ctb + process.kf + process.dr )
process.schedule = cms.Schedule( process.tt )

# create options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing( 'analysis' )
# specify input MC
Samples = [
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/00b3d04b-4c7b-4506-8d82-9538fb21ee19.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0390df7b-7c2a-45d0-9bdb-e13f6565f65a.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0919ec05-0bdc-4f85-9a21-3a167117ea5e.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/09e5d64c-d0de-442c-943d-620927afc59c.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0ce69121-29ff-4f67-babf-c3372a273ce6.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0eccd3e4-cbe0-403e-a574-750ec48804fc.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1118e869-7f5d-4d34-98ae-43d09b52c437.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1967aa5c-4fb4-4039-ac75-2f6fcdc8d0b0.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/19fb8c07-8b8c-418f-aa6f-8bedc4f8c8c5.root',
  '/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1a364293-a1ab-42a0-bbda-afc6d404cb2e.root'
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
  #skipEvents = cms.untracked.uint32( 1 ),
  noEventSort = cms.untracked.bool( True ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' ),
)
process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )
process.MessageLogger.cerr.enableStatistics = False
process.TFileService = cms.Service( "TFileService", fileName = cms.string( "Hist.root" ) )
