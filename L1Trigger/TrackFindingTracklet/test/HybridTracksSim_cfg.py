################################################################################################
# This script runs DTC + displaced tracklet + TM & DR & KF simulation
# allowing to identify problems quickly during developement.
# This script is a specialized and light-weight version of L1TrackNtupleMaker_cfg.py
# To run execute do
# cmsRun L1Trigger/TrackFindingTracklet/test/HybridTracksSim_cfg.py
# where the arguments take default values if you don't specify them. You can change defaults below.
#################################################################################################

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

# load code that associates stubs with mctruth
process.load( 'SimTracker.TrackTriggerAssociation.StubAssociator_cff' )
# load code that analyzes mc truth
process.load( 'L1Trigger.TrackTrigger.AnalyzerMC_cff' )
# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
# load code that analyzes DTCStubs
process.load( 'L1Trigger.TrackerDTC.Analyzer_cff' )
# L1 tracking => hybrid emulation 
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
#--- Load code that analyzes hybrid emulation 
process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
# load code that fits hybrid tracks
process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
from L1Trigger.TrackFindingTracklet.Customize_cff import *
simConfig( process )

process.AnalyzerTracklet.InputTag = ( "l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks" )

# build schedule
process.mc       = cms.Sequence( process.StubAssociator          + process.AnalyzerMC       )
process.dtc      = cms.Sequence( process.ProducerDTC             + process.AnalyzerDTC      )
process.tracklet = cms.Sequence( process.L1TExtendedHybridTracks + process.AnalyzerTracklet )
process.sim      = cms.Sequence( process.ProducerSim             + process.AnalyzerSim      )
process.tt       = cms.Path( process.mc + process.dtc + process.tracklet + process.sim )
process.schedule = cms.Schedule( process.tt )

# create options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing( 'analysis' )
# specify input MC
#from MCsamples.Scripts.getCMSdata_cfi import *
#from MCsamples.Scripts.getCMSlocaldata_cfi import *
#from MCsamples.RelVal_1260_D88.PU200_TTbar_14TeV_cfi import *
#inputMC = getCMSdataFromCards()
Samples = [
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/0f0bcfd3-dafe-4dda-8d39-9765f6eae68e.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/30f8edb9-b6a9-4596-a0ab-fccb07611910.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/33c9fb5a-93ea-4a29-abad-ec75900e7c55.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/3be594c4-9067-4ee7-b2c7-39f8894328e2.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/4f8965fd-fda6-414c-bc3e-92598ba7b251.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/5ec71aac-41ad-4366-8229-f4b09aabd6a5.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/612d219f-fc9c-4c05-b5f9-70649771f569.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/6d36a59f-b984-493d-a750-dfeb6b2515c0.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/74dcce3e-e636-4502-8dd9-64eecbb5e31d.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/7e10bf8d-6b5d-4d60-bf44-3f7e6cc831a2.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/82f7377e-b0a1-4fee-8ce4-4b45cbe23182.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/8eb00500-1562-4ef0-bec0-d95ce1ce3c02.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/a07eea4e-d869-41e9-8514-f1bce62fdbd2.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/abce35f4-2c95-4edd-aad3-1010eac9a426.root",
  ##"/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/e06e0c3d-b694-4a43-87f9-4c2a9bebaa17.root"
  #'/store/relval/CMSSW_15_1_0_pre5/RelValDoubleMuFlatPt1To100/GEN-SIM-DIGI-RAW/150X_mcRun4_realistic_v1_RV269_Run4D110_noPU-v1/2590000/076c038a-eb20-4be9-8bff-04af5917c436.root'
  "file:/data/store/CMSSW_15_1_0_pre5_RelValTTbar_14TeV_TuneCP5_GEN-SIM-DIGI-RAW_PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2_2590000.root"
  #"file:/data/store/CMSSW_15_1_0_pre5_RelValDoubleMuFlatPt1To100_GEN-SIM-DIGI-RAW_150X_mcRun4_realistic_v1_RV269_Run4D110_noPU-v1_2590000.root"
]
#Samples = ["/store/trimmed.root"]
options.register( 'inputMC', Samples, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed" )
# specify number of events to process.
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
process.TFileService = cms.Service( "TFileService", fileName = cms.string( "Hist.root" ) )

if ( False ):
  process.out = cms.OutputModule (
    "PoolOutputModule",
    fileName = cms.untracked.string("TrackletTracks.root"),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = process.RAWSIMEventContent.outputCommands,
  )
  process.out.outputCommands.append('drop  *')
  process.out.outputCommands.append('keep  *_TTStubsFromPhase2TrackerDigis_ClusterAccepted_*')
  process.out.outputCommands.append('keep  *_TTStubsFromPhase2TrackerDigis_StubAccepted_*')
  process.out.outputCommands.append('keep  *_Cleaner_*_*')
  process.out.outputCommands.append('keep  *_StubAssociator_*_*')
  process.out.outputCommands.append('keep  *_ProducerDTC_*_*')
  process.out.outputCommands.append('keep  *_l1tTTTracksFromTrackletEmulation_*_*')
  process.FEVToutput_step = cms.EndPath( process.out )
  process.schedule.append( process.FEVToutput_step )
