# this compares event by event the output of the C++ emulation with the ModelSim simulation of the firmware
import FWCore.ParameterSet.Config as cms

process = cms.Process( "Demo" )
process.load( 'Configuration.Geometry.GeometryExtended2026D49Reco_cff' )
process.load( 'Configuration.Geometry.GeometryExtended2026D49_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'Configuration.StandardSequences.L1TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag( process.GlobalTag, 'auto:phase2_realistic', '' )

# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.ProducerED_cff' )
# cosutmize TT algorithm
from L1Trigger.TrackerDTC.Customize_cff import *
producerUseTMTT(process)
analyzerUseTMTT(process)
#--- Load code that produces tfp Stubs
process.load( 'L1Trigger.TrackerTFP.Producer_cff' )
#--- Load code that demonstrates tfp Stubs
process.load( 'L1Trigger.TrackerTFP.Demonstrator_cff' )

# build schedule
process.tt = cms.Sequence (  process.TrackerDTCProducer
                           + process.TrackerTFPProducerGP
                           + process.TrackerTFPProducerHT
                           + process.TrackerTFPProducerMHT
                           + process.TrackerTFPProducerSF
                           + process.TrackerTFPProducerSFout
                           + process.TrackerTFPProducerKFin
                           + process.TrackerTFPProducerKF
                          )
process.demo = cms.Path( process.tt + process.TrackerTFPDemonstrator )
process.schedule = cms.Schedule( process.demo )

# create options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing( 'analysis' )
# specify input MC
Samples = {
  #'/store/relval/CMSSW_11_2_0_pre6/RelValSingleMuFlatPt1p5To8/GEN-SIM-DIGI-RAW/112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/20000/4E15C795-152F-A040-8AF4-5AF5F97EB996.root'
  #'/store/relval/CMSSW_11_2_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/20000/05D370F9-B42A-9C40-A1D8-BDD261A43A0D.root'
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/0074C44A-BBE2-6849-965D-CB73FE0C0E6C.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/00BE971B-A866-B34C-9EE3-48EF12C78C8D.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/011CEE57-5477-AE4D-A91F-4853259170CE.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/043F7E57-B485-7344-93A8-CE1C0BF92AF0.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/05C7E3C6-19A5-BB4C-84AF-A6F6D0DEFAE2.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/087810CE-EBBB-CF47-A9D7-23C96F1E77BE.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/13BBEBBC-BB36-8A46-8E7E-AB69FFB63CB5.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/15865F8B-C3F5-594E-8E21-F4A330A96FA9.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/187113F3-C0B0-514C-AFFD-28E7DCA3E524.root',
  '/store/relval/CMSSW_11_2_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/1A0446C6-74F8-9749-BC4A-1D3A296CA4BA.root'
}
options.register( 'inputMC', Samples, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed" )
# specify number of events to process.
options.register( 'Events',100,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Number of Events to analyze" )
options.parseArguments()

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )
process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring( options.inputMC ),
  #skipEvents = cms.untracked.uint32( 914 ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)
process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )