import FWCore.ParameterSet.Config as cms

process = cms.Process("TestDTCStubAnalyzer")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff') 
process.load('Configuration.Geometry.GeometryExtendedRun4D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load('L1Trigger.TrackerDTC.AnalyzerDTCStubs_cff')

Sample = ["/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV_TuneCP5/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_RV269_Run4D110_PU-v2/2590000/0f0bcfd3-dafe-4dda-8d39-9765f6eae68e.root",
          "/store/relval/CMSSW_15_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2590000/cb282e49-ab69-4449-8d93-c465ecffe707.root",
          "/store/relval/CMSSW_15_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2590000/ddee2a64-84a8-4ca7-b4f1-5df00da1b669.root",
          "/store/relval/CMSSW_15_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2590000/37bac9a7-5c2c-4d29-9bc9-2ed5d429bdda.root"]

process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring(Sample)
)

process.run = cms.Path(process.ProducerDTC)

process.TFileService = cms.Service( "TFileService", fileName = cms.string( "OT_DTC_Stubs.root" ) )