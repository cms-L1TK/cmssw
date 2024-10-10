################################################################################################
# This script runs DTC + prompt tracklet + KF interface + new KF emulator with analyzer for each step
# allowing to identify problems quickly during developement.
# This script is a specialized and light-weight version of L1TrackNtupleMaker_cfg.py
# To run execute do
# cmsRun L1Trigger/TrackFindingTracklet/test/HybridTracksNewKF_cfg.py
# where the arguments take default values if you don't specify them. You can change defaults below.
#################################################################################################

import FWCore.ParameterSet.Config as cms

process = cms.Process( "Demo" )
process.load( 'FWCore.MessageService.MessageLogger_cfi' )
process.load( 'Configuration.EventContent.EventContent_cff' )
process.load( 'Configuration.Geometry.GeometryExtended2026D98Reco_cff' ) 
process.load( 'Configuration.Geometry.GeometryExtended2026D98_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# load code that associates stubs with mctruth
process.load( 'SimTracker.TrackTriggerAssociation.StubAssociator_cff' )
# load code that analyzes TPs
process.load( 'SimTracker.TrackTriggerAssociation.Analyzer_cff' )
# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
# load code that analyzes DTCStubs
process.load( 'L1Trigger.TrackerDTC.Analyzer_cff' )
# L1 tracking => hybrid emulation 
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
from L1Trigger.TrackFindingTracklet.Customize_cff import *
fwConfig( process )
process.l1tTTTracksFromTrackletEmulation.readMoreMcTruth = False
#--- Load code that analyzes hybrid emulation 
process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
# load code that fits hybrid tracks
process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )

# build schedule
process.mc       = cms.Sequence( process.StubAssociator  + process.AnalyzerSA       )
process.dtc      = cms.Sequence( process.ProducerDTC     + process.AnalyzerDTC      )
process.tracklet = cms.Sequence( process.L1THybridTracks + process.AnalyzerTracklet )
process.tm       = cms.Sequence( process.ProducerTM      + process.AnalyzerTM       )
process.dr       = cms.Sequence( process.ProducerDR      + process.AnalyzerDR       )
process.kf       = cms.Sequence( process.ProducerKF      + process.AnalyzerKF       )
process.tq       = cms.Sequence( process.ProducerTQ                                 )
process.tfp      = cms.Sequence( process.ProducerTFP     + process.AnalyzerTFP      )
process.tt       = cms.Path( process.mc + process.dtc + process.tracklet + process.tm + process.dr + process.kf + process.tq + process.tfp )
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
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/00b68219-8585-406f-88d0-84da05a13280.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/03249c40-1c3e-41b2-925b-f8e3c14e78fc.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/1320c08c-5d74-4390-bbd6-51b29c09985f.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/219cfad7-744d-41c7-807e-e06acc5adac8.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/38cdeb1a-2eff-4734-8bc1-926883b7d1ff.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/6cd2e735-60ce-4033-813f-9a97b2aba581.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/8286cc60-6061-47ef-b919-e8ebabe0bec7.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/8f3bf9a1-43d0-4a03-8b07-e1ea38ce53fe.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/9154e90c-9cd0-4c28-8f27-9a1be278d5bd.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/9510ed6d-491e-4e62-a37f-7e020070a269.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/ac71e3c7-6c1d-4c34-bb0a-091dcbd02081.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/b4373015-0c29-4ccd-93d2-a744f37e58b0.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/b4d49071-2547-40f4-9324-bda51f310d95.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/bd46f473-1f38-4aaa-aef8-5c28e85e0d05.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/c07aeb69-60bf-4c2e-8657-ed1273c646c2.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/c33b6ffe-a5ec-4231-a49a-b23ded34cce7.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/e71d2389-bf11-4dd2-b6af-db478e04912e.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/e9e5f3b1-6a35-46b2-979f-43de94e23309.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/f8ced16d-dc52-4c53-a699-c0052aa269c9.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/08059b66-ac7c-4a79-be7a-5e94d6582319.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/083e6dea-2bdc-4882-8765-706be93d7d48.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/3071112d-efdd-4e7a-94f6-035337d0e652.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/3b8f3b90-a916-4b14-98f8-c93d303b3832.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/4f0be00b-b8bf-4d35-b797-b102e7456a78.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/5923733b-3761-456c-bb45-b36c9273fa5f.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/765189f0-b2ba-4fea-9e34-6797d376d231.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/7f5bcab5-2e42-4475-b154-2cc3f3d04d7b.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/909636f8-9794-45bd-aed5-730383e16182.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/afa4a500-aac5-49ea-8cd7-c93982ff5bb1.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/c430343a-27e1-46cd-9102-2eb03dc04333.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/c4332089-c1a5-408b-8d88-76b851ca8f3c.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/d0005024-a266-4b14-8165-a285baff534f.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/d32f4499-d80c-4435-84bf-2ff57d3c8ab3.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/d6c2ba2f-ccd4-4c23-ae05-43e56684675a.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/ed47842b-9db9-40ca-88c3-896b5b6242e0.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/f065339b-4deb-42e5-bb9f-d3714e5b12b1.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/fb019890-2e17-4bc5-aac2-040be5424e8e.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/ff7e5009-7b38-4fe9-bc00-9cb895f1db08.root',
  #'/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/ffc26379-9571-4929-bac4-1e882a6caef9.root'
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/0b2b0b0b-f312-48a8-9d46-ccbadc69bbfd.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/0c3cb20d-8556-450d-b4f0-e5c754818f74.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/0eafa2b4-711a-43ec-be1c-7e564c294a9a.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/1450b1bb-171e-495e-a767-68e2796d95c2.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/15498564-9cf0-4219-aab7-f97b3484b122.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/1838a806-316b-4f53-9d22-5b3856019623.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/1a34eb87-b9a3-47fb-b945-57e6f775fcac.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/1add5b2e-19cb-4581-956d-271907d03b72.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/1bed1837-ef65-4e07-a2ac-13c705b20fc1.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/1d057884-72bd-4353-8375-ec4616c00a33.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/1f75fe73-12d7-41c7-80b0-19b988571f15.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/20357610-cd50-4e51-b2f1-bac317762b91.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/2464c1ff-8de5-45a1-9376-2ac27c529190.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/28f2ac9b-b0a6-44a1-b10e-32ea9f59b611.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/29b8a5ec-f1e6-41f5-9eb8-228a5faa92da.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/2bda740e-2589-4b17-8b06-acdbb7c281af.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/2cffac22-a285-439a-ab14-3ec3a58fd15f.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/2d6fa796-3236-4fbf-a1b0-f5fa3db0686a.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/2dbb8df5-bc9a-492f-9fc7-99991bca1c21.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/3161d0a6-efc9-409b-88ca-c64a7131b97a.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/32f4fdb1-ea7c-4cda-ba98-e26975a52d44.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/3645d756-7fa3-4d0a-8c05-572f466ac938.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/375aa37c-fc24-4b3d-8e96-6a3cd9963789.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/390491c5-0e8d-4548-9077-f79d595d3a02.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/3a9968e3-9403-4c51-9030-0e6fd74788e8.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/3e05598b-4ee8-4fe7-a67f-2c2477f6c470.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/3e0ad616-c2de-4d1a-9b6e-b065fe0861f3.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/43d6b9da-7e64-4dd3-8c28-251f4fbdefb8.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/4517ad52-5a8c-4c7a-96b0-febf9530a9d4.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/4e29cfae-4d75-44c7-8e7e-34b4c4c90f00.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/4e627aac-0008-4f50-bb92-987f07dc0da9.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/4fa6929d-1f92-404d-8a5e-67391427c70b.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/54f1df94-6ee0-4b21-85da-853953ceddf9.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/5665b92b-eef9-4f4e-9f23-32de71e79aa9.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/56e1df5a-3a09-412f-ad5d-98bc23daa6e5.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/5b5c9cef-6a2f-4b88-a442-8d4805985ec4.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/5f3e43bd-3c7f-4f67-8639-7e500c4be0fe.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/5fd5749b-0f44-4213-b393-3ea283dc1176.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/609ddb77-6a09-451d-bd2c-c079629e63d6.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/60e42382-95dd-41eb-838b-e30662a898e0.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/6388751b-83d7-4cf6-8bb3-ab32d333a776.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/69ddd3b1-eb25-41a8-8dc0-d564ee70672f.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/6d56ec0d-666c-42a9-95b7-a23aa7c50f5a.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/709c9606-dd42-44d5-b220-c8345a75ad6e.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/7223614f-d8f6-4245-8055-12e64532ec2d.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/735659c5-68c0-4654-95e8-e9b6c545d1bb.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/768b6e5e-07cf-48da-bdc4-64f382c0a39c.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/7968ec81-dfe8-4831-b042-c763c05b84cf.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/7ad6be04-acaa-4b73-98ea-5c463e89dc6f.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/7c8f2946-7b05-4004-803f-5273e646837f.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/7d2d33c2-c5b1-445e-b738-8bccd2af5c68.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/7f51a425-cf98-41dd-b8e8-d638940ac4d9.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/81bd28bb-0a7f-4db7-a34b-8657ced22ff2.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/82117053-66a6-463e-ac78-13bbcd7c56f2.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/858edf87-673d-4830-b432-1d22ccfde2aa.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/88200c11-c54d-4941-a8b5-850f7a458ac3.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/8a732918-ac36-4c09-8847-cbfc50a80f4a.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/91b55b13-e0f4-482e-a56a-27ebd66d8534.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/96fb3724-6857-4b2a-a6df-3415aa290d43.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/a012f608-a91d-4353-87cd-4f8bb5bfb8a7.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/a6832206-6762-49dd-8841-e4040cbe7cd7.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/a7195537-3e8a-4a80-9636-ec42df6709b3.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/a8b3b3da-0083-4579-90d6-b775c8a87040.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/a9975bdc-5d9d-4c99-a01a-1ff36d7635bf.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/adc048ac-044e-4199-b120-c62969434a22.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/ae7b84cf-044f-4d9d-8ffd-03e9f0713f91.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/aef41e6d-4299-42da-ab04-3343d65abbf8.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/b2f3d6a7-188e-48fb-93f0-7dbafdf0481c.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/b416cb61-2287-4a1a-b736-92c2ed041fb6.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/b776a639-6cdd-41ce-8698-51b9c9167378.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/ba7ee146-ea40-450a-b4f7-74727c089eb1.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/bb307357-f3a3-47d1-b951-bd633b5ec055.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/bb67feda-4a7f-4469-a382-7ecd9fe77ddf.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/bbbdd79c-9123-4617-8c12-4e859fdd823f.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/bdd8fc4f-4183-4f81-9bf4-92f81dc910af.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/c3da9637-cd80-4fe7-b9fe-c0790bbc1a3b.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/c5f424f5-d6f2-4665-8cb7-05427f980c37.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/c708603a-e563-4c14-81cd-5e6149824ee6.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/c8a4afc6-c07b-434b-8751-1ea90b0b082c.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/ca938992-5547-41a0-8790-6b80a11feada.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/cb6e15d9-5796-4aa9-9ac0-8e51ecf47f39.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/d3f7a4bd-4108-4046-8bbf-4c6739b772fa.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/d698fbcc-2865-4950-b913-377579637985.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/d91923ed-cdfa-4df3-ade1-8fde83f8a09c.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/da3595cc-e3d5-40f8-bb4b-1601216bc3a1.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/db4ee0f9-ade1-4329-82fa-a78997dec695.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/dc02f812-7434-49b1-b820-7de4f90acb75.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/dd241ccd-7488-47cb-b049-140e013aa71d.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/e46876ab-ee6e-4af7-bd19-3c8986334b33.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/e54cfbf6-6bd7-4441-96f9-e332b2ba040f.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/e67662b5-3cd5-42cd-8222-e6797b4264b2.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/eacfd61d-6a31-450d-bbd2-4f4533083c35.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/eb3f386a-5c0d-4f0c-a388-e2b40f99ba4e.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/ecc968c1-9afe-4a28-9abc-413948e54641.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/f2f8e734-ef96-4ff3-bc06-8e808e9a5419.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/f82eca78-0458-4d1e-886a-20570188d76b.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/f8c1d660-a1de-4f19-b46f-c63717032647.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/faeba23e-e477-4bd0-8c8b-77e5fd7e73d0.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/fc01f938-5334-4208-bde6-cc60c9830477.root',
  '/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/fefae469-22d6-455e-a222-144e8e63043e.root'
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
  #skipEvents = cms.untracked.uint32( 353+20+110+90 ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)
process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )
process.MessageLogger.cerr.enableStatistics = False
process.TFileService = cms.Service( "TFileService", fileName = cms.string( "Hist.root" ) )

if ( False ):
  process.out = cms.OutputModule (
    "PoolOutputModule",
    fileName = cms.untracked.string("L1Tracks.root"),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring('drop *', 'keep *_TTTrack*_*_*', 'keep *_TTStub*_*_*' )
  )
  process.FEVToutput_step = cms.EndPath( process.out )
  process.schedule.append( process.FEVToutput_step )
