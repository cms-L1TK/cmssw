# this compares event by event the output of the C++ emulation with the ModelSim simulation of the firmware
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

# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
# cosutmize TT algorithm
from L1Trigger.TrackerDTC.Customize_cff import *
producerUseTMTT(process)
analyzerUseTMTT(process)
#--- Load code that produces tfp Stubs
process.load( 'L1Trigger.TrackerTFP.Producer_cff' )
from L1Trigger.TrackerTFP.Customize_cff import *
setupUseTMTT( process )
#--- Load code that demonstrates tfp Stubs
process.load( 'L1Trigger.TrackerTFP.Demonstrator_cff' )

# build schedule
process.tt = cms.Sequence (  process.TrackerDTCProducer
                           + process.TrackerTFPProducerGP
                           + process.TrackerTFPProducerHT
                           + process.TrackerTFPProducerMHT
                           + process.TrackerTFPProducerZHT
                           + process.TrackerTFPProducerTTFound
                           + process.TrackerTFPProducerKFin
                          )
process.demo = cms.Path( process.tt + process.TrackerTFPDemonstrator )
process.schedule = cms.Schedule( process.demo )

# create options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing( 'analysis' )
# specify input MC
Samples = [
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/00b3d04b-4c7b-4506-8d82-9538fb21ee19.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/35d2ec8f-be56-4f30-ada0-4d35a561f1ea.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/677fe665-1c93-4e29-876d-d39eb5540184.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/95b4f689-b24b-4a02-8432-6b0e6e9db822.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/d0777d76-db01-4bc3-8833-0e297b5095c7.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0390df7b-7c2a-45d0-9bdb-e13f6565f65a.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/38ec0a8d-1140-4bd2-b69f-94ff8577c421.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/6ad38b6e-1108-4525-b497-0811e2e4b826.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/9853e109-a978-4832-8f1f-fc311b409939.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/d3a22ce7-ca19-4cca-9290-9622acd5753c.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0919ec05-0bdc-4f85-9a21-3a167117ea5e.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/3b9fb9b9-ec7b-45c4-a26c-ee3403461173.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/6ca457b1-4415-414c-a87c-342ed56cc6a7.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/9a2172a4-70bb-49b2-82f6-cc97517d496c.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/d895c18f-845d-4a11-9a0c-1c3f98d316f2.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/09e5d64c-d0de-442c-943d-620927afc59c.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/3bb5aab6-b3ee-4d19-8ec2-7aa776171176.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/6d465243-cf6e-42e4-b14d-c944636af5ba.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/9f750704-3eeb-4072-8d25-68cdc0235c90.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/d8a67776-9942-4489-b042-686bdb5b6196.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0ce69121-29ff-4f67-babf-c3372a273ce6.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/3f664e50-afea-421a-ae4e-a54a4cf6b2fc.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/714c4692-2a55-4a96-b41d-b9829766e8e7.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/a0b8d8f6-b405-461f-b9e6-6f5cf32f5675.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/df2db8fc-deab-47df-8dd6-a6e7fc8c32d5.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/0eccd3e4-cbe0-403e-a574-750ec48804fc.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/429dab4d-9f42-41b1-9798-a32b3f8c87d0.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/71e1ab5a-818e-4a02-aafe-78d0ef02826a.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/a465ff4f-f985-4191-b7cd-7b794b86c41e.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/e0dbeedb-d615-4839-ace3-dda55cf41e72.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1118e869-7f5d-4d34-98ae-43d09b52c437.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/470e682a-a48a-4928-bb0c-1d71c673819e.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/74443259-d9ef-43d5-bc50-85c18d19b2b6.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/aaaf4370-de26-4037-87b8-967c6e4e5fee.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/e25cc87c-fe12-4252-9fb3-15808b0315b3.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1967aa5c-4fb4-4039-ac75-2f6fcdc8d0b0.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/4e94d058-2744-48a5-860c-2678df096976.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/7b432181-b6fe-4003-90b8-7746c51c6693.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/ac95b939-6ed5-4af7-9ac5-c1210f9c8bd8.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/e6fab173-bdbd-4cfc-a5b4-0a517a497d9a.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/19fb8c07-8b8c-418f-aa6f-8bedc4f8c8c5.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/5398b6ba-fa68-4ccd-9876-9323d33eb3fc.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/7c54d6fc-9188-4e50-866c-abd260d97c47.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/ad27659a-dbd7-426f-af99-9cc6c794e442.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/e77dffdb-1db3-4a23-9869-1900b9b7809f.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1a364293-a1ab-42a0-bbda-afc6d404cb2e.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/54c0512c-3099-4b95-91c7-cf7ed17d8c71.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/7d3ebe57-2d6a-465f-a8b3-4067a773313c.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/af0c26df-52ea-44f1-a7ef-91c24f97b71f.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/e94070fc-09a4-4f2a-83c7-063d67b0bb7c.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1e1d5661-929b-4e8e-a123-ed3e579f14fc.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/54db435c-b726-4d94-9dff-b7e7399f7ff6.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/80be925f-0091-4ed8-bff8-c3bef37cba0f.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/b4500789-f4b3-4f4e-9eb2-4ad105217928.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/eb3e94b8-004a-4c88-a69f-59827f0e2cba.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/1f8ab4bb-1357-4621-a9ce-a61e09c913da.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/573a6963-97f0-431b-8738-298a5ca51855.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/82660447-bd8d-4e0c-b4a3-d3843bdbd3a7.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/b8c918a7-050f-44d5-8325-174770371510.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/eeea8311-1611-4d25-9e17-a3fc6e4a5324.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/2660d57d-101c-4903-8113-ee4e2d26aa78.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/5a344a2f-8976-4b14-8e32-f727f521873f.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/8c17163d-7d77-4c29-9147-7f175b4e9f4d.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/c49898e8-5bf4-4a8f-a6f8-7f6650621886.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/ef78cb26-294c-4557-a134-54fd01f6a7c0.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/28b3eb19-e1c9-49c2-a948-1656d514e43e.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/5ed94184-f81f-4df3-81fd-9e4a8f09bcbb.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/8d1d2408-218b-4588-aa57-7e641ef5d9d6.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/c4fa2f2d-c4b4-4eb9-83cf-24b6f34365fd.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/f1f78e78-7a0d-474a-ba4f-a31881f8c385.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/2997186c-8593-4247-990d-1cd968fb4d39.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/6090563a-ef31-44e5-a92f-45f07fbbd6ba.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/8f929a6a-b28f-471d-93f2-2a4d6a330409.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/c53bc90b-c441-4cf9-b9ce-ba88f122763e.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/f3a5b46e-4748-4c68-b6af-8318f1252e46.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/2b966c6d-cf7f-42ec-ac93-53ad49a8f04d.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/63cca44d-f2da-491b-ae8b-03e78d523152.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/8fc59b58-1e55-4524-8f5f-8919464869d6.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/ccca391d-bc50-4dc7-87c2-f8eff54be7f1.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/f3e4b5f1-d436-4f71-8f64-29cfdbc2678d.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/2ca6b4d1-685b-4044-8019-f3edee31bb0b.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/647cfbc3-837a-4076-a2cb-912f96d711d9.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/92de7daf-3271-4bd4-ab42-36cb362488e0.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/ccf7bb93-e520-4cf9-ac15-8c12f837ca92.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/f9d8acb9-9a32-4abf-a39c-c9f36e2dc167.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/33527dfd-f0fc-4e08-b56c-09ecaa1b9a93.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/64aa5786-d901-4607-9cb4-e269724761c5.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/93180bab-e6d1-4a1d-bcf8-6cbde082be22.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/ce82089c-1a1a-440c-b7ae-4b7692208bce.root',
'/store/relval/CMSSW_12_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v2_2026D88PU200-v1/2590000/ff6b1688-35a1-4dc9-8325-17abb035cb2c.root'
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
  skipEvents = cms.untracked.uint32( 100+182+123+111 ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)
process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )
