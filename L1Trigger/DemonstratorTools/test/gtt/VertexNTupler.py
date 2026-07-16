# this compares event by event the output of the C++ emulation with the ModelSim simulation of the firmware
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load( 'FWCore.MessageService.MessageLogger_cfi' )
process.load( 'Configuration.EventContent.EventContent_cff' )
process.load( 'Configuration.Geometry.GeometryExtendedRun4D110Reco_cff' ) 
process.load( 'Configuration.Geometry.GeometryExtendedRun4D110_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load( 'L1Trigger.TrackerDTC.DTC_cff' )
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
process.load( 'L1Trigger.TrackFindingTracklet.Demonstrator_cff' )
from L1Trigger.TrackFindingTracklet.Customize_cff import *
#reducedConfig( process )
fwConfig( process )

# build schedule
process.TrackProcessorEmulation = cms.Sequence (  process.ProducerDTC
                                                + process.L1THybridTracks
                                                + process.ProducerTM
                                                + process.ProducerDR
                                                + process.ProducerKF
                                                + process.ProducerTQ
                                                + process.ProducerTFP
                                                )

def get_data_from_file(filepath):
    samples = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # Remove whitespace and newline characters
                line = line.strip()
                # Skip empty lines
                if line:
                    samples.append(line)
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
        return []
    except Exception as e:
        print(f"Error reading file: {e}")
        return []
    
    return samples

Samples = get_data_from_file("/home/hep/am2023/new-investigation/CMSSW_15_1_0_pre4/src/Files.txt")[0]
#Samples = ["/store/mc/Phase2Spring24DIGIRECOMiniAOD/DisplacedSUSY_stopToBottom_M-800_500mm_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2810000/c63cb0cc-b447-4220-8f72-8005a07aa8a8.root"]

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring( Samples )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("vertices.root")
)

process.vertexAnalyzer = cms.EDAnalyzer("SimpleVertexAnalyzer",
    SimVertexInputTag = cms.InputTag("g4SimHits",""),
    TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
    MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
    MyProcess = cms.int32(1),
    DebugMode = cms.bool(False),
    L1Tk_nPar = cms.int32(4), # use 4 or 5-parameter L1 tracking?
    L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
    TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
    TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
    TP_minPt = cms.double(1.9),       # only save TPs with pt > X GeV
    TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
    TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
)

process.p = cms.Path(process.vertexAnalyzer)