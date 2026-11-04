####################################################
# Analyze TTClusters & TTStubs                     #
# (P2 objects used for L1 tracking)                #
# using either AnalyzerClusterStub.cc,             #
# or AnalyzerCluster.cc.                           #
# or AnalyzerStub.cc.                              #
# (Edit this file to decide which to run)          #
####################################################

################################################################################
# import of general framework
################################################################################
import FWCore.ParameterSet.Config as cms

process = cms.Process('AnalyzerClusterStub')

GEOMETRY = "D121"

################################################################################
# global tag
################################################################################
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

################################################################################
# load the specific tracker geometry
################################################################################
process.load('Configuration.Geometry.GeometryExtendedRun4' + GEOMETRY + 'Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4' + GEOMETRY +'_cff')

################################################################################
# load the magnetic field
################################################################################
process.load('Configuration.StandardSequences.MagneticField_cff')

################################################################################
# define the source and maximum number of events to generate and simulate
################################################################################
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_20_0_0_patch1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/150X_mcRun4_realistic_v1_Run4_D121_pp_noPU_STD_Reference_20260903_192421_RV19-v1/2830000/0c820980-7257-49d8-b9b4-1fcae28a2f2a.root')
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

################################################################################
# load the analyzers
################################################################################
# Old analyzer for clusters & stubs (rates, origin of stubs etc.)
process.AnalyzerClusterStub = cms.EDAnalyzer("AnalyzerClusterStub",
  DebugMode = cms.bool(True),
  TTClusterInputTag = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive"),
  TTStubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
  TTClusterAssocInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterInclusive"),
  TTStubAssocInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
  TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
  TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),   
)

# New analyzer for clusters (position resolution w.r.t. SimHits etc.)
process.AnalyzerCluster = cms.EDAnalyzer("AnalyzerCluster",
  trackingParticles= cms.InputTag("mix","MergedTrackTruth"),
  ttclusters = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive"),
  links = cms.InputTag("simSiPixelDigis", "Tracker"),
  simhitsbarrel = cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
  simhitsendcap = cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"),
  simtracks = cms.InputTag("g4SimHits"),
  ECasRings = cms.bool(True),
  SimTrackMinPt = cms.double(2.)
)

# New analyzer for stub (bend resolution etc.)
process.AnalyzerStub = cms.EDAnalyzer("AnalyzerStub",
  trackingParticles= cms.InputTag("mix","MergedTrackTruth"),
  ttstubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
  ttstubassoc = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
  ECasRings = cms.bool(True)
)

#################################################################################################
# define output file and message logger
#################################################################################################
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('file:analysis.root'),
    closeFileFast = cms.untracked.bool(True)
)

#################################################################################################
# define the final path to be fed to cmsRun
#################################################################################################

#-- Choose which analyzer to run

# This is the old one for both clusters & stubs.
#process.p = cms.Path( process.AnalyzerClusterStub )

# This is the new one for clusters alone (includes resolution plots)
#process.p = cms.Path( process.AnalyzerCluster )

# This is the new one for stubs alone (includes resolution plots)
process.p = cms.Path( process.AnalyzerStub )



