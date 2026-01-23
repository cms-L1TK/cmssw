import FWCore.ParameterSet.Config as cms

l1VertexAnalyzer = cms.EDAnalyzer(
    "L1VertexAnalyzer",
    vertices = cms.InputTag("l1tVertexProducer", "L1Vertices"),
    debug = cms.untracked.bool(False),
)
