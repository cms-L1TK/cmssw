import FWCore.ParameterSet.Config as cms

l1TrackZ0ResolutionAnalyzer = cms.EDAnalyzer(
    "L1TrackZ0ResolutionAnalyzer",
    vertices        = cms.InputTag("l1tVertexFinder", "L1Vertices"),
    selectedTracks  = cms.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelected"),
    associatedTracks= cms.InputTag("l1tTrackVertexAssociationProducer", "Level1TTTracksSelectedAssociated"),
    debug           = cms.untracked.bool(False),
)
