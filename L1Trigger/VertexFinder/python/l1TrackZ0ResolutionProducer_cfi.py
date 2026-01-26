import FWCore.ParameterSet.Config as cms

l1TrackZ0NtupleProducer = cms.EDAnalyzer(
    "L1TrackZ0NtupleProducer",
    vertices=cms.InputTag("l1tVertexFinder", "L1Vertices"),
    associatedTracks=cms.InputTag("l1tTrackVertexAssociationProducer", "Level1TTTracksSelectedAssociated"),

    # untracked
    debug=cms.untracked.bool(False),
    storeDz=cms.untracked.bool(True),
)
