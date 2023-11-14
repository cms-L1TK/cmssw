import FWCore.ParameterSet.Config as cms
from L1Trigger.TrackTrigger.TrackQualityParams_cfi import *

TrackFindingTrackletProducer_params = cms.PSet (

  InputTag              = cms.InputTag( "l1tTTTracksFromTrackletEmulation", "Level1TTTracks"), #
  InputTagDTC           = cms.InputTag( "TrackerDTCProducer", "StubAccepted"),              #
  LabelTBout            = cms.string  ( "TrackFindingTrackletProducerTBout" ),              #
  LabelDRin             = cms.string  ( "TrackFindingTrackletProducerDRin"  ),              #
  LabelDR               = cms.string  ( "TrackFindingTrackletProducerDR"    ),              #
  LabelCTB              = cms.string  ( "TrackFindingTrackletProducerKFin"  ),              #
  LabelKF               = cms.string  ( "TrackFindingTrackletProducerKF"    ),              #
  LabelTT               = cms.string  ( "TrackFindingTrackletProducerTT"    ),              #
  LabelAS               = cms.string  ( "TrackFindingTrackletProducerAS"    ),              #
  LabelKFout            = cms.string  ( "TrackFindingTrackletProducerKFout" ),              #
  BranchStubsAccepted   = cms.string  ( "StubAccepted"  ),  # branch for prodcut with passed stubs
  BranchTracksAccepted  = cms.string  ( "TrackAccepted" ),  # branch for prodcut with passed tracks
  BranchStubsTruncated  = cms.string  ( "StubLost"      ),  # branch for prodcut with lost stubs
  BranchTracksTruncated = cms.string  ( "TracksLost"    ),  # branch for prodcut with lost tracks
  EnableTruncation      = cms.bool    ( True  ),                                            # enable emulation of truncation for TBout, KF, KFin, lost stubs are filled in BranchLost
  PrintKFDebug          = cms.bool    ( True ),                                            # print end job internal unused MSB
  UseTTStubResiduals    = cms.bool    ( False ),                                            # stub residuals are recalculated from seed parameter and TTStub position

)
