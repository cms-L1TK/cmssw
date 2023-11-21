import FWCore.ParameterSet.Config as cms

TrackerTFPProducer_params = cms.PSet (

  LabelDTC = cms.string( "TrackerDTCProducer"    ),      #
  LabelPP  = cms.string( "TrackerTFPProducerPP"  ),      #
  LabelGP  = cms.string( "TrackerTFPProducerGP"  ),      #
  LabelHT  = cms.string( "TrackerTFPProducerHT"  ),      #
  LabelCTB = cms.string( "TrackerTFPProducerCTB" ),      #
  LabelKF  = cms.string( "TrackerTFPProducerKF"  ),      #
  LabelDR  = cms.string( "TrackerTFPProducerDR"  ),      #
  LabelTFP = cms.string( "TrackerTFPProducerTFP" ),      #
  BranchStubsAccepted   = cms.string( "StubAccepted"  ),  # branch for prodcut with passed stubs
  BranchTracksAccepted  = cms.string( "TrackAccepted" ),  # branch for prodcut with passed tracks
  BranchStubsTruncated  = cms.string( "StubLost"      ),  # branch for prodcut with lost stubs
  BranchTracksTruncated = cms.string( "TracksLost"    ),  # branch for prodcut with lost tracks
  EnableTruncation = cms.bool  ( True ),                # enable emulation of truncation, lost stubs are filled in BranchLost
  PrintKFDebug     = cms.bool  ( False )                 # print end job internal unused MSB

)