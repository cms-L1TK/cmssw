# configuration for hybrid track reconstruction chain emulating EDProducer

import FWCore.ParameterSet.Config as cms

TrackFindingTrackletProducer_params = cms.PSet (

  BDT              = cms.FileInPath( "L1Trigger/TrackFindingTracklet/data/HYBRID_NEW_KF_TQ_XGBOOST_vl.json" ), #
  InputTagTTDTC    = cms.InputTag  ( "ProducerDTC", "StubAccepted" ),
  InputTagTracklet = cms.InputTag  ( "l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks" ), # Tracklet TTTracks
  InputLabelTFP    = cms.string( "ProducerTQ"      ),                  #
  InputLabelTQ     = cms.string( "ProducerKF"      ),                  #
  InputLabelKF     = cms.string( "ProducerDR"      ),                  #
  InputLabelDR     = cms.string( "ProducerTM"      ),                  #
  InputLabelASP    = cms.string( "l1tTTTracksFromTrackletEmulation" ), #
  InputLabelTM     = cms.string( "l1tTTTracksFromTrackletEmulation" ), #
  InputLabelPP     = cms.string( "ProducerDTC"     ),                  #
  BranchStubs      = cms.string( "StubAccepted"    ),                  #
  BranchAllStubs   = cms.string( "AllStubAccepted" ),                  #
  BranchTracks     = cms.string( "TrackAccepted"   ),                  #
  BranchTTTracks   = cms.string( "TTTrackAccepted" ),                  #
  BranchTruncated  = cms.string( "Truncated"       ),                  #
  PrintKFDebug     = cms.bool  ( False ),                              # print end job internal unused MSB

)
