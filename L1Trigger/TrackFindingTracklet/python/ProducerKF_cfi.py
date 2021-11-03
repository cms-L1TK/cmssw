import FWCore.ParameterSet.Config as cms

TrackFindingTrackletProducerKF_params = cms.PSet (

  InputTag             = cms.InputTag( "TrackletTracksFromTrackletEmulation", "Level1TTTracks"), #
  InputTagDTC          = cms.InputTag( "TrackerDTCProducer", "StubAccepted"),                    #
  LabelKFin            = cms.string  ( "TrackFindingTrackletProducerKFin"  ),                    #
  LabelKF              = cms.string  ( "TrackFindingTrackletProducerKF"    ),                    #
  LabelTT              = cms.string  ( "TrackFindingTrackletProducerTT"    ),                    #
  LabelAS              = cms.string  ( "TrackFindingTrackletProducerAS"    ),                    #
  LabelKFout           = cms.string  ( "TrackFindingTrackletProducerKFout" ),                    #
  BranchAcceptedStubs  = cms.string  ( "StubAccepted"  ),                                        #
  BranchAcceptedTracks = cms.string  ( "TrackAccepted" ),                                        #
  BranchLostStubs      = cms.string  ( "StubLost"      ),                                        #
  BranchLostTracks     = cms.string  ( "TrackLost"     ),                                        #
  CheckHistory         = cms.bool    ( False ),                                                  # checks if input sample production is configured as current process
  EnableTruncation     = cms.bool    ( True  )                                                   # enable emulation of truncation, lost stubs are filled in BranchLost

  SummerChain         = cms.bool  ( False ),                                                    # reduce l1 tracking to summer chain configuration
  SummerChainChannels = cms.vint32( 0, 1, 25, 2, 26, 4, 28, 5, 29, 6, 30, 7, 31, 8, 32, 9, 33 ) # map of used tfp channels in summer chain config
)