import FWCore.ParameterSet.Config as cms
from L1Trigger.TrackTrigger.TrackQualityParams_cfi import *

TrackTriggerDataFormats_params = cms.PSet (

  UseHybrid = cms.bool( True ),

  TrackQualityPSet = cms.PSet(TrackQualityParams),

  ZHoughTransform = cms.PSet (

    NumBinsZT  = cms.int32( 2 ),
    NumBinsCot = cms.int32( 2 ),
    NumStages  = cms.int32( 5 )

  ),

  KalmanFilter = cms.PSet (

    RangeFactor = cms.double( 2.0 ) # search window of each track parameter in initial uncertainties

  ),

  DuplicateRemoval = cms.PSet (
    WidthPhi0    = cms.int32( 12 ), # number of bist used for phi0
    WidthQoverPt = cms.int32( 15 ), # number of bist used for qOverPt
    WidthCot     = cms.int32( 16 ), # number of bist used for cot(theta)
    WidthZ0      = cms.int32( 12 )  # number of bist used for z0
  ),

)
