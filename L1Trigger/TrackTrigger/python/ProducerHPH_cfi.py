import FWCore.ParameterSet.Config as cms
from L1Trigger.TrackFindingTMTT.TMTrackProducer_Defaults_cfi import TMTrackProducer_params

# ParameterSet used by HitPatternHelper

HitPatternHelper_params = cms.PSet (

  hphDebug   = cms.bool(False),   
  useNewKF   = cms.bool(False),
  oldKFPSet  = cms.PSet(TMTrackProducer_params.EtaSectors),
  deltaTanL  = cms.double(0.125),
  cotNbins   = cms.int32(512),
  z0Nbins    = cms.int32(512)
  
)
