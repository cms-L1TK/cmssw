import FWCore.ParameterSet.Config as cms

# ParameterSet used by HitPatternHelper

HitPatternHelper_params = cms.PSet (

  hphDebug   = cms.bool(False),   
  useNewKF   = cms.bool(False),
  chosenRofZ = cms.double(50.),
  deltaTanL  = cms.double(0.125),
  etaRegions = cms.vdouble(-2.4,-2.08,-1.68,-1.26,-0.90,-0.62,-0.41,-0.20,0.0,0.20,0.41,0.62,0.90,1.26,1.68,2.08,2.4)

)
