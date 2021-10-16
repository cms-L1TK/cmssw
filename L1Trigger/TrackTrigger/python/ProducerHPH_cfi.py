import FWCore.ParameterSet.Config as cms

HitPatternHelper_params = cms.PSet (

  hphDebug   = cms.bool(False),   
  useNewKF   = cms.bool(False),
  chosenRofZ = cms.double(50.0),
  deltaTanL = cms.double(0.125)

)
