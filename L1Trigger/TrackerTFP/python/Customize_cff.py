# function to alter TrackTriggerSetup to provide TMTT configuration

import FWCore.ParameterSet.Config as cms

def setupUseTMTT(process):
  process.TrackTriggerSetup.TrackingParticle.MinPt = 3.0
  process.TrackTriggerSetup.TMTT.WidthR   = 11
  process.TrackTriggerSetup.TMTT.WidthPhi = 14
  process.TrackTriggerSetup.TMTT.WidthZ   = 13
  process.TrackTriggerSetup.TrackFinding.MinPt        =  3.0
  process.TrackTriggerSetup.TrackFinding.MinPtCand    =  3.0
  process.TrackTriggerSetup.TrackFinding.MaxEta       =  2.4
  process.TrackTriggerSetup.TrackFinding.ChosenRofPhi = 67.24
  process.TrackTriggerSetup.TrackFinding.NumLayers    =  8
  process.TrackTriggerSetup.HoughTransform.MinLayers = 5
  process.TrackTriggerSetup.CleanTrackBuilder.MaxStubs = 4
  process.TrackTriggerSetup.KalmanFilter.RangeFactor     =  3.0
  process.TrackTriggerSetup.KalmanFilter.NumWorker       =  4
  process.TrackTriggerSetup.KalmanFilter.MaxLayers       =  8
  process.TrackTriggerSetup.KalmanFilter.MaxSeedingLayer =  3
  process.TrackTriggerSetup.KalmanFilter.MaxGaps         =  2
  process.TrackTriggerSetup.KalmanFilter.BaseShift       = 0
  process.TrackTriggerSetup.KalmanFilter.ShiftChi20      = 0
  process.TrackTriggerSetup.KalmanFilter.ShiftChi21      = 0
  process.TrackTriggerSetup.KalmanFilter.CutChi2         = 2.0
  
  process.TrackTriggerKalmanFilterFormats.BaseShiftx0           =  -4
  process.TrackTriggerKalmanFilterFormats.BaseShiftx1           = -10
  process.TrackTriggerKalmanFilterFormats.BaseShiftx2           =  -5
  process.TrackTriggerKalmanFilterFormats.BaseShiftx3           =  -6
  process.TrackTriggerKalmanFilterFormats.BaseShiftr0           =  -9
  process.TrackTriggerKalmanFilterFormats.BaseShiftr1           =  -3
  process.TrackTriggerKalmanFilterFormats.BaseShiftS00          =  -6
  process.TrackTriggerKalmanFilterFormats.BaseShiftS01          = -12
  process.TrackTriggerKalmanFilterFormats.BaseShiftS12          =  -2
  process.TrackTriggerKalmanFilterFormats.BaseShiftS13          =  -3
  process.TrackTriggerKalmanFilterFormats.BaseShiftK00          =  -5
  process.TrackTriggerKalmanFilterFormats.BaseShiftK10          = -12
  process.TrackTriggerKalmanFilterFormats.BaseShiftK21          = -12
  process.TrackTriggerKalmanFilterFormats.BaseShiftK31          = -12
  process.TrackTriggerKalmanFilterFormats.BaseShiftR00          =  -8
  process.TrackTriggerKalmanFilterFormats.BaseShiftR11          =   3
  process.TrackTriggerKalmanFilterFormats.BaseShiftInvR00Approx = -36
  process.TrackTriggerKalmanFilterFormats.BaseShiftInvR11Approx = -47
  process.TrackTriggerKalmanFilterFormats.BaseShiftInvR00Cor    = -15
  process.TrackTriggerKalmanFilterFormats.BaseShiftInvR11Cor    = -15
  process.TrackTriggerKalmanFilterFormats.BaseShiftInvR00       = -27
  process.TrackTriggerKalmanFilterFormats.BaseShiftInvR11       = -38
  process.TrackTriggerKalmanFilterFormats.BaseShiftC00          =  11
  process.TrackTriggerKalmanFilterFormats.BaseShiftC01          =   4
  process.TrackTriggerKalmanFilterFormats.BaseShiftC11          =  -2
  process.TrackTriggerKalmanFilterFormats.BaseShiftC22          =   8
  process.TrackTriggerKalmanFilterFormats.BaseShiftC23          =   8
  process.TrackTriggerKalmanFilterFormats.BaseShiftC33          =   7
  process.TrackTriggerKalmanFilterFormats.BaseShiftr02          =  -2
  process.TrackTriggerKalmanFilterFormats.BaseShiftr12          =  10
  process.TrackTriggerKalmanFilterFormats.BaseShiftchi20        = -13
  process.TrackTriggerKalmanFilterFormats.BaseShiftchi21        = -12

  return process
