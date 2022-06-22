import FWCore.ParameterSet.Config as cms

def setupUseTMTT(process):
    process.TrackTriggerKalmanFilterFormats.UseHybrid = False
    process.TrackTriggerSetup.TrackingParticle.MinPt = 3.0
    process.TrackTriggerSetup.Firmware.MaxdPhi = 0.01
    process.TrackTriggerSetup.TrackFinding.MinPt = 3.0
    process.TrackTriggerSetup.TrackFinding.MaxEta = 2.4
    process.TrackTriggerSetup.TrackFinding.ChosenRofPhi = 67.24
    process.TrackTriggerSetup.KalmanFilterIn.ShiftRangePhi = 0
    process.TrackTriggerSetup.KalmanFilterIn.ShiftRangeZ   = 0
    process.TrackTriggerSetup.KalmanFilter.NumWorker = 2
    return process