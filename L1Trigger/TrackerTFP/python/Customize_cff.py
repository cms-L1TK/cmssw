import FWCore.ParameterSet.Config as cms

def setupUseTMTT(process):
    process.TrackTriggerSetup.TrackingParticle.MinPt = 3.0
    process.TrackTriggerSetup.TrackFinding.MinPt = 3.0
    process.TrackTriggerSetup.TrackFinding.MaxEta = 2.4
    process.TrackTriggerSetup.TrackFinding.ChosenRofPhi = 67.24
    process.TrackTriggerSetup.KalmanFilterIn.ShiftRangePhi = 0
    process.TrackTriggerSetup.KalmanFilterIn.ShiftRangeZ   = 0
    process.TrackTriggerSetup.KalmanFilter.NumWorker = 2
    process.TrackTriggerSetup.KalmanFilter.MaxLayers = 7
    process.TrackTriggerSetup.KalmanFilter.ShiftInitialC00 = -3+13
    process.TrackTriggerSetup.KalmanFilter.ShiftInitialC11 = -3+13
    process.TrackTriggerSetup.KalmanFilter.ShiftInitialC22 = -4+14
    process.TrackTriggerSetup.KalmanFilter.ShiftInitialC33 = -4+14
    process.TrackTriggerKalmanFilterFormats.BaseShiftx0           =  -1
    process.TrackTriggerKalmanFilterFormats.BaseShiftx1           =  -8
    process.TrackTriggerKalmanFilterFormats.BaseShiftx2           =   0
    process.TrackTriggerKalmanFilterFormats.BaseShiftx3           =  -2
    process.TrackTriggerKalmanFilterFormats.BaseShiftv0           =  -4
    process.TrackTriggerKalmanFilterFormats.BaseShiftv1           =  10
    process.TrackTriggerKalmanFilterFormats.BaseShiftr0           =  -7
    process.TrackTriggerKalmanFilterFormats.BaseShiftr1           =   0
    process.TrackTriggerKalmanFilterFormats.BaseShiftS00          =   0
    process.TrackTriggerKalmanFilterFormats.BaseShiftS01          =  -7
    process.TrackTriggerKalmanFilterFormats.BaseShiftS12          =   7
    process.TrackTriggerKalmanFilterFormats.BaseShiftS13          =   5
    process.TrackTriggerKalmanFilterFormats.BaseShiftK00          = -16
    process.TrackTriggerKalmanFilterFormats.BaseShiftK10          = -23
    process.TrackTriggerKalmanFilterFormats.BaseShiftK21          = -21
    process.TrackTriggerKalmanFilterFormats.BaseShiftK31          = -23
    process.TrackTriggerKalmanFilterFormats.BaseShiftR00          =  -4
    process.TrackTriggerKalmanFilterFormats.BaseShiftR11          =   9
    process.TrackTriggerKalmanFilterFormats.BaseShiftInvR00Approx = -27
    process.TrackTriggerKalmanFilterFormats.BaseShiftInvR11Approx = -40
    process.TrackTriggerKalmanFilterFormats.BaseShiftInvR00Cor    = -15
    process.TrackTriggerKalmanFilterFormats.BaseShiftInvR11Cor    = -15
    process.TrackTriggerKalmanFilterFormats.BaseShiftInvR00       = -25
    process.TrackTriggerKalmanFilterFormats.BaseShiftInvR11       = -37
    process.TrackTriggerKalmanFilterFormats.BaseShiftC00          =   5
    process.TrackTriggerKalmanFilterFormats.BaseShiftC01          =  -2
    process.TrackTriggerKalmanFilterFormats.BaseShiftC11          =  -7
    process.TrackTriggerKalmanFilterFormats.BaseShiftC22          =   6
    process.TrackTriggerKalmanFilterFormats.BaseShiftC23          =   4
    process.TrackTriggerKalmanFilterFormats.BaseShiftC33          =   4
    return process