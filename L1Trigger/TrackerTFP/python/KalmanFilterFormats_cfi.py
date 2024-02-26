# configuration of internal KF variable bases which can be shifted by powers of 2 w.r.t. KF output track parameter
# TrackerTFPProducer_params.PrintKFDebug printouts unused MSB for each variable, so that one could consider decreasing the basseshift by that amount
# numerical instabillity (negative C00, C11, C22, C33) requires smaller baseshifts of related variables (rx, Sxx, Kxx, Rxx, invRxx)
# if a variable overflows an Exception will be thrown and the corresponding baseshift needs to be increased.
import FWCore.ParameterSet.Config as cms

TrackTriggerKalmanFilterFormats_params = cms.PSet (
  
  BaseShiftx0           = cms.vint32(  -5+1,  -5+1,  -5+1,   -5+1,   -5+1,   -5+1,   -5+1,  -5+1 ),
  BaseShiftx1           = cms.vint32( -10+1, -10+1, -10+1,  -10+1,  -10+1,  -10+1,  -10+1, -10+1 ),
  BaseShiftx2           = cms.vint32(  -5+1,  -5+1,  -5+1,   -5+1,   -5+1,   -5+1,   -5+1,  -5+1 ),
  BaseShiftx3           = cms.vint32(  -5+2,  -5+2,  -5+2,   -5+2,   -5+2,   -5+2,   -5+2,  -5+2 ),
  BaseShiftr0           = cms.vint32(  -9+2,  -9+2,  -9+2,   -9+2,   -9+2,   -9+2,   -9+2,  -9+2 ),
  BaseShiftr1           = cms.vint32(  -3+1,  -3+1,  -3+1,   -3+1,   -3+1,   -3+1,   -3+1,  -3+1 ),

  BaseShiftS00          = cms.vint32(   0,     0,     0+2,    1+2,    2+3,     3+3,   3+2,   1+2 ),
  BaseShiftS01          = cms.vint32(  -6,    -6,    -5+2,   -4+2,   -4+3,    -4+3,  -5+3,  -6+2 ),
  BaseShiftS12          = cms.vint32(   6,     6,     5+1,    6+1,    7+2,     9+2,   9+1,   7+1 ),
  BaseShiftS13          = cms.vint32(   5,     5,     5+1,    6+1,    6+2,     6+2,   6+1,   4+2 ),

  BaseShiftK00          = cms.vint32(  -8,    -8,    -7-1,   -8-1,   -8-1,   -9-1,  -10+0, -11+0 ),
  BaseShiftK10          = cms.vint32( -13,   -13,   -12-1,  -13-1,  -14+0,  -15+0,  -16-1, -17-1 ),
  BaseShiftK21          = cms.vint32( -13,   -13,   -12-1,  -13-1,  -13-1,  -14-1,  -14-1, -15-1 ),
  BaseShiftK31          = cms.vint32( -13,   -13,   -12-1,  -14+0,  -14+0,  -15+0,  -16-1, -17-1 ),

  BaseShiftR00          = cms.vint32(  -7,    -7,    -7+2,   -6+3,   -5+5,   -3+3,   -3+3,  -4+3 ),
  BaseShiftR11          = cms.vint32(   4,     4,     3+1,    4+2,    6+3,    9+1,    9+1,   7+2 ),
  BaseShiftInvR00Approx = cms.vint32( -16,   -16,   -17-2,  -16-1,  -15-5,  -17-3,  -17-3, -18-3 ),
  BaseShiftInvR11Approx = cms.vint32( -26,   -26,   -26-2,  -27-2,  -27-2,  -30+0,  -34-1, -34-2 ),
  BaseShiftInvR00Cor    = cms.vint32( -15,   -15,   -15,    -15,    -15,    -15,    -15,   -15   ),
  BaseShiftInvR11Cor    = cms.vint32( -15,   -15,   -15,    -15,    -15,    -15,    -15,   -15   ),
  BaseShiftInvR00       = cms.vint32( -16,   -16,   -26-2,  -25-3,  -25-4,  -26-3,  -27-3, -27-3 ),
  BaseShiftInvR11       = cms.vint32( -26,   -26,   -35-2,  -36-2,  -37-1,  -38-1,  -43-1, -43-2 ),
  
  BaseShiftC00          = cms.vint32(   7,     8+1,   6+2,    5+5,    9+2,    9+1,    6+2,   4+2 ),
  BaseShiftC01          = cms.vint32(   1,     3+1,   0+2,   -1+3,    2+2,    1+2,   -2+3,  -3+2 ),
  BaseShiftC11          = cms.vint32(  -4,    -2+1,  -5+2,   -6+2,   -5+2,   -6+3,   -8+3,  -9+3 ),
  BaseShiftC22          = cms.vint32(   8,     8+0,   6+0,    7+2,   10+1,   10+0,    7+0,   5+1 ),
  BaseShiftC23          = cms.vint32(   7,     8+1,   5+1,    4+2,    6+2,    7+0,    5+0,   4-1 ),
  BaseShiftC33          = cms.vint32(   7,     7+2,   4+2,    3+3,    5+2,    5+1,    4+1,   3+1 ),

  BaseShiftr02          = cms.vint32(  -7,    -7,    -4,     -3,     -3+2,   -3+2,   -2+1,  -2-1 ),
  BaseShiftr12          = cms.vint32(   1,     1,     7+2,    8,      9+1,   10+1,    9,     9+1 ),
  BaseShiftchi20        = cms.vint32( -13,   -13,   -13+3,  -12+2,  -11+2,  -11+2,  -10,   -10-2 ),
  BaseShiftchi21        = cms.vint32( -14,   -14,   -13+3,  -12+0,  -11+1,  -11+1,  -10-4, -10-4 )

)