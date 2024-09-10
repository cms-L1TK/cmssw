# configuration of internal KF variable bases which can be shifted by powers of 2 w.r.t. KF output track parameter
# TrackerTFPProducer_params.PrintKFDebug printouts unused MSB for each variable, so that one could consider decreasing the basseshift by that amount
# numerical instabillity (negative C00, C11, C22, C33) requires smaller baseshifts of related variables (rx, Sxx, Kxx, Rxx, invRxx)
# if a variable overflows an Exception will be thrown and the corresponding baseshift needs to be increased.
import FWCore.ParameterSet.Config as cms

TrackTriggerKalmanFilterFormats_params = cms.PSet (

  EnableIntegerEmulation = cms.bool( True ),

  WidthR00 = cms.int32( 20 ),
  WidthR11 = cms.int32( 20 ),
  
  WidthC00 = cms.int32( 20 ),
  WidthC01 = cms.int32( 20 ),
  WidthC11 = cms.int32( 20 ),
  WidthC22 = cms.int32( 20 ),
  WidthC23 = cms.int32( 20 ),
  WidthC33 = cms.int32( 20 ),

  BaseShiftx0           = cms.int32(   0 ),
  BaseShiftx1           = cms.int32(  -6 ),
  BaseShiftx2           = cms.int32(   0 ),
  BaseShiftx3           = cms.int32(  -1 ),
  BaseShiftr0           = cms.int32(  -7 ),
  BaseShiftr1           = cms.int32(   0 ),

  BaseShiftS00          = cms.int32(  -2 ),
  BaseShiftS01          = cms.int32( -10 ),
  BaseShiftS12          = cms.int32(   1 ),
  BaseShiftS13          = cms.int32(  -1 ),

  BaseShiftK00          = cms.int32(  -7 ),
  BaseShiftK10          = cms.int32( -13 ),
  BaseShiftK21          = cms.int32( -13 ),
  BaseShiftK31          = cms.int32( -13 ),

  BaseShiftR00          = cms.int32(  -3 ),
  BaseShiftR11          = cms.int32(   6 ),
  BaseShiftInvR00Approx = cms.int32( -41 ),
  BaseShiftInvR11Approx = cms.int32( -50 ),
  BaseShiftInvR00Cor    = cms.int32( -15 ),
  BaseShiftInvR11Cor    = cms.int32( -15 ),
  BaseShiftInvR00       = cms.int32( -32 ),
  BaseShiftInvR11       = cms.int32( -41 ),
  
  BaseShiftC00         = cms.int32(   8  ),
  BaseShiftC01         = cms.int32(   3 ),
  BaseShiftC11         = cms.int32(  -4 ),
  BaseShiftC22         = cms.int32(   5 ),
  BaseShiftC23         = cms.int32(   6 ),
  BaseShiftC33         = cms.int32(   5 ),

  BaseShiftr02          = cms.int32(  -2 ),
  BaseShiftr12          = cms.int32(  10 ),
  BaseShiftchi20        = cms.int32( -10 ),
  BaseShiftchi21        = cms.int32( -10 )

)