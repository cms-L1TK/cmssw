# defines PSet to assign tracklet tracks and stubs to output channel based on their Pt or seed type as well as DTC stubs to input channel
import FWCore.ParameterSet.Config as cms

ChannelAssignment_params = cms.PSet (

  WidthLayerId = cms.int32(4),
  WidthStubId = cms.int32(10),
  WidthPSTilt = cms.int32(1),
  PtBoundaries = cms.vdouble( 3.0, 5.0, 8.0, 12.0, 24.0 ), # positive pt Boundaries in GeV (symmetric negatives are assumed), first boundary is pt cut, last boundary is infinity, defining ot bins used by DR
  DRinDepthMemory = cms.int32(32),
  DRNumComparisonModules = cms.int32(16),
  DRMinIdenticalStubs = cms.int32(3),
  KFinDepthMemory = cms.int32(128),

  SeedTypes = cms.vstring( "L1L2", "L2L3", "L3L4", "L5L6", "D1D2", "D3D4", "L1D1", "L2D1" ), # seed types used in tracklet algorithm (position gives int value)

  SeedTypesSeedLayers = cms.PSet (      # seeding layers of seed types using default layer id [barrel: 1-6, discs: 11-15]
    L1L2 = cms.vint32(  1,  2 ),
    L2L3 = cms.vint32(  2,  3 ),
    L3L4 = cms.vint32(  3,  4 ),
    L5L6 = cms.vint32(  5,  6 ),
    D1D2 = cms.vint32( 11, 12 ),
    D3D4 = cms.vint32( 13, 14 ),
    L1D1 = cms.vint32(  1, 11 ),
    L2D1 = cms.vint32(  2, 11 )
  ),
  SeedTypesProjectionLayers = cms.PSet ( # layers a seed types can project to using default layer id [barrel: 1-6, discs: 11-15]
    L1L2 = cms.vint32(  3,  4,  5,  6, 11, 12, 13, 14 ),
    L2L3 = cms.vint32(  1,  4,  5,  6, 11, 12, 13, 14 ),
    L3L4 = cms.vint32(  1,  2,  5,  6, 11, 12 ),
    L5L6 = cms.vint32(  1,  2,  3,  4 ),
    D1D2 = cms.vint32(  1,  2, 13, 14, 15 ),
    D3D4 = cms.vint32(  1, 11, 12, 15 ),
    L1D1 = cms.vint32( 12, 13, 14, 15 ),
    L2D1 = cms.vint32(  1, 12, 13, 14 )
  ),

  IRChannelsIn = cms.vint32( range(0, 48) ) # vector of DTC id indexed by connected IR module id (from order in processingmodules.dat)

)
