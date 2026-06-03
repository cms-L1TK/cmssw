# functions to alter configurations

import FWCore.ParameterSet.Config as cms

# configures track finding s/w to behave as track finding f/w
def fwConfig(process):
  process.l1tTTTracksFromTrackletEmulation.Fakefit = True
  process.l1tTTTracksFromTrackletEmulation.RemovalType = ""
  process.l1tTTTracksFromTrackletEmulation.DoMultipleMatches = False
  process.l1tTTTracksFromTrackletEmulation.StoreTrackBuilderOutput = True
  process.TrackFindingTrackletSetup.DR.UseDTCStubs = False
  process.TrackFindingTrackletSetup.DR.UseTTStubs = False

# configures track finding s/w to behave as a subchain of processing steps
def reducedConfig(process):
  #fwConfig(process) disabled for now as new KF not working with reduced config
  process.TrackFindingTrackletSetup.Firmware.FreqBEHigh = 240 # Frequency of DTC & KF (determines truncation)
  process.TrackFindingTrackletSetup.KF.NumWorker = 1
  process.TrackFindingTrackletSetup.SeedTypes = cms.vstring( "L5L6" )
  process.TrackFindingTrackletSetup.SeedTypesSeedLayers = cms.PSet( L5L6 = cms.vint32( 5,  6 ) )
  process.TrackFindingTrackletSetup.SeedTypesProjectionLayers = cms.PSet( L5L6 = cms.vint32(  1,  2,  3,  4 ) )
  # this are tt::Setup::dtcId in order as in process.l1tTTTracksFromTrackletEmulation.processingModulesFile translated by 
  # reverssing naming logic described in L1FPGATrackProducer
  # TO DO: Eliminate cfg param IRChannelsIn by taking this info from Tracklet wiring map.
  process.TrackFindingTrackletSetup.IRChannelsIn = cms.vint32( 0, 1, 25, 2, 26, 4, 5, 29, 6, 30, 7, 31, 8, 9, 33 )
  process.l1tTTTracksFromTrackletEmulation.Reduced = True
  process.l1tTTTracksFromTrackletEmulation.DoMultipleMatches = False
  process.l1tTTTracksFromTrackletEmulation.memoryModulesFile = 'L1Trigger/TrackFindingTracklet/data/memorymodules_reduced.dat'
  process.l1tTTTracksFromTrackletEmulation.processingModulesFile = 'L1Trigger/TrackFindingTracklet/data/processingmodules_reduced.dat'
  process.l1tTTTracksFromTrackletEmulation.wiresFile = 'L1Trigger/TrackFindingTracklet/data/wires_reduced.dat'

# configures displaced track finding followed by Track Processing sim 
def simConfig(process):
  process.l1tTTTracksFromExtendedTrackletEmulation.Fakefit = True
  process.l1tTTTracksFromExtendedTrackletEmulation.RemovalType = ""
  process.l1tTTTracksFromExtendedTrackletEmulation.DoMultipleMatches = False
  process.TrackFindingTrackletSetup.TB.SeedTypes = ( "L1L2", "L2L3", "L3L4", "L5L6", "D1D2", "D3D4", "L1D1", "L2D1", "L2L3L4", "L4L5L6", "L2L3D1", "D1D2L2" )
  process.TrackFindingTrackletSetup.TB.SeedTypesSeedLayers = cms.PSet (
      L1L2   = cms.vint32(  1,  2     ),
      L2L3   = cms.vint32(  2,  3     ),
      L3L4   = cms.vint32(  3,  4     ),
      L5L6   = cms.vint32(  5,  6     ),
      D1D2   = cms.vint32( 11, 12     ),
      D3D4   = cms.vint32( 13, 14     ),
      L1D1   = cms.vint32(  1, 11     ),
      L2D1   = cms.vint32(  2, 11     ),
      L2L3L4 = cms.vint32(  2,  3,  4 ),
      L4L5L6 = cms.vint32(  4,  5,  6 ),
      L2L3D1 = cms.vint32(  2,  3, 11 ),
      D1D2L2 = cms.vint32( 11, 12,  2 )
  )
  process.TrackFindingTrackletSetup.TB.SeedTypesProjectionLayers = cms.PSet (
      L1L2   = cms.vint32(  3,  4,  5,  6, 11, 12, 13, 14 ),
      L2L3   = cms.vint32(  1,  4,  5,  6, 11, 12, 13, 14 ),
      L3L4   = cms.vint32(  1,  2,  5,  6, 11, 12 ),
      L5L6   = cms.vint32(  1,  2,  3,  4 ),
      D1D2   = cms.vint32(  1,  2, 13, 14, 15 ),
      D3D4   = cms.vint32(  1, 11, 12, 15 ),
      L1D1   = cms.vint32( 12, 13, 14, 15 ),
      L2D1   = cms.vint32(  1, 12, 13, 14 ),
      L2L3L4 = cms.vint32(  1,  5,  6, 11, 12, 13 ),
      L4L5L6 = cms.vint32(  1,  2,  3 ),
      L2L3D1 = cms.vint32(  1, 12, 13, 14 ),
      D1D2L2 = cms.vint32(  1, 13, 14 )
  )
  process.TrackFindingTrackletSetup.TM.MuxOrder = ( "L1L2", "L2L3", "L1D1", "D1D2", "D3D4", "L2D1", "L2L3D1", "D1D2L2", "L3L4", "L2L3L4", "L5L6", "L4L5L6" )

# configures pure tracklet algorithm (as opposed to Hybrid algorithm)
def trackletConfig(process):
  process.l1tTTTracksFromTrackletEmulation.fitPatternFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/fitpattern.txt') 
