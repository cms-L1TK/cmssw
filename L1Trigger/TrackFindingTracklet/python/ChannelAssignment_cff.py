import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackFindingTracklet.ChannelAssignment_cfi import ChannelAssignment_params

ChannelAssignment = cms.ESProducer("trackFindingTracklet::ProducerChannelAssignment", ChannelAssignment_params)