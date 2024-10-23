# ESProducer to provide and calculate dataformats used by Kalman Filter emulator enabling tuning of bit widths

import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackerTFP.KalmanFilterFormats_cfi import TrackTriggerKalmanFilterFormats_params

TrackTriggerKalmanFilterFormats = cms.ESProducer("trackerTFP::ProducerFormatsKF", TrackTriggerKalmanFilterFormats_params)