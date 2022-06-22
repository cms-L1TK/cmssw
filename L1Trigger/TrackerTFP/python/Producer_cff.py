# Produce L1 tracks with TMTT C++ emulation
import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackTrigger.Setup_cff import TrackTriggerSetup
from L1Trigger.TrackerTFP.Producer_cfi import TrackerTFPProducer_params
from L1Trigger.TrackerTFP.DataFormats_cff import TrackTriggerDataFormats
from L1Trigger.TrackerTFP.LayerEncoding_cff import TrackTriggerLayerEncoding
from L1Trigger.TrackerTFP.KalmanFilterFormats_cff import TrackTriggerKalmanFilterFormats

TrackerTFPProducerGP = cms.EDProducer( 'trackerTFP::ProducerGP', TrackerTFPProducer_params )
TrackerTFPProducerHT = cms.EDProducer( 'trackerTFP::ProducerHT', TrackerTFPProducer_params )
TrackerTFPProducerMHT = cms.EDProducer( 'trackerTFP::ProducerMHT', TrackerTFPProducer_params )
TrackerTFPProducerZHT = cms.EDProducer( 'trackerTFP::ProducerZHT', TrackerTFPProducer_params )
TrackerTFPProducerTTFound = cms.EDProducer( 'trackerTFP::ProducerTTTrackFound', TrackerTFPProducer_params )
TrackerTFPProducerKFin = cms.EDProducer( 'trackerTFP::ProducerKFin', TrackerTFPProducer_params )
TrackerTFPProducerKF = cms.EDProducer( 'trackerTFP::ProducerKF', TrackerTFPProducer_params )
TrackerTFPProducerTTFitted = cms.EDProducer( 'trackerTFP::ProducerTTTrackFitted', TrackerTFPProducer_params )
TrackerTFPProducerAS = cms.EDProducer( 'trackerTFP::ProducerAS', TrackerTFPProducer_params )