# Produce L1 tracks with TMTT C++ emulation
import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackTrigger.Setup_cff import TrackTriggerSetup
from L1Trigger.TrackerTFP.Producer_cfi import TrackerTFPProducer_params
from L1Trigger.TrackerTFP.DataFormats_cff import TrackTriggerDataFormats
from L1Trigger.TrackerTFP.LayerEncoding_cff import TrackTriggerLayerEncoding
from L1Trigger.TrackerTFP.KalmanFilterFormats_cff import TrackTriggerKalmanFilterFormats

TrackerTFPProducerPP = cms.EDProducer( 'trackerTFP::ProducerPP', TrackerTFPProducer_params )
TrackerTFPProducerGP = cms.EDProducer( 'trackerTFP::ProducerGP', TrackerTFPProducer_params )
TrackerTFPProducerHT = cms.EDProducer( 'trackerTFP::ProducerHT', TrackerTFPProducer_params )
TrackerTFPProducerCTB = cms.EDProducer( 'trackerTFP::ProducerCTB', TrackerTFPProducer_params )
TrackerTFPProducerKF = cms.EDProducer( 'trackerTFP::ProducerKF', TrackerTFPProducer_params )
TrackerTFPProducerDR = cms.EDProducer( 'trackerTFP::ProducerDR', TrackerTFPProducer_params )
TrackerTFPProducerTFP = cms.EDProducer( 'trackerTFP::ProducerTFP', TrackerTFPProducer_params )