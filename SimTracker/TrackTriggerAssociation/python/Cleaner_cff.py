import FWCore.ParameterSet.Config as cms

from SimTracker.TrackTriggerAssociation.Cleaner_cfi import Cleaner_params

Cleaner = cms.EDProducer( 'tt::Cleaner', Cleaner_params )
