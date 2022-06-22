import FWCore.ParameterSet.Config as cms

from SimTracker.TrackTriggerAssociation.Clean_cfi import Clean_params

CleanTP = cms.EDProducer( 'tt::CleanTP', Clean_params )