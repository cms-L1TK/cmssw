import FWCore.ParameterSet.Config as cms

from SimTracker.TrackTriggerAssociation.Clean_cfi import Clean_params

CleanTV = cms.EDProducer( 'tt::CleanTV', Clean_params )