import FWCore.ParameterSet.Config as cms

from SimTracker.TrackTriggerAssociation.Clean_cfi import Clean_params

CleanAssoc = cms.EDProducer( 'tt::CleanAssoc', Clean_params )