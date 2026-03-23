import FWCore.ParameterSet.Config as cms
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import TTStubAlgorithm_official_Phase2TrackerDigi_ as _TTStubAlgorithm

TTStubAlgorithm_official_Phase2TrackerDigi_ = _TTStubAlgorithm.clone(
        cosmics = cms.bool(True),
)

# Set the preferred hit matching algorithms.
# We prefer the global geometry algorithm for now in order not to break
# anything. Override with process.TTStubAlgorithm_PSimHit_ = ...,
# etc. in your configuration.
TTStubAlgorithm_Phase2TrackerDigi_ = cms.ESPrefer("TTStubAlgorithm_official_Phase2TrackerDigi_")
