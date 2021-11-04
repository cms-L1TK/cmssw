#ifndef L1Trigger_TrackFindingTracklet_TrackBuilderChannel_h
#define L1Trigger_TrackFindingTracklet_TrackBuilderChannel_h

#include "FWCore/Framework/interface/data_default_record_trait.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackBuilderChannelRcd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"

#include <vector>

namespace trackFindingTracklet {

  /*! \class  trackFindingTracklet::TrackBuilderChannel
   *  \brief  Class to assign tracklet tracks ans stubs to channel
   *          based on their Pt or seed type
   *  \author Thomas Schuh
   *  \date   2020, Nov; updated 2021 Oct
   */
  class TrackBuilderChannel {
  public:
    TrackBuilderChannel() {}
    TrackBuilderChannel(const edm::ParameterSet& iConfig, const tt::Setup* setup);
    ~TrackBuilderChannel() {}
    // sets channelId of given TTTrackRef, return false if track outside pt range
    bool channelId(const TTTrackRef& ttTrackRef, int& channelId);
    // number of used channels
    int numChannels() const { return numChannels_; }
    // sets layerId of given TTStubRef and TTTrackRef, returns false if seeed stub
    bool layerId(const TTTrackRef& ttTrackRef, const TTStubRef& ttStubRef, int& layerId);
    // max number layers a sedd type may project to
    int maxNumProjectionLayers() const { return maxNumProjectionLayers_; }

  private:
    // helper class to store configurations
    const tt::Setup* setup_;
    // reduce l1 tracking to summer chain configuration
    bool summerChain_;
    // use tracklet seed type as channel id if False, binned track pt used if True
    bool useDuplicateRemoval_;
    // pt Boundaries in GeV, last boundary is infinity
    std::vector<double> boundaries_;
    // seed type names
    std::vector<std::string> seedTypeNames_;
    // number of used seed types in tracklet algorithm
    int numSeedTypes_;
    // number of used channels
    int numChannels_;
    // max number layers a sedd type may project to
    int maxNumProjectionLayers_;
    // seeding layers of seed types using default layer id [barrel: 1-6, discs: 11-15]
    std::vector<std::vector<int>> seedTypesSeedLayers_;
    // layers a seed types can project to using default layer id [barrel: 1-6, discs: 11-15]
    std::vector<std::vector<int>> seedTypesProjectionLayers_;
  };

}  // namespace trackFindingTracklet

EVENTSETUP_DATA_DEFAULT_RECORD(trackFindingTracklet::TrackBuilderChannel, trackFindingTracklet::TrackBuilderChannelRcd);

#endif