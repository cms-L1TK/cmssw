#include "L1Trigger/TrackFindingTracklet/interface/ChannelAssignment.h"

#include <vector>

using namespace std;
using namespace edm;
using namespace tt;

namespace trklet {

  ChannelAssignment::ChannelAssignment(const edm::ParameterSet& iConfig, const Setup* setup)
      : setup_(setup),
        useDuplicateRemoval_(iConfig.getParameter<bool>("UseDuplicateRemoval")),
        boundaries_(iConfig.getParameter<vector<double>>("PtBoundaries")),
        seedTypeNames_(iConfig.getParameter<vector<string>>("SeedTypes")),
        numSeedTypes_(seedTypeNames_.size()),
        numChannels_(useDuplicateRemoval_ ? 2 * boundaries_.size() : numSeedTypes_),
        maxNumProjectionLayers_(iConfig.getParameter<int>("MaxNumProjectionLayers")),
        channelEncoding_(iConfig.getParameter<int>("IRChannelsIn")) {
    const ParameterSet& pSetSeedTypesSeedLayers = iConfig.getParameter<ParameterSet>("SeedTypesSeedLayers");
    const ParameterSet& pSetSeedTypesProjectionLayers = iConfig.getParameter<ParameterSet>("SeedTypesProjectionLayers");
    seedTypesSeedLayers_.reserve(numSeedTypes_);
    seedTypesProjectionLayers_.reserve(numSeedTypes_);
    for (const string& s : seedTypeNames_) {
      seedTypesSeedLayers_.emplace_back(pSetSeedTypesSeedLayers.getParameter<vector<int>>(s));
      seedTypesProjectionLayers_.emplace_back(pSetSeedTypesProjectionLayers.getParameter<vector<int>>(s));
    }
  }

  // sets channelId of given TTTrackRef, return false if track outside pt range
  bool ChannelAssignment::channelId(const TTTrackRef& ttTrackRef, int& channelId) {
    if (!useDuplicateRemoval_) {
      const int seedType = ttTrackRef->trackSeedType();
      if (seedType >= numSeedTypes_) {
        cms::Exception exception("logic_error");
        exception << "TTTracks form seed type" << seedType << " not in supported list: (";
        for (const auto& s : seedTypeNames_)
          exception << s << " ";
        exception << ").";
        exception.addContext("trackFindingTracklet:ChannelAssignment:channelId");
        throw exception;
      }
      channelId = ttTrackRef->phiSector() * numSeedTypes_ + seedType;
      return true;
    }
    const double pt = 2. * setup_->invPtToDphi() / abs(rInv);
    channelId = -1;
    for (double boundary : boundaries_) {
      if (pt < boundary)
        break;
      else
        channelId++;
    }
    if (channelId == -1)
      return false;
    channelId = rInv < 0. ? channelId : numChannels_ - channelId - 1;
    channelId += phiSector * numChannels_;
    return true;
  }

  // sets layerId of given TTStubRef and TTTrackRef, returns false if seeed stub
  bool ChannelAssignment::layerId(const TTTrackRef& ttTrackRef, const TTStubRef& ttStubRef, int& layerId) {
    layerId = -1;
    if (seedType < 0 || seedType >= numSeedTypes_) {
      cms::Exception exception("logic_error");
      exception.addContext("trackFindingTracklet::ChannelAssignment::layerId");
      exception << "TTTracks with with seed type " << seedType << " not supported.";
      throw exception;
    }
    const int layer = setup_->layerId(ttStubRef);
    const vector<int>& seedingLayers = seedTypesSeedLayers_[seedType];
    if (find(seedingLayers.begin(), seedingLayers.end(), layer) != seedingLayers.end())
      return false;
    const vector<int>& projectingLayers = seedTypesProjectionLayers_[seedType];
    const auto pos = find(projectingLayers.begin(), projectingLayers.end(), layer);
    if (pos == projectingLayers.end()) {
      const string& name = seedTypeNames_[seedType];
      cms::Exception exception("logic_error");
      exception.addContext("trackFindingTracklet::ChannelAssignment::layerId");
      exception << "TTStub from layer " << layer << " (barrel: 1-6; discs: 11-15) from seed type " << name
                << " not supported.";
      throw exception;
    }
    layerId = distance(projectingLayers.begin(), pos);
    return true;
  }

  // sets layerId of given TTStubRef and TTTrackRef, returns false if seeed stub
  bool TrackBuilderChannel::layerId(const TTTrackRef& ttTrackRef, const TTStubRef& ttStubRef, int& layerId) const {
    return this->layerId(ttTrackRef->trackSeedType(), ttStubRef, layerId);
  }

  // return tracklet layerId (barrel: [0-5], endcap: [6-10]) for given TTStubRef
  int TrackBuilderChannel::trackletLayerId(const TTStubRef& ttStubRef) const {
    static constexpr int offsetBarrel = 1;
    static constexpr int offsetDisks = 5;
    return setup_->layerId(ttStubRef) - (setup_->barrel(ttStubRef) ? offsetBarrel : offsetDisks);
  }

}  // namespace trklet