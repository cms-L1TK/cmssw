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
        channelEncoding_(iConfig.getParameter<vector<int>>("IRChannelsIn")) {
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
  bool ChannelAssignment::channelId(const TTTrackRef& ttTrackRef, int& channelId) const {
    const int phiSector = ttTrackRef->phiSector();
    bool valid(true);
    channelId = -1;
    if (!useDuplicateRemoval_) {
      const int seedType = ttTrackRef->trackSeedType();
      checkSeedType(seedType);
      channelId = seedType;
    } else {
      const double rInv = ttTrackRef->rInv();
      const double pt = 2. * setup_->invPtToDphi() / abs(rInv);
      for (double boundary : boundaries_) {
        if (pt < boundary)
          break;
        else
          channelId++;
      }
      if (channelId == -1)
        valid = false;
      channelId = rInv < 0. ? channelId : numChannels_ - channelId - 1;
    }
    channelId += phiSector * numChannels_;
    return valid;
  }

  // sets layerId of given TTStubRef and seedType, returns false if seeed stub
  bool ChannelAssignment::layerId(int seedType, const TTStubRef& ttStubRef, int& layerId) const {
    layerId = -1;
    checkSeedType(seedType);
    const int layer = setup_->layerId(ttStubRef);
    const vector<int>& seedingLayers = seedTypesSeedLayers_[seedType];
    if (find(seedingLayers.begin(), seedingLayers.end(), layer) != seedingLayers.end())
      return false;
    const vector<int>& projectingLayers = seedTypesProjectionLayers_[seedType];
    const auto pos = find(projectingLayers.begin(), projectingLayers.end(), layer);
    if (pos == projectingLayers.end()) {
      const string& name = seedTypeNames_[seedType];
      cms::Exception exception("logic_error");
      exception.addContext("ChannelAssignment::ChannelAssignment::layerId");
      exception << "TTStub from layer " << layer << " (barrel: 1-6; discs: 11-15) from seed type " << name
                << " not supported.";
      throw exception;
    }
    layerId = distance(projectingLayers.begin(), pos);
    return true;
  }

  // return tracklet layerId (barrel: [0-5], endcap: [6-10]) for given TTStubRef
  int ChannelAssignment::trackletLayerId(const TTStubRef& ttStubRef) const {
    static constexpr int offsetBarrel = 1;
    static constexpr int offsetDisks = 5;
    return setup_->layerId(ttStubRef) - (setup_->barrel(ttStubRef) ? offsetBarrel : offsetDisks);
  }

  // checks is seedType is supported
  void ChannelAssignment::checkSeedType(int seedType) const {
    if (seedType >= 0 && seedType < numSeedTypes_)
      return;
    cms::Exception exception("logic_error");
    exception << "TTTracks form seed type" << seedType << " not in supported list: (";
    for (const auto& s : seedTypeNames_)
      exception << s << " ";
    exception << ").";
    exception.addContext("ChannelAssignment:checkSeedType:channelId");
    throw exception;
  }

}  // namespace trklet