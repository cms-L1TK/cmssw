#include "L1Trigger/TrackFindingTracklet/interface/TrackBuilderChannel.h"

#include <vector>

using namespace std;
using namespace edm;

namespace trackFindingTracklet {

  TrackBuilderChannel::TrackBuilderChannel(const edm::ParameterSet& iConfig, const Setup* setup)
      : setup_(setup),
        summerChain_(iConfig.getParameter<bool>("SummerChain")),
        useDuplicateRemoval_(iConfig.getParameter<bool>("UseDuplicateRemoval")),
        boundaries_(iConfig.getParameter<vector<double>>("PtBoundaries")),
        seedTypeNames_(iConfig.getParameter<vector<string>>("SeedTypes")),
        numSeedTypes_(seedTypeNames_.size()),
        maxNumProjectionLayers_(iConfig.getParameter<int>("MaxNumProjectionLayers"))
  {
    if (summerChain_)
      numChannels_ = 1;
    else if (useDuplicateRemoval_)
      numChannels_ = 2 * boundaries_.size();
    else 
      numChannels_ = numSeedTypes_;
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
  bool TrackBuilderChannel::channelId(const TTTrackRef& ttTrackRef, int& channelId) {
    if (summerChain_) {
      if (ttTrackRef->trackSeedType() != 0) {
        cms::Exception exception("logic_error");
        exception << "TTTracks form seed type L1L2 not supported in summer chain configuration.";
        exception.addContext("trackFindingTracklet:TrackBuilderChannel:channelId");
        throw exception;
      }
      channelId = ttTrackRef->phiSector();
      return true;
    }
    if (!useDuplicateRemoval_) {
      channelId = ttTrackRef->phiSector() * numSeedTypes_ + ttTrackRef->trackSeedType();
      return true;
    }
    const double pt = ttTrackRef->momentum().perp();
    channelId = -1;
    for (double boundary : boundaries_) {
      if (pt < boundary)
        break;
      else
        channelId++;
    }
    if (channelId == -1)
      return false;
    channelId = ttTrackRef->rInv() < 0. ? channelId : numChannels_ - channelId - 1;
    channelId += ttTrackRef->phiSector() * numChannels_;
    return true;
  }

  // sets layerId of given TTStubRef and TTTrackRef, returns false if seeed stub
  bool TrackBuilderChannel::layerId(const TTTrackRef& ttTrackRef, const TTStubRef& ttStubRef, int& layerId) {
    layerId = -1;
    const int seedType = ttTrackRef->trackSeedType();
    if (seedType < 0 || seedType >= numSeedTypes_) {
      cms::Exception exception("logic_error");
      exception.addContext("trackFindingTracklet::TrackBuilderChannel::layerId");
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
      exception.addContext("trackFindingTracklet::TrackBuilderChannel::layerId");
      exception << "TTStub from layer " << layer << " (barrel: 1-6; discs: 11-15) from seed type " << name << " not supported.";
      throw exception;
    }
    layerId = distance(projectingLayers.begin(), pos);
     return true;
   }

}  // namespace trackFindingTracklet