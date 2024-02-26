#include "L1Trigger/TrackFindingTracklet/interface/DR.h"

#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace edm;
using namespace tt;
using namespace trackerTFP;

namespace trklet {

  DR::DR(const ParameterSet& iConfig,
         const Setup* setup,
         const DataFormats* dataFormats,
         const ChannelAssignment* channelAssignment,
         int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        channelAssignment_(channelAssignment),
        region_(region) {}

  // read in and organize input tracks and stubs
  void DR::consume(const StreamsTrack& streamsTrack, const StreamsStub& streamsStub) {
    auto nonNullTrack = [](int& sum, const FrameTrack& frame) { return sum += (frame.first.isNonnull() ? 1 : 0); };
    auto nonNullStub = [](int& sum, const FrameStub& frame) { return sum += (frame.first.isNonnull() ? 1 : 0); };
    // count tracks and stubs and reserve corresponding vectors
    int sizeStubs(0);
    const int offset = region_ * setup_->numLayers();
    const StreamTrack& streamTrack = streamsTrack[region_];
    input_.reserve(streamTrack.size());
    const int sizeTracks = accumulate(streamTrack.begin(), streamTrack.end(), 0, nonNullTrack);
    for (int layer = 0; layer < setup_->numLayers(); layer++) {
      const StreamStub& streamStub = streamsStub[offset + layer];
      sizeStubs += accumulate(streamStub.begin(), streamStub.end(), 0, nonNullStub);
    }
    tracks_.reserve(sizeTracks);
    stubs_.reserve(sizeStubs);
    // transform input data into handy structs
    for (int frame = 0; frame < (int)streamTrack.size(); frame++) {
      const FrameTrack& frameTrack = streamTrack[frame];
      if (frameTrack.first.isNull()) {
        input_.push_back(nullptr);
        continue;
      }
      vector<Stub*> stubs;
      stubs.reserve(setup_->numLayers());
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        const FrameStub& frameStub = streamsStub[offset + layer][frame];
        if (frameStub.first.isNull())
          continue;
        double dZ, dPhi, z, phi;
        TTBV ttBV = frameStub.second;
        dataFormats_->format(Variable::dZ, Process::ctb).extract(ttBV, dZ);
        dataFormats_->format(Variable::dPhi, Process::ctb).extract(ttBV, dPhi);;
        dataFormats_->format(Variable::z, Process::ctb).extract(ttBV, z);
        dataFormats_->format(Variable::phi, Process::ctb).extract(ttBV, phi);
        ttBV >>= dataFormats_->format(Variable::r, Process::ctb).width();
        const TTBV stubId(ttBV, channelAssignment_->widthStubId(), 0, true);
        const TTBV ctb(frameStub.second, dataFormats_->width(Variable::dZ, Process::ctb) + dataFormats_->width(Variable::dPhi, Process::ctb) + dataFormats_->width(Variable::z, Process::ctb) + dataFormats_->width(Variable::phi, Process::ctb) + dataFormats_->width(Variable::r, Process::ctb), 0);
        const FrameStub fs(frameStub.first, "1" + ctb.str());
        stubs_.emplace_back(fs, stubId.val(), layer, phi, z, dPhi, dZ);
        stubs.push_back(&stubs_.back());
      }
      tracks_.emplace_back(frameTrack, stubs);
      input_.push_back(&tracks_.back());
    }
    // remove all gaps between end and last track
    for (auto it = input_.end(); it != input_.begin();)
      it = (*--it) ? input_.begin() : input_.erase(it);
  }

  // fill output products
  void DR::produce(StreamsStub& accpetedStubs,
                   StreamsTrack& acceptedTracks,
                   StreamsStub& lostStubs,
                   StreamsTrack& lostTracks) {
    const int offset = region_ * setup_->numLayers();
    // remove duplicated tracks, no merge of stubs, one stub per layer expected
    vector<Track*> cms(channelAssignment_->numComparisonModules(), nullptr);
    vector<Track*>& tracks = input_;
    for (Track*& track : tracks) {
      if (!track)
        // gaps propagate trough chain and appear in output stream
        continue;
      for (Track*& trackCM : cms) {
        if (!trackCM) {
          // tracks used in CMs don't propagate trough chain and do not appear in output stream unaltered
          trackCM = track;
          track = nullptr;
          break;
        }
        if (equalEnough(track, trackCM)) {
          // tracks compared in CMs propagate trough chain and appear in output stream as gap if identified as duplicate or unaltered elsewise
          if (better(track, trackCM))
            trackCM = track;
          track = nullptr;
          break;
        }
      }
    }
    // remove first number of CMs nullptr
    const int gaps = min((int)tracks.size(), channelAssignment_->numComparisonModules());
    tracks.erase(tracks.begin(), next(tracks.begin(), gaps));
    // add cms tracks
    tracks.insert(tracks.end(), cms.begin(), cms.end());
    // remove all gaps between end and last track
    for (auto it = tracks.end(); it != tracks.begin();)
      it = (*--it) ? tracks.begin() : tracks.erase(it);
    // store output
    StreamTrack& streamTrack = acceptedTracks[region_];
    streamTrack.reserve(tracks.size());
    for (int layer = 0; layer < setup_->numLayers(); layer++)
      accpetedStubs[offset + layer].reserve(tracks.size());
    for (Track* track : tracks) {
      if (!track) {
        streamTrack.emplace_back(FrameTrack());
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          accpetedStubs[offset + layer].emplace_back(FrameStub());
        continue;
      }
      streamTrack.push_back(track->frame_);
      TTBV hitPattern(0, setup_->numLayers());
      for (Stub* stub : track->stubs_) {
        hitPattern.set(stub->channel_);
        accpetedStubs[offset + stub->channel_].push_back(stub->frame_);
      }
      for (int layer : hitPattern.ids(false))
        accpetedStubs[offset + layer].emplace_back(FrameStub());
    }
  }

  // compares two tracks, returns true if those are considered duplicates
  bool DR::equalEnough(Track* t0, Track* t1) const {
    int same(0);
    for (int layer = 0; layer < setup_->numLayers(); layer++) {
      auto onLayer = [layer](Stub* stub) { return stub->channel_ == layer; };
      const auto s0 = find_if(t0->stubs_.begin(), t0->stubs_.end(), onLayer);
      const auto s1 = find_if(t1->stubs_.begin(), t1->stubs_.end(), onLayer);
      if (s0 != t0->stubs_.end() && s1 != t1->stubs_.end() && **s0 == **s1)
        same++;
    }
    return same >= channelAssignment_->minIdenticalStubs();
  }

  bool DR::better(Track* lhs, Track* rhs) const {
    if (lhs->nConsistentStubs_ > rhs->nConsistentStubs_)
      return lhs;
    else if (lhs->nConsistentStubs_ == rhs->nConsistentStubs_) {
      if (lhs->chi2_ < rhs->chi2_)
        return lhs;
    }
    return rhs;
  }

}  // namespace trklet