#include "L1Trigger/TrackFindingTracklet/interface/DuplicateRemoval.h"

#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace edm;
using namespace tt;

namespace trklet {

  DuplicateRemoval::DuplicateRemoval(const ParameterSet& iConfig,
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
  void DuplicateRemoval::consume(const StreamsTrack& streamsTrack, const StreamsStub& streamsStub) {
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
        // create handy stub
        const StubTM stubTM(frameStub, dataFormats_);
        const StubDR stubDR(stubTM);
        stubs_.emplace_back(stubDR.frame(), stubTM.stubId(), layer);
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
  void DuplicateRemoval::produce(StreamsTrack& streamsTrack, StreamsStub& streamsStub) {
    const int offset = region_ * setup_->numLayers();
    // remove duplicated tracks, no merge of stubs, one stub per layer expected
    vector<Track*> cms(channelAssignment_->numComparisonModules(), nullptr);
    for (Track*& track : input_) {
      if (!track)
        // gaps propagate through chain and appear in output stream
        continue;
      for (Track*& trackCM : cms) {
        if (!trackCM) {
          // tracks used in CMs propagate through chain and do appear in output stream unaltered
          trackCM = track;
          break;
        }
        if (equalEnough(track, trackCM)) {
          // tracks compared in CMs propagate through chain and appear in output stream as gap if identified as duplicate or unaltered elsewise
          track = nullptr;
          break;
        }
      }
    }
    // remove all gaps between end and last track
    for (auto it = input_.end(); it != input_.begin();)
      it = (*--it) ? input_.begin() : input_.erase(it);
    // store output
    StreamTrack& streamTrack = streamsTrack[region_];
    streamTrack.reserve(input_.size());
    for (int layer = 0; layer < setup_->numLayers(); layer++)
      streamsStub[offset + layer].reserve(input_.size());
    for (Track* track : input_) {
      if (!track) {
        streamTrack.emplace_back(FrameTrack());
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          streamsStub[offset + layer].emplace_back(FrameStub());
        continue;
      }
      streamTrack.push_back(track->frame_);
      TTBV hitPattern(0, setup_->numLayers());
      for (Stub* stub : track->stubs_) {
        hitPattern.set(stub->layer_);
        streamsStub[offset + stub->layer_].push_back(stub->frame_);
      }
      for (int layer : hitPattern.ids(false))
        streamsStub[offset + layer].emplace_back(FrameStub());
    }
  }

  // compares two tracks, returns true if those are considered duplicates
  bool DuplicateRemoval::equalEnough(Track* t0, Track* t1) const {
    int same(0);
    for (int layer = 0; layer < setup_->numLayers(); layer++) {
      auto onLayer = [layer](Stub* stub) { return stub->layer_ == layer; };
      const auto s0 = find_if(t0->stubs_.begin(), t0->stubs_.end(), onLayer);
      const auto s1 = find_if(t1->stubs_.begin(), t1->stubs_.end(), onLayer);
      if (s0 != t0->stubs_.end() && s1 != t1->stubs_.end() && **s0 == **s1)
        same++;
    }
    return same >= channelAssignment_->minIdenticalStubs();
  }

}  // namespace trklet
