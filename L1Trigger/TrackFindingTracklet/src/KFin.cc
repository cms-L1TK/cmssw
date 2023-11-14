#include "L1Trigger/TrackFindingTracklet/interface/KFin.h"

#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace edm;
using namespace tt;
//using namespace trackerTFP;

namespace trklet {

  KFin::KFin(const ParameterSet& iConfig, const Setup* setup, const ChannelAssignment* channelAssignment, int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        channelAssignment_(channelAssignment),
        region_(region),
        input_(channelAssignment_->numNodesDR()) {}

  // read in and organize input tracks and stubs
  void KFin::consume(const StreamsTrack& streamsTrack, const StreamsStub& streamsStub) {
    const int offsetTrack = region_ * channelAssignment_->numNodesDR();
    auto valid = [](int& sum, const FrameTrack& frame) { return sum += (frame.first.isNonnull() ? 1 : 0); };
    // count tracks and stubs and reserve corresponding vectors
    int sizeTracks(0);
    for (int channel = 0; channel < channelAssignment_->numNodesDR(); channel++) {
      const int streamTrackId = offsetTrack + channel;
      const StreamTrack& streamTrack = streamsTrack[streamTrackId];
      input_[channel].reserve(streamTrack.size());
      sizeTracks += accumulate(streamTrack.begin(), streamTrack.end(), 0, valid);
    }
    tracks_.reserve(sizeTracks);
    // transform input data into handy structs
    for (int channel = 0; channel < channelAssignment_->numNodesDR(); channel++) {
      vector<Track*>& input = input_[channel];
      const int streamTrackId = offsetTrack + channel;
      const int offsetStub = streamTrackId * setup_->numLayers();
      const StreamTrack& streamTrack = streamsTrack[streamTrackId];
      for (int frame = 0; frame < (int)streamTrack.size(); frame++) {
        const FrameTrack& frameTrack = streamTrack[frame];
        if (frameTrack.first.isNull()) {
          input.push_back(nullptr);
          continue;
        }
        vector<FrameStub> stubs;
        stubs.reserve(setup_->numLayers());
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          stubs.push_back(streamsStub[offsetStub + layer][frame]);
        tracks_.emplace_back(frameTrack, stubs);
        input.push_back(&tracks_.back());
      }
      // remove all gaps between end and last track
      for (auto it = input.end(); it != input.begin();)
        it = (*--it) ? input.begin() : input.erase(it);
    }
  }

  // fill output products
  void KFin::produce(StreamsStub& accpetedStubs,
                     StreamsTrack& acceptedTracks,
                     StreamsStub& lostStubs,
                     StreamsTrack& lostTracks) {
    // merge number of nodes DR to number of Nodes KF and store result
    static const int nMux = channelAssignment_->numNodesDR() / setup_->kfNumWorker();
    const int offsetTrack = region_ * setup_->kfNumWorker();
    for (int nodeKF = 0; nodeKF < setup_->kfNumWorker(); nodeKF++) {
      const int offset = nodeKF * nMux;
      deque<Track*> accepted;
      deque<Track*> lost;
      vector<deque<Track*>> stacks(nMux);
      vector<deque<Track*>> inputs(nMux);
      for (int channel = 0; channel < nMux; channel++) {
        const vector<Track*>& input = input_[offset + channel];
        inputs[channel] = deque<Track*>(input.begin(), input.end());
      }
      // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
      while (!all_of(inputs.begin(), inputs.end(), [](const deque<Track*>& tracks) { return tracks.empty(); }) or
             !all_of(stacks.begin(), stacks.end(), [](const deque<Track*>& tracks) { return tracks.empty(); })) {
        // fill input fifos
        for (int channel = 0; channel < nMux; channel++) {
          deque<Track*>& stack = stacks[channel];
          Track* track = pop_front(inputs[channel]);
          if (track)
            stack.push_back(track);
        }
        // merge input fifos to one stream, prioritizing higher input channel over lower channel
        bool nothingToRoute(true);
        for (int channel = nMux - 1; channel >= 0; channel--) {
          Track* track = pop_front(stacks[channel]);
          if (track) {
            nothingToRoute = false;
            accepted.push_back(track);
            break;
          }
        }
        if (nothingToRoute)
          accepted.push_back(nullptr);
      }
      // truncate if desired
      if (enableTruncation_ && (int)accepted.size() > setup_->numFrames()) {
        const auto limit = next(accepted.begin(), setup_->numFrames());
        copy_if(limit, accepted.end(), back_inserter(lost), [](const Track* track) { return track; });
        accepted.erase(limit, accepted.end());
      }
      // remove all gaps between end and last track
      for (auto it = accepted.end(); it != accepted.begin();)
        it = (*--it) ? accepted.begin() : accepted.erase(it);
      // fill products StreamsStub& accpetedStubs, StreamsTrack& acceptedTracks, StreamsStub& lostStubs, StreamsTrack& lostTracks
      const int channelTrack = offsetTrack + nodeKF;
      const int offsetStub = channelTrack * setup_->numLayers();
      // fill lost tracks and stubs without gaps
      lostTracks[channelTrack].reserve(lost.size());
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        lostStubs[offsetStub + layer].reserve(lost.size());
      for (Track* track : lost) {
        lostTracks[channelTrack].emplace_back(track->frame_);
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          lostStubs[offsetStub + layer].emplace_back(track->stubs_[layer]);
      }
      // fill accepted tracks and stubs with gaps
      acceptedTracks[channelTrack].reserve(accepted.size());
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        accpetedStubs[offsetStub + layer].reserve(accepted.size());
      for (Track* track : accepted) {
        if (!track) {  // fill gap
          acceptedTracks[channelTrack].emplace_back(FrameTrack());
          for (int layer = 0; layer < setup_->numLayers(); layer++)
            accpetedStubs[offsetStub + layer].emplace_back(FrameStub());
          continue;
        }
        acceptedTracks[channelTrack].emplace_back(track->frame_);
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          accpetedStubs[offsetStub + layer].emplace_back(track->stubs_[layer]);
      }
    }
  }

  // remove and return first element of deque, returns nullptr if empty
  template <class T>
  T* KFin::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

}  // namespace trklet
