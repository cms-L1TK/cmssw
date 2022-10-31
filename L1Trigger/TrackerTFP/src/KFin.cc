#include "L1Trigger/TrackerTFP/interface/KFin.h"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <deque>
#include <vector>
#include <cmath>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  KFin::KFin(const ParameterSet& iConfig, const Setup* setup, const DataFormats* dataFormats, const LayerEncoding* layerEncoding, int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        layerEncoding_(layerEncoding),
        region_(region),
        input_(dataFormats_->numChannel(Process::zht)) {}

  // read in and organize input product (fill vector input_)
  void KFin::consume(const StreamsStub& streams) {
    const int offset = region_ * dataFormats_->numChannel(Process::zht);
    auto validFrame = [](int& sum, const FrameStub& frame) { return sum += frame.first.isNonnull() ? 1 : 0; };
    int nStubsZHT(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::zht); channel++) {
      const StreamStub& stream = streams[channel + offset];
      input_.reserve(stream.size());
      nStubsZHT += accumulate(stream.begin(), stream.end(), 0, validFrame);
    }
    stubsZHT_.reserve(nStubsZHT);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::zht); channel++) {
      for (const FrameStub& frame : streams[channel + offset]) {
        // Store input stubs in vector, so rest of KFin algo can work with pointers to them (saves CPU)
        StubZHT* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubsZHT_.emplace_back(frame, dataFormats_, channel);
          stub = &stubsZHT_.back();
        }
        input_[channel].push_back(stub);
      }
    }
    stubsKFin_.reserve(nStubsZHT);
  }

  // read in and remove TTTrackRefs from argument
  void KFin::consume(vector<TTTrackRef>& ttTrackRefs) {
    // count all tracks of this processing region
    int trackId(-1);
    auto trackIdChange = [&trackId](int& sum, const StubZHT& stub) {
      if (stub.trackId() == trackId)
        return sum;
      trackId = stub.trackId();
      return ++sum;
    };
    const int nTracks = accumulate(stubsZHT_.begin(), stubsZHT_.end(), 0, trackIdChange);
    tracksKFin_.reserve(nTracks);
    ttTrackRefs_.reserve(nTracks);
    // move them
    const auto limit = next(ttTrackRefs.end(), -nTracks);
    move(limit, ttTrackRefs.end(), back_inserter(ttTrackRefs_));
    ttTrackRefs.erase(limit, ttTrackRefs.end());
  }

  // fill output products
  void KFin::produce(StreamsTrack& acceptedTracks, StreamsStub& acceptedStubs, StreamsTrack& lostTracks, StreamsStub& lostStubs) {
    const int offsetOut = region_ * dataFormats_->numChannel(Process::kfin);
    auto putS = [](const deque<StubKFin*>& stubs, StreamStub& stream) {
      stream.reserve(stubs.size());
      for (StubKFin* stub : stubs)
        stream.emplace_back(stub ? stub->frame() : FrameStub());
    };
    auto putT = [](const deque<TrackKFin*>& tracks, StreamTrack& stream) {
      stream.reserve(tracks.size());
      for (TrackKFin* track : tracks)
        stream.emplace_back(track ? track->frame() : FrameTrack());
    };
    // loop over worker
    for (int channelOut = 0; channelOut < dataFormats_->numChannel(Process::kfin); channelOut++) {
      const int offsetIn = channelOut * setup_->kfinNumMuxedChannel();
      deque<TrackKFin*> tracks;
      deque<TrackKFin*> tracksLost;
      vector<deque<StubKFin*>> stubs(setup_->numLayers());
      vector<deque<StubKFin*>> stubsLost(setup_->numLayers());
      for (int channelIn = 0; channelIn < setup_->kfinNumMuxedChannel(); channelIn++) {
        const vector<StubZHT*> input = input_[channelIn + offsetIn];
        // identify tracks in input container
        int id;
        auto different = [&id](StubZHT* stub) { return !stub || id != stub->trackId(); };
        for (auto it = input.begin(); it != input.end();) {
          const auto start = find_if(it, input.end(), [](StubZHT* stub){ return stub; });
          id = (*start)->trackId();
          const auto end = find_if(start, input.end(), different);
          vector<StubZHT*> track(start, end);
          track.erase(remove(track.begin(), track.end(), nullptr), track.end());
          if (track.empty())
            break;
          // convert
          const StubZHT& stub = *track.front();
          const TTBV& maybePattern = layerEncoding_->maybePattern(stub.zT());
          const vector<int>& layerEncoding = layerEncoding_->layerEncoding(stub.zT());
          tracksKFin_.emplace_back(pop_back(ttTrackRefs_), dataFormats_, maybePattern, stub.inv2R(), stub.phiT(), stub.zT());
          vector<vector<StubKFin*>> layerStubs(setup_->numLayers());
          for (vector<StubKFin*>& layer : layerStubs)
            layer.reserve(track.size());
          for (StubZHT* stub : track) {
            const int layerId = setup_->layerId(stub->ttStubRef());
            const auto it = find(layerEncoding.begin(), layerEncoding.end(), layerId);
            const int kfLayerId = min((int)distance(layerEncoding.begin(), it), setup_->numLayers() - 1);
            stubsKFin_.emplace_back(*stub, kfLayerId);
            layerStubs[kfLayerId].push_back(&stubsKFin_.back());
          }
          auto maxSize = [](int& max, const vector<StubKFin*>& layer){ return max = std::max(max, (int)layer.size()); };
          int length = accumulate(layerStubs.begin(), layerStubs.end(), 0, maxSize);
          for (vector<StubKFin*>& layer : layerStubs)
            layer.resize(length, nullptr);
          // apply truncation
          if (length > setup_->kfinMaxStubsPerLayer()) {
            tracksLost.push_back(&tracksKFin_.back());
            tracksLost.insert(tracksLost.end(), length - setup_->kfinMaxStubsPerLayer() - 1, nullptr);
            length = setup_->kfinMaxStubsPerLayer();
            for (int layerId = 0; layerId < setup_->numLayers(); layerId++) {
              vector<StubKFin*>& layer = layerStubs[layerId];
              copy(next(layer.begin(), length), layer.end(), back_inserter(stubsLost[layerId]));
              layer.resize(length);
            }
          }
          tracks.push_back(&tracksKFin_.back());
          tracks.insert(tracks.end(), length - 1, nullptr);
          for (int layerId = 0; layerId < setup_->numLayers(); layerId++) {
            vector<StubKFin*>& layer = layerStubs[layerId];
            copy(layer.begin(), layer.end(), back_inserter(stubs[layerId]));
          }
          // set begin of next track
          it = end;
        }
      }
      // apply truncation
      if (enableTruncation_ && (int)tracks.size() > setup_->numFrames()) {
        const auto limitT = next(tracks.begin(), setup_->numFrames());
        copy(limitT, tracks.end(), back_inserter(tracksLost));
        tracks.erase(limitT, tracks.end());
        for (int layerId = 0; layerId < setup_->numLayers(); layerId++) {
          deque<StubKFin*>& layer = stubs[layerId];
          const auto limitS = next(layer.begin(), setup_->numFrames());
          copy(limitS, layer.end(), back_inserter(stubsLost[layerId]));
          layer.erase(limitS, layer.end());
        }
      }
      // put data on ed products
      putT(tracks, acceptedTracks[channelOut + offsetOut]);
      putT(tracksLost, lostTracks[channelOut + offsetOut]);
      const int offsetStub = (channelOut + offsetOut) * setup_->numLayers();
      for (int layerId = 0; layerId < setup_->numLayers(); layerId++) {
        putS(stubs[layerId], acceptedStubs[layerId + offsetStub]);
        putS(stubsLost[layerId], lostStubs[layerId + offsetStub]);
      }
    }
  }

  // remove and return last element of vector, returns nullRef if empty
  TTTrackRef KFin::pop_back(std::vector<TTTrackRef>& ttTrackRefs) const {
    TTTrackRef ttTrackRef;
    if (!ttTrackRefs.empty()) {
      ttTrackRef = ttTrackRefs.back();
      ttTrackRefs.pop_back();
    }
    return ttTrackRef;
  }

}  // namespace trackerTFP