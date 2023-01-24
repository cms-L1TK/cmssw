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
      input_[channel].reserve(stream.size());
      nStubsZHT += accumulate(stream.begin(), stream.end(), 0, validFrame);
    }
    stubsZHT_.reserve(nStubsZHT);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::zht); channel++) {
      vector<stringstream> sss(8);
      vector<StubZHT*>& stubs = input_[channel];
      for (const FrameStub& frame : streams[channel + offset]) {
        // Store input stubs in vector, so rest of KFin algo can work with pointers to them (saves CPU)
        StubZHT* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubsZHT_.emplace_back(frame, dataFormats_, channel);
          stub = &stubsZHT_.back();
          sss[0] << setw(7) << stub->newTrk() << " ";
          sss[1] << setw(7) << dataFormats_->format(Variable::r, Process::zht).integer(stub->r()) << " ";
          sss[2] << setw(7) << dataFormats_->format(Variable::phi, Process::zht).integer(stub->phi()) << " ";
          sss[3] << setw(7) << dataFormats_->format(Variable::z, Process::zht).integer(stub->z()) << " ";
          sss[4] << setw(7) << stub->layer() << " ";
          sss[5] << setw(7) << stub->module() << " ";
          sss[6] << setw(7) << stub->phiT() << " ";
          sss[7] << setw(7) << stub->zT() << " ";
        } else
          for (stringstream& ss : sss)
            ss << setw(7) << "x" << " ";
        stubs.push_back(stub);
      }
      if (region_ == 7 && channel == 11) {
        cout << "KFin consume" << endl;
        for (const stringstream& ss: sss)
          cout << ss.str() << endl;
        //throw cms::Exception("...");
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
    // prepare input data
    vector<deque<TrackKFin*>> inputTracks(dataFormats_->numChannel(Process::zht));
    vector<deque<StubKFin*>> inputStubs(dataFormats_->numChannel(Process::zht));
    prepare(inputTracks, inputStubs);
    // loop over worker
    for (int channelOut = 0; channelOut < dataFormats_->numChannel(Process::kfin); channelOut++) {
      const int indexTracks = region_ * dataFormats_->numChannel(Process::kfin) + channelOut;
      const int offsetStubs = indexTracks * setup_->numLayers();
      deque<TrackKFin*> tracks;
      vector<deque<StubKFin*>> stubs(setup_->numLayers());
      vector<deque<TrackKFin*>> streamsTracks;
      vector<deque<StubKFin*>> streamsStubs;
      streamsTracks.reserve(setup_->kfinNumMuxedChannel());
      streamsStubs.reserve(setup_->kfinNumMuxedChannel());
      for (int channelIn = 0; channelIn < setup_->kfinNumMuxedChannel(); channelIn++) {
        const int index = channelIn * dataFormats_->numChannel(Process::kfin) + channelOut;
        streamsTracks.push_back(inputTracks[index]);
        streamsStubs.push_back(inputStubs[index]);
      }
      vector<stringstream> sssT(7);
      // route
      route(streamsTracks, tracks);
      /*for (TrackKFin* track : tracks) {
        if (!track) {
          for (stringstream& ss : sssT)
            ss << "x ";
          continue;
        }
        sssT[0] << setw(7) << track->trackId() / setup_->kfinMaxTracks() << " ";
        sssT[1] << setw(7) << track->trackId() % setup_->kfinMaxTracks() << " ";
        sssT[2] << setw(7) << track->hitPattern() << " ";
        sssT[3] << setw(7) << track->lmap().val() << " ";
        sssT[4] << setw(7) << track->maybePattern() << " ";
        sssT[5] << setw(7) << dataFormats_->format(Variable::phiT, Process::kfin).integer(track->phiT()) << " ";
        sssT[6] << setw(7) << dataFormats_->format(Variable::zT, Process::kfin).integer(track->zT()) << " ";
      }*/
      /*if (region_ == 1) {
        for (const stringstream& ss : sssT)
          cout << ss.str() << endl;
        throw cms::Exception("...");
      }*/
      route(streamsStubs, stubs);
      // aplly truncation
      if (enableTruncation_) {
        if ((int)tracks.size() > setup_->numFrames())
          tracks.resize(setup_->numFrames());
        for (deque<StubKFin*>& stream : stubs)
          if ((int)stream.size() > setup_->numFrames())
            stream.resize(setup_->numFrames());
      }
      vector<stringstream> sssS(8);
      for (StubKFin* stub : stubs[1]) {
        if (!stub) {
          for (stringstream& ss : sssS)
            ss << "x ";
          continue;
        }
        sssS[0] << setw(7) << stub->trackId() / setup_->kfinMaxTracks() << " ";
        sssS[1] << setw(7) << stub->trackId() % setup_->kfinMaxTracks() << " ";
        sssS[2] << setw(7) << stub->stubId() << " ";
        sssS[3] << setw(7) << dataFormats_->format(Variable::r, Process::kfin).integer(stub->r()) << " ";
        sssS[4] << setw(7) << dataFormats_->format(Variable::phi, Process::kfin).integer(stub->phi()) << " ";
        sssS[5] << setw(7) << dataFormats_->format(Variable::z, Process::kfin).integer(stub->z()) << " ";
        sssS[6] << setw(7) << dataFormats_->format(Variable::phi, Process::kfin).integer(stub->dPhi()) << " ";
        sssS[7] << setw(7) << dataFormats_->format(Variable::z, Process::kfin).integer(stub->dZ()) << " ";
      }
      /*if (region_ == 8 && channelOut == 0) {
        for (const stringstream& ss : sssS)
          cout << ss.str() << endl;
        throw cms::Exception("...");
      }*/
      // cycle event, remove all gaps
      tracks.erase(remove(tracks.begin(), tracks.end(), nullptr), tracks.end());
      for (deque<StubKFin*>& stream : stubs)
        stream.erase(remove(stream.begin(), stream.end(), nullptr), stream.end());
      // sort according to track id
      sort(tracks.begin(), tracks.end(), [](TrackKFin* lhs, TrackKFin* rhs){ return lhs->trackId() < rhs->trackId(); });
      for (deque<StubKFin*>& stream : stubs) {
        // sort according to stub id on layer
        sort(stream.begin(), stream.end(), [](StubKFin* lhs, StubKFin* rhs){ return lhs->stubId() < rhs->stubId(); });
        // sort according to track id
        stable_sort(stream.begin(), stream.end(), [](StubKFin* lhs, StubKFin* rhs){ return lhs->trackId() < rhs->trackId(); });
      }
      // add all gaps
      const int size = accumulate(tracks.begin(), tracks.end(), 0, [](int& sum, TrackKFin* track){ return sum += track->size(); });
      for (int frame = 0; frame < size;) {
        const int trackId = tracks[frame]->trackId();
        const int length = tracks[frame]->size();
        tracks.insert(next(tracks.begin(), frame + 1), length - 1, nullptr);
        for (int layer = 0; layer < setup_->numLayers(); layer++) {
          deque<StubKFin*>& stream = stubs[layer];
          if (frame > (int)stream.size()) {
            stream.insert(stream.end(), length, nullptr);
            continue;
          }
          const auto begin = next(stream.begin(), frame);
          const auto end = find_if(begin, stream.end(), [trackId](StubKFin* stub){ return stub->trackId() != trackId; });
          stream.insert(end, length - distance(begin, end), nullptr);
        }
        frame += length;
      }
      // remove all invalidated tracks due to truncation
      vector<stringstream> sssT0(4);
      for (TrackKFin* track : tracks) {
        if (!track) {
          for (stringstream& ss : sssT0)
            ss << setw(7) << "x" << " ";
          continue;
        }
        sssT0[0] << setw(7) << dataFormats_->format(Variable::inv2R, Process::kfin).integer(track->inv2R()) << " ";
        sssT0[1] << setw(7) << dataFormats_->format(Variable::phiT, Process::kfin).integer(track->phiT()) << " ";
        sssT0[2] << setw(7) << dataFormats_->format(Variable::zT, Process::kfin).integer(track->zT()) << " ";
        sssT0[3] << setw(7) << track->maybePattern().val() << " ";
      }
      if (region_ == 7 && channelOut == 1) {
        cout << "T0" << endl;
        for (const stringstream& ss : sssT0)
          cout << ss.str() << endl;
        //throw cms::Exception("...");
      }
      for (int frame = 0; frame < size;) {
        const int trackId = tracks[frame]->trackId();
        const int size = tracks[frame]->size();
        TTBV hitPattern(0, setup_->numLayers());
        for (const deque<StubKFin*>& stream : stubs) {
          for (StubKFin* stub : stream) {
            if (!stub || stub->trackId() != trackId)
              continue;
            hitPattern.set(stub->layer());
            break;
          }
        }
        if (hitPattern.count() < setup_->kfMinLayers()) {
          tracks[frame] = nullptr;
          for (int layer = 0; layer < setup_->numLayers(); layer++)
            for (int sFrame = frame; sFrame < frame + size; sFrame++)
              stubs[layer][sFrame] = nullptr;
        }
        frame += size;
      }
      // truncate if desired
      if (enableTruncation_ && size > setup_->numFrames()) {
        tracks.resize(setup_->numFrames());
        for (deque<StubKFin*> stream : stubs)
          stream.resize(setup_->numFrames());
      }
      // store
      vector<stringstream> sssT1(4);
      for (TrackKFin* track : tracks) {
        if (!track) {
          for (stringstream& ss : sssT1)
            ss << setw(7) << "x" << " ";
          continue;
        }
        sssT1[0] << setw(7) << dataFormats_->format(Variable::inv2R, Process::kfin).integer(track->inv2R()) << " ";
        sssT1[1] << setw(7) << dataFormats_->format(Variable::phiT, Process::kfin).integer(track->phiT()) << " ";
        sssT1[2] << setw(7) << dataFormats_->format(Variable::zT, Process::kfin).integer(track->zT()) << " ";
        sssT1[3] << setw(7) << track->maybePattern().val() << " ";
      }
      if (region_ == 7 && channelOut == 1) {
        cout << "T1" << endl;
        for (const stringstream& ss : sssT1)
          cout << ss.str() << endl;
        throw cms::Exception("...");
      }
      acceptedTracks[indexTracks].reserve(size);
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        acceptedStubs[offsetStubs + layer].reserve(size);
      transform(tracks.begin(), tracks.end(), back_inserter(acceptedTracks[indexTracks]), [](TrackKFin* track){ return (track ? track->frame() : FrameTrack()); });
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        transform(stubs[layer].begin(), stubs[layer].end(), back_inserter(acceptedStubs[offsetStubs + layer]), [](StubKFin* stub){ return (stub ? stub->frame() : FrameStub()); });
    }
  }

  //
  void KFin::route(vector<deque<StubKFin*>>& input, vector<deque<StubKFin*>>& outputs) const {
    for (int channelOut = 0; channelOut < (int)outputs.size(); channelOut++) {
      deque<StubKFin*>& output = outputs[channelOut];
      vector<deque<StubKFin*>> inputs(input);
      for (deque<StubKFin*>& stream : inputs) {
        for (StubKFin*& stub : stream)
          if (stub && stub->layer() != channelOut)
            stub = nullptr;
        for (auto it = stream.end(); it != stream.begin();)
          it = (*--it) ? stream.begin() : stream.erase(it);
      }
      vector<deque<StubKFin*>> stacks(input.size());
      // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
      while (!all_of(inputs.begin(), inputs.end(), [](const deque<StubKFin*>& stubs) { return stubs.empty(); }) or
              !all_of(stacks.begin(), stacks.end(), [](const deque<StubKFin*>& stubs) { return stubs.empty(); })) {
        // fill input fifos
        for (int channel = 0; channel < (int)input.size(); channel++) {
          deque<StubKFin*>& stack = stacks[channel];
          StubKFin* stub = pop_front(inputs[channel]);
          if (stub) {
            if (enableTruncation_ && (int)stack.size() == setup_->kfinDepthMemory() - 1)
              pop_front(stack);
            stack.push_back(stub);
          }
        }
        // merge input fifos to one stream, prioritizing higher input channel over lower channel
        bool nothingToRoute(true);
        for (int channel = (int)input.size() - 1; channel >= 0; channel--) {
          StubKFin* stub = pop_front(stacks[channel]);
          if (stub) {
            nothingToRoute = false;
            output.push_back(stub);
            break;
          }
        }
        if (nothingToRoute)
          output.push_back(nullptr);
      }
    }
  }

  //
  void KFin::route(vector<deque<TrackKFin*>>& inputs, deque<TrackKFin*>& output) const {
    vector<deque<TrackKFin*>> stacks(inputs.size());
    // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
    while (!all_of(inputs.begin(), inputs.end(), [](const deque<TrackKFin*>& tracks) { return tracks.empty(); }) or
            !all_of(stacks.begin(), stacks.end(), [](const deque<TrackKFin*>& tracks) { return tracks.empty(); })) {
      // fill input fifos
      for (int channel = 0; channel < (int)inputs.size(); channel++) {
        deque<TrackKFin*>& stack = stacks[channel];
        TrackKFin* track = pop_front(inputs[channel]);
        if (track) {
          if (enableTruncation_ && (int)stack.size() == setup_->kfinDepthMemory() - 1)
            pop_front(stack);
          stack.push_back(track);
        }
      }
      // merge input fifos to one stream, prioritizing higher input channel over lower channel
      bool nothingToRoute(true);
      for (int channel = (int)inputs.size() - 1; channel >= 0; channel--) {
        TrackKFin* track = pop_front(stacks[channel]);
        if (track) {
          nothingToRoute = false;
          output.push_back(track);
          break;
        }
      }
      if (nothingToRoute)
        output.push_back(nullptr);
    }
  }

  // peprare input data, associate stub collections with TTTrackRefs
  void KFin::prepare(vector<deque<TrackKFin*>>& streamsTracks, vector<deque<StubKFin*>>& streamsStubs) {
    for (int channel = 0; channel < dataFormats_->numChannel(Process::zht); channel++) {
      const int offset = channel / dataFormats_->numChannel(Process::kfin) * setup_->kfinMaxTracks();
      const vector<StubZHT*>& input = input_[channel];
      deque<TrackKFin*>& tracks = streamsTracks[channel];
      deque<StubKFin*>& stubs = streamsStubs[channel];
      // identify tracks in input container
      auto newTrk = [](StubZHT* stub) { return stub && stub->newTrk(); };
      const auto begin = find_if(input.begin(), input.end(), [](StubZHT* stub){ return stub; });
      const int startGap = distance(input.begin(), begin);
      tracks.insert(tracks.end(), startGap, nullptr);
      stubs.insert(stubs.end(), startGap, nullptr);
      int sumSize = startGap - setup_->zhtMinLayers();
      int trackId(0);
      for (auto it = begin; it != input.end();) {
        auto end = find_if(it, input.end(), newTrk);
        if (end != input.end())
          end = next(end, 1);
        vector<StubZHT*> track(it, end);
        if (all_of(track.begin(), track.end(), [](StubZHT* s){ return !s; }))
          break;
        // set begin of next track
        it = end;
        // restore clock accurancy
        sumSize += (int)track.size();
        const int deltaStubs = sumSize - (int)stubs.size();
        if (deltaStubs > 0)
          stubs.insert(stubs.end(), deltaStubs, nullptr);
        const int deltaTracks = sumSize - (int)tracks.size();
        if (deltaTracks > 0)
          tracks.insert(tracks.end(), deltaTracks, nullptr);
        track.erase(remove(track.begin(), track.end(), nullptr), track.end());
        // convert
        const StubZHT& stubZHT = *track.front();
        TTBV hitPattern(0, setup_->numLayers());
        vector<int> layerCounts(setup_->numLayers(), 0);
        const auto start = stubsKFin_.end();
        const vector<int>& le = layerEncoding_->layerEncoding(stubZHT.zT());
        for (StubZHT* stub : track) {
          const int layerId = setup_->layerId(stub->ttStubRef());
          const auto it = find(le.begin(), le.end(), layerId);
          const int kfLayerId = min((int)distance(le.begin(), it), setup_->numLayers() - 1);
          hitPattern.set(kfLayerId);
          if (layerCounts[kfLayerId] < setup_->kfinMaxStubsPerLayer())
            stubsKFin_.emplace_back(*stub, kfLayerId, offset + trackId, layerCounts[kfLayerId]++);
        }
        const TTBV& maybePattern = layerEncoding_->maybePattern(stubZHT.zT());
        tracksKFin_.emplace_back(pop_back(ttTrackRefs_), dataFormats_, maybePattern, hitPattern, layerCounts, stubZHT.inv2R(), stubZHT.phiT(), stubZHT.zT(), offset + trackId);
        // pattern reco
        if (hitPattern.count() < setup_->kfMinLayers())
          continue;
        // fill streams
        if (trackId++ >= setup_->kfinMaxTracks())
          continue;
        tracks.push_back(&tracksKFin_.back());
        for (auto sit = start; sit != stubsKFin_.end(); sit++)
          stubs.push_back(&*sit);
      }
      vector<stringstream> sss(7);
      for (StubKFin* stub : stubs) {
          if (!stub) {
            for (stringstream& ss : sss)
              ss << setw(7) << "x" << " ";
            continue;
          }
          sss[0] << setw(7) << stub->trackId() << " ";
          sss[1] << setw(7) << stub->layer() << " ";
          sss[2] << setw(7) << dataFormats_->format(Variable::r, Process::zht).integer(stub->r()) << " ";
          sss[3] << setw(7) << dataFormats_->format(Variable::phi, Process::zht).integer(stub->phi()) << " ";
          sss[4] << setw(7) << dataFormats_->format(Variable::z, Process::zht).integer(stub->z()) << " ";
          sss[5] << setw(7) << dataFormats_->format(Variable::phi, Process::zht).integer(stub->dPhi()) << " ";
          sss[6] << setw(7) << dataFormats_->format(Variable::z, Process::zht).integer(stub->dZ()) << " ";
      }
      /*if (region_ == 0 && channel == 2*1+0) {
        for (const stringstream& ss : sss)
          cout << ss.str() << endl;
        throw cms::Exception("....");
      }*/
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

}  // namespace trackerTFP