#include "L1Trigger/TrackerTFP/interface/SeedFilter.h"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <deque>
#include <vector>
#include <set>
#include <utility>
#include <cmath>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  SeedFilter::SeedFilter(const ParameterSet& iConfig, const Setup* setup, const DataFormats* dataFormats, const LayerEncoding* layerEncoding, int region) :
    enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
    setup_(setup),
    dataFormats_(dataFormats),
    layerEncoding_(layerEncoding),
    region_(region),
    cot_(dataFormats_->format(Variable::cot, Process::sf)),
    zT_(dataFormats_->format(Variable::zT, Process::sf)),
    z_(dataFormats_->format(Variable::z, Process::sf)),
    input_(dataFormats_->numChannel(Process::mht)) {}

  // read in and organize input product
  void SeedFilter::consume(const StreamsStub& streams) {
    auto valid = [](int& sum, const FrameStub& frame){ return sum += (frame.first.isNonnull() ? 1 : 0); };
    int nStubs(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[region_ * dataFormats_->numChannel(Process::mht) + channel];
      nStubs += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubsMHT_.reserve(nStubs);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[region_ * dataFormats_->numChannel(Process::mht) + channel];
      vector<StubMHT*>& stubs = input_[channel];
      stubs.reserve(stream.size());
      for (const FrameStub& frame : stream) {
        StubMHT* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubsMHT_.emplace_back(frame, dataFormats_);
          stub = &stubsMHT_.back();
        }
        stubs.push_back(stub);
      }
    }
  }

  // fill output products
  void SeedFilter::produce(StreamsStub& accepted, StreamsStub& lost) {
    auto kill = [this](const vector<StubSF*>& stubs) {
      TTBV hitPattern(0, setup_->numLayers());
      set<int> psLayers;
      TTBV psPattern(0, setup_->numLayers());
      for (StubSF* stub : stubs) {
        hitPattern.set(stub->layer());
        if (setup_->psModule(stub->ttStubRef()))
          psLayers.insert(stub->layer());
      }
      // at least 4 layers
      if (hitPattern.count() < setup_->kfMinLayers())
        return true;
      // at least 2 ps stubs (two stubs in first 3 layers?)
      if ((int)psLayers.size() < 2)
        return true;
      return false;
    };
    auto lessSkippedLayers = [this](const vector<StubSF*>& lhs, const vector<StubSF*>& rhs) {
      auto numSkippedLayers = [this](const vector<StubSF*>& stubs) {
        TTBV hitPattern(0, setup_->numLayers());
        for (StubSF* stub : stubs)
          hitPattern.set(stub->layer());
        return hitPattern.count(0, hitPattern.pmEncode(), false);
      };
      return numSkippedLayers(lhs) < numSkippedLayers(rhs);
    };
    auto smallerChi2 = [this](const vector<StubSF*>& lhs, const vector<StubSF*>& rhs) {
      auto chi2 = [this](const vector<StubSF*>& stubs) {
        double d(0.);
        for (StubSF* stub : stubs)
          d += stub->chi2();
        return d / ((int)stubs.size() - 2);
      };
      return chi2(lhs) < chi2(rhs);
    };
    auto moreLayers = [this](const vector<StubSF*>& lhs, const vector<StubSF*>& rhs) {
      auto numLayers = [this](const vector<StubSF*>& stubs) {
        TTBV hitPattern(0, setup_->numLayers());
        for (StubSF* stub : stubs)
          hitPattern.set(stub->layer());
        return hitPattern.count();
      };
      return numLayers(lhs) > numLayers(rhs);
    };
    const int numLayers = 16;
    static const vector<pair<int, int>> seedingLayers = {{1, 2}, {1, 3}, {2, 3}, {1, 11}, {2, 11}, {1, 12}, {2, 12}, {11, 12}};
    const double dT = setup_->chosenRofPhi() - setup_->chosenRofZ();
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      StreamStub& acceptedStream = accepted[region_ * dataFormats_->numChannel(Process::mht) + channel];
      StreamStub& lostStream = lost[region_ * dataFormats_->numChannel(Process::mht) + channel];
      deque<StubSF*> output;
      vector<StubMHT*>& stubs = input_[channel];
      stubs.erase(remove(stubs.begin(), stubs.end(), nullptr), stubs.end());
      set<int> ids;
      for (StubMHT* stub : stubs)
        if (stub)
          ids.insert(stub->trackId());
      vector<vector<StubMHT*>> tracks(ids.size());
      int i(0);
      for (auto it = stubs.begin(); it != stubs.end();) {
        const auto start = it;
        const int id = (*it)->trackId();
        auto different = [id](StubMHT* stub){ return id != stub->trackId(); };
        it = find_if(it, stubs.end(), different);
        vector<StubMHT*>& track = tracks[i++];
        track.reserve(distance(start, it));
        copy(start, it, back_inserter(track));
      }
      for (const vector<StubMHT*>& track : tracks) {
        const int eta = track.front()->sectorEta();
        const double dZT = (sinh(setup_->boundarieEta(eta + 1)) - sinh(setup_->boundarieEta(eta))) * setup_->chosenRofZ();
        vector<vector<StubMHT*>> layerStubs(numLayers);
        for (vector<StubMHT*>& layer : layerStubs)
          layer.reserve(track.size());
        for (StubMHT* stub : track)
          layerStubs[layerId(stub)].push_back(stub);
        deque<pair<StubMHT*, StubMHT*>> seeds;
        for (const pair<int, int>& seedingLayer : seedingLayers)
          for (StubMHT* inner : layerStubs[seedingLayer.first])
            if (setup_->psModule(inner->ttStubRef()))
              for (StubMHT* outer : layerStubs[seedingLayer.second])
                if (setup_->psModule(outer->ttStubRef()))
                  if (outer->r() > inner->r())
                    seeds.emplace_back(inner, outer);
        vector<vector<StubSF*>> seedTracks;
        seedTracks.reserve(seeds.size());
        int pos(0);
        for (const pair<StubMHT*, StubMHT*>& seed : seeds) {
          const double r1 = seed.first->r() + dT;
          const double r2 = seed.second->r() + dT;
          const double z1 = seed.first->z();
          const double z2 = seed.second->z();
          const double r = (r2 + r1) / 2.;
          const double z = (z2 + z1) / 2.;
          const double dr = (r2 - r1);
          const double dz = (z2 - z1);
          const double cot = cot_.digi(dz / dr);
          const double cotGlobal = cot + setup_->sectorCot(eta);
          const double zT = zT_.digi(z - cot * r);
          const double z0 = zT - cot * setup_->chosenRofZ();
          if (abs(z0) > setup_->beamWindowZ() + cot_.base() * setup_->chosenRofZ())
            continue;
          if (abs(zT) > dZT / 2. + zT_.base() / 2.)
            continue;
          //if (!cot_.inRange(cot))
            //continue;
          const vector<int>& layerEncoding = layerEncoding_->layerEncoding(eta, zT_.toUnsigned(zT_.integer(zT)), cot_.toUnsigned(cot_.integer(cot)));
          vector<StubSF*> stubsSF;
          stubsSF.reserve(track.size());
          for (StubMHT* stub : track) {
            const double sr = stub->r() + dT;
            const double sz = stub->z();
            const double chi = sz - (zT + sr * cot);
            const double dZ = zT_.base() + cot_.base() * abs(sr) + setup_->dZ(stub->ttStubRef(), cotGlobal);
            if (abs(chi) < dZ / 2.) {
              int layer = distance(layerEncoding.begin(), find(layerEncoding.begin(), layerEncoding.end(), layerId(stub)));
              if (layer >= setup_->numLayers())
                layer = setup_->numLayers() - 1;
              const double v = setup_->v1(stub->ttStubRef(), cotGlobal);
              const double cov = pow(cot_.base() * stub->r(), 2) + pow(zT_.base(), 2);
              const double chi2 = pow(chi, 2) / (v + cov);
              stubsSF_.emplace_back(*stub, layer, zT_.integer(zT), cot_.integer(cot), chi2);
              stubsSF.push_back(&stubsSF_.back());
            }
            pos++;
          }
          if (pos >= setup_->numFrames())
            break;
          set<int> ids;
          for (StubSF* stub : stubsSF)
            ids.insert(stub->layer());
          if ((int)ids.size() >= setup_->sfMinLayers())
            seedTracks.push_back(stubsSF);
        }
        if (seedTracks.empty())
          continue;
        seedTracks.erase(remove_if(seedTracks.begin(), seedTracks.end(), kill), seedTracks.end());
        stable_sort(seedTracks.begin(), seedTracks.end(), lessSkippedLayers);
        stable_sort(seedTracks.begin(), seedTracks.end(), smallerChi2);
        stable_sort(seedTracks.begin(), seedTracks.end(), moreLayers);
        const vector<StubSF*>& stubs = seedTracks.front();
        copy(stubs.begin(), stubs.end(), back_inserter(output));
      }
      if (enableTruncation_ && ((int)output.size() > setup_->numFrames())) {
        const auto limit = next(output.begin(), setup_->numFrames());
        lostStream.reserve(distance(limit, output.end()));
        transform(limit, output.end(), back_inserter(lostStream), [](StubSF* stub){ return stub->frame(); });
        output.erase(limit, output.end());
      }
      acceptedStream.reserve(output.size());
      transform(output.begin(), output.end(), back_inserter(acceptedStream), [](StubSF* stub){ return stub->frame(); });
    }
  }

  // stub layerId (barrel: 1-6, endcap: 11-15)
  int SeedFilter::layerId(StubMHT* stub) const {
    int layer = stub->layer();
    bool barrel = setup_->barrel(stub->ttStubRef());
    int lay = -1;
    if (layer == 0)
      lay = 1;
    else if (layer == 1)
      lay = 2;
    else if (layer == 2)
      lay = barrel ? 6 : 11;
    else if (layer == 3)
      lay = barrel ? 5 : 12;
    else if (layer == 4)
      lay = barrel ? 4 : 13;
    else if (layer == 5)
      lay = 14;
    else if (layer == 6)
      lay = barrel ? 3 : 15;
    return lay;
  }

  // remove and return first element of deque, returns nullptr if empty
  template<class T>
  T* SeedFilter::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

} // namespace trackerTFP