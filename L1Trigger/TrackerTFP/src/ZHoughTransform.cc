#include "L1Trigger/TrackerTFP/interface/ZHoughTransform.h"

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

  ZHoughTransform::ZHoughTransform(const ParameterSet& iConfig,
                                   const Setup* setup,
                                   const DataFormats* dataFormats,
                                   int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        region_(region),
        input_(dataFormats->numChannel(Process::mht))  {}

  // read in and organize input product (fill vector input_)
  void ZHoughTransform::consume(const StreamsStub& streams) {
    auto valid = [](int sum, const FrameStub& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    const int offset = region_ * dataFormats_->numChannel(Process::mht);
    int nStubsMHT(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      nStubsMHT += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubs_.reserve(nStubsMHT);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      vector<StubMHT*>& stubs = input_[channel];
      stubs.reserve(stream.size());
      // Store input stubs in vector, so rest of ZHT algo can work with pointers to them (saves CPU)
      for (const FrameStub& frame : stream) {
        StubMHT* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubs_.emplace_back(frame, dataFormats_);
          stub = &stubs_.back();
        }
        stubs.push_back(stub);
      }
    }
  }

  // fill output products
  void ZHoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) const {
    // loop over worker
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      // input data
      const vector<StubMHT*>& input = input_[channel];
      // output data
      StreamStub& stream = accepted[region_ * dataFormats_->numChannel(Process::mht) + channel];
      stream.reserve(input.size());
      // identify tracks in input container
      int id;
      auto different = [&id](StubMHT* stub) { return !stub || id != stub->trackId(); };
      for (auto it = input.begin(); it != input.end();) {
        if (!*it) {
          stream.emplace_back();
          it++;
          continue;
        }
        id = (*it)->trackId();
        const auto start = it;
        const auto end = find_if(start, input.end(), different);
        const vector<StubMHT*> track(start, end);
        // run single track through r-z hough trasnform ans store result
        produce(track, stream);
        // set begin of next track
        it = end;
      }
    }
  }

  // run single track through r-z hough transform ans store result
  void ZHoughTransform::produce(const vector<StubMHT*>& track, StreamStub& accepted) const {
    static const double baseCot = dataFormats_->base(Variable::cot, Process::zht) * 2.;
    static const double baseZT = dataFormats_->base(Variable::zT, Process::zht) * 2.;
    int maxLayer(-1);
    int maxLayerPS(-1);
    pair<int, int> limitsCot(setup_->zhtNumBinsCot() / 2, -1);
    pair<int, int> limitsZT(setup_->zhtNumBinsZT() / 2, -1);
    for (int um = 0; um < setup_->zhtNumBinsCot(); um++) {
      const int im = um - setup_->zhtNumBinsCot() / 2.;
      const double m = (im + .5) * baseCot;
      for (int uc = 0; uc < setup_->zhtNumBinsZT(); uc++) {
        const int ic = uc - setup_->zhtNumBinsZT() / 2.;
        const double c = (ic + .5) * baseZT;
        TTBV hitPattern(0, setup_->numLayers());
        TTBV hitPatternPS(0, setup_->numLayers());
        for (StubMHT* stub : track) {
          const double r = stub->r() + setup_->chosenRofPhi() - setup_->chosenRofZ();
          const double chi = stub->z() - (r * m + c);
          const double dChi = stub->dZ();
          if (abs(2. * chi) < baseZT + abs(r) * baseCot + dChi) {
            const int layer = stub->layer();
            hitPattern.set(layer);
            if (stub->ps())
              hitPatternPS.set(layer);
          }
        }
        if (hitPattern.count() < setup_->zhtMinLayers() || hitPatternPS.count() < setup_->zhtMinLayersPS())
          continue;
        if (hitPatternPS.count() < maxLayerPS || hitPattern.count() < maxLayer)
          continue;
        if (hitPatternPS.count() == maxLayerPS && hitPattern.count() == maxLayer) {
          limitsCot = make_pair(min(limitsCot.first, im), max(limitsCot.second, im));
          limitsZT = make_pair(min(limitsZT.first, ic), max(limitsZT.second, ic));
        } else {
          maxLayer = hitPattern.count();
          maxLayerPS = hitPatternPS.count();
          limitsCot = make_pair(im, im);
          limitsZT = make_pair(ic, ic);
        }
      }
    }
    const int im = (limitsCot.first + limitsCot.second);
    const int ic = (limitsZT.first + limitsZT.second);
    const double m = (im + .5) * baseCot / 2.;
    const double c = (ic + .5) * baseZT / 2.;
    vector<StubZHT> stubs;
    stubs.reserve(track.size());
    TTBV hitPattern(0, setup_->numLayers());
    TTBV hitPatternPS(0, setup_->numLayers());
    for (StubMHT* stub : track) {
      const double r = stub->r() + setup_->chosenRofPhi() - setup_->chosenRofZ();
      const double z = stub->z() - (r * m + c);
      const double dZ = stub->dZ();
      if (abs(2. * z) < baseZT * 2. + abs(r) * baseCot * 2. + dZ) {
        stubs.emplace_back(*stub, z, im, ic);
        const int layer = stub->layer();
        hitPattern.set(layer);
        if (stub->ps())
          hitPatternPS.set(layer);
      }
    }
    if (hitPattern.count() < setup_->zhtMinLayers() || hitPatternPS.count() < setup_->zhtMinLayersPS())
      stubs.clear();
    for (const StubZHT& stub : stubs)
      accepted.emplace_back(stub.frame());
    accepted.insert(accepted.end(), track.size() - stubs.size(), FrameStub());
  }

}  // namespace trackerTFP
