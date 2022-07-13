#include "L1Trigger/TrackerTFP/interface/MiniHoughTransform.h"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <deque>
#include <vector>
#include <utility>
#include <cmath>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  MiniHoughTransform::MiniHoughTransform(const ParameterSet& iConfig,
                                         const Setup* setup,
                                         const DataFormats* dataFormats,
                                         int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        region_(region),
        input_(dataFormats_->numChannel(Process::ht)) {}

  // read in and organize input product (fill vector input_)
  void MiniHoughTransform::consume(const StreamsStub& streams) {
    auto valid = [](int sum, const FrameStub& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    int nStubsHT(0);
    for (int binInv2R = 0; binInv2R < dataFormats_->numChannel(Process::ht); binInv2R++) {
      const StreamStub& stream = streams[region_ * dataFormats_->numChannel(Process::ht) + binInv2R];
      nStubsHT += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubs_.reserve(nStubsHT);
    for (int binInv2R = 0; binInv2R < dataFormats_->numChannel(Process::ht); binInv2R++) {
      const int inv2R = dataFormats_->format(Variable::inv2R, Process::ht).toSigned(binInv2R);
      const StreamStub& stream = streams[region_ * dataFormats_->numChannel(Process::ht) + binInv2R];
      vector<StubHT*>& stubs = input_[binInv2R];
      stubs.reserve(stream.size());
      // Store input stubs in vector, so rest of MHT algo can work with pointers to them (saves CPU)
      for (const FrameStub& frame : stream) {
        StubHT* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubs_.emplace_back(frame, dataFormats_, inv2R);
          stub = &stubs_.back();
        }
        stubs.push_back(stub);
      }
    }
  }

  // fill output products
  void MiniHoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) const {
    // loop over worker
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      // input data
      const vector<StubHT*>& input = input_[channel];
      // output data
      StreamStub& stream = accepted[region_ * dataFormats_->numChannel(Process::mht) + channel];
      stream.reserve(input.size());
      // identify tracks in input container
      int id;
      auto different = [&id](StubHT* stub) { return id != stub->trackId(); };
      for (auto it = input.begin(); it != input.end();) {
        id = (*it)->trackId();
        const auto start = it;
        const auto end = find_if(start, input.end(), different);
        const vector<StubHT*> track(start, end);
        // run single track through r-z hough trasnform ans store result
        produce(track, stream);
        // set begin of next track
        it = end;
      }
    }
  }

  // run single track through r-z hough trasnform ans store result
  void MiniHoughTransform::produce(const vector<StubHT*>& track, StreamStub& accepted) const {
    static const double baseInv2R = dataFormats_->base(Variable::inv2R, Process::mht) * 2.;
    static const double basePhiT = dataFormats_->base(Variable::phiT, Process::mht) * 2.;
    int maxLayer(-1);
    pair<int, int> limitsInv2R(setup_->mhtNumBinsInv2R() / 2, -1);
    pair<int, int> limitsPhiT(setup_->mhtNumBinsPhiT() / 2, -1);
    for (int um = 0; um < setup_->mhtNumBinsInv2R(); um++) {
      const int im = um - setup_->mhtNumBinsInv2R() / 2.;
      const double m = (im + .5) * baseInv2R;
      for (int uc = 0; uc < setup_->mhtNumBinsPhiT(); uc++) {
        const int ic = uc - setup_->mhtNumBinsPhiT() / 2.;
        const double c = (ic + .5) * basePhiT;
        TTBV hitPattern(0, setup_->numLayers());
        TTBV hitPatternPS(0, setup_->numLayers());
        for (StubHT* stub : track) {
          const double r = stub->r();
          const double chi = stub->phi() - (r * m + c);
          const double dChi = stub->dPhi();
          if (abs(2. * chi) < basePhiT + abs(r) * baseInv2R + dChi) {
            const int layer = stub->layer();
            hitPattern.set(layer);
            if (stub->ps())
              hitPatternPS.set(layer);
          }
        }
        if (hitPattern.count() < setup_->mhtMinLayers() || hitPatternPS.count() < setup_->mhtMinLayersPS())
          continue;
        if (hitPattern.count() < maxLayer)
          continue;
        if (hitPattern.count() == maxLayer) {
          limitsInv2R = make_pair(min(limitsInv2R.first, im), max(limitsInv2R.second, im));
          limitsPhiT = make_pair(min(limitsPhiT.first, ic), max(limitsPhiT.second, ic));
        } else {
          maxLayer = hitPattern.count();
          limitsInv2R = make_pair(im, im);
          limitsPhiT = make_pair(ic, ic);
        }
      }
    }
    const double im = (limitsInv2R.first + limitsInv2R.second);
    const double ic = (limitsPhiT.first + limitsPhiT.second);
    const double m = (im + .5) * baseInv2R / 2.;
    const double c = (ic + .5) * basePhiT / 2.;
    vector<StubMHT> stubs;
    stubs.reserve(track.size());
    TTBV hitPattern(0, setup_->numLayers());
    TTBV hitPatternPS(0, setup_->numLayers());
    for (StubHT* stub : track) {
      const double r = stub->r();
      const double phi = stub->phi() - (r * m + c);
      const double dPhi = stub->dPhi();
      if (abs(2. * phi) < basePhiT * 2. + abs(r) * baseInv2R * 2. + dPhi) {
        stubs.emplace_back(*stub, phi, im, ic);
        const int layer = stub->layer();
        hitPattern.set(layer);
        if (stub->ps())
          hitPatternPS.set(layer);
      }
    };
    if (hitPattern.count() < setup_->mhtMinLayers() || hitPatternPS.count() < setup_->mhtMinLayersPS())
      stubs.clear();
    for (const StubMHT& stub : stubs)
      accepted.emplace_back(stub.frame());
    accepted.insert(accepted.end(), track.size() - stubs.size(), FrameStub());
  }

}  // namespace trackerTFP
