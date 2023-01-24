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
        input_(dataFormats_->numChannel(Process::ht)),
        basePhi_(dataFormats_->base(Variable::phi, Process::mht)),
        baseInv2R_(dataFormats_->base(Variable::inv2R, Process::mht)),
        basePhiT_(dataFormats_->base(Variable::phiT, Process::mht)) {
    baseInv2R_ = dataFormats_->base(Variable::inv2R, Process::ht) / setup_->mhtNumBinsInv2R();
    basePhiT_ = dataFormats_->base(Variable::phiT, Process::ht) / setup_->mhtNumBinsPhiT();
  }

  // read in and organize input product (fill vector input_)
  void MiniHoughTransform::consume(const StreamsStub& streams) {
    auto valid = [](int& sum, const FrameStub& frame) { return sum += (frame.first.isNonnull() ? 1 : 0); };
    const int offset = region_ * dataFormats_->numChannel(Process::ht);
    int nStubsHT(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::ht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      nStubsHT += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubsHT_.reserve(nStubsHT);
    stubsMHT_.reserve(nStubsHT);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::ht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      vector<StubHT*>& stubs = input_[channel];
      stubs.reserve(stream.size());
      // Store input stubs in vector, so rest of MHT algo can work with pointers to them (saves CPU)
      int n(0);
      for (const FrameStub& frame : stream) {
        StubHT* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubsHT_.emplace_back(frame, dataFormats_, channel);
          stub = &stubsHT_.back();
        }
        if (n++ < setup_->numFramesIO())
          stubs.push_back(stub);
      }
    }
  }

  // fill output products
  void MiniHoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) {
    const int offset = region_ * dataFormats_->numChannel(Process::mht);
    // loop over worker
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      // input data
      const vector<StubHT*>& input = input_[channel];
      // output data
      StreamStub& stream = accepted[offset + channel];
      StreamStub& streamLost = lost[offset + channel];
      // identify tracks in input container
      int id;
      auto different = [&id](StubHT* stub) { return id != stub->trackId(); };
      auto last = [](StubMHT* stub) { return stub && stub->newTrk(); };
      auto invalid = [](StubMHT* stub){ return !stub || !stub->valid(); };
      deque<StubMHT*> stubsArray;
      int sumSize = -setup_->mhtMinLayers();
      for (auto it = input.begin(); it != input.end();) {
        id = (*it)->trackId();
        const auto start = it;
        const auto end = find_if(start, input.end(), different);
        const vector<StubHT*> track(start, end);
        // restore clock accurancy
        sumSize += (int)track.size();
        const int delta = sumSize - (int)stubsArray.size();
        if (delta > 0)
          stubsArray.insert(stubsArray.end(), delta, nullptr);
        // run single track through r-phi hough trasnform and store result
        array(track, stubsArray);
        // set begin of next track
        it = end;
      }
      deque<StubMHT*> stubsCell;
      sumSize = -setup_->mhtMinLayers();
      for (auto it = stubsArray.begin(); it != stubsArray.end();) {
        const auto start = it;
        auto end = find_if(start, stubsArray.end(), last);
        if (end != stubsArray.end())
          end++;
        vector<StubMHT*> track(start, end);
        // restore clock accurancy
        sumSize += (int)track.size();
        const int delta = sumSize - (int)stubsCell.size();
        if (delta > 0)
          stubsCell.insert(stubsCell.end(), delta, nullptr);
        track.erase(remove_if(track.begin(), track.end(), invalid), track.end());
        if (!track.empty())
          track.back()->setNewTrk();
        // run single track through final r-phi cell and store result
        cell(track, stubsCell);
        // set begin of next track
        it = end;
      }
      // aplly truncation
      if (enableTruncation_ && (int)stubsCell.size() >= setup_->numFrames()) {
        const auto limit = next(stubsCell.begin(), setup_->numFrames());
        auto valid = [](int& sum, StubMHT* stub) { return sum += stub ? 1 : 0; };
        const int nLost = accumulate(limit, stubsCell.end(), 0, valid);
        streamLost.reserve(nLost);
        for (auto it = limit; it != stubsCell.end(); it++)
          if (*it)
            streamLost.emplace_back((*it)->frame());
        stubsCell.erase(limit, stubsCell.end());
      }
      // cosmetics -- remove gaps at the end of stream
      for (auto it = stubsCell.end(); it != stubsCell.begin();)
        it = (*--it) == nullptr ? stubsCell.erase(it) : stubsCell.begin();
      if (!stubsCell.empty())
        stubsCell.back()->setNewTrk();
      // store final tracks
      stream.reserve(stubsCell.size());
      auto toFrame = [](StubMHT* stub) { return stub ? stub->frame() : FrameStub(); };
      transform(stubsCell.begin(), stubsCell.end(), back_inserter(stream), toFrame);
    }
  }

  // run single track through r-phi hough transform and store result
  void MiniHoughTransform::array(const vector<StubHT*>& track, deque<StubMHT*>& stream) {
    bool valid(false);
    for (int um = 0; um < setup_->mhtNumBinsInv2R(); um++) {
      const int im = um - setup_->mhtNumBinsInv2R() / 2.;
      const double m = (im + .5) * baseInv2R_;
      vector<TTBV> hitPatterns(setup_->mhtNumBinsPhiT(), TTBV(0, setup_->numLayers()));
      for (StubHT* stub : track) {
        const double chi = digi(stub->phi() - stub->r() * m, basePhi_);
        const double dChi = digi((stub->dPhi() + abs(stub->r()) * baseInv2R_) / 2., basePhi_);
        const double chiMin = chi - dChi;
        const double chiMax = chi + dChi;
        int cMin = floor(chiMin / basePhiT_ + 1.e-12);
        int cMax = floor(chiMax / basePhiT_ + 1.e-12);
        if (cMin > setup_->mhtNumBinsPhiT() / 2 - 1 || cMax < - setup_->mhtNumBinsPhiT() / 2)
          continue;
        cMin = max(cMin, - setup_->mhtNumBinsPhiT() / 2) + setup_->mhtNumBinsPhiT() / 2;
        cMax = min(cMax, setup_->mhtNumBinsPhiT() / 2 - 1) + setup_->mhtNumBinsPhiT() / 2;
        for (int c = cMin; c <= cMax; c++) {
          stub->mark(um, c);
          hitPatterns[c].set(stub->layer());
        }
      }
      for (int uc = 0; uc < setup_->mhtNumBinsPhiT(); uc++) {
        if (hitPatterns[uc].count() < setup_->mhtMinLayers())
          continue;
        valid = true;
        for (StubHT* stub : track)
          stub->hit(um, uc);
      }
    }
    if (!valid)
      return;
    for (StubHT* stub : track) {
      stubsMHT_.emplace_back(*stub, stub == track.back());
      stream.push_back(&stubsMHT_.back());
    }
  }

  // run single track through final r-phi cell and store result
  void MiniHoughTransform::cell(const vector<StubMHT*>& track, deque<StubMHT*>& stream) const {
    TTBV hitPattern(0, setup_->numLayers());
    for (StubMHT* stub : track)
      hitPattern.set(stub->layer());
    if (hitPattern.count() < setup_->mhtMinLayers())
      return;
    stream.insert(stream.end(), track.begin(), track.end());
  }

}  // namespace trackerTFP
