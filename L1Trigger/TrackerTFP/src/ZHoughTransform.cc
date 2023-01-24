#include "L1Trigger/TrackerTFP/interface/ZHoughTransform.h"

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

  ZHoughTransform::ZHoughTransform(const ParameterSet& iConfig,
                                   const Setup* setup,
                                   const DataFormats* dataFormats,
                                   int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        region_(region),
        input_(dataFormats->numChannel(Process::mht)),
        baseZ_(dataFormats_->base(Variable::z, Process::zht)),
        validZ0_(setup_->zhtNumBinsCot(), TTBV(0, setup_->zhtNumBinsZT())) {
    const Format<Variable::zT, Process::gp> gp(setup);
    const double dCot = (gp.base() + 2. * setup->beamWindowZ()) / setup->chosenRofZ();
    baseCot_ = dataFormats_->base(Variable::z, Process::dtc) / dataFormats_->base(Variable::r, Process::dtc);
    const int shiftCot = ceil(log2(dCot / setup->zhtNumBinsCot() / baseCot_));
    baseCot_ *= pow(2., shiftCot);
    baseZT_ = gp.base() / setup->zhtNumBinsZT();
    // prepare cut on z0 boundaries
    const double limitZ0 = setup_->beamWindowZ() + (baseZT_ + setup_->chosenRofZ() * baseCot_) / 2.;
    for (int um = 0; um < setup_->zhtNumBinsCot(); um++) {
      TTBV& ttbv = validZ0_[um];
      const int im = um - setup_->zhtNumBinsCot() / 2.;
      const double m = (im + .5) * baseCot_;
      for (int uc = setup_->zhtNumBinsZT() - 1; uc > -1; uc--) {
        const int ic = uc - setup_->zhtNumBinsZT() / 2.;
        const double c = (ic + .5) * baseZT_;
        const double z0 = c - setup_->chosenRofZ() * m;
        if (abs(z0) < limitZ0)
          ttbv.set(uc);
      }
    }
  }

  // read in and organize input product (fill vector input_)
  void ZHoughTransform::consume(const StreamsStub& streams) {
    auto valid = [](int sum, const FrameStub& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    const int offset = region_ * dataFormats_->numChannel(Process::mht);
    int nStubsMHT(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      nStubsMHT += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubsMHT_.reserve(nStubsMHT);
    stubsZHT_.reserve(nStubsMHT);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      vector<StubMHT*>& stubs = input_[channel];
      stubs.reserve(stream.size());
      // Store input stubs in vector, so rest of ZHT algo can work with pointers to them (saves CPU)
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
  void ZHoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) {
    const int offset = region_ * dataFormats_->numChannel(Process::zht);
    // loop over worker
    for (int channel = 0; channel < dataFormats_->numChannel(Process::zht); channel++) {
      // input data
      const vector<StubMHT*>& input = input_[channel];
      // output data
      StreamStub& stream = accepted[offset + channel];
      StreamStub& streamLost = lost[offset + channel];
      // identify tracks in input container
      int id;
      auto different = [&id](StubMHT* stub) { return !stub || id != stub->trackId(); };
      auto last = [](StubZHT* stub) { return stub && stub->newTrk(); };
      auto invalid = [](StubZHT* stub){ return !stub || !stub->valid(); };
      deque<StubZHT*> stubsArray;
      int sumSize = -setup_->zhtMinLayers();
      for (auto it = input.begin(); it != input.end();) {
        const auto start = find_if(it, input.end(), [](StubMHT* stub){ return stub; });
        id = (*start)->trackId();
        const auto end = find_if(start, input.end(), different);
        vector<StubMHT*> track(start, end);
        // restore clock accurancy
        sumSize += (int)track.size();
        const int delta = sumSize - (int)stubsArray.size();
        if (delta > 0)
          stubsArray.insert(stubsArray.end(), delta, nullptr);
        track.erase(remove(track.begin(), track.end(), nullptr), track.end());
        // run single track through r-z hough trasnform and store result
        array(track, stubsArray);
        // set begin of next track
        it = end;
      }
      deque<StubZHT*> stubsCell;
      sumSize = -setup_->mhtMinLayers();
      for (auto it = stubsArray.begin(); it != stubsArray.end();) {
        const auto start = it;
        auto end = find_if(start, stubsArray.end(), last);
        if (end != stubsArray.end())
          end++;
        vector<StubZHT*> track(start, end);
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
        auto valid = [](int& sum, StubZHT* stub) { return sum += stub ? 1 : 0; };
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
      if (region_ == 7 && channel == 11) {
        vector<stringstream> sss(4);
        for (StubZHT* stub : stubsCell) {
          if (!stub) {
            for (stringstream& ss : sss)
              ss << setw(7) << "x" << " ";
            continue;
          }
          sss[0] << setw(7) << stub->newTrk() << " ";
          sss[1] << setw(7) << dataFormats_->format(Variable::r, Process::zht).integer(stub->r()) << " ";
          sss[2] << setw(7) << dataFormats_->format(Variable::phi, Process::zht).integer(stub->phi()) << " ";
          sss[3] << setw(7) << dataFormats_->format(Variable::z, Process::zht).integer(stub->z()) << " ";
        }
        cout << "ZHT" << endl;
        for (const stringstream& ss : sss)
          cout << ss.str() << endl;
      }
      stream.reserve(stubsCell.size());
      auto toFrame = [](StubZHT* stub) { return stub ? stub->frame() : FrameStub(); };
      transform(stubsCell.begin(), stubsCell.end(), back_inserter(stream), toFrame);
    }
  }

  // run single track through r-z hough transform and store result
  void ZHoughTransform::array(const vector<StubMHT*>& track, deque<StubZHT*>& stream) {
    if (track.empty())
      return;
    const double limitZT = baseZT_ * setup_->zhtNumBinsZT() / 2.;
    bool valid(false);
    for (int um = 0; um < setup_->zhtNumBinsCot(); um++) {
      const TTBV& validZ0 = validZ0_[um];
      const int im = um - setup_->zhtNumBinsCot() / 2.;
      const double m = (im + .5) * baseCot_;
      vector<pair<TTBV, TTBV>> hitPatterns(setup_->zhtNumBinsZT(), {TTBV(0, setup_->numLayers()), TTBV(0, setup_->numLayers())});
      for (StubMHT* stub : track) {
        const double r = stub->r() + setup_->chosenRofPhi() - setup_->chosenRofZ();
        const double chi = digi(stub->z() - r * m, baseZ_);
        const double dChi = digi((stub->dZ() + abs(r) * baseCot_) / 2., baseZ_);
        const double zTMin = chi - dChi;
        const double zTMax = chi + dChi;
        // cut on sector boundaries
        if (zTMin >= limitZT || zTMax < -limitZT)
          continue;
        int cMin = floor(max(zTMin, -limitZT) / baseZT_ + 1.e-12);
        int cMax = floor(min(zTMax, limitZT) / baseZT_ + 1.e-12);
        cMin = max(cMin, - setup_->zhtNumBinsZT() / 2) + setup_->zhtNumBinsZT() / 2;
        cMax = min(cMax, setup_->zhtNumBinsZT() / 2 - 1) + setup_->zhtNumBinsZT() / 2;
        for (int c = cMin; c <= cMax; c++) {
          // cut on z0 boundaries
          if (!validZ0[c])
            continue;
          stub->mark(um, c);
          hitPatterns[c].first.set(stub->layer());
          if (stub->ps())
            hitPatterns[c].second.set(stub->layer());
        }
      }
      for (int uc = 0; uc < setup_->zhtNumBinsZT(); uc++) {
        const pair<TTBV, TTBV>& hitPattern = hitPatterns[uc];
        if (hitPattern.first.count() < setup_->zhtMinLayers() || hitPattern.second.count() < setup_->zhtMinLayersPS())
          continue;
        valid = true;
        for (StubMHT* stub : track)
          stub->hit(um, uc);
      }
    }
    if (!valid)
      return;
    for (StubMHT* stub : track) {
      stubsZHT_.emplace_back(*stub, stub == track.back());
      stream.push_back(&stubsZHT_.back());
    }
  }

  // run single track through final r-z cell and store result
  void ZHoughTransform::cell(const vector<StubZHT*>& track, deque<StubZHT*>& stream) const {
    pair<TTBV, TTBV> hitPattern(TTBV(0, setup_->numLayers()), TTBV(0, setup_->numLayers()));
    for (StubZHT* stub : track) {
      hitPattern.first.set(stub->layer());
      if (stub->ps())
        hitPattern.second.set(stub->layer());
    }
    if (hitPattern.first.count() < setup_->zhtMinLayers())
      return;
    if (hitPattern.second.count() < setup_->zhtMinLayersPS())
      return;
    stream.insert(stream.end(), track.begin(), track.end());
  }

}  // namespace trackerTFP
