#include "L1Trigger/TrackerTFP/interface/HoughTransform.h"

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

  HoughTransform::HoughTransform(const ParameterSet& iConfig,
                                 const Setup* setup,
                                 const DataFormats* dataFormats,
                                 int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        inv2R_(dataFormats_->format(Variable::inv2R, Process::ht)),
        phiT_(dataFormats_->format(Variable::phiT, Process::ht)),
        region_(region),
        input_(setup_->htNumBinsInv2R(), vector<vector<StubGP*>>(dataFormats_->numChannel(Process::gp))) {}

  // read in and organize input product (fill vector input_)
  void HoughTransform::consume(const StreamsStub& streams) {
    const int offset = region_ * dataFormats_->numChannel(Process::gp);
    auto validFrame = [](int sum, const FrameStub& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    int nStubsGP(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::gp); channel++) {
      const StreamStub& stream = streams[channel + offset];
      const int n = accumulate(stream.begin(), stream.end(), 0, validFrame);
      for (vector<vector<StubGP*>>& binInv2R : input_)
        binInv2R[channel].reserve(n);
      nStubsGP += n;
    }
    stubsGP_.reserve(nStubsGP);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::gp); channel++) {
      for (const FrameStub& frame : streams[channel + offset]) {
        // Store input stubs in vector, so rest of HT algo can work with pointers to them (saves CPU)
        StubGP* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubsGP_.emplace_back(frame, dataFormats_, channel);
          stub = &stubsGP_.back();
        }
        for (int binInv2R = 0; binInv2R < setup_->htNumBinsInv2R(); binInv2R++)
          input_[binInv2R][channel].push_back(stub && stub->inInv2RBin(binInv2R) ? stub : nullptr);
      }
    }
    auto validStub = [](int& sum, StubGP* stub) { return sum += stub ? 1 : 0; };
    int nStubsHT(0);
    for (const vector<vector<StubGP*>>& binInv2R : input_)
      for (const vector<StubGP*>& sector : binInv2R)
        nStubsHT += 2 * accumulate(sector.begin(), sector.end(), 0, validStub);
    stubsHT_.reserve(nStubsHT);
  }

  // fill output products
  void HoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) {
    const int offset = region_ * dataFormats_->numChannel(Process::ht);
    auto put = [](const deque<StubHT*>& stubs, StreamStub& stream) {
      stream.reserve(stubs.size());
      for (StubHT* stub : stubs)
        stream.emplace_back(stub ? stub->frame() : FrameStub());
    };
    for (int binInv2R = 0; binInv2R < setup_->htNumBinsInv2R(); binInv2R++) {
      const int channel = binInv2R + offset;
      const int inv2R = inv2R_.toSigned(binInv2R);
      deque<StubHT*> stubsAccepted;
      deque<StubHT*> stubsLost;
      for (vector<StubGP*>& sectorStubs : input_[binInv2R]) {
        deque<StubHT*> sector;
        // associate stubs with inv2R and phiT bins
        fillIn(inv2R, sectorStubs, sector, stubsLost);
        // ht collects all stubs before readout starts -> remove all gaps
        sector.erase(remove(sector.begin(), sector.end(), nullptr), sector.end());
        // identify tracks
        readOut(sector, stubsAccepted);
      }
      // apply truncation
      if (enableTruncation_ && (int)stubsAccepted.size() > setup_->numFrames()) {
        const auto limit = next(stubsAccepted.begin(), setup_->numFrames());
        copy_if(limit, stubsAccepted.end(), back_inserter(stubsLost), [](StubHT* stub) { return stub; });
        stubsAccepted.erase(limit, stubsAccepted.end());
        // cosmetics -- remove gaps at the end of stream
        for (auto it = stubsAccepted.end(); it != stubsAccepted.begin();)
          it = (*--it) == nullptr ? stubsAccepted.erase(it) : stubsAccepted.begin();
        stubsAccepted.back()->newTrk();
      }
      // store found tracks
      put(stubsAccepted, accepted[channel]);
      // store lost stubs
      put(stubsLost, lost[channel]);
    }
  }

  // associate stubs with phiT bins in this inv2R column
  void HoughTransform::fillIn(int inv2R, vector<StubGP*>& input, deque<StubHT*>& accepted, deque<StubHT*>& lost) {
    // fifo, used to store stubs which belongs to a second possible track
    deque<StubHT*> stack;
    reverse(input.begin(), input.end());
    // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
    while (!input.empty() || !stack.empty()) {
      StubHT* stubHT = nullptr;
      StubGP* stubGP = pop_back(input);
      if (stubGP) {
        const double phiT = stubGP->phi() - inv2R_.floating(inv2R) * stubGP->r();
        const int major = phiT_.integer(phiT);
        if (major >= -setup_->htNumBinsPhiT() / 2 && major < setup_->htNumBinsPhiT() / 2) {
          // major candidate has pt > threshold (3 GeV)
          // stubHT records which HT bin this stub is added to
          stubsHT_.emplace_back(*stubGP, inv2R, major);
          stubHT = &stubsHT_.back();
        }
        const double chi = phiT - phiT_.floating(major);
        if (abs(stubGP->r() * inv2R_.base()) + 2. * abs(chi) >= phiT_.base()) {
          // stub belongs to two candidates
          const int minor = chi >= 0. ? major + 1 : major - 1;
          if (minor >= -setup_->htNumBinsPhiT() / 2 && minor < setup_->htNumBinsPhiT() / 2) {
            // second (minor) candidate has pt > threshold (3 GeV)
            stubsHT_.emplace_back(*stubGP, inv2R, minor);
            if (enableTruncation_ && (int)stack.size() == setup_->htDepthMemory() - 1)
              // buffer overflow
              lost.push_back(pop_front(stack));
            // store minor stub in fifo
            stack.push_back(&stubsHT_.back());
          }
        }
      }
      // take a minor stub if no major stub available
      accepted.push_back(stubHT ? stubHT : pop_front(stack));
    }
    // truncate to many input stubs
    if (!enableTruncation_ || (int)accepted.size() <= setup_->numFrames())
      return;
    const auto limit = next(accepted.begin(), setup_->numFrames());
    copy_if(limit, accepted.end(), back_inserter(lost), [](StubHT* stub) { return stub; });
    accepted.erase(limit, accepted.end());
  }

  // identify tracks
  void HoughTransform::readOut(const deque<StubHT*>& input, deque<StubHT*>& output) const {
    // used to recognise in which order tracks are found
    TTBV trkFoundPhiTs(0, setup_->htNumBinsPhiT());
    // hitPattern for all possible tracks, used to find tracks
    vector<TTBV> patternHits(setup_->htNumBinsPhiT(), TTBV(0, setup_->numLayers()));
<<<<<<< HEAD
    // found unsigned phiTs, ordered in time
    vector<int> binsPhiT;
    // stub container for all possible tracks
    vector<vector<StubHT*>> tracks(setup_->htNumBinsPhiT());
    for (int binPhiT = 0; binPhiT < setup_->htNumBinsPhiT(); binPhiT++) {
      const int phiT = phiT_.toSigned(binPhiT);
      auto samePhiT = [phiT](int sum, StubHT* stub) { return sum + (stub->phiT() == phiT); };
      const int numAccepted = accumulate(acceptedSector.begin(), acceptedSector.end(), 0, samePhiT);
      const int numLost = accumulate(lostSector.begin(), lostSector.end(), 0, samePhiT);
      tracks[binPhiT].reserve(numAccepted + numLost);
    }
    for (StubHT* stub : acceptedSector) {
      const int binPhiT = phiT_.toUnsigned(stub->phiT());
=======
    // found phiTs, ordered in time
    vector<int> phiTs;
    phiTs.reserve(setup_->htNumBinsPhiT());
    for (StubHT* stub : input) {
      const int binPhiT = stub->phiTlocal() + setup_->htNumBinsPhiT() / 2;
>>>>>>> ebfddccae6a (KFin emulation added, mini hts and kf updated.)
      TTBV& pattern = patternHits[binPhiT];
      pattern.set(stub->layer());
      if (trkFoundPhiTs[binPhiT] || pattern.count() < setup_->htMinLayers())
        continue;
      // first time track found
      trkFoundPhiTs.set(binPhiT);
      phiTs.push_back(stub->phiTlocal());
    }
    // read out found tracks ordered as found
    for (int phiT : phiTs) {
      copy_if(input.begin(), input.end(), back_inserter(output), [phiT](StubHT* stub){ return stub->phiTlocal() == phiT; });
      output.back()->newTrk();
    }
  }

  // remove and return first element of deque, returns nullptr if empty
  template <class T>
  T* HoughTransform::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

<<<<<<< HEAD
}  // namespace trackerTFP
=======
  // remove and return last element of vector, returns nullptr if empty
  template <class T>
  T* HoughTransform::pop_back(vector<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.back();
      ts.pop_back();
    }
    return t;
  }

}  // namespace trackerTFP
>>>>>>> ebfddccae6a (KFin emulation added, mini hts and kf updated.)
