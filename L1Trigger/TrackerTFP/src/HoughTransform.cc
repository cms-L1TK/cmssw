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
                                 deque<StubHT>& stubs)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        inv2R_(&dataFormats_->format(Variable::inv2R, Process::ht)),
        phiT_(&dataFormats_->format(Variable::phiT, Process::ht)),
        phi_(&dataFormats_->format(Variable::phi, Process::ht)),
        stubsHT_(stubs) {}

  // fill output products
  void HoughTransform::produce(const vector<vector<StubGP*>>& streamsIn,
                               vector<deque<StubHT*>>& streamsOut,
                               vector<deque<StubHT*>>& streamsTrunc) {
    static const int numChannelIn = dataFormats_->numChannel(Process::gp);
    static const int numChannelOut = dataFormats_->numChannel(Process::ht);
    static const int chan = setup_->kfNumWorker();
    static const int mux = numChannelOut / chan;
    for (int channelOut = 0; channelOut < numChannelOut; channelOut++) {
      const int inv2Ru = mux * (channelOut % chan) + channelOut / chan;
      const int inv2R = inv2R_->toSigned(inv2Ru);
      deque<StubHT*>& output = streamsOut[channelOut];
      deque<StubHT*>& truncated = streamsTrunc[channelOut];
      for (int channelIn = numChannelIn - 1; channelIn >= 0; channelIn--) {
        const vector<StubGP*>& input = streamsIn[channelIn];
        deque<StubHT*> stubs;
        // associate stubs with inv2R and phiT bins
        fillIn(inv2R, channelIn, input, stubs, truncated);
        // apply truncation
        if (enableTruncation_ && (int)stubs.size() > setup_->numFrames()) {
          const auto limit = next(stubs.begin(), setup_->numFrames());
          copy_if(limit, stubs.end(), back_inserter(truncated), [](StubHT* stub) { return stub; });
          stubs.erase(limit, stubs.end());
        }
        // ht collects all stubs before readout starts -> remove all gaps
        stubs.erase(remove(stubs.begin(), stubs.end(), nullptr), stubs.end());
        // identify tracks
        readOut(stubs, output);
      }
      // apply truncation
      if (enableTruncation_ && (int)output.size() > setup_->numFrames()) {
        const auto limit = next(output.begin(), setup_->numFrames());
        copy_if(limit, output.end(), back_inserter(truncated), [](StubHT* stub) { return stub; });
        output.erase(limit, output.end());
      }
    }
  }

  // associate stubs with phiT bins in this inv2R column
  void HoughTransform::fillIn(
      int inv2R, int sector, const vector<StubGP*>& input, deque<StubHT*>& output, deque<StubHT*>& truncated) {
    auto inv2RrangeCheck = [inv2R](StubGP* stub) {
      return (stub && stub->inv2RMin() <= inv2R && stub->inv2RMax() >= inv2R) ? stub : nullptr;
    };
    auto convert = [this, sector](StubGP* stub, int inv2R, int phiTht) {
      static const DataFormat& gp = dataFormats_->format(Variable::phiT, Process::gp);
      const double r = stub->r();
      const double phi = stub->phi() - (inv2R_->floating(inv2R) * r + phiT_->floating(phiTht));
      const double z = stub->z();
      const TTBV& layer = stub->layer();
      const int phiTgp = gp.toSigned(sector % setup_->gpNumBinsPhiT());
      const int phiT = phiT_->integer(gp.floating(phiTgp) + phiT_->floating(phiTht));
      const int zT = sector / setup_->gpNumBinsPhiT() - setup_->gpNumBinsZT() / 2;
      stubsHT_.emplace_back(*stub, r, phi, z, layer, phiT, zT);
      return &stubsHT_.back();
    };
    // Latency of ht fifo firmware
    static constexpr int latency = 1;
    // static delay container
    deque<StubHT*> delay(latency, nullptr);
    // fifo, used to store stubs which belongs to a second possible track
    deque<StubHT*> stack;
    // stream of incroming stubs
    deque<StubGP*> stream;
    transform(input.begin(), input.end(), back_inserter(stream), inv2RrangeCheck);
    // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
    while (!stream.empty() || !stack.empty() ||
           !all_of(delay.begin(), delay.end(), [](const StubHT* stub) { return !stub; })) {
      StubHT* stubHT = nullptr;
      StubGP* stubGP = pop_front(stream);
      if (stubGP) {
        const double phiT = stubGP->phi() - inv2R_->floating(inv2R) * stubGP->r();
        const int major = phiT_->integer(phiT);
        if (major >= -setup_->htNumBinsPhiT() / 2 && major < setup_->htNumBinsPhiT() / 2) {
          // major candidate has pt > threshold (3 GeV)
          // stubHT records which HT bin this stub is added to
          stubHT = convert(stubGP, inv2R, major);
        }
        const double chi = phi_->digi(phiT - phiT_->floating(major));
        if (abs(stubGP->r() * inv2R_->base()) + 2. * abs(chi) >= phiT_->base()) {
          // stub belongs to two candidates
          const int minor = chi >= 0. ? major + 1 : major - 1;
          if (minor >= -setup_->htNumBinsPhiT() / 2 && minor < setup_->htNumBinsPhiT() / 2) {
            // second (minor) candidate has pt > threshold (3 GeV)
            StubHT* stub = convert(stubGP, inv2R, minor);
            delay.push_back(stub);
          }
        }
      }
      // add nullptr to delay pipe if stub didn't fill any cell
      if ((int)delay.size() == latency)
        delay.push_back(nullptr);
      // take fifo latency into account (read before write)
      StubHT* stub = pop_front(delay);
      if (stub) {
        if (enableTruncation_ && (int)stack.size() == setup_->htDepthMemory() - 1)
          // buffer overflow
          truncated.push_back(pop_front(stack));
        // store minor stub in fifo
        stack.push_back(stub);
      }
      // take a minor stub if no major stub available
      output.push_back(stubHT ? stubHT : pop_front(stack));
    }
  }

  // identify tracks
  void HoughTransform::readOut(const deque<StubHT*>& input, deque<StubHT*>& output) const {
    auto toBinPhiT = [this](StubHT* stub) {
      static const DataFormat& gp = dataFormats_->format(Variable::phiT, Process::gp);
      const double phiT = phiT_->floating(stub->phiT());
      const double local = phiT - gp.digi(phiT);
      return phiT_->integer(local) + setup_->htNumBinsPhiT() / 2;
    };
    auto toLayerId = [this](StubHT* stub) {
      static const DataFormat& layer = dataFormats_->format(Variable::layer, Process::ctb);
      return stub->layer().val(layer.width());
    };
    // used to recognise in which order tracks are found
    TTBV trkFoundPhiTs(0, setup_->htNumBinsPhiT());
    // hitPattern for all possible tracks, used to find tracks
    vector<TTBV> patternHits(setup_->htNumBinsPhiT(), TTBV(0, setup_->numLayers()));
    // found phiTs, ordered in time
    vector<int> phiTs;
    phiTs.reserve(setup_->htNumBinsPhiT());
    for (StubHT* stub : input) {
      const int binPhiT = toBinPhiT(stub);
      const int layerId = toLayerId(stub);
      TTBV& pattern = patternHits[binPhiT];
      pattern.set(layerId);
      if (trkFoundPhiTs[binPhiT] || pattern.count() < setup_->htMinLayers())
        continue;
      // first time track found
      trkFoundPhiTs.set(binPhiT);
      phiTs.push_back(binPhiT);
    }
    // read out found tracks ordered as found
    for (int phiT : phiTs) {
      auto samePhiT = [phiT, toBinPhiT, this](StubHT* stub) { return toBinPhiT(stub) == phiT; };
      // read out stubs in reverse order to emulate f/w (backtracking linked list)
      copy_if(input.rbegin(), input.rend(), back_inserter(output), samePhiT);
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

}  // namespace trackerTFP