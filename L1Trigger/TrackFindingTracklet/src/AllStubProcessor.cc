#include "L1Trigger/TrackFindingTracklet/interface/AllStubProcessor.h"

#include <vector>
#include <deque>
#include <set>
#include <numeric>
#include <algorithm>

namespace trklet {

  AllStubProcessor::AllStubProcessor(const Setup* setup, const DataFormats* dataFormats, int region)
      : setup_(setup), dataFormats_(dataFormats), region_(region), streams_(setup_->aspNumInput()) {}

  // read in and organize input tracks and stubs
  void AllStubProcessor::consume(const tt::StreamsStub& streamsStub) {
    // count stubs
    int num(0);
    const int offset = region_ * setup_->aspNumInput();
    for (int chan = 0; chan < setup_->aspNumInput(); chan++)
      num += streamsStub[offset + chan].size();
    stubs_.reserve(num);
    for (int chan = 0; chan < setup_->aspNumInput(); chan++) {
      int stubId = std::pow(2, setup_->aspWidthStubId()) * (chan / setup_->aspNumChannelPerLayer());
      std::deque<StubASP*>& stream = streams_[chan];
      for (const tt::FrameStub& frame : streamsStub[offset + chan]) {
        if (frame.first.isNull()) {
          stream.push_back(nullptr);
          continue;
        }
        // parse bits and convert stub
        const trackerDTC::SensorModule* sm = setup_->sensorModule(frame.first);
        const trackerDTC::SensorModule::Type type = sm->type();
        const int widthBend = setup_->stubWidthBend(type);
        const int widthAlpha = setup_->stubWidthAlpha(type);
        const int widthPhi = setup_->stubWidthPhi(type);
        const int widthZ = setup_->stubWidthZ(type);
        const int widthR = setup_->stubWidthR(type);
        const double baseAlpha = setup_->stubBasePhi(type);
        const double basePhi = setup_->stubBasePhi(type);
        const double baseZ = setup_->stubBaseZ(type);
        const double baseR = setup_->stubBaseR(type);
        TTBV ttBV = frame.second;
        //const int bend = ttBV.val(widthBend, 0, true);
        ttBV >>= widthBend;
        TTBV alpha;
        if (type == trackerDTC::SensorModule::Disk2S)
          alpha = TTBV(ttBV, widthAlpha, 0, true);
        ttBV >>= widthAlpha;
        // stub phi w.r.t. to nonant edge in rad
        double phi = ttBV.val(widthPhi, 0, true) * basePhi;
        ttBV >>= widthPhi;
        // stub z w.r.t. an offset in cm
        double z = ttBV.val(widthZ, 0, true) * baseZ;
        ttBV >>= widthZ;
        // stub r w.r.t. an offset in cm
        double r = ttBV.val(widthR, 0, sm->barrel()) * baseR;
        // stub r in cm
        if (type == trackerDTC::SensorModule::Disk2S) {
          const double row = alpha.val(baseAlpha);
          r = setup_->stubDiskR(sm->layerIndex(), r);
          r = std::sqrt(std::pow(r, 2) + std::pow(row * setup_->cbcPitch(), 2));
        } else
          r += sm->offsetR() + setup_->regChosenRofPhi();
        // stub z in cm
        z += sm->offsetZ();
        // stub phi in rad
        phi = tt::deltaPhi(phi - setup_->stubRangePhi() / 2. + region_ * setup_->regRangePhiT());
        stubs_.emplace_back(frame.first, dataFormats_, r, phi, z, stubId++);
        stream.push_back(&stubs_.back());
      }
    }
  }

  // fill output products
  void AllStubProcessor::produce(tt::StreamsStub& streamsStub) {
    // mux streams
    int idOut(0);
    for (const std::vector<int>& idsIn : setup_->aspMux()) {
      tt::StreamStub& stream = streamsStub[idOut];
      // count stubs
      int num(0);
      for (int idIn : idsIn)
        num += streams_[idIn].size();
      stream.reserve(num);
      // mux sstubs
      for (int idIn : idsIn)
        for (StubASP* stub : streams_[idIn])
          stream.emplace_back(stub->frame());
      // truncate stubs
      if (setup_->enableTruncation() && num > setup_->numFrames())
        stream.resize(setup_->numFrames());
    }
  }

}  // namespace trklet
