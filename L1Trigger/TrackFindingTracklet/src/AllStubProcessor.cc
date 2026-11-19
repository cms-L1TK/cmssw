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
    const int offsetR = region_ * setup_->aspNumInput();
    for (int chan = 0; chan < setup_->aspNumInput(); chan++)
      num += streamsStub[offsetR + chan].size();
    stubs_.reserve(num);
    for (int layer = 0; layer < setup_->aspNumLayers(); layer++) {
      const int offsetL = layer * setup_->aspNumChannelPerLayer();
      for (int chan = 0; chan < setup_->aspNumChannelPerLayer(); chan++) {
        const int index = offsetL + chan;
        int stubId(chan * std::pow(2, setup_->aspWidthStubId()));
        std::deque<StubASP*>& stream = streams_[index];
        for (const tt::FrameStub& frame : streamsStub[offsetR + index]) {
          if (frame.first.isNull()) {
            stream.push_back(nullptr);
            continue;
          }
          // parse bits and convert stub
          const trackerDTC::SensorModule* sm = setup_->sensorModule(frame.first);
          const trackerDTC::SensorModule::Type type = sm->type();
          const int widthBend = setup_->aspWidthBend(type);
          const int widthAlpha = setup_->aspWidthAlpha(type);
          const int widthPhi = setup_->aspWidthPhi(type);
          const int widthZ = setup_->aspWidthZ(type);
          const int widthR = setup_->aspWidthR(type);
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
          // stub phi in rad
          double phi = ttBV.val(widthPhi) * basePhi;
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
          } else if (sm->barrel())
            r += setup_->stubLayerR(sm->layerIndex());
          // stub phi wrt nonant center in rad
          phi -= setup_->stubRangePhi() / 2.;
          // stub z in cm
          z += sm->offsetZ();
          // store stub
          stubs_.emplace_back(frame.first, dataFormats_, r, phi, z, stubId++);
          stream.push_back(&stubs_.back());
        }
      }
    }
  }

  // fill output products
  void AllStubProcessor::produce(tt::StreamsStub& streamsStub) {
    // mux streams
    int idOut(region_ * setup_->aspMux().size());
    for (const std::vector<int>& idsIn : setup_->aspMux()) {
      tt::StreamStub& stream = streamsStub[idOut++];
      // count stubs
      int num(0);
      for (int idIn : idsIn)
        num += streams_[idIn].size();
      stream.reserve(num);
      // mux stubs
      for (int idIn : idsIn)
        for (StubASP* stub : streams_[idIn])
          stream.emplace_back(stub->frame());
      // truncate stubs
      if (setup_->enableTruncation() && num > setup_->numFrames())
        stream.resize(setup_->numFrames());
    }
  }

}  // namespace trklet
