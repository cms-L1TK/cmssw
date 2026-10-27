#include "L1Trigger/TrackFindingTracklet/interface/DuplicateRemoval.h"

#include <vector>
#include <numeric>
#include <algorithm>

namespace trklet {

  DuplicateRemoval::DuplicateRemoval(const Setup* setup, const DataFormats* dataFormats, int region, const TTDTC& ttDTC)
      : setup_(setup), dataFormats_(dataFormats), region_(region) {
    if (setup_->drUseDTCStubs()) {
      // prep dtc stub container
      int size(0);
      for (int channel : ttDTC.tfpChannels())
        for (const tt::FrameStub& frame : ttDTC.stream(region, channel))
          if (frame.first.isNonnull())
            size++;
      dtc_.reserve(size);
      // fill dtc stub container
      for (int channel : ttDTC.tfpChannels())
        for (const tt::FrameStub& frame : ttDTC.stream(region, channel))
          if (frame.first.isNonnull())
            dtc_.push_back(frame);
    }
  }

  // read in and organize input tracks and stubs
  void DuplicateRemoval::consume(const tt::StreamsTrack& streamsTrack, const tt::StreamsStub& streamsStub) {
    auto nonNullTrack = [](int sum, const tt::FrameTrack& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    auto nonNullStub = [](int sum, const tt::FrameStub& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    // count tracks and stubs and reserve corresponding vectors
    int sizeStubs(0);
    const int offset = region_ * setup_->tmNumLayers();
    const tt::StreamTrack& streamTrack = streamsTrack[region_];
    stream_.reserve(streamTrack.size());
    const int sizeTracks = std::accumulate(streamTrack.begin(), streamTrack.end(), 0, nonNullTrack);
    for (int layer = 0; layer < setup_->tmNumLayers(); layer++) {
      const tt::StreamStub& streamStub = streamsStub[offset + layer];
      sizeStubs += std::accumulate(streamStub.begin(), streamStub.end(), 0, nonNullStub);
    }
    tracks_.reserve(sizeTracks);
    stubs_.reserve(sizeStubs);
    // store tracks and stubs
    for (int frame = 0; frame < static_cast<int>(streamTrack.size()); frame++) {
      const tt::FrameTrack& frameTrack = streamTrack[frame];
      if (frameTrack.first.isNull()) {
        stream_.push_back(nullptr);
        continue;
      }
      tracks_.emplace_back();
      Track& track = tracks_.back();
      track.stubs_ = std::vector<Stub*>(setup_->tmNumLayers(), nullptr);
      // parse track bits
      TTBV ttBV(frameTrack.second);
      track.cot_ = ttBV.val(setup_->tbWidthCot(), 0, true) * setup_->tbBaseCot();
      ttBV >>= setup_->tbWidthCot();
      track.z0_ = ttBV.val(setup_->tbWidthZ0(), 0, true) * setup_->tbBaseZ0();
      ttBV >>= setup_->tbWidthZ0();
      track.phi0_ = ttBV.val(setup_->tbWidthPhi0(), 0, true) * setup_->tbBasePhi0();
      ttBV >>= setup_->tbWidthPhi0();
      track.inv2R_ = ttBV.val(setup_->tbWidthInv2R(), 0, true) * setup_->tbBaseInv2R();
      ttBV >>= setup_->tbWidthInv2R();
      const int seedType = ttBV.val(setup_->tbWidthSeedType());
      // parse stubs
      const std::vector<int>& layersSeed = setup_->tbSeedLayers(seedType);
      for (int layer = 0; layer < setup_->tmNumLayers(); layer++) {
        const tt::FrameStub& frameStub = streamsStub[offset + layer][frame];
        if (frameStub.first.isNull())
          continue;
        stubs_.emplace_back(frameStub.first);
        Stub& stub = stubs_.back();
        stub.sm_ = setup_->sensorModule(frameStub.first);
        const auto itS = std::find(layersSeed.begin(), layersSeed.end(), stub.sm_->layerId());
        if (itS == layersSeed.end()) {
          const int widthR = setup_->tbWidthR(stub.sm_->type());
          const int widthPhi = setup_->tbWidthPhi(stub.sm_->type());
          const int widthRZ = setup_->tbWidthRZ(stub.sm_->type());
          const double baseR = setup_->tbBaseR(stub.sm_->type());
          const double basePhi = setup_->tbBasePhi(stub.sm_->type());
          const double baseRZ = setup_->tbBaseRZ(stub.sm_->type());
          TTBV ttBV(frameStub.second);
          stub.z_ = ttBV.val(widthRZ, 0, true) * baseRZ;
          ttBV >>= widthRZ;
          stub.phi_ = ttBV.val(widthPhi, 0, true) * basePhi;
          ttBV >>= widthPhi;
          stub.r_ = ttBV.val(widthR, 0, stub.sm_->barrel()) * baseR;
          ttBV >>= widthR;
          stub.stubId_ = ttBV.val(setup_->tbWidthStubId());
          stub.seed_ = false;
        } else {
          stub.stubId_ = TTBV(frameStub.second).val(setup_->tbWidthStubId());
          stub.seed_ = true;
        }
        track.stubs_[layer] = &stubs_.back();
      }
      // store output track
      track.trackDR_ = TrackDR(frameTrack.first, dataFormats_, seedType);
      stream_.push_back(&tracks_.back());
    }
    // remove all gaps between end and last track
    for (auto it = stream_.end(); it != stream_.begin();)
      it = (*--it) ? stream_.begin() : stream_.erase(it);
  }

  // fill output products
  void DuplicateRemoval::produce(tt::StreamsTrack& streamsTrack, tt::StreamsStub& streamsStub) {
    // remove duplicated tracks, no merge of stubs, one stub per layer expected
    algo();
    // unify stubs
    unify();
    // calc stub position and uncertainties
    pos();
    // replace stubs with DTC stubs
    if (setup_->drUseDTCStubs())
      dtc();
    // replace stubs with TT stubs
    if (setup_->drUseTTStubs())
      tt();
    // remove HO corrections from all stubs
    if (setup_->kfUseSimulation())
      sim();
    // calc stub uncertainties
    delta();
    // base transformation
    redigi();
    // store output
    store(streamsTrack, streamsStub);
  }

  // remove duplicated tracks, no merge of stubs, one stub per layer expected
  void DuplicateRemoval::algo() {
    std::vector<Track*> cms(setup_->drNumComparisonModules(), nullptr);
    for (Track*& track : stream_) {
      if (!track)
        // gaps propagate through chain and appear in output stream
        continue;
      for (Track*& trackCM : cms) {
        if (!trackCM) {
          // tracks used in CMs propagate through chain and do appear in output stream unaltered
          trackCM = track;
          break;
        }
        if (equalEnough(track, trackCM)) {
          // tracks compared in CMs propagate through chain and appear in output stream as gap if identified as duplicate or unaltered elsewise
          track = nullptr;
          break;
        }
      }
    }
    // remove all gaps between end and last track
    for (auto it = stream_.end(); it != stream_.begin();)
      it = (*--it) ? stream_.begin() : stream_.erase(it);
  }

  // unify stubs
  void DuplicateRemoval::unify() {
    for (Track* track : stream_) {
      if (!track)
        continue;
      track->phi0_ -= tt::digiR(setup_->stubRangePhi() / 2., setup_->tbBasePhi0());
      track->inv2R_ *= -1.;
      for (Stub* stub : track->stubs_) {
        if (!stub)
          continue;
        if (stub->seed_) {
          if (stub->sm_->barrel()) {
            stub->r_ = tt::digi(setup_->stubLayerR(stub->sm_->layerIndex()), setup_->tbBaseR());
            const double z = tt::digi(track->z0_ + stub->r_ * track->cot_, setup_->tbBaseZ0());
            if (std::abs(z) > tt::digi(setup_->tbMinZ(), setup_->tbBaseZ0()) && stub->sm_->layerIndex() == 0)
              stub->r_ = tt::digi(setup_->tbInnerRadius(), setup_->tbBaseR());
          } else {
            const double z = tt::digi(setup_->stubDiskZ(stub->sm_->layerIndex()), setup_->tbBaseZ());
            const double cot = tt::digi(std::abs(track->cot_), setup_->drBaseCot());
            const double invCot = tt::digi(1. / cot, setup_->drBaseInvCot());
            const double side = track->cot_ < 0. ? -1. : 1.;
            stub->r_ = tt::digi((z - side * track->z0_) * invCot, setup_->tbBaseR());
          }
        } else {
          if (stub->sm_->barrel()) {
            stub->r_ += tt::digi(setup_->stubLayerR(stub->sm_->layerIndex()), setup_->tbBaseR());
          } else
            stub->z_ = tt::digi(-stub->z_ * track->cot_, setup_->tbBaseZ());
          if (stub->sm_->type() == trackerDTC::SensorModule::Disk2S)
            stub->r_ = tt::digi(setup_->stubDiskR(stub->sm_->layerIndex(), stub->r_), setup_->tbBaseR());
        }
      }
    }
  }

  // calc stub position
  void DuplicateRemoval::pos() {
    // unify stubs
    for (Track* track : stream_) {
      if (!track)
        continue;
      for (Stub* stub : track->stubs_) {
        if (!stub)
          continue;
        stub->phi_ += track->phi0_ + stub->r_ * track->inv2R_;
        stub->z_ += track->z0_ + stub->r_ * track->cot_;
      }
    }
  }

  // replace stubs with DTC stubs
  void DuplicateRemoval::dtc() {
    for (Track* track : stream_) {
      if (!track)
        continue;
      for (Stub* stub : track->stubs_) {
        if (!stub)
          continue;
        auto via = [stub](const tt::FrameStub& fs) { return fs.first == stub->ttStubRef_; };
        const tt::FrameStub& fs = *std::find_if(dtc_.begin(), dtc_.end(), via);
        const GlobalPoint gp = setup_->stubPosDTC(fs, region_);
        stub->r_ = gp.perp();
        stub->phi_ = tt::deltaPhi(gp.phi() - region_ * setup_->regRangePhiT());
        stub->z_ = gp.z();
        // apply ho corrections
        correct(track, stub);
      }
    }
  }

  // replace stubs with TT stubs
  void DuplicateRemoval::tt() {
    for (Track* track : stream_) {
      if (!track)
        continue;
      for (Stub* stub : track->stubs_) {
        if (!stub)
          continue;
        const GlobalPoint gp = setup_->stubPosTT(stub->ttStubRef_);
        stub->r_ = gp.perp();
        stub->phi_ = tt::deltaPhi(gp.phi() - region_ * setup_->regRangePhiT());
        stub->z_ = gp.z();
        // apply ho corrections
        correct(track, stub);
      }
    }
  }

  // remove HO corrections from all stubs
  void DuplicateRemoval::sim() {
    for (Track& track : tracks_) {
      for (Stub*& stub : track.stubs_) {
        if (!stub)
          continue;
        correct(track, stub, false);
      }
    }
  }

  // add (or remove) HO corrections from stub
  void DuplicateRemoval::correct(const Track& track, Stub* stub, bool add) {
    if (track.inv2R_ == 0.)
      return;
    const double trackPhiL = track.phi0_ + track.inv2R_ * stub->r_;
    const double trackZL = track.z0_ + track.cot_ * stub->r_;
    const double trackPhi = track.phi0_ + std::asin(stub->r_ * track.inv2R_);
    const double trackZ = track.z0_ + track.cot_ * std::asin(stub->r_ * track.inv2R_) / track.inv2R_;
    if (add) {
      stub->phi_ += trackPhiL - trackPhi;
      stub->z_ += trackZL - trackZ;
    } else {
      stub->phi_ -= trackPhiL - trackPhi;
      stub->z_ -= trackZL - trackZ;
    }
  }

  // calc stub uncertainties
  void DuplicateRemoval::delta() {
    for (Track* track : stream_) {
      if (!track)
        continue;
      const double cot = tt::digi(std::abs(track->cot_), setup_->drBaseCot());
      const double inv2R = tt::digi(std::abs(track->inv2R_), setup_->drBaseInv2R());
      const double scat =
          tt::digi(setup_->smScattering() * inv2R * .5, dataFormats_->base(Variable::dPhi, Process::dr));
      for (Stub* stub : track->stubs_) {
        if (!stub)
          continue;
        const double r = tt::digi(stub->r_, setup_->drBaseR());
        const double res = tt::digi(stub->sm_->pitchRow() / r * .5, dataFormats_->base(Variable::dPhi, Process::dr));
        const double dT =
            stub->sm_->pitchCol() * (setup_->smTiltApproxSlope() * cot + setup_->smTiltApproxIntercept()) / 2.;
        const double dB = stub->sm_->pitchCol() / 2.;
        const double dD = stub->sm_->pitchCol() * cot / 2.;
        const int id = stub->sm_->layerIndexCombined();
        const bool t = id < setup_->sysNumBarrelLayerPS() &&
                       std::abs(stub->z_) > tt::digi(setup_->tbTiltedLayerLimitsZ(id), setup_->tbBaseZ());
        stub->dPhi_ = scat + res;
        stub->dZ_ = dD;
        if (stub->sm_->barrel())
          stub->dZ_ = t ? dT : dB;
        stub->dPhi_ += dataFormats_->base(Variable::dPhi, Process::dr);
        stub->dZ_ += dataFormats_->base(Variable::dZ, Process::dr);
      }
    }
  }

  // base transformation
  void DuplicateRemoval::redigi() {
    for (Track* track : stream_) {
      if (!track)
        continue;
      TTBV hitPattern(0, setup_->drNumLayers());
      for (Stub* stub : track->stubs_) {
        if (!stub)
          continue;
        const bool rounded = !stub->seed_ && stub->sm_->type() == trackerDTC::SensorModule::Type::DiskPS;
        stub->r_ = rounded ? tt::digiR(stub->r_, setup_->tbBaseR()) : tt::digi(stub->r_, setup_->tbBaseR());
        stub->phi_ = tt::digi(stub->phi_, setup_->tbBasePhi());
        stub->z_ = tt::digi(stub->z_, setup_->tbBaseZ());
        stub->dPhi_ = tt::digi(stub->dPhi_, dataFormats_->base(Variable::dPhi, Process::dr));
        stub->dZ_ = tt::digi(stub->dZ_, dataFormats_->base(Variable::dZ, Process::dr));
        // base transform into high precision format
        stub->r_ =
            tt::redigi(stub->r_, setup_->tbBaseR(), dataFormats_->base(Variable::r, Process::dr), setup_->widthDSPbu());
        stub->phi_ = tt::redigi(
            stub->phi_, setup_->tbBasePhi(), dataFormats_->base(Variable::phi, Process::dr), setup_->widthDSPbu());
        stub->z_ =
            tt::redigi(stub->z_, setup_->tbBaseZ(), dataFormats_->base(Variable::z, Process::dr), setup_->widthDSPbu());
      }
    }
  }

  // store output
  void DuplicateRemoval::store(tt::StreamsTrack& streamsTrack, tt::StreamsStub& streamsStub) const {
    const int offset = region_ * setup_->drNumLayers();
    tt::StreamTrack& streamTrack = streamsTrack[region_];
    streamTrack.reserve(stream_.size());
    for (int layer = 0; layer < setup_->drNumLayers(); layer++)
      streamsStub[offset + layer].reserve(stream_.size());
    for (Track* track : stream_) {
      if (!track) {
        streamTrack.emplace_back();
        for (int layer = 0; layer < setup_->drNumLayers(); layer++)
          streamsStub[offset + layer].emplace_back();
        continue;
      }
      streamTrack.push_back(track->trackDR_.frame());
      int id(0);
      for (Stub* s : track->stubs_) {
        if (!s) {
          streamsStub[offset + id++].emplace_back();
          continue;
        }
        const StubDR stubDR(s->ttStubRef_, dataFormats_, s->r_, s->phi_, s->z_, s->dPhi_, s->dZ_);
        streamsStub[offset + id++].push_back(stubDR.frame());
      }
    }
  }

  // compares two tracks, returns true if those are considered duplicates
  bool DuplicateRemoval::equalEnough(Track* t0, Track* t1) const {
    int same(0);
    for (int layer = 0; layer < setup_->tmNumLayers(); layer++) {
      const Stub* s0 = t0->stubs_[layer];
      const Stub* s1 = t1->stubs_[layer];
      if (s0 && s1 && s0->stubId_ == s1->stubId_)
        same++;
    }
    return same >= setup_->drMinIdenticalStubs();
  }

}  // namespace trklet
