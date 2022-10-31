#include "L1Trigger/TrackFindingTracklet/interface/KFin.h"

#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace edm;
using namespace tt;
using namespace trackerTFP;

namespace trklet {

  KFin::KFin(const ParameterSet& iConfig,
             const Setup* setup,
             const DataFormats* dataFormats,
             const LayerEncoding* layerEncoding,
             const ChannelAssignment* channelAssignment,
             int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        layerEncoding_(layerEncoding),
        channelAssignment_(channelAssignment),
        region_(region),
        input_(channelAssignment_->numNodesDR()) {}

  // read in and organize input tracks and stubs
  void KFin::consume(const StreamsTrack& streamsTrack, const StreamsStub& streamsStub) {
    const int offsetTrack = region_ * channelAssignment_->numNodesDR();
    auto nonNullTrack = [](int& sum, const FrameTrack& frame) { return sum += (frame.first.isNonnull() ? 1 : 0); };
    auto nonNullStub = [](int& sum, const FrameStub& frame) { return sum += (frame.first.isNonnull() ? 1 : 0); };
    // count tracks and stubs and reserve corresponding vectors
    int sizeTracks(0);
    int sizeStubs(0);
    for (int channel = 0; channel < channelAssignment_->numNodesDR(); channel++) {
      const int streamTrackId = offsetTrack + channel;
      const int offsetStub = streamTrackId * setup_->numLayers();
      const StreamTrack& streamTrack = streamsTrack[streamTrackId];
      input_[channel].reserve(streamTrack.size());
      sizeTracks += accumulate(streamTrack.begin(), streamTrack.end(), 0, nonNullTrack);
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        const StreamStub& streamStub = streamsStub[offsetStub + layer];
        sizeStubs += accumulate(streamStub.begin(), streamStub.end(), 0, nonNullStub);
      }
    }
    tracks_.reserve(sizeTracks);
    stubs_.reserve(sizeStubs);
    // transform input data into handy structs
    for (int channel = 0; channel < channelAssignment_->numNodesDR(); channel++) {
      vector<Track*>& input = input_[channel];
      const int streamTrackId = offsetTrack + channel;
      const int offsetStub = streamTrackId * setup_->numLayers();
      const StreamTrack& streamTrack = streamsTrack[streamTrackId];
      for (int frame = 0; frame < (int)streamTrack.size(); frame++) {
        const FrameTrack& frameTrack = streamTrack[frame];
        if (frameTrack.first.isNull()) {
          input.push_back(nullptr);
          continue;
        }
        //convert track parameter
        const double r2Inv = digi(-ttTrackRef->rInv() / 2., baseUinv2R_);
        const double phi0U =
            digi(tt::deltaPhi(ttTrackRef->phi() - region_ * setup_->baseRegion() + setup_->hybridRangePhi() / 2.),
                 baseUphiT_);
        const double phi0S = digi(phi0U - setup_->hybridRangePhi() / 2., baseUphiT_);
        const double cot = digi(ttTrackRef->tanL(), baseUcot_);
        const double z0 = digi(ttTrackRef->z0(), baseUzT_);
        const double phiT = digi(phi0S + r2Inv * digi(setup_->chosenRofPhi(), baseUr_), baseUphiT_);
        const double zT = digi(z0 + cot * digi(setup_->chosenRofZ(), baseUr_), baseUzT_);
        // kill tracks outside of fiducial range
        if (abs(phiT) > setup_->baseRegion() / 2. || abs(zT) > setup_->hybridMaxCot() * setup_->chosenRofZ() ||
            abs(z0) > setup_->beamWindowZ()) {
          input.push_back(nullptr);
          continue;
        }
        // convert stubs
        vector<Stub*> stubs;
        stubs.reserve(setup_->numLayers());
        for (int layer = 0; layer < setup_->numLayers(); layer++) {
          const FrameStub& frameStub = streamsStub[offsetStub + layer][frame];
          if (frameStub.first.isNull())
            continue;
          // parse residuals from tt::Frame and take r and layerId from tt::TTStubRef
          const bool barrel = setup_->barrel(ttStubRef);
          const int layerIdTracklet = setup_->trackletLayerId(ttStubRef);
          const double basePhi = barrel ? settings_->kphi1() : settings_->kphi(layerIdTracklet);
          const double baseRZ = barrel ? settings_->kz(layerIdTracklet) : settings_->kz();
          const int widthRZ = barrel ? settings_->zresidbits() : settings_->rresidbits();
          TTBV hw(frameStub.second);
          const TTBV hwRZ(hw, widthRZ, 0, true);
          hw >>= widthRZ;
          const TTBV hwPhi(hw, settings_->phiresidbits(), 0, true);
          hw >>= settings_->phiresidbits();
          const double r = digi(setup_->stubR(hw, ttStubRef) - setup_->chosenRofPhi(), baseUr_);
          double phi = hwPhi.val(basePhi);
          if (basePhi > baseUphi_)
            phi += baseUphi_ / 2.;
          const double z = digi(hwRZ.val(baseRZ) * (barrel ? 1. : -cot), baseUz_);
          // determine module type
          bool psTilt;
          if (barrel) {
            const double posZ = (r + digi(setup_->chosenRofPhi(), baseUr_)) * cot + z0 + z;
            const int indexLayerId = setup_->indexLayerId(ttStubRef);
            const double limit = setup_->limitsTiltedZ(indexLayerId);
            psTilt = abs(posZ) < limit;
          } else
            psTilt = setup_->psModule(ttStubRef);
          if (useTTStubResiduals_) {
            const GlobalPoint gp = setup_->stubPos(ttStubRef);
            const double ttR = r;
            const double ttZ = gp.z() - (z0 + (ttR + setup_->chosenRofPhi()) * cot);
            stubs_.emplace_back(ttStubRef, layerId, ttR, phi, ttZ, psTilt);
          } else
            stubs_.emplace_back(ttStubRef, layerId, r, phi, z, psTilt);
          stubs.push_back(&stubs_.back());
        }
        // create fake seed stubs, since TrackBuilder doesn't output these stubs, required by the KF.
        for (int layerId : channelAssignment_->seedingLayers(channel)) {
          const vector<TTStubRef>& ttStubRefs = ttTrackRef->getStubRefs();
          auto sameLayer = [this, layerId](const TTStubRef& ttStubRef) {
            return setup_->layerId(ttStubRef) == layerId;
          };
          const TTStubRef& ttStubRef = *find_if(ttStubRefs.begin(), ttStubRefs.end(), sameLayer);
          const bool barrel = setup_->barrel(ttStubRef);
          double r;
          if (barrel)
            r = digi(setup_->hybridLayerR(layerId - setup_->offsetLayerId()) - setup_->chosenRofPhi(), baseUr_);
          else {
            r = (z0 +
                 digi(setup_->hybridDiskZ(layerId - setup_->offsetLayerId() - setup_->offsetLayerDisks()), baseUzT_)) *
                digi(1. / digi(abs(cot), baseCot), baseInvCot_);
            r = digi(r - digi(setup_->chosenRofPhi(), baseUr_), baseUr_);
          }
          static constexpr double phi = 0.;
          static constexpr double z = 0.;
          // determine module type
          bool psTilt;
          if (barrel) {
            const double posZ =
                digi(digi(setup_->hybridLayerR(layerId - setup_->offsetLayerId()), baseUr_) * cot + z0, baseUz_);
            const int indexLayerId = setup_->indexLayerId(ttStubRef);
            const double limit = digi(setup_->limitsTiltedZ(indexLayerId), baseUz_);
            psTilt = abs(posZ) < limit;
          } else
            psTilt = true;
          const GlobalPoint gp = setup_->stubPos(ttStubRef);
          const double ttR = gp.perp() - setup_->chosenRofPhi();
          const double ttZ = gp.z() - (z0 + (ttR + setup_->chosenRofPhi()) * cot);
          if (useTTStubResiduals_)
            stubs_.emplace_back(ttStubRef, layerId, ttR, phi, ttZ, psTilt);
          else
            stubs_.emplace_back(ttStubRef, layerId, r, phi, z, psTilt);
          stubs.push_back(&stubs_.back());
        }
        const bool valid = frame < setup_->numFrames() ? true : enableTruncation_;
        tracks_.emplace_back(ttTrackRef, valid, r2Inv, phiT, cot, zT, stubs);
        input.push_back(&tracks_.back());
      }
      // remove all gaps between end and last track
      for (auto it = input.end(); it != input.begin();)
        it = (*--it) ? input.begin() : input.erase(it);
    }
  }

  // fill output products
  void KFin::produce(StreamsStub& accpetedStubs,
                     StreamsTrack& acceptedTracks,
                     StreamsStub& lostStubs,
                     StreamsTrack& lostTracks) {
    // calculate stub uncertainties
    static constexpr int usedMSBpitchOverRaddr = 1;
    static const double baseRlut =
        dataFormats_->base(Variable::r, Process::kfin) *
        pow(2, dataFormats_->width(Variable::r, Process::zht) - setup_->widthAddrBRAM18() + usedMSBpitchOverRaddr);
    static const double baseRinvR = dataFormats_->base(Variable::r, Process::kfin) *
                                    pow(2, dataFormats_->width(Variable::r, Process::zht) - setup_->widthAddrBRAM18());
    static const double basePhi =
        dataFormats_->base(Variable::inv2R, Process::kfin) * dataFormats_->base(Variable::r, Process::kfin);
    static const double baseInvR =
        pow(2.,
            ceil(log2(dataFormats_->base(Variable::r, Process::kfin) / setup_->tbInnerRadius())) -
                setup_->widthDSPbu()) /
        dataFormats_->base(Variable::r, Process::kfin);
    static const double maxCot = sinh(setup_->maxEta()) + setup_->beamWindowZ() / setup_->chosenRofZ();
    static constexpr int usedMSBCotLutaddr = 3;
    static const double baseCotLut = pow(2., ceil(log2(maxCot)) - setup_->widthAddrBRAM18() + usedMSBCotLutaddr);
    // base transform into high precision TMTT format
    for (Track& track : tracks_) {
      track.inv2R_ = redigi(track.inv2R_, baseUinv2R_, baseHinv2R_, setup_->widthDSPbu());
      track.phiT_ = redigi(track.phiT_, baseUphiT_, baseHphiT_, setup_->widthDSPbu());
      track.cot_ = redigi(track.cot_, baseUcot_, baseHcot_, setup_->widthDSPbu());
      track.zT_ = redigi(track.zT_, baseUzT_, baseHzT_, setup_->widthDSPbu());
      for (Stub* stub : track.stubs_) {
        stub->r_ = redigi(stub->r_, baseUr_, baseHr_, setup_->widthDSPbu());
        stub->phi_ = redigi(stub->phi_, baseUphi_, baseHphi_, setup_->widthDSPbu());
        stub->z_ = redigi(stub->z_, baseUz_, baseHz_, setup_->widthDSPbu());
      }
    }
    // find sector
    for (Track& track : tracks_) {
      const int sectorPhi = track.phiT_ < 0. ? 0 : 1;
      track.phiT_ -= (sectorPhi - .5) * setup_->baseSector();
      const int sectorEta = floor(track.zT_ / dataFormats_->format(Variable::zT, Process::zht).range()) + setup_->numSectorsEta() / 2;
      if (sectorEta > setup_->numSectorsEta() - 1 || sectorEta < 0) {
        track.valid_ = false;
        continue;
      }
      const double sectorCot = (sectorEta - setup_->numSectorsEta() / 2 + .5) * dataFormats_->format(Variable::cot, Process::zht).range();
      const double sectorZT = setup_->chosenRofZ() * sectorCot;
      track.cot_ = track.cot_ - digi(sectorCot, baseHcot_);
      track.zT_ = track.zT_ - digi(sectorZT, baseHzT_);
      track.sector_ = sectorPhi * setup_->numSectorsEta() + sectorEta;
    }
    // base transform into TMTT format
    for (Track& track : tracks_) {
      if (!track.valid_)
        continue;
      // store track parameter shifts
      const double dinv2R = digi(track.inv2R_ - digi(track.inv2R_, baseLinv2R_), baseHinv2R_);
      const double dphiT = digi(track.phiT_ - digi(track.phiT_, baseLphiT_), baseHphiT_);
      const double dcot = digi(track.cot_ - digi(track.cot_, baseLcot_), baseHcot_);
      const double dzT = digi(track.zT_ - digi(track.zT_, baseLzT_), baseHzT_);
      // shift track parameter;
      track.inv2R_ = digi(track.inv2R_, baseLinv2R_);
      track.phiT_ = digi(track.phiT_, baseLphiT_);
      track.cot_ = digi(track.cot_, baseLcot_);
      track.zT_ = digi(track.zT_, baseLzT_);
      // range checks
      if (!dataFormats_->format(Variable::inv2R, Process::kfin).inRange(track.inv2R_, true))
        track.valid_ = false;
      if (!dataFormats_->format(Variable::phiT, Process::kfin).inRange(track.phiT_, true))
        track.valid_ = false;
      if (!dataFormats_->format(Variable::cot, Process::kfin).inRange(track.cot_, true))
        track.valid_ = false;
      if (!dataFormats_->format(Variable::zT, Process::kfin).inRange(track.zT_, true))
        track.valid_ = false;
      if (!track.valid_)
        continue;
      // adjust stub residuals by track parameter shifts
      for (Stub* stub : track.stubs_) {
        const double dphi = digi(dphiT + stub->r_ * dinv2R, baseHphi_);
        const double r = stub->r_ + digi(setup_->chosenRofPhi() - setup_->chosenRofZ(), baseHr_);
        const double dz = digi(dzT + r * dcot, baseHz_);
        stub->phi_ = digi(stub->phi_ + dphi, baseLphi_);
        stub->z_ = digi(stub->z_ + dz, baseLz_);
        // range checks
        if (!dataFormats_->format(Variable::phi, Process::kfin).inRange(stub->phi_))
          stub->valid_ = false;
        if (!dataFormats_->format(Variable::z, Process::kfin).inRange(stub->z_))
          stub->valid_ = false;
      }
    }
    // encode layer id
    for (Track& track : tracks_) {
      if (!track.valid_)
        continue;
      // lookup maybe layers
      track.maybe_ = layerEncoding_->maybePattern(track.zT_);
      // lookup layerEncoding
      const vector<int>& layerEncoding = layerEncoding_->layerEncoding(track.zT_);
      for (Stub* stub : track.stubs_) {
        if (!stub->valid_)
          continue;
        // replace layerId by encoded layerId
        const auto it = find(layerEncoding.begin(), layerEncoding.end(), stub->layer_);
        stub->layer_ = min((int)distance(layerEncoding.begin(), it), setup_->numLayers() - 1);
        // kill stubs from layers which can't be crossed by track
        if (it == layerEncoding.end())
          stub->valid_ = false;
      }
    }
    // kill tracks with not enough layer
    for (Track& track : tracks_) {
      if (!track.valid_)
        continue;
      TTBV hits(0, setup_->numLayers());
      for (const Stub* stub : track.stubs_)
        if (stub->valid_)
          hits.set(stub->layer_);
      if (hits.count() < setup_->kfMinLayers())
        track.valid_ = false;
    }
    // calculate stub uncertainties
    for (Track& track : tracks_) {
      if (!track.valid_)
        continue;
      const double zTSector = dataFormats_->format(Variable::cot, Process::gp).digi(track.zT_);
      const double cotSector = zTSector / setup_->chosenRofZ();
      const double inv2R = abs(track.inv2R_);
      for (Stub* stub : track.stubs_) {
        const bool barrel = setup_->barrel(stub->ttStubRef_);
        const bool ps = barrel ? setup_->psModule(stub->ttStubRef_) : stub->psTilt_;
        const bool tilt = barrel ? (ps && !stub->psTilt_) : false;
        const double length = ps ? setup_->pitchColPS() : setup_->pitchCol2S();
        const double pitch = ps ? setup_->pitchRowPS() : setup_->pitchRow2S();
        const double pitchOverR = digi(pitch / (digi(stub->r_, baseR) + setup_->chosenRofPhi()), basePhi);
        const double r = digi(stub->r_, baseRinvR) + setup_->chosenRofPhi();
        const double sumdz = track.zT_ + stub->z_;
        const double dZ = digi(sumdz - digi(setup_->chosenRofZ(), baseLr_) * track.cot_, baseLcot_ * baseLr_);
        const double sumcot = track.cot_ + digi(cotSector, baseHcot_);
        const double cot = digi(abs(dZ * digi(1. / r, baseInvR) + sumcot), baseCotLut);
        double lengthZ = length;
        double lengthR = 0.;
        if (!barrel) {
          lengthZ = length * cot;
          lengthR = length;
        } else if (tilt) {
          lengthZ = length * abs(setup_->tiltApproxSlope() * cot + setup_->tiltApproxIntercept());
          lengthR = setup_->tiltUncertaintyR();
        }
        const double scat = digi(setup_->scattering(), baseR);
        stub->dZ_ = lengthZ + baseZ;
        stub->dPhi_ = (scat + digi(lengthR, baseR)) * inv2R + pitchOverR;
        stub->dPhi_ = digi(stub->dPhi_, basePhi) + basePhi;
      }
    }
    // fill products StreamsStub& accpetedStubs, StreamsTrack& acceptedTracks, StreamsStub& lostStubs, StreamsTrack& lostTracks
    auto frameTrack = [this](Track* track) {
      const TTBV maybe(track->maybe_);
      const TTBV inv2R(dataFormats_->format(Variable::inv2R, Process::kfin).ttBV(track->inv2R_));
      const TTBV phiT(dataFormats_->format(Variable::phiT, Process::kfin).ttBV(track->phiT_));
      const TTBV cot(dataFormats_->format(Variable::cot, Process::kfin).ttBV(track->cot_));
      const TTBV zT(dataFormats_->format(Variable::zT, Process::kfin).ttBV(track->zT_));
      return FrameTrack(track->ttTrackRef_, "1" + maybe.str() + inv2R.str() + phiT.str() + cot.str() + zT.str());
    };
    auto frameStub = [this](Track* track, int layer) {
      auto equal = [layer](Stub* stub) { return stub->channel_ == layer; };
      const auto it = find_if(track->stubs_.begin(), track->stubs_.end(), equal);
      if (it == track->stubs_.end())
        return FrameStub();
      Stub* stub = *it;
      const TTBV r(dataFormats_->format(Variable::r, Process::kfin).ttBV(stub->r_));
      const TTBV phi(dataFormats_->format(Variable::phi, Process::kfin).ttBV(stub->phi_));
      const TTBV z(dataFormats_->format(Variable::z, Process::kfin).ttBV(stub->z_));
      const TTBV dPhi(dataFormats_->format(Variable::dPhi, Process::kfin).ttBV(stub->dPhi_));
      const TTBV dZ(dataFormats_->format(Variable::dZ, Process::kfin).ttBV(stub->dZ_));
      return FrameStub(stub->ttStubRef_, Frame("1" + r.str() + phi.str() + z.str() + dPhi.str() + dZ.str()));
    };
    // merge number of nodes DR to number of Nodes KF and store result
    static const int nMux = channelAssignment_->numNodesDR() / setup_->kfNumWorker();
    const int offsetTrack = region_ * setup_->kfNumWorker();
    for (int nodeKF = 0; nodeKF < setup_->kfNumWorker(); nodeKF++) {
      const int offset = nodeKF * nMux;
      deque<Track*> accepted;
      deque<Track*> lost;
      vector<deque<Track*>> stacks(nMux);
      vector<deque<Track*>> inputs(nMux);
      for (int channel = 0; channel < nMux; channel++) {
        const vector<Track*>& input = input_[offset + channel];
        inputs[channel] = deque<Track*>(input.begin(), input.end());
      }
      // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
      while (!all_of(inputs.begin(), inputs.end(), [](const deque<Track*>& tracks) { return tracks.empty(); }) or
             !all_of(stacks.begin(), stacks.end(), [](const deque<Track*>& tracks) { return tracks.empty(); })) {
        // fill input fifos
        for (int channel = 0; channel < nMux; channel++) {
          deque<Track*>& stack = stacks[channel];
          Track* track = pop_front(inputs[channel]);
          if (track)
            stack.push_back(track);
        }
        // merge input fifos to one stream, prioritizing higher input channel over lower channel
        bool nothingToRoute(true);
        for (int channel = nMux - 1; channel >= 0; channel--) {
          Track* track = pop_front(stacks[channel]);
          if (track) {
            nothingToRoute = false;
            accepted.push_back(track);
            break;
          }
        }
        if (nothingToRoute)
          accepted.push_back(nullptr);
      }
      // truncate if desired
      if (enableTruncation_ && (int)accepted.size() > setup_->numFrames()) {
        const auto limit = next(accepted.begin(), setup_->numFrames());
        copy_if(limit, accepted.end(), back_inserter(lost), [](const Track* track) { return track; });
        accepted.erase(limit, accepted.end());
      }
      // remove all gaps between end and last track
      for (auto it = accepted.end(); it != accepted.begin();)
        it = (*--it) ? accepted.begin() : accepted.erase(it);
      // fill products StreamsStub& accpetedStubs, StreamsTrack& acceptedTracks, StreamsStub& lostStubs, StreamsTrack& lostTracks
      const int channelTrack = offsetTrack + nodeKF;
      const int offsetStub = channelTrack * setup_->numLayers();
      // fill lost tracks and stubs without gaps
      lostTracks[channelTrack].reserve(lost.size());
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        lostStubs[offsetStub + layer].reserve(lost.size());
      for (Track* track : lost) {
        lostTracks[channelTrack].emplace_back(frameTrack(track));
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          lostStubs[offsetStub + layer].emplace_back(frameStub(track, layer));
      }
      // fill accepted tracks and stubs with gaps
      acceptedTracks[channelTrack].reserve(accepted.size());
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        accpetedStubs[offsetStub + layer].reserve(accepted.size());
      for (Track* track : accepted) {
        if (!track) {  // fill gap
          acceptedTracks[channelTrack].emplace_back(FrameTrack());
          for (int layer = 0; layer < setup_->numLayers(); layer++)
            accpetedStubs[offsetStub + layer].emplace_back(FrameStub());
          continue;
        }
        acceptedTracks[channelTrack].emplace_back(frameTrack(track));
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          accpetedStubs[offsetStub + layer].emplace_back(frameStub(track, layer));
      }
    }
  }

  // remove and return first element of deque, returns nullptr if empty
  template <class T>
  T* KFin::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

}  // namespace trklet
