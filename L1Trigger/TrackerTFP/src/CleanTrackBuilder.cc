#include "L1Trigger/TrackerTFP/interface/CleanTrackBuilder.h"

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

  CleanTrackBuilder::CleanTrackBuilder(const ParameterSet& iConfig,
                                       const Setup* setup,
                                       const DataFormats* dataFormats,
                                       vector<StubCTB>& stubs,
                                       vector<TrackCTB>& tracks)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        stubsCTB_(stubs),
        tracksCTB_(tracks) {
    stubs_.reserve(stubs.capacity());
    tracks_.reserve(tracks.capacity());
  }

  // fill output products
  void CleanTrackBuilder::produce(const vector<vector<StubHT*>>& streamsIn,
                                  vector<deque<TrackCTB*>>& regionTracks,
                                  vector<vector<deque<StubCTB*>>>& regionStubs) {
    static int region = 0;
    static const int numChannelIn = dataFormats_->numChannel(Process::ht);
    static const int numChannelOut = dataFormats_->numChannel(Process::ctb);
    static const int numChannel = numChannelIn / numChannelOut;
    static const int numLayers = setup_->numLayers();
    // loop over worker
    for (int channelOut = 0; channelOut < numChannelOut; channelOut++) {
      // clean input tracks
      const int offsetIn = channelOut * numChannel;
      vector<deque<Track*>> streamsT(numChannel);
      vector<deque<Stub*>> streamsS(numChannel);
      for (int cin = 0; cin < numChannel; cin++)
        cleanStream(streamsIn[cin + offsetIn], streamsT[cin], streamsS[cin], offsetIn + cin);
      // route
      deque<Track*> tracks;
      vector<deque<Stub*>> stubs(numLayers);
      route(streamsT, tracks);
      route(streamsS, stubs);
      // sort
      sort(tracks, stubs);
      // convert
      deque<TrackCTB*>& channelTracks = regionTracks[channelOut];
      vector<deque<StubCTB*>>& channelStubs = regionStubs[channelOut];
      convert(tracks, stubs, channelTracks, channelStubs);
    }
    region++;
  }

  //
  void CleanTrackBuilder::cleanStream(const vector<StubHT*>& input,
                                      deque<Track*>& tracks,
                                      deque<Stub*>& stubs,
                                      int channelId) {
    static const DataFormat dfInv2R = dataFormats_->format(Variable::inv2R, Process::ht);
    static const int chan = setup_->kfNumWorker();
    static const int mux = setup_->htNumBinsInv2R() / chan;
    const int inv2Ru = chan * (channelId % mux) + channelId / mux;
    const double inv2R = dfInv2R.floating(dfInv2R.toSigned(inv2Ru));
    const int offset = (channelId % mux) * setup_->ctbMaxTracks();
    int trackId = offset;
    // identify tracks in input container
    int id;
    auto toTrkId = [this](StubHT* stub) {
      static const DataFormat& phiT = dataFormats_->format(Variable::phiT, Process::ht);
      static const DataFormat& zT = dataFormats_->format(Variable::zT, Process::ht);
      return (phiT.ttBV(stub->phiT()) + zT.ttBV(stub->zT())).val();
    };
    auto different = [&id, toTrkId](StubHT* stub) { return id != toTrkId(stub); };
    int delta = -setup_->htMinLayers();
    int old = 0;
    for (auto it = input.begin(); it != input.end();) {
      id = toTrkId(*it);
      const auto start = it;
      const auto end = find_if(start, input.end(), different);
      const vector<StubHT*> track(start, end);
      // restore clock accurancy
      delta += (int)track.size() - old;
      old = track.size();
      if (delta > 0) {
        stubs.insert(stubs.end(), delta, nullptr);
        tracks.insert(tracks.end(), delta, nullptr);
        delta = 0;
      }
      // run single track through r-phi and r-z hough transform
      cleanTrack(track, tracks, stubs, inv2R, trackId++);
      if (trackId - offset == setup_->ctbMaxTracks())
        break;
      // set begin of next track
      it = end;
    }
  }

  // run single track through r-phi and r-z hough transform
  void CleanTrackBuilder::cleanTrack(
      const vector<StubHT*>& track, deque<Track*>& tracks, deque<Stub*>& stubs, double inv2R, int trackId) {
    static const int numBinsInv2R = setup_->ctbNumBinsInv2R();
    static const int numBinsPhiT = setup_->ctbNumBinsPhiT();
    static const int numBinsZ0 = setup_->ctbNumBinsZ0();
    static const int numBinsZT = setup_->ctbNumBinsZT();
    static const double baseInv2R = dataFormats_->base(Variable::inv2R, Process::ht) / numBinsInv2R;
    static const double basePhiT = dataFormats_->base(Variable::phiT, Process::ht) / numBinsPhiT;
    //static const double baseZ0 = dataFormats_->range(Variable::zT, Process::dr) / numBinsZ0;
    static const double baseZ0 = 2. * setup_->beamWindowZ() / numBinsZ0;
    static const double baseZT = dataFormats_->base(Variable::zT, Process::ht) / numBinsZT;
    auto toLayerId = [](StubHT* stub) { return stub->layer().val(3); };
    auto toDPhi = [this, inv2R](StubHT* stub) {
      static const DataFormat df = dataFormats_->format(Variable::dPhi, Process::ctb);
      const bool barrel = stub->layer()[5];
      const bool ps = stub->layer()[4];
      const bool tilt = stub->layer()[3];
      const double pitchRow = ps ? setup_->pitchRowPS() : setup_->pitchRow2S();
      const double pitchCol = ps ? setup_->pitchColPS() : setup_->pitchCol2S();
      const double pitchColR = barrel ? (tilt ? setup_->tiltUncertaintyR() : 0.0) : pitchCol;
      const double r = stub->r() + setup_->chosenRofPhi();
      const double dPhi = pitchRow / r + (setup_->scattering() + pitchColR) * abs(inv2R) + df.base();
      return df.digi(dPhi);
    };
    auto toDZ = [this](StubHT* stub) {
      static const DataFormat dfZT = dataFormats_->format(Variable::zT, Process::ht);
      static const DataFormat dfDZ = dataFormats_->format(Variable::dZ, Process::ctb);
      static const double m = setup_->tiltApproxSlope();
      static const double c = setup_->tiltApproxIntercept();
      const bool barrel = stub->layer()[5];
      const bool ps = stub->layer()[4];
      const bool tilt = stub->layer()[3];
      const double pitchCol = ps ? setup_->pitchColPS() : setup_->pitchCol2S();
      const double zT = dfZT.floating(stub->zT());
      const double cot = abs(zT) / setup_->chosenRofZ();
      const double sigmaZ = barrel ? (tilt ? m * cot + c : pitchCol) : cot * pitchCol;
      const double dZ = sigmaZ + dfDZ.base();
      return dfDZ.digi(dZ);
    };
    vector<Stub*> tStubs;
    tStubs.reserve(track.size());
    vector<TTBV> hitPatternPhi(numBinsInv2R * numBinsPhiT, TTBV(0, setup_->numLayers()));
    vector<TTBV> hitPatternZ(numBinsZ0 * numBinsZT, TTBV(0, setup_->numLayers()));
    TTBV tracksPhi(0, numBinsInv2R * numBinsPhiT);
    TTBV tracksZ(0, numBinsZ0 * numBinsZT);
    // identify finer tracks each stub is consistent with
    for (StubHT* stub : track) {
      const int layerId = toLayerId(stub);
      const double dPhi = toDPhi(stub);
      const double dZ = toDZ(stub);
      // r - phi HT
      auto phiT = [stub](double inv2R, double dPhi) { return inv2R * stub->r() + stub->phi() + dPhi / 2.; };
      TTBV hitsPhi(0, numBinsInv2R * numBinsPhiT);
      for (int binInv2R = 0; binInv2R < numBinsInv2R; binInv2R++) {
        const int offset = binInv2R * numBinsPhiT;
        const double inv2RMin = (binInv2R - numBinsInv2R / 2.) * baseInv2R;
        const double inv2RMax = inv2RMin + baseInv2R;
        const auto phiTs = {phiT(inv2RMin, dPhi), phiT(inv2RMax, dPhi), phiT(inv2RMin, -dPhi), phiT(inv2RMax, -dPhi)};
        const int binPhiTMin = floor(*min_element(phiTs.begin(), phiTs.end()) / basePhiT + 1.e-12) + numBinsPhiT / 2;
        const int binPhiTMax = floor(*max_element(phiTs.begin(), phiTs.end()) / basePhiT + 1.e-12) + numBinsPhiT / 2;
        for (int binPhiT = 0; binPhiT < numBinsPhiT; binPhiT++)
          if (binPhiT >= binPhiTMin && binPhiT <= binPhiTMax)
            hitsPhi.set(offset + binPhiT);
      }
      // check for tracks on the fly
      for (int phi : hitsPhi.ids()) {
        hitPatternPhi[phi].set(layerId);
        if (hitPatternPhi[phi].count() >= setup_->ctbMinLayers())
          tracksPhi.set(phi);
      }
      // r - z HT
      auto z = [this, stub](double zT, double z0) {
        return z0 + (stub->r() + setup_->chosenRofPhi()) / setup_->chosenRofZ() * (zT - z0);
      };
      static const DataFormat& dfZ = dataFormats_->format(Variable::z, Process::ctb);
      TTBV hitsZ(0, numBinsZ0 * numBinsZT);
      const double z0 = dfZ.digi(stub->z() - dZ / 2.);
      const double z1 = dfZ.digi(stub->z() + dZ / 2.);
      for (int binZ0 = 0; binZ0 < numBinsZ0; binZ0++) {
        const int offset = binZ0 * numBinsZT;
        const double z0Min = (binZ0 - numBinsZ0 / 2.) * baseZ0;
        const double z0Max = z0Min + baseZ0;
        for (int binZT = 0; binZT < numBinsZT; binZT++) {
          const double zTMin = (binZT - numBinsZT / 2.) * baseZT;
          const double zTMax = zTMin + baseZT;
          const auto zs = {z(zTMin, z0Min), z(zTMax, z0Min), z(zTMin, z0Max), z(zTMax, z0Max)};
          const double zMin = dfZ.digi(*min_element(zs.begin(), zs.end()));
          const double zMax = dfZ.digi(*max_element(zs.begin(), zs.end()));
          if (z1 >= zMin && z0 <= zMax)
            hitsZ.set(offset + binZT);
        }
      }
      // check for tracks on the fly
      for (int z : hitsZ.ids()) {
        hitPatternZ[z].set(layerId);
        if (hitPatternZ[z].count() >= setup_->ctbMinLayers())
          tracksZ.set(z);
      }
      // store stubs consistent finer tracks
      stubs_.emplace_back(stub, trackId, hitsPhi, hitsZ, layerId, dPhi, dZ);
      tStubs.push_back(&stubs_.back());
    }
    // clean
    tracks.insert(tracks.end(), tStubs.size() - 1, nullptr);
    tracks_.emplace_back(setup_, trackId, tracksPhi, tracksZ, tStubs, inv2R);
    tracks.push_back(&tracks_.back());
    stubs.insert(stubs.end(), tStubs.begin(), tStubs.end());
  }

  //
  void CleanTrackBuilder::Stub::update(const TTBV& phi, const TTBV& z, vector<int>& ids, int max) {
    auto consistent = [](const TTBV& stub, const TTBV& track) {
      for (int id : track.ids())
        if (stub[id])
          return true;
      return false;
    };
    if (consistent(hitsPhi_, phi) && consistent(hitsZ_, z) && ids[layerId_] < max)
      stubId_ = ids[layerId_]++;
    else
      valid_ = false;
  }

  // construct Track from Stubs
  CleanTrackBuilder::Track::Track(const Setup* setup,
                                  int trackId,
                                  const TTBV& hitsPhi,
                                  const TTBV& hitsZ,
                                  const std::vector<Stub*>& stubs,
                                  double inv2R)
      : valid_(true), stubs_(stubs), trackId_(trackId), hitsPhi_(hitsPhi), hitsZ_(hitsZ), inv2R_(inv2R) {
    vector<int> stubIds(setup->numLayers(), 0);
    for (Stub* stub : stubs_)
      stub->update(hitsPhi_, hitsZ_, stubIds, setup->ctbMaxStubs());
    const int nLayer =
        accumulate(stubIds.begin(), stubIds.end(), 0, [](int& sum, int i) { return sum += (i > 0 ? 1 : 0); });
    if (nLayer < setup->ctbMinLayers())
      valid_ = false;
    size_ = *max_element(stubIds.begin(), stubIds.end());
  }

  //
  void CleanTrackBuilder::route(vector<deque<Stub*>>& input, vector<deque<Stub*>>& outputs) const {
    for (int channelOut = 0; channelOut < (int)outputs.size(); channelOut++) {
      deque<Stub*>& output = outputs[channelOut];
      vector<deque<Stub*>> inputs(input);
      for (deque<Stub*>& stream : inputs) {
        for (Stub*& stub : stream)
          if (stub && (!stub->valid_ || stub->layerId_ != channelOut))
            stub = nullptr;
        for (auto it = stream.end(); it != stream.begin();)
          it = (*--it) ? stream.begin() : stream.erase(it);
      }
      vector<deque<Stub*>> stacks(input.size());
      // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
      while (!all_of(inputs.begin(), inputs.end(), [](const deque<Stub*>& stubs) { return stubs.empty(); }) or
             !all_of(stacks.begin(), stacks.end(), [](const deque<Stub*>& stubs) { return stubs.empty(); })) {
        // fill input fifos
        for (int channel = 0; channel < (int)input.size(); channel++) {
          deque<Stub*>& stack = stacks[channel];
          Stub* stub = pop_front(inputs[channel]);
          if (stub) {
            if (enableTruncation_ && (int)stack.size() == setup_->ctbDepthMemory() - 1)
              pop_front(stack);
            stack.push_back(stub);
          }
        }
        // merge input fifos to one stream, prioritizing higher input channel over lower channel
        bool nothingToRoute(true);
        //for (int channel = (int)input.size() - 1; channel >= 0; channel--) {
        for (int channel = 0; channel < (int)input.size(); channel++) {
          Stub* stub = pop_front(stacks[channel]);
          if (stub) {
            nothingToRoute = false;
            output.push_back(stub);
            break;
          }
        }
        if (nothingToRoute)
          output.push_back(nullptr);
      }
    }
  }

  //
  void CleanTrackBuilder::route(vector<deque<Track*>>& inputs, deque<Track*>& output) const {
    vector<deque<Track*>> stacks(inputs.size());
    // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
    while (!all_of(inputs.begin(), inputs.end(), [](const deque<Track*>& tracks) { return tracks.empty(); }) or
           !all_of(stacks.begin(), stacks.end(), [](const deque<Track*>& tracks) { return tracks.empty(); })) {
      // fill input fifos
      for (int channel = 0; channel < (int)inputs.size(); channel++) {
        deque<Track*>& stack = stacks[channel];
        Track* track = pop_front(inputs[channel]);
        if (track && track->valid_) {
          if (enableTruncation_ && (int)stack.size() == setup_->ctbDepthMemory() - 1)
            pop_front(stack);
          stack.push_back(track);
        }
      }
      // merge input fifos to one stream, prioritizing higher input channel over lower channel
      bool nothingToRoute(true);
      //for (int channel = (int)inputs.size() - 1; channel >= 0; channel--) {
      for (int channel = 0; channel < (int)inputs.size(); channel++) {
        Track* track = pop_front(stacks[channel]);
        if (track) {
          nothingToRoute = false;
          output.push_back(track);
          break;
        }
      }
      if (nothingToRoute)
        output.push_back(nullptr);
    }
  }

  // sort
  void CleanTrackBuilder::sort(deque<Track*>& tracks, vector<deque<Stub*>>& stubs) const {
    // aplly truncation
    if (enableTruncation_) {
      if ((int)tracks.size() > setup_->numFrames())
        tracks.resize(setup_->numFrames());
      for (deque<Stub*>& stream : stubs)
        if ((int)stream.size() > setup_->numFrames())
          stream.resize(setup_->numFrames());
    }
    // cycle event, remove all gaps
    tracks.erase(remove(tracks.begin(), tracks.end(), nullptr), tracks.end());
    for (deque<Stub*>& stream : stubs)
      stream.erase(remove(stream.begin(), stream.end(), nullptr), stream.end());
    // prepare sort according to track id arrival order
    vector<int> trackIds;
    trackIds.reserve(tracks.size());
    transform(tracks.begin(), tracks.end(), back_inserter(trackIds), [](Track* track) { return track->trackId_; });
    auto cleaned = [&trackIds](Stub* stub) {
      return find(trackIds.begin(), trackIds.end(), stub->trackId_) == trackIds.end();
    };
    auto order = [&trackIds](auto lhs, auto rhs) {
      const auto l = find(trackIds.begin(), trackIds.end(), lhs->trackId_);
      const auto r = find(trackIds.begin(), trackIds.end(), rhs->trackId_);
      return distance(r, l) < 0;
    };
    for (deque<Stub*>& stream : stubs) {
      // remove stubs from removed tracks
      stream.erase(remove_if(stream.begin(), stream.end(), cleaned), stream.end());
      // sort according to stub id on layer
      stable_sort(stream.begin(), stream.end(), [](Stub* lhs, Stub* rhs) { return lhs->stubId_ < rhs->stubId_; });
      // sort according to track id arrival order
      stable_sort(stream.begin(), stream.end(), order);
    }
    // add all gaps
    const int size =
        accumulate(tracks.begin(), tracks.end(), 0, [](int& sum, Track* track) { return sum += track->size_; });
    for (int frame = 0; frame < size;) {
      const int trackId = tracks[frame]->trackId_;
      const int length = tracks[frame]->size_;
      tracks.insert(next(tracks.begin(), frame + 1), length - 1, nullptr);
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        deque<Stub*>& stream = stubs[layer];
        if (frame >= (int)stream.size()) {
          stream.insert(stream.end(), length, nullptr);
          continue;
        }
        const auto begin = next(stream.begin(), frame);
        const auto end = find_if(begin, stream.end(), [trackId](Stub* stub) { return stub->trackId_ != trackId; });
        stream.insert(end, length - distance(begin, end), nullptr);
      }
      frame += length;
    }
  }

  //
  void CleanTrackBuilder::convert(const deque<Track*>& iTracks,
                                  const vector<deque<Stub*>>& iStubs,
                                  deque<TrackCTB*>& oTracks,
                                  vector<deque<StubCTB*>>& oStubs) {
    for (int iFrame = 0; iFrame < (int)iTracks.size();) {
      Track* track = iTracks[iFrame];
      if (!track) {
        oTracks.push_back(nullptr);
        for (deque<StubCTB*>& stubs : oStubs)
          stubs.push_back(nullptr);
        iFrame++;
        continue;
      }
      StubHT* s = nullptr;
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        for (int n = 0; n < track->size_; n++) {
          Stub* stub = iStubs[layer][iFrame + n];
          if (!stub) {
            oStubs[layer].push_back(nullptr);
            continue;
          }
          s = stub->stubHT_;
          const double r = s->r();
          const double phi = s->phi();
          const double z = s->z();
          const double dPhi = stub->dPhi_;
          const double dZ = stub->dZ_;
          stubsCTB_.emplace_back(*s, r, phi, z, dPhi, dZ);
          oStubs[layer].push_back(&stubsCTB_.back());
        }
      }
      const double inv2R = track->inv2R_;
      const double phiT = dataFormats_->format(Variable::phiT, Process::ctb).floating(s->phiT());
      const double zT = dataFormats_->format(Variable::zT, Process::ctb).floating(s->zT());
      tracksCTB_.emplace_back(TTTrackRef(), dataFormats_, inv2R, phiT, zT);
      oTracks.push_back(&tracksCTB_.back());
      oTracks.insert(oTracks.end(), track->size_ - 1, nullptr);
      iFrame += track->size_;
    }
  }

  // remove and return first element of deque, returns nullptr if empty
  template <class T>
  T* CleanTrackBuilder::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

  void CleanTrackBuilder::put(TrackCTB* track,
                              const vector<vector<StubCTB*>>& stubs,
                              int region,
                              TTTracks& ttTracks) const {
    static const double dPhi = dataFormats_->format(Variable::phiT, Process::ctb).range();
    const double invR = -track->inv2R() * 2.;
    const double phi0 = deltaPhi(track->phiT() - track->inv2R() * setup_->chosenRofPhi() + region * dPhi);
    const double zT = track->zT();
    const double cot = zT / setup_->chosenRofZ();
    TTBV hits(0, setup_->numLayers());
    double chi2phi(0.);
    double chi2z(0.);
    const int nStubs = accumulate(
        stubs.begin(), stubs.end(), 0, [](int& sum, const vector<StubCTB*>& layer) { return sum += layer.size(); });
    vector<TTStubRef> ttStubRefs;
    ttStubRefs.reserve(nStubs);
    for (int layer = 0; layer < setup_->numLayers(); layer++) {
      for (StubCTB* stub : stubs[layer]) {
        hits.set(layer);
        chi2phi += pow(stub->phi(), 2) / pow(stub->dPhi(), 2);
        chi2z += pow(stub->z(), 2) / pow(stub->dZ(), 2);
        ttStubRefs.push_back(stub->frame().first);
      }
    }
    static constexpr int nPar = 4;
    static constexpr double d0 = 0.;
    static constexpr double z0 = 0;
    static constexpr double trkMVA1 = 0.;
    static constexpr double trkMVA2 = 0.;
    static constexpr double trkMVA3 = 0.;
    const int hitPattern = hits.val();
    const double bField = setup_->bField();
    TTTrack<Ref_Phase2TrackerDigi_> ttTrack(
        invR, phi0, cot, z0, d0, chi2phi, chi2z, trkMVA1, trkMVA2, trkMVA3, hitPattern, nPar, bField);
    ttTrack.setStubRefs(ttStubRefs);
    ttTrack.setPhiSector(region);
    ttTracks.emplace_back(ttTrack);
  }

}  // namespace trackerTFP
