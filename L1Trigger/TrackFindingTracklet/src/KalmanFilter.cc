#include "L1Trigger/TrackFindingTracklet/interface/KalmanFilter.h"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <deque>
#include <vector>
#include <set>
#include <cmath>

using namespace std;
using namespace edm;
using namespace tt;
using namespace trackerTFP;

namespace trklet {

  KalmanFilter::KalmanFilter(const ParameterSet& iConfig,
                             const Setup* setup,
                             const DataFormats* dataFormats,
                             const LayerEncoding* layerEncoding,
                             KalmanFilterFormats* kalmanFilterFormats,
                             int region,
                             TTTracks& ttTracks)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        use5ParameterFit_(iConfig.getParameter<bool>("Use5ParameterFit")),
        setup_(setup),
        dataFormats_(dataFormats),
        layerEncoding_(layerEncoding),
        kalmanFilterFormats_(kalmanFilterFormats),
        region_(region),
        ttTracks_(ttTracks),
        layer_(0),
        x0_(&kalmanFilterFormats_->format(VariableKF::x0)),
        x1_(&kalmanFilterFormats_->format(VariableKF::x1)),
        x2_(&kalmanFilterFormats_->format(VariableKF::x2)),
        x3_(&kalmanFilterFormats_->format(VariableKF::x3)),
        H00_(&kalmanFilterFormats_->format(VariableKF::H00)),
        H12_(&kalmanFilterFormats_->format(VariableKF::H12)),
        m0_(&kalmanFilterFormats_->format(VariableKF::m0)),
        m1_(&kalmanFilterFormats_->format(VariableKF::m1)),
        v0_(&kalmanFilterFormats_->format(VariableKF::v0)),
        v1_(&kalmanFilterFormats_->format(VariableKF::v1)),
        r0_(&kalmanFilterFormats_->format(VariableKF::r0)),
        r1_(&kalmanFilterFormats_->format(VariableKF::r1)),
        S00_(&kalmanFilterFormats_->format(VariableKF::S00)),
        S01_(&kalmanFilterFormats_->format(VariableKF::S01)),
        S12_(&kalmanFilterFormats_->format(VariableKF::S12)),
        S13_(&kalmanFilterFormats_->format(VariableKF::S13)),
        K00_(&kalmanFilterFormats_->format(VariableKF::K00)),
        K10_(&kalmanFilterFormats_->format(VariableKF::K10)),
        K21_(&kalmanFilterFormats_->format(VariableKF::K21)),
        K31_(&kalmanFilterFormats_->format(VariableKF::K31)),
        R00_(&kalmanFilterFormats_->format(VariableKF::R00)),
        R11_(&kalmanFilterFormats_->format(VariableKF::R11)),
        R00Rough_(&kalmanFilterFormats_->format(VariableKF::R00Rough)),
        R11Rough_(&kalmanFilterFormats_->format(VariableKF::R11Rough)),
        invR00Approx_(&kalmanFilterFormats_->format(VariableKF::invR00Approx)),
        invR11Approx_(&kalmanFilterFormats_->format(VariableKF::invR11Approx)),
        invR00Cor_(&kalmanFilterFormats_->format(VariableKF::invR00Cor)),
        invR11Cor_(&kalmanFilterFormats_->format(VariableKF::invR11Cor)),
        invR00_(&kalmanFilterFormats_->format(VariableKF::invR00)),
        invR11_(&kalmanFilterFormats_->format(VariableKF::invR11)),
        C00_(&kalmanFilterFormats_->format(VariableKF::C00)),
        C01_(&kalmanFilterFormats_->format(VariableKF::C01)),
        C11_(&kalmanFilterFormats_->format(VariableKF::C11)),
        C22_(&kalmanFilterFormats_->format(VariableKF::C22)),
        C23_(&kalmanFilterFormats_->format(VariableKF::C23)),
        C33_(&kalmanFilterFormats_->format(VariableKF::C33)) {}

  // read in and organize input tracks and stubs
  void KalmanFilter::consume(const StreamsTrack& streamsTrack, const StreamsStub& streamsStub) {
    static const int numLayers = setup_->numLayers();
    const int offset = region_ * numLayers;
    const StreamTrack& streamTrack = streamsTrack[region_];
    const int numTracks = accumulate(streamTrack.begin(), streamTrack.end(), 0, [](int& sum, const FrameTrack& f) {
      return sum += (f.first.isNull() ? 0 : 1);
    });
    int numStubs(0);
    for (int layer = 0; layer < numLayers; layer++) {
      const StreamStub& streamStub = streamsStub[offset + layer];
      numStubs += accumulate(streamStub.begin(), streamStub.end(), 0, [](int& sum, const FrameStub& f) {
        return sum += (f.first.isNull() ? 0 : 1);
      });
    }
    tracks_.reserve(numTracks);
    stubs_.reserve(numStubs);
    int trackId(0);
    for (int frame = 0; frame < (int)streamTrack.size(); frame++) {
      const FrameTrack& frameTrack = streamTrack[frame];
      if (frameTrack.first.isNull()) {
        stream_.push_back(nullptr);
        continue;
      }
      tracks_.emplace_back(frameTrack, dataFormats_);
      TrackCTB* track = &tracks_.back();
      vector<StubCTB*> stubs(numLayers, nullptr);
      for (int layer = 0; layer < numLayers; layer++) {
        const FrameStub& frameStub = streamsStub[offset + layer][frame];
        if (frameStub.first.isNull())
          continue;
        stubs_.emplace_back(frameStub, dataFormats_);
        stubs[layer] = &stubs_.back();
      }
      const TTBV& maybePattern = layerEncoding_->maybePattern(track->zT());
      states_.emplace_back(kalmanFilterFormats_, track, stubs, maybePattern, trackId++);
      stream_.push_back(&states_.back());
    }
  }

  // fill output products
  void KalmanFilter::produce(StreamsStub& streamsStub,
                             StreamsTrack& streamsTrack,
                             int& numAcceptedStates,
                             int& numLostStates) {
    // 5 parameter fit simulation
    if (use5ParameterFit_) {
      // Propagate state to each layer in turn, updating it with all viable stub combinations there, using KF maths
      for (layer_ = 0; layer_ < setup_->numLayers(); layer_++)
        addLayer();
    } else {  // 4 parameter fit emulation
      // seed building
      for (layer_ = 0; layer_ < setup_->kfMaxSeedingLayer(); layer_++)
        addSeedLayer();
      // calulcate seed parameter
      calcSeeds();
      // Propagate state to each layer in turn, updating it with all viable stub combinations there, using KF maths
      for (layer_ = setup_->kfNumSeedStubs(); layer_ < setup_->numLayers(); layer_++)
        addLayer();
    }
    // apply final cuts
    finalize();
    // count total number of final states
    const int nStates =
        accumulate(stream_.begin(), stream_.end(), 0, [](int& sum, State* state) { return sum += (state ? 1 : 0); });
    // apply truncation
    if (enableTruncation_ && (int)stream_.size() > setup_->numFramesHigh())
      stream_.resize(setup_->numFramesHigh());
    // cycle event, remove gaps
    stream_.erase(remove(stream_.begin(), stream_.end(), nullptr), stream_.end());
    // store number of states which got taken into account
    numAcceptedStates += (int)stream_.size();
    // store number of states which got not taken into account due to truncation
    numLostStates += nStates - (int)stream_.size();
    // best track per candidate selection
    accumulator();
    // Transform States into output products
    conv(streamsStub, streamsTrack);
  }

  // calculates the helix params & their cov. matrix from a pair of stubs
  void KalmanFilter::calcSeeds() {
    auto update = [this](State* s) {
      m0_->updateRangeActual(s->m0());
      m1_->updateRangeActual(s->m1());
      v0_->updateRangeActual(s->v0());
      v1_->updateRangeActual(s->v1());
      H00_->updateRangeActual(s->H00());
      H12_->updateRangeActual(s->H12());
    };
    for (State*& state : stream_) {
      if (!state)
        continue;
      State* s1 = state->parent();
      State* s0 = s1->parent();
      update(s0);
      update(s1);
      static const double rangeInvdH = 1. / setup_->kfMinSeedDeltaR();
      static const double rangeInvdH2 = rangeInvdH * rangeInvdH;
      static const int widthInvdH = setup_->widthDSPbu();
      static const int widthInvdH2 = setup_->widthDSPbu();
      static const int widthHv = setup_->widthDSPab();
      static const int widthH2v = setup_->widthDSPau();
      static const double baseH = H00_->base();
      static const int baseDiffInvdH = ceil(log2(rangeInvdH * pow(2., -widthInvdH) * baseH));
      static const int baseDiffInvdH2 = ceil(log2(rangeInvdH2 * pow(2., -widthInvdH2) * baseH * baseH));
      static const int baseDiffHv0 = 1 + v0_->width() + H00_->width() - widthHv;
      static const int baseDiffHv1 = 1 + v1_->width() + H12_->width() - widthHv;
      static const int baseDiffH2v0 = 1 + v0_->width() - 1 + 2 * H00_->width() - widthH2v;
      static const int baseDiffH2v1 = 1 + v1_->width() - 1 + 2 * H12_->width() - widthH2v;
      static const double baseH2 = baseH * baseH;
      static const double baseInvdH = pow(2., baseDiffInvdH) / baseH;
      static const double baseInvdH2 = pow(2., baseDiffInvdH2) / baseH2;
      static const double baseHm0 = baseH * m0_->base();
      static const double baseHm1 = baseH * m1_->base();
      static const double baseHv0 = baseH * v0_->base() * pow(2, baseDiffHv0);
      static const double baseHv1 = baseH * v1_->base() * pow(2, baseDiffHv1);
      static const double baseH2v0 = baseH2 * v0_->base() * pow(2, baseDiffH2v0);
      static const double baseH2v1 = baseH2 * v1_->base() * pow(2, baseDiffH2v1);
      const double dH = floor((s1->H00() - s0->H00() + 1.e-11) / baseH) * baseH;
      const double invdH = digi(1.0 / dH, baseInvdH);
      const double invdH2 = digi(1.0 / dH / dH, baseInvdH2);
      const double H02 = digi(s0->H00() * s0->H00(), baseH2);
      const double H12 = digi(s1->H00() * s1->H00(), baseH2);
      const double H22 = digi(s0->H12() * s0->H12(), baseH2);
      const double H32 = digi(s1->H12() * s1->H12(), baseH2);
      const double H1m0 = digi(s1->H00() * s0->m0(), baseHm0);
      const double H0m1 = digi(s0->H00() * s1->m0(), baseHm0);
      const double H3m2 = digi(s1->H12() * s0->m1(), baseHm1);
      const double H2m3 = digi(s0->H12() * s1->m1(), baseHm1);
      const double H1v0 = digi(s1->H00() * s0->v0(), baseHv0);
      const double H0v1 = digi(s0->H00() * s1->v0(), baseHv0);
      const double H3v2 = digi(s1->H12() * s0->v1(), baseHv1);
      const double H2v3 = digi(s0->H12() * s1->v1(), baseHv1);
      const double H12v0 = digi(H12 * s0->v0(), baseH2v0);
      const double H02v1 = digi(H02 * s1->v0(), baseH2v0);
      const double H32v2 = digi(H32 * s0->v1(), baseH2v1);
      const double H22v3 = digi(H22 * s1->v1(), baseH2v1);
      const double x0 = x0_->digi((s1->m0() - s0->m0()) * invdH);
      const double x2 = x2_->digi((s1->m1() - s0->m1()) * invdH);
      const double x1 = x1_->digi((H1m0 - H0m1) * invdH);
      const double x3 = x3_->digi((H3m2 - H2m3) * invdH);
      const double C00 = C00_->digi((s1->v0() + s0->v0()) * invdH2);
      const double C22 = C22_->digi((s1->v1() + s0->v1()) * invdH2);
      const double C01 = C01_->digi(-(H1v0 + H0v1) * invdH2);
      const double C23 = C23_->digi(-(H3v2 + H2v3) * invdH2);
      const double C11 = C11_->digi((H12v0 + H02v1) * invdH2);
      const double C33 = C33_->digi((H32v2 + H22v3) * invdH2);
      // create updated state
      states_.emplace_back(State(s1, {x0, x1, x2, x3, 0., 0., 0., C00, C11, C22, C33, C01, C23, 0., 0., 0.}));
      state = &states_.back();
      x0_->updateRangeActual(x0);
      x1_->updateRangeActual(x1);
      x2_->updateRangeActual(x2);
      x3_->updateRangeActual(x3);
      C00_->updateRangeActual(C00);
      C01_->updateRangeActual(C01);
      C11_->updateRangeActual(C11);
      C22_->updateRangeActual(C22);
      C23_->updateRangeActual(C23);
      C33_->updateRangeActual(C33);
    }
  }

  // apply final cuts
  void KalmanFilter::finalize() {
    for (State*& state : stream_) {
      if (!state)
        continue;
      TrackCTB* track = state->track();
      const double inv2R = track->inv2R() + state->x0();
      const double phiT = track->phiT() + state->x1();
      const double cot = track->zT() / setup_->chosenRofZ() + +state->x2();
      const double zT = track->zT() + state->x3();
      const double z0 = zT - setup_->chosenRofZ() * cot;
      // pt cut
      const bool validX0 = dataFormats_->format(Variable::inv2R, Process::dr).inRange(inv2R);
      // cut on phi sector boundaries
      const bool validX1 = dataFormats_->format(Variable::phiT, Process::dr).inRange(phiT);
      // cot cut
      const bool validX2 = dataFormats_->format(Variable::cot, Process::dr).inRange(cot);
      // zT cut
      const bool validX3 = dataFormats_->format(Variable::zT, Process::dr).inRange(zT);
      // z0 cut
      const bool invaldiZ0 = abs(z0) > setup_->beamWindowZ();
      // stub residual cut
      State* s = state;
      TTBV hitPattern(0, setup_->numLayers());
      TTBV hitPatternPS(0, setup_->numLayers());
      while ((s = s->parent())) {
        StubCTB* stub = s->stub();
        const double r = stub->r();
        double phi = stub->phi() - (state->x1() + r * state->x0());
        if (use5ParameterFit_)
          phi -= state->x0() / (r + setup_->chosenRofPhi());
        const double rz = r + H00_->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
        const double z = stub->z() - (state->x3() + rz * state->x2());
        const bool validPhi = dataFormats_->format(Variable::phi, Process::kf).inRange(phi);
        const bool validZ = dataFormats_->format(Variable::z, Process::kf).inRange(z);
        if (!validPhi || !validZ)
          continue;
        hitPattern.set(s->layer());
        if (setup_->psModule(stub->frame().first))
          hitPatternPS.set(s->layer());
      }
      // layer cut
      bool invalidLayers = hitPattern.count() < setup_->kfMinLayers();
      bool invalidLayersPS = hitPatternPS.count() < setup_->kfMinLayersPS();
      // apply
      if (invalidLayers || invalidLayersPS || !validX0 || !validX1 || !validX2 || !validX3 || invaldiZ0)
        state = nullptr;
    }
  }

  // Transform States into output products
  void KalmanFilter::conv(StreamsStub& streamsStub, StreamsTrack& streamsTrack) {
    const int offset = region_ * setup_->numLayers();
    StreamTrack& streamTrack = streamsTrack[region_];
    streamTrack.reserve(stream_.size());
    for (int layer = 0; layer < setup_->numLayers(); layer++)
      streamsStub[offset + layer].reserve(stream_.size());
    for (State* state : stream_) {
      TrackCTB* track = state->track();
      const double inv2R = track->inv2R() + state->x0();
      const double phiT = track->phiT() + state->x1();
      const double cot = track->zT() / setup_->chosenRofZ() + state->x2();
      const double zT = track->zT() + state->x3();
      const TTTrackRef& ttTrackRef = track->frame().first;
      const Track trackDR(ttTrackRef, dataFormats_, Process::dr, inv2R, phiT, cot, zT);
      streamTrack.emplace_back(trackDR.frame());
      TTBV hitPattern(0, setup_->numLayers());
      State* s = state;
      while ((s = s->parent())) {
        StubCTB* stub = s->stub();
        const double r = stub->r();
        double phi = stub->phi() - (state->x1() + r * state->x0());
        if (use5ParameterFit_)
          phi -= state->x0() / (r + setup_->chosenRofPhi());
        const double rz = r + H00_->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
        const double z = stub->z() - (state->x3() + rz * state->x2());
        const double dPhi = stub->dPhi();
        const double dZ = stub->dZ();
        const bool validPhi = dataFormats_->format(Variable::phi, Process::kf).inRange(phi);
        const bool validZ = dataFormats_->format(Variable::z, Process::kf).inRange(z);
        if (!validPhi || !validZ)
          continue;
        const StubKF stubKF(*stub, r, phi, z, dPhi, dZ);
        streamsStub[offset + s->layer()].emplace_back(stubKF.frame());
        hitPattern.set(s->layer());
      }
      for (int layer : hitPattern.ids(false))
        streamsStub[offset + layer].emplace_back(FrameStub());
      // store d0 in copied TTTracks
      if (use5ParameterFit_) {
        ttTracks_.emplace_back(ttTrackRef->rInv(),
                               ttTrackRef->phi(),
                               ttTrackRef->tanL(),
                               ttTrackRef->z0(),
                               state->x4(),
                               ttTrackRef->chi2XY(),
                               ttTrackRef->chi2Z(),
                               ttTrackRef->trkMVA1(),
                               ttTrackRef->trkMVA2(),
                               ttTrackRef->trkMVA3(),
                               ttTrackRef->hitPattern(),
                               5,
                               setup_->bField());
        ttTracks_.back().setPhiSector(ttTrackRef->phiSector());
        ttTracks_.back().setEtaSector(ttTrackRef->etaSector());
        ttTracks_.back().setTrackSeedType(ttTrackRef->trackSeedType());
        ttTracks_.back().setStubPtConsistency(ttTrackRef->stubPtConsistency());
        ttTracks_.back().setStubRefs(ttTrackRef->getStubRefs());
      }
    }
  }

  // adds a layer to states
  void KalmanFilter::addLayer() {
    // Latency of KF Associator block firmware
    static constexpr int latency = 5;
    // dynamic state container for clock accurate emulation
    deque<State*> streamOutput;
    // Memory stack used to handle combinatorics
    deque<State*> stack;
    // static delay container
    deque<State*> delay(latency, nullptr);
    // each trip corresponds to a f/w clock tick
    // done if no states to process left, taking as much time as needed
    while (!stream_.empty() || !stack.empty() ||
           !all_of(delay.begin(), delay.end(), [](const State* state) { return state == nullptr; })) {
      State* state = pop_front(stream_);
      // Process a combinatoric state if no (non-combinatoric?) state available
      if (!state)
        state = pop_front(stack);
      streamOutput.push_back(state);
      // The remainder of the code in this loop deals with combinatoric states.
      if (state)
        state = state->comb(states_, layer_);
      delay.push_back(state);
      state = pop_front(delay);
      if (state)
        stack.push_back(state);
    }
    stream_ = streamOutput;
    // Update state with next stub using KF maths
    for (State*& state : stream_)
      if (state)
        update(state);
  }

  // adds a layer to states to build seeds
  void KalmanFilter::addSeedLayer() {
    // Latency of KF Associator block firmware
    static constexpr int latency = 5;
    // dynamic state container for clock accurate emulation
    deque<State*> streamOutput;
    // Memory stack used to handle combinatorics
    deque<State*> stack;
    // static delay container
    deque<State*> delay(latency, nullptr);
    // each trip corresponds to a f/w clock tick
    // done if no states to process left, taking as much time as needed
    while (!stream_.empty() || !stack.empty() ||
           !all_of(delay.begin(), delay.end(), [](const State* state) { return state == nullptr; })) {
      State* state = pop_front(stream_);
      // Process a combinatoric state if no (non-combinatoric?) state available
      if (!state)
        state = pop_front(stack);
      streamOutput.push_back(state);
      // The remainder of the code in this loop deals with combinatoric states.
      if (state)
        state = state->combSeed(states_, layer_);
      delay.push_back(state);
      state = pop_front(delay);
      if (state)
        stack.push_back(state);
    }
    stream_ = streamOutput;
    // Update state with next stub using KF maths
    for (State*& state : stream_)
      if (state)
        state = state->update(states_, layer_);
  }

  // best state selection
  void KalmanFilter::accumulator() {
    // prepare arrival order
    vector<int> trackIds;
    trackIds.reserve(stream_.size());
    for (State* state : stream_) {
      const int trackId = state->trackId();
      if (find_if(trackIds.begin(), trackIds.end(), [trackId](int id) { return id == trackId; }) == trackIds.end())
        trackIds.push_back(trackId);
    }
    // sort in number of skipped layers
    auto numSkippedLayers = [](State* state) {
      const TTBV& hitPattern = state->hitPattern();
      TTBV pattern = state->maybePattern();
      pattern |= hitPattern;
      return pattern.count(0, hitPattern.pmEncode(true), false);
    };
    auto lessSkippedLayers = [numSkippedLayers](State* lhs, State* rhs) {
      return numSkippedLayers(lhs) < numSkippedLayers(rhs);
    };
    stable_sort(stream_.begin(), stream_.end(), lessSkippedLayers);
    // sort in number of consistent stubs
    auto isConsistent = [this](State* state, StubCTB* stub) {
      double phi = stub->phi() - (state->x1() + stub->r() * state->x0());
      if (use5ParameterFit_)
        phi -= state->x4() / (stub->r() + setup_->chosenRofPhi());
      const double rz = stub->r() + H00_->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
      const double z = stub->z() - (state->x3() + rz * state->x2());
      return m0_->digi(abs(phi)) - 1.e-12 < stub->dPhi() / 2. && m1_->digi(abs(z)) - 1.e-12 < stub->dZ() / 2.;
    };
    auto numConsistentLayers = [isConsistent](State* state) {
      int num(0);
      State* s = state;
      while ((s = s->parent()))
        if (isConsistent(state, s->stub()))
          num++;
      return num;
    };
    auto moreConsistentLayers = [numConsistentLayers](State* lhs, State* rhs) {
      return numConsistentLayers(lhs) > numConsistentLayers(rhs);
    };
    stable_sort(stream_.begin(), stream_.end(), moreConsistentLayers);
    // sort in number of consistent ps stubs
    auto isConsistentPS = [this](State* state, StubCTB* stub) {
      if (!setup_->psModule(stub->frame().first))
        return false;
      double phi = stub->phi() - (state->x1() + stub->r() * state->x0());
      if (use5ParameterFit_)
        phi -= state->x4() / (stub->r() + setup_->chosenRofPhi());
      const double rz = stub->r() + H00_->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
      const double z = stub->z() - (state->x3() + rz * state->x2());
      return m0_->digi(abs(phi)) - 1.e-12 < stub->dPhi() / 2. && m1_->digi(abs(z)) - 1.e-12 < stub->dZ() / 2.;
    };
    auto numConsistentLayersPS = [isConsistentPS](State* state) {
      int num(0);
      State* s = state;
      while ((s = s->parent()))
        if (isConsistentPS(state, s->stub()))
          num++;
      return num;
    };
    auto moreConsistentLayersPS = [numConsistentLayersPS](State* lhs, State* rhs) {
      return numConsistentLayersPS(lhs) > numConsistentLayersPS(rhs);
    };
    stable_sort(stream_.begin(), stream_.end(), moreConsistentLayersPS);
    // sort in track id as arrived
    auto order = [&trackIds](auto lhs, auto rhs) {
      const auto l = find(trackIds.begin(), trackIds.end(), lhs->trackId());
      const auto r = find(trackIds.begin(), trackIds.end(), rhs->trackId());
      return distance(r, l) < 0;
    };
    stable_sort(stream_.begin(), stream_.end(), order);
    // keep first state (best due to previous sorts) per track id
    stream_.erase(
        unique(stream_.begin(), stream_.end(), [](State* lhs, State* rhs) { return lhs->trackId() == rhs->trackId(); }),
        stream_.end());
  }

  // updates state
  void KalmanFilter::update4(State*& state) {
    if (state->layer() != layer_)
      return;
    // All variable names & equations come from Fruhwirth KF paper http://dx.doi.org/10.1016/0168-9002%2887%2990887-4", where F taken as unit matrix. Stub uncertainties projected onto (phi,z), assuming no correlations between r-phi & r-z planes.
    // stub phi residual wrt input helix
    const double m0 = state->m0();
    // stub z residual wrt input helix
    const double m1 = state->m1();
    // stub projected phi uncertainty squared);
    const double v0 = state->v0();
    // stub projected z uncertainty squared
    const double v1 = state->v1();
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofPhi
    const double H00 = state->H00();
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofZ
    const double H12 = state->H12();
    m0_->updateRangeActual(m0);
    m1_->updateRangeActual(m1);
    v0_->updateRangeActual(v0);
    v1_->updateRangeActual(v1);
    H00_->updateRangeActual(H00);
    H12_->updateRangeActual(H12);
    // helix inv2R wrt input helix
    double x0 = state->x0();
    // helix phi at radius ChosenRofPhi wrt input helix
    double x1 = state->x1();
    // helix cot(Theta) wrt input helix
    double x2 = state->x2();
    // helix z at radius chosenRofZ wrt input helix
    double x3 = state->x3();
    // cov. matrix
    double C00 = state->C00();
    double C01 = state->C01();
    double C11 = state->C11();
    double C22 = state->C22();
    double C23 = state->C23();
    double C33 = state->C33();
    // stub phi residual wrt current state
    const double r0C = x1_->digi(m0 - x1);
    const double r0 = r0_->digi(r0C - x0 * H00);
    // stub z residual wrt current state
    const double r1C = x3_->digi(m1 - x3);
    const double r1 = r1_->digi(r1C - x2 * H12);
    // matrix S = H*C
    const double S00 = S00_->digi(C01 + H00 * C00);
    const double S01 = S01_->digi(C11 + H00 * C01);
    const double S12 = S12_->digi(C23 + H12 * C22);
    const double S13 = S13_->digi(C33 + H12 * C23);
    // Cov. matrix of predicted residuals R = V+HCHt = C+H*St
    const double R00 = R00_->digi(v0 + S01 + H00 * S00);
    const double R11 = R11_->digi(v1 + S13 + H12 * S12);
    // improved dynamic cancelling
    //const int msb0 = R00_->width();
    const int msb0 = max(0, (int)ceil(log2(R00 / R00_->base())));
    const int shift0 = R00_->width() - msb0;
    const double shiftedR00 = R00 * pow(2., shift0);
    const double R00Rough = R00Rough_->digi(shiftedR00);
    const double invR00Approx = invR00Approx_->digi(1. / R00Rough);
    const double invR00Cor = invR00Cor_->digi(2. - invR00Approx * shiftedR00);
    const double invR00 = invR00_->digi(invR00Approx * invR00Cor);
    //const int msb1 = R11_->width();
    const int msb1 = max(0, (int)ceil(log2(R11 / R11_->base())));
    const int shift1 = R11_->width() - msb1;
    const double shiftedR11 = R11 * pow(2., shift1);
    const double R11Rough = R11Rough_->digi(shiftedR11);
    const double invR11Approx = invR11Approx_->digi(1. / R11Rough);
    const double invR11Cor = invR11Cor_->digi(2. - invR11Approx * shiftedR11);
    const double invR11 = invR11_->digi(invR11Approx * invR11Cor);
    // Kalman gain matrix K = S*R(inv)
    const double shiftedS00 = S00_->digi(S00 * pow(2., shift0));
    const double shiftedS01 = S01_->digi(S01 * pow(2., shift0));
    const double shiftedS12 = S12_->digi(S12 * pow(2., shift1));
    const double shiftedS13 = S13_->digi(S13 * pow(2., shift1));
    const double K00 = K00_->digi(shiftedS00 * invR00);
    const double K10 = K10_->digi(shiftedS01 * invR00);
    const double K21 = K21_->digi(shiftedS12 * invR11);
    const double K31 = K31_->digi(shiftedS13 * invR11);
    // Updated helix params, their cov. matrix
    x0 = x0_->digi(x0 + r0 * K00);
    x1 = x1_->digi(x1 + r0 * K10);
    x2 = x2_->digi(x2 + r1 * K21);
    x3 = x3_->digi(x3 + r1 * K31);
    C00 = C00_->digi(C00 - S00 * K00);
    C01 = C01_->digi(C01 - S01 * K00);
    C11 = C11_->digi(C11 - S01 * K10);
    C22 = C22_->digi(C22 - S12 * K21);
    C23 = C23_->digi(C23 - S13 * K21);
    C33 = C33_->digi(C33 - S13 * K31);
    // update variable ranges to tune variable granularity
    r0_->updateRangeActual(r0);
    r1_->updateRangeActual(r1);
    S00_->updateRangeActual(S00);
    S01_->updateRangeActual(S01);
    S12_->updateRangeActual(S12);
    S13_->updateRangeActual(S13);
    R00_->updateRangeActual(R00);
    R11_->updateRangeActual(R11);
    R00Rough_->updateRangeActual(R00Rough);
    invR00Approx_->updateRangeActual(invR00Approx);
    invR00Cor_->updateRangeActual(invR00Cor);
    invR00_->updateRangeActual(invR00);
    R11Rough_->updateRangeActual(R11Rough);
    invR11Approx_->updateRangeActual(invR11Approx);
    invR11Cor_->updateRangeActual(invR11Cor);
    invR11_->updateRangeActual(invR11);
    K00_->updateRangeActual(K00);
    K10_->updateRangeActual(K10);
    K21_->updateRangeActual(K21);
    K31_->updateRangeActual(K31);
    // create updated state
    states_.emplace_back(State(state, {x0, x1, x2, x3, 0., C00, C11, C22, C33, C01, C23, 0., 0., 0.}));
    state = &states_.back();
    x0_->updateRangeActual(x0);
    x1_->updateRangeActual(x1);
    x2_->updateRangeActual(x2);
    x3_->updateRangeActual(x3);
    C00_->updateRangeActual(C00);
    C01_->updateRangeActual(C01);
    C11_->updateRangeActual(C11);
    C22_->updateRangeActual(C22);
    C23_->updateRangeActual(C23);
    C33_->updateRangeActual(C33);
  }

  // updates state
  void KalmanFilter::update5(State*& state) {
    if (state->layer() != layer_)
      return;
    const double m0 = state->m0();
    const double m1 = state->m1();
    const double v0 = state->v0();
    const double v1 = state->v1();
    const double H00 = state->H00();
    const double H12 = state->H12();
    const double H04 = state->H04();
    double x0 = state->x0();
    double x1 = state->x1();
    double x2 = state->x2();
    double x3 = state->x3();
    double x4 = state->x4();
    double C00 = state->C00();
    double C01 = state->C01();
    double C11 = state->C11();
    double C22 = state->C22();
    double C23 = state->C23();
    double C33 = state->C33();
    double C44 = state->C44();
    double C40 = state->C40();
    double C41 = state->C41();
    const double r0 = m0 - x1 - x0 * H00 - x4 / H04;
    const double r1 = m1 - x3 - x2 * H12;
    const double S00 = C01 + H00 * C00 + C40 / H04;
    const double S01 = C11 + H00 * C01 + C41 / H04;
    const double S12 = C23 + H12 * C22;
    const double S13 = C33 + H12 * C23;
    const double S04 = C41 + H00 * C40 + C44 / H04;
    const double R00 = v0 + S01 + H00 * S00 + S04 / H04;
    const double R11 = v1 + S13 + H12 * S12;
    const double K00 = S00 / R00;
    const double K10 = S01 / R00;
    const double K21 = S12 / R11;
    const double K31 = S13 / R11;
    const double K40 = S04 / R00;
    x0 += r0 * K00;
    x1 += r0 * K10;
    x2 += r1 * K21;
    x3 += r1 * K31;
    x4 += r0 * K40;
    C00 -= S00 * K00;
    C01 -= S01 * K00;
    C11 -= S01 * K10;
    C22 -= S12 * K21;
    C23 -= S13 * K21;
    C33 -= S13 * K31;
    C44 -= S04 * K40;
    C40 -= S04 * K00;
    C41 -= S04 * K10;
    states_.emplace_back(State(state, {x0, x1, x2, x3, x4, C00, C11, C22, C33, C01, C23, C44, C40, C41}));
    state = &states_.back();
  }

  // remove and return first element of deque, returns nullptr if empty
  template <class T>
  T* KalmanFilter::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

}  // namespace trklet
