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
                             int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        layerEncoding_(layerEncoding),
        kalmanFilterFormats_(kalmanFilterFormats),
        region_(region),
        layer_(0),
        x0_(setup_->numLayers()),
        x1_(setup_->numLayers()),
        x2_(setup_->numLayers()),
        x3_(setup_->numLayers()),
        H00_(setup_->numLayers()),
        H12_(setup_->numLayers()),
        m0_(setup_->numLayers()),
        m1_(setup_->numLayers()),
        v0_(setup_->numLayers()),
        v1_(setup_->numLayers()),
        r0_(setup_->numLayers()),
        r1_(setup_->numLayers()),
        S00_(setup_->numLayers()),
        S01_(setup_->numLayers()),
        S12_(setup_->numLayers()),
        S13_(setup_->numLayers()),
        K00_(setup_->numLayers()),
        K10_(setup_->numLayers()),
        K21_(setup_->numLayers()),
        K31_(setup_->numLayers()),
        R00_(setup_->numLayers()),
        R11_(setup_->numLayers()),
        R00Rough_(setup_->numLayers()),
        R11Rough_(setup_->numLayers()),
        invR00Approx_(setup_->numLayers()),
        invR11Approx_(setup_->numLayers()),
        invR00Cor_(setup_->numLayers()),
        invR11Cor_(setup_->numLayers()),
        invR00_(setup_->numLayers()),
        invR11_(setup_->numLayers()),
        C00_(setup_->numLayers()),
        C01_(setup_->numLayers()),
        C11_(setup_->numLayers()),
        C22_(setup_->numLayers()),
        C23_(setup_->numLayers()),
        C33_(setup_->numLayers()),
        r02_(setup_->numLayers()),
        r12_(setup_->numLayers()),
        chi20_(setup_->numLayers()),
        chi21_(setup_->numLayers()) {
    for (int layer = 0; layer < setup_->numLayers(); layer++) {
      x0_[layer] = &kalmanFilterFormats_->format(VariableKF::x0, layer);
      x1_[layer] = &kalmanFilterFormats_->format(VariableKF::x1, layer);
      x2_[layer] = &kalmanFilterFormats_->format(VariableKF::x2, layer);
      x3_[layer] = &kalmanFilterFormats_->format(VariableKF::x3, layer);
      H00_[layer] = &kalmanFilterFormats_->format(VariableKF::H00, layer);
      H12_[layer] = &kalmanFilterFormats_->format(VariableKF::H12, layer);
      m0_[layer] = &kalmanFilterFormats_->format(VariableKF::m0, layer);
      m1_[layer] = &kalmanFilterFormats_->format(VariableKF::m1, layer);
      v0_[layer] = &kalmanFilterFormats_->format(VariableKF::v0, layer);
      v1_[layer] = &kalmanFilterFormats_->format(VariableKF::v1, layer);
      r0_[layer] = &kalmanFilterFormats_->format(VariableKF::r0, layer);
      r1_[layer] = &kalmanFilterFormats_->format(VariableKF::r1, layer);
      S00_[layer] = &kalmanFilterFormats_->format(VariableKF::S00, layer);
      S01_[layer] = &kalmanFilterFormats_->format(VariableKF::S01, layer);
      S12_[layer] = &kalmanFilterFormats_->format(VariableKF::S12, layer);
      S13_[layer] = &kalmanFilterFormats_->format(VariableKF::S13, layer);
      K00_[layer] = &kalmanFilterFormats_->format(VariableKF::K00, layer);
      K10_[layer] = &kalmanFilterFormats_->format(VariableKF::K10, layer);
      K21_[layer] = &kalmanFilterFormats_->format(VariableKF::K21, layer);
      K31_[layer] = &kalmanFilterFormats_->format(VariableKF::K31, layer);
      R00_[layer] = &kalmanFilterFormats_->format(VariableKF::R00, layer);
      R11_[layer] = &kalmanFilterFormats_->format(VariableKF::R11, layer);
      R00Rough_[layer] = &kalmanFilterFormats_->format(VariableKF::R00Rough, layer);
      R11Rough_[layer] = &kalmanFilterFormats_->format(VariableKF::R11Rough, layer);
      invR00Approx_[layer] = &kalmanFilterFormats_->format(VariableKF::invR00Approx, layer);
      invR11Approx_[layer] = &kalmanFilterFormats_->format(VariableKF::invR11Approx, layer);
      invR00Cor_[layer] = &kalmanFilterFormats_->format(VariableKF::invR00Cor, layer);
      invR11Cor_[layer] = &kalmanFilterFormats_->format(VariableKF::invR11Cor, layer);
      invR00_[layer] = &kalmanFilterFormats_->format(VariableKF::invR00, layer);
      invR11_[layer] = &kalmanFilterFormats_->format(VariableKF::invR11, layer);
      C00_[layer] = &kalmanFilterFormats_->format(VariableKF::C00, layer);
      C01_[layer] = &kalmanFilterFormats_->format(VariableKF::C01, layer);
      C11_[layer] = &kalmanFilterFormats_->format(VariableKF::C11, layer);
      C22_[layer] = &kalmanFilterFormats_->format(VariableKF::C22, layer);
      C23_[layer] = &kalmanFilterFormats_->format(VariableKF::C23, layer);
      C33_[layer] = &kalmanFilterFormats_->format(VariableKF::C33, layer);
      r02_[layer] = &kalmanFilterFormats_->format(VariableKF::r02, layer);
      r12_[layer] = &kalmanFilterFormats_->format(VariableKF::r12, layer);
      chi20_[layer] = &kalmanFilterFormats_->format(VariableKF::chi20, layer);
      chi21_[layer] = &kalmanFilterFormats_->format(VariableKF::chi21, layer);
    }
  }

  // read in and organize input tracks and stubs
  void KalmanFilter::consume(const StreamsTrack& streamsTrack, const StreamsStub& streamsStub) {
    static const int numLayers = setup_->numLayers();
    const int offset = region_ * numLayers;
    const StreamTrack& streamTrack = streamsTrack[region_];
    const int numTracks = accumulate(streamTrack.begin(), streamTrack.end(), 0, [](int& sum, const FrameTrack& f){ return sum += (f.first.isNull() ? 0 : 1); });
    int numStubs(0);
    for (int layer = 0; layer < numLayers; layer++) {
      const StreamStub& streamStub = streamsStub[offset + layer];
      numStubs += accumulate(streamStub.begin(), streamStub.end(), 0, [](int& sum, const FrameStub& f){ return sum += (f.first.isNull() ? 0 : 1); });
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
    // Propagate state to each layer in turn, updating it with all viable stub combinations there, using KF maths
    for (layer_ = 0; layer_ < setup_->numLayers(); layer_++)
      addLayer();
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

  // apply final cuts
  void KalmanFilter::finalize() {
    for (State*& state : stream_) {
      if (!state)
        continue;
      // layer cut
      bool invalidLayers = state->hitPattern().count() < setup_->kfMinLayers();
      // pt cut
      const bool invalidX0 =
          abs(state->x0() + state->track()->inv2R()) >
          setup_->invPtToDphi() / setup_->minPt() + dataFormats_->format(Variable::inv2R, Process::ht).base();
      // cut on phi sector boundaries
      const bool invalidX1 =
          abs(state->x1() + state->track()->phiT()) > dataFormats_->format(Variable::phiT, Process::gp).range() / 2.;
      // z0 cut
      static const DataFormat& dfZT = dataFormats_->format(Variable::zT, Process::kf);
      const double z0 = dfZT.digi(state->x3() - H12_[0]->digi(setup_->chosenRofZ()) * state->x2());
      const bool invaldiZ0 = abs(z0) > dfZT.digi(setup_->beamWindowZ());
      // stub residual cut
      State* s = state;
      TTBV hits(0, setup_->numLayers());
      while ((s = s->parent())) {
        StubCTB* stub = s->stub();
        const double r = stub->r();
        const double phi = stub->phi() - (state->x1() + r * state->x0());
        const double rz = r + H00_[0]->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
        const double z = stub->z() - (state->x3() + rz * state->x2());
        if (dataFormats_->format(Variable::phi, Process::kf).inRange(phi) &&
            dataFormats_->format(Variable::z, Process::kf).inRange(z))
          hits.set(s->layer());
      }
      if (hits.count() < setup_->kfMinLayers())
        invalidLayers = true;
      // apply
      if (invalidLayers || invalidX0 || invalidX1 || invaldiZ0)
        state = nullptr;
    }
  }

  // Transform States into output products
  void KalmanFilter::conv(StreamsStub& streamsStub, StreamsTrack& streamsTrack) const {
    const int offset = region_ * setup_->numLayers();
    StreamTrack& streamTrack = streamsTrack[region_];
    streamTrack.reserve(stream_.size());
    for (int layer = 0; layer < setup_->numLayers(); layer++)
      streamsStub[offset + layer].reserve(stream_.size());
    for (State* state : stream_) {
      TrackCTB* track = state->track();
      const double inv2R = track->inv2R() + state->x0();
      const double phiT = track->phiT() + state->x1();
      const double cot = state->x2();
      const double zT = track->zT() + state->x3();
      const TTTrackRef& ttTrackRef = track->frame().first;
      const TTBV hwinv2R(dataFormats_->format(Variable::inv2R, Process::kf).ttBV(inv2R));
      const TTBV hwphiT(dataFormats_->format(Variable::phiT, Process::kf).ttBV(phiT));
      const TTBV hwcot(dataFormats_->format(Variable::cot, Process::kf).ttBV(cot));
      const TTBV hwzT(dataFormats_->format(Variable::zT, Process::kf).ttBV(zT));
      streamTrack.emplace_back(ttTrackRef, "1" + hwinv2R.str() + hwphiT.str() + hwcot.str() + hwzT.str() + "0");
      State* s = state;
      while ((s = s->parent())) {
        StubCTB* stub = s->stub();
        const double r = stub->r();
        const double phi = stub->phi() - (state->x1() + r * state->x0());
        const double rz = r + H00_[0]->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
        const double z = stub->z() - (state->x3() + rz * state->x2());
        const double dPhi = stub->dPhi();
        const double dZ = stub->dZ();
        if (!dataFormats_->format(Variable::phi, Process::kf).inRange(phi) || !dataFormats_->format(Variable::z, Process::kf).inRange(z))
          continue;
        const TTStubRef& ttStubRef = stub->frame().first;
        const TTBV hwr(dataFormats_->format(Variable::r, Process::kf).ttBV(r));
        const TTBV hwphi(dataFormats_->format(Variable::phi, Process::kf).ttBV(phi));
        const TTBV hwz(dataFormats_->format(Variable::z, Process::kf).ttBV(z));
        const TTBV hwdPhi(dataFormats_->format(Variable::dPhi, Process::ctb).ttBV(dPhi));
        const TTBV hwdZ(dataFormats_->format(Variable::dZ, Process::ctb).ttBV(dZ));
        streamsStub[offset + s->layer()].emplace_back(ttStubRef, "1" + hwr.str() + hwphi.str() + hwz.str() + hwdPhi.str() + hwdZ.str());
      }
      for (int layer : state->hitPattern().ids(false))
        streamsStub[offset + layer].emplace_back(FrameStub());
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
      if (state && !state->isDone())
        update(state);
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
      const double phi = stub->phi() - (state->x1() + stub->r() * state->x0());
      const double rz = stub->r() + H00_[0]->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
      const double z = stub->z() - (state->x3() + rz * state->x2());
      return m0_[0]->digi(abs(phi)) - 1.e-12 < stub->dPhi() / 2. && m1_[0]->digi(abs(z)) - 1.e-12 < stub->dZ() / 2.;
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
  void KalmanFilter::update(State*& state) {
    static const int shifChi20 = setup_->kfShiftChi20();
    static const int shifChi21 = setup_->kfShiftChi21();
    //static const double chi2cut = pow(2.0, setup_->kfPowCutChi2());
    static const double chi2cut = 9.e9;
    if (state->isSkip()) {
      if (state->trackPattern()[layer_ + 1])
        state = state->unskip(states_, layer_ + 1);
      return;
    }
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
    m0_[layer_]->updateRangeActual(m0);
    m1_[layer_]->updateRangeActual(m1);
    v0_[layer_]->updateRangeActual(v0);
    v1_[layer_]->updateRangeActual(v1);
    H00_[layer_]->updateRangeActual(H00);
    H12_[layer_]->updateRangeActual(H12);
    // helix inv2R wrt input helix
    double x0 = 0.;
    // helix phi at radius ChosenRofPhi wrt input helix
    double x1 = 0.;
    // helix cot(Theta) wrt input helix
    double x2 = 0.;
    // helix z at radius chosenRofZ wrt input helix
    double x3 = 0.;
    // cov. matrix
    double C00 = 0.;
    double C01 = 0.;
    double C11 = 0.;
    double C22 = 0.;
    double C23 = 0.;
    double C33 = 0.;
    // chi2s
    double chi20 = 0.;
    double chi21 = 0.;
    stringstream ss00;
    if (state->hitPattern().count() == 1) {
      const State* p = state->parent();
      const double dH00 = H00 - p->H00();
      const double dH002 = dH00 * dH00;
      const double dH12 = H12 - p->H12();
      const double dH122 = dH12 * dH12;
      x0 = (m0 - p->m0()) / dH00;
      x1 = (H00 * p->m0() - m0 * p->H00()) / dH00;
      x2 = (m1 - p->m1()) / dH12;
      x3 = (H12 * p->m1() - m1 * p->H12()) / dH12;
      C00 = (v0 + p->v0()) / dH002;
      C01 = -(H00 * p->v0() + p->H00() * v0) / dH002;
      C11 = (H00 * H00 * p->v0() + p->H00() * p->H00() * v0) / dH002;
      C22 = (v1 + p->v1()) / dH122;
      C23 = -(H12 * p->v1() + p->H12() * v1) / dH122;
      C33 = (H12 * H12 * p->v1() + p->H12() * p->H12() * v1) / dH122;
    } else if (state->hitPattern().count() > 1) {
      x0 = state->x0();
      x1 = state->x1();
      x2 = state->x2();
      x3 = state->x3();
      C00 = state->C00();
      C01 = state->C01();
      C11 = state->C11();
      C22 = state->C22();
      C23 = state->C23();
      C33 = state->C33();
      chi20 = state->chi20();
      chi21 = state->chi21();
      // stub phi residual wrt current state
      const double r0C = x1_[layer_]->digi(m0 - x1);
      const double r0 = r0_[layer_]->digi(r0C - x0 * H00);
      // stub z residual wrt current state
      const double r1C = x3_[layer_]->digi(m1 - x3);
      const double r1 = r1_[layer_]->digi(r1C - x2 * H12);
      // matrix S = H*C
      const double S00 = S00_[layer_]->digi(C01 + H00 * C00);
      const double S01 = S01_[layer_]->digi(C11 + H00 * C01);
      const double S12 = S12_[layer_]->digi(C23 + H12 * C22);
      const double S13 = S13_[layer_]->digi(C33 + H12 * C23);
      // Cov. matrix of predicted residuals R = V+HCHt = C+H*St
      const double R00C = S01_[layer_]->digi(v0 + S01);
      const double R00 = R00_[layer_]->digi(R00C + H00 * S00);
      const double R11C = S13_[layer_]->digi(v1 + S13);
      const double R11 = R11_[layer_]->digi(R11C + H12 * S12);
      // improved dynamic cancelling
      const int msb0 = R00_[layer_]->width();
      //const int msb0 = max(0, (int)ceil(log2(R00 / R00_[layer_]->base()) - 1.e-12));
      const double R00Rough = R00Rough_[layer_]->digi(R00 * pow(2., R00_[layer_]->width() - msb0));
      const double invR00Approx = invR00Approx_[layer_]->digi(1. / R00Rough);
      const double invR00Cor = invR00Cor_[layer_]->digi(2. - invR00Approx * R00 * pow(2., R00_[layer_]->width() - msb0));
      const double invR00 = invR00_[layer_]->digi(invR00Approx * invR00Cor);
      const int msb1 = R11_[layer_]->width();
      //const int msb1 = max(0, (int)ceil(log2(R11 / R11_[layer_]->base()) - 1.e-12));
      const double R11Rough = R11Rough_[layer_]->digi(R11 * pow(2., R11_[layer_]->width() - msb1));
      const double invR11Approx = invR11Approx_[layer_]->digi(1. / R11Rough);
      const double invR11Cor = invR11Cor_[layer_]->digi(2. - invR11Approx * R11 * pow(2., R11_[layer_]->width() - msb1));
      const double invR11 = invR11_[layer_]->digi(invR11Approx * invR11Cor);
      // Kalman gain matrix K = S*R(inv)
      const double K00 = K00_[layer_]->digi(S00 * invR00 * pow(2., R00_[layer_]->width() - msb0));
      const double K10 = K10_[layer_]->digi(S01 * invR00 * pow(2., R00_[layer_]->width() - msb0));
      const double K21 = K21_[layer_]->digi(S12 * invR11 * pow(2., R11_[layer_]->width() - msb1));
      const double K31 = K31_[layer_]->digi(S13 * invR11 * pow(2., R11_[layer_]->width() - msb1));
      /*ss00 << std::fixed << std::showpoint << std::setprecision(12);
      ss00 << msb0 << " " << R00 * invR00 * pow(2., R00_[layer_]->width() - msb0) << endl;
      ss00 << C00 << " " << C00 / C00_[layer_]->base() << endl;
      ss00 << R00 << " " << (v0 + S01 + H00 * S00) << " " << (R00) / R00_[layer_]->base() << " " << (v0 + S01 + H00 * S00) / R00_[layer_]->base() << endl;
      ss00 << invR00 * pow(2., R00_[layer_]->width() - msb0) << " " << (1. / R00) << " " << (invR00 * pow(2., R00_[layer_]->width() - msb0)) / invR00_[layer_]->base() << " " << (1. / R00) / invR00_[layer_]->base() << endl;
      ss00 << K00 << " " << K00_[layer_]->digi((C01 + H00 * C00) / R00) << " " << K00_[layer_]->integer(K00) << " " << K00_[layer_]->integer((C01 + H00 * C00) / R00) << " " << K00_[layer_]->integer(S00 / R00) << " " << K00_[layer_]->integer((C01 + H00 * C00) * invR00 * pow(2., R00_[layer_]->width() - msb0)) << endl;
      ss00 << S00 << " " << (C01 + H00 * C00) << " " << S00 / S00_[layer_]->base() << " " << (C01 + H00 * C00) / S00_[layer_]->base() << endl;
      ss00 << C00_[layer_]->digi(C00 - S00 * K00) << " " << (C00 - (C01 + H00 * C00) * (C01 + H00 * C00) / R00) << " " << (C00 - S00 * K00) / C00_[layer_]->base() << " " << (C00 - (C01 + H00 * C00) * (C01 + H00 * C00) / R00) / C00_[layer_]->base() << endl;
      State* p = state;
      while ((p = p->parent())) {
        if (p->hitPattern().count() == 2)
          ss00 << "z0" << p->x3() - setup_->chosenRofZ() * p->x2() << endl;
        if (p->hitPattern().count() == 1)
          ss00 << "dH00 " << p->H00() - p->parent()->H00() << endl;
      }*/
      // Updated helix params, their cov. matrix & chi2s
      x0 = x0_[layer_]->digi(x0 + r0 * K00);
      x1 = x1_[layer_]->digi(x1 + r0 * K10);
      x2 = x2_[layer_]->digi(x2 + r1 * K21);
      x3 = x3_[layer_]->digi(x3 + r1 * K31);
      C00 = C00_[layer_]->digi(C00 - S00 * K00);
      C01 = C01_[layer_]->digi(C01 - S01 * K00);
      C11 = C11_[layer_]->digi(C11 - S01 * K10);
      C22 = C22_[layer_]->digi(C22 - S12 * K21);
      C23 = C23_[layer_]->digi(C23 - S13 * K21);
      C33 = C33_[layer_]->digi(C33 - S13 * K31);
      // squared residuals
      const double r02 = r02_[layer_]->digi(r0 * r0);
      const double r12 = r12_[layer_]->digi(r1 * r1);
      chi20 = chi20_[layer_]->digi(chi20 + r02 * invR00 * pow(2., shifChi20));
      chi21 = chi21_[layer_]->digi(chi21 + r12 * invR11 * pow(2., shifChi21));
      // update variable ranges to tune variable granularity
      r0_[layer_]->updateRangeActual(r0);
      r1_[layer_]->updateRangeActual(r1);
      S00_[layer_]->updateRangeActual(S00);
      S01_[layer_]->updateRangeActual(S01);
      S12_[layer_]->updateRangeActual(S12);
      S13_[layer_]->updateRangeActual(S13);
      R00_[layer_]->updateRangeActual(R00);
      R11_[layer_]->updateRangeActual(R11);
      R00Rough_[layer_]->updateRangeActual(R00Rough);
      invR00Approx_[layer_]->updateRangeActual(invR00Approx);
      invR00Cor_[layer_]->updateRangeActual(invR00Cor);
      invR00_[layer_]->updateRangeActual(invR00);
      R11Rough_[layer_]->updateRangeActual(R11Rough);
      invR11Approx_[layer_]->updateRangeActual(invR11Approx);
      invR11Cor_[layer_]->updateRangeActual(invR11Cor);
      invR11_[layer_]->updateRangeActual(invR11);
      K00_[layer_]->updateRangeActual(K00);
      K10_[layer_]->updateRangeActual(K10);
      K21_[layer_]->updateRangeActual(K21);
      K31_[layer_]->updateRangeActual(K31);
      //r02_[layer_]->updateRangeActual(r02);
      //r12_[layer_]->updateRangeActual(r12);
    }
    // cut on eta sector boundaries
    //const bool invalidX3 = abs(x3) > dataFormats_->format(Variable::zT, Process::gp).base() / 2.;
    // cut on triple found inv2R window
    const bool invalidX0 = abs(x0) > 1.5 * dataFormats_->format(Variable::inv2R, Process::ht).base();
    // cut on triple found phiT window
    const bool invalidX1 = abs(x1) > 1.5 * dataFormats_->format(Variable::phiT, Process::ht).base();
    // cot cut
    const bool invalidX2 = abs(x2) > dataFormats_->format(Variable::cot, Process::gp).base() / 2.;
    // z0 cut
    const bool invalidZ0 = abs(x3 - x2 * setup_->chosenRofZ()) > setup_->beamWindowZ();
    // chi2 cut
    const double dof = state->hitPattern().count() - 1;
    const double chi2 = dof > 0 ? (chi20 + chi21) / 2. / dof : 0.;
    const bool validChi2 = chi2 < chi2cut;
    //if (invalidX0 || invalidX1 || invalidX2 || invalidZ0 || !validChi2) {
    if (invalidX0 || invalidX1 || invalidX2 || invalidZ0) {
      state = nullptr;
      return;
    }
    //if (C00 < 0)
      //cout << ss00.str();
    // create updated state
    states_.emplace_back(State(state, {x0, x1, x2, x3, chi20, chi21, C00, C11, C22, C33, C01, C23}));
    state = &states_.back();
    x0_[layer_]->updateRangeActual(x0);
    x1_[layer_]->updateRangeActual(x1);
    x2_[layer_]->updateRangeActual(x2);
    x3_[layer_]->updateRangeActual(x3);
    C00_[layer_]->updateRangeActual(C00);
    C01_[layer_]->updateRangeActual(C01);
    C11_[layer_]->updateRangeActual(C11);
    C22_[layer_]->updateRangeActual(C22);
    C23_[layer_]->updateRangeActual(C23);
    C33_[layer_]->updateRangeActual(C33);
    chi20_[layer_]->updateRangeActual(chi20);
    chi21_[layer_]->updateRangeActual(chi21);
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
