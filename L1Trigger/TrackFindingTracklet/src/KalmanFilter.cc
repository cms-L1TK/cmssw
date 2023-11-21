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
                             vector<TrackKF>& tracks,
                             vector<StubKF>& stubs)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        layerEncoding_(layerEncoding),
        kalmanFilterFormats_(kalmanFilterFormats),
        tracks_(tracks),
        stubs_(stubs),
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
        C33_(&kalmanFilterFormats_->format(VariableKF::C33)) {
    C00_->updateRangeActual(pow(dataFormats_->base(Variable::inv2R, Process::ht), 2) *
                                pow(2, setup_->kfShiftInitialC00()) -
                            .5 * C00_->base());
    C11_->updateRangeActual(pow(dataFormats_->base(Variable::phiT, Process::ht), 2) *
                                pow(2, setup_->kfShiftInitialC11()) -
                            .5 * C11_->base());
    C22_->updateRangeActual(pow(dataFormats_->base(Variable::cot, Process::gp), 2) *
                                pow(2, setup_->kfShiftInitialC22()) -
                            .5 * C22_->base());
    C33_->updateRangeActual(pow(dataFormats_->base(Variable::zT, Process::gp), 2) *
                                pow(2, setup_->kfShiftInitialC33()) -
                            .5 * C33_->base());
  }

  // fill output products
  void KalmanFilter::produce(const vector<vector<TrackCTB*>>& tracksIn,
                             const vector<vector<StubCTB*>>& stubsIn,
                             vector<vector<TrackKF*>>& tracksOut,
                             vector<vector<vector<StubKF*>>>& stubsOut,
                             int& numAcceptedStates,
                             int& numLostStates) {
    static const int numChannel = dataFormats_->numChannel(Process::kf);
    static const int numLayers = setup_->numLayers();
    for (int channel = 0; channel < numChannel; channel++) {
      //static int debug = 0;
      const int offsetL = channel * numLayers;
      const vector<TrackCTB*>& tracksChannel = tracksIn[channel];
      // proto state creation
      deque<State*> stream;
      int trackId(0);
      for (int frame = 0; frame < (int)tracksChannel.size();) {
        TrackCTB* track = tracksChannel[frame];
        if (!track) {
          frame++;
          continue;
        }
        const auto begin = next(tracksChannel.begin(), frame);
        const auto end = find_if(begin + 1, tracksChannel.end(), [](TrackCTB* track) { return track; });
        const int size = distance(begin, end);
        vector<vector<StubCTB*>> stubs(numLayers);
        for (vector<StubCTB*>& layer : stubs)
          layer.reserve(size);
        for (int layer = 0; layer < numLayers; layer++) {
          const vector<StubCTB*>& layerAll = stubsIn[layer + offsetL];
          vector<StubCTB*>& layerTrack = stubs[layer];
          for (int frameS = 0; frameS < size; frameS++) {
            StubCTB* stub = layerAll[frameS + frame];
            if (!stub)
              break;
            layerTrack.push_back(stub);
          }
        }
        const TTBV& maybePattern = layerEncoding_->maybePattern(track->zT());
        states_.emplace_back(kalmanFilterFormats_, track, stubs, maybePattern, trackId++);
        stream.insert(stream.end(), size - 1, nullptr);
        stream.push_back(&states_.back());
        frame += size;
      }
      // Propagate state to each layer in turn, updating it with all viable stub combinations there, using KF maths
      for (layer_ = 0; layer_ < numLayers; layer_++)
        addLayer(stream);
      // apply final cuts
      finalize(stream);
      // count total number of final states
      const int nStates =
          accumulate(stream.begin(), stream.end(), 0, [](int& sum, State* state) { return sum += (state ? 1 : 0); });
      // apply truncation
      if (enableTruncation_ && (int)stream.size() > setup_->numFrames())
        stream.resize(setup_->numFrames());
      // cycle event, remove gaps
      stream.erase(remove(stream.begin(), stream.end(), nullptr), stream.end());
      // store number of states which got taken into account
      numAcceptedStates += (int)stream.size();
      // store number of states which got not taken into account due to truncation
      numLostStates += nStates - (int)stream.size();
      // best track per candidate selection
      accumulator(stream);
      // Transform States into Tracks
      vector<TrackKF*>& tracks = tracksOut[channel];
      vector<vector<StubKF*>>& stubs = stubsOut[channel];
      conv(stream, tracks, stubs);
    }
  }

  // apply final cuts
  void KalmanFilter::finalize(deque<State*>& stream) {
    for (State*& state : stream) {
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
      const double z0 = dfZT.digi(state->x3() - H12_->digi(setup_->chosenRofZ()) * state->x2());
      const bool invaldiZ0 = abs(z0) > dfZT.digi(setup_->beamWindowZ());
      // stub residual cut
      State* s = state;
      TTBV hits(0, setup_->numLayers());
      while ((s = s->parent())) {
        StubCTB* stub = s->stub();
        const double r = stub->r();
        const double phi = stub->phi() - (state->x1() + r * state->x0());
        const double rz = r + H00_->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
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

  // Transform States into Tracks
  void KalmanFilter::conv(const deque<State*>& states, vector<TrackKF*>& tracks, vector<vector<StubKF*>>& stubs) {
    static const DataFormat& dfInv2R = dataFormats_->format(Variable::inv2R, Process::ht);
    static const DataFormat& dfPhiT = dataFormats_->format(Variable::phiT, Process::ht);
    const int nTracks =
        accumulate(states.begin(), states.end(), 0, [](int& sum, State* s) { return sum += (s ? 1 : 0); });
    tracks.reserve(nTracks);
    for (vector<StubKF*>& layer : stubs)
      layer.reserve(nTracks);
    for (State* state : states) {
      State* s = state;
      while ((s = s->parent())) {
        StubCTB* stub = s->stub();
        const double r = stub->r();
        const double phi = stub->phi() - (state->x1() + r * state->x0());
        const double rz = r + H00_->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
        const double z = stub->z() - (state->x3() + rz * state->x2());
        const double dPhi = stub->dPhi();
        const double dZ = stub->dZ();
        if (dataFormats_->format(Variable::phi, Process::kf).inRange(phi) &&
            dataFormats_->format(Variable::z, Process::kf).inRange(z)) {
          stubs_.emplace_back(*stub, r, phi, z, dPhi, dZ);
          stubs[s->layer()].push_back(&stubs_.back());
        } else
          stubs[s->layer()].push_back(nullptr);
      }
      for (int layer : state->hitPattern().ids(false))
        stubs[layer].push_back(nullptr);
      TrackCTB* track = state->track();
      const double inv2R = track->inv2R() + state->x0();
      const double phiT = track->phiT() + state->x1();
      const double cot = state->x2();
      const double zT = track->zT() + state->x3();
      const bool inInv2R = dfInv2R.integer(inv2R) == dfInv2R.integer(track->inv2R());
      const bool inPhiT = dfPhiT.integer(phiT) == dfPhiT.integer(track->phiT());
      const TTBV match(inInv2R && inPhiT, 1);
      tracks_.emplace_back(*track, inv2R, phiT, cot, zT, match);
      tracks.push_back(&tracks_.back());
    }
  }

  // adds a layer to states
  void KalmanFilter::addLayer(deque<State*>& stream) {
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
    while (!stream.empty() || !stack.empty() ||
           !all_of(delay.begin(), delay.end(), [](const State* state) { return state == nullptr; })) {
      State* state = pop_front(stream);
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
    stream = streamOutput;
    // Update state with next stub using KF maths
    for (State*& state : stream)
      if (state && !state->isDone())
        update(state);
  }

  // best state selection
  void KalmanFilter::accumulator(deque<State*>& stream) {
    // prepare arrival order
    vector<int> trackIds;
    trackIds.reserve(stream.size());
    for (State* state : stream) {
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
    stable_sort(stream.begin(), stream.end(), lessSkippedLayers);
    // sort in number of consistent stubs
    auto isConsistent = [this](State* state, StubCTB* stub) {
      const double phi = stub->phi() - (state->x1() + stub->r() * state->x0());
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
    stable_sort(stream.begin(), stream.end(), moreConsistentLayers);
    // sort in track id as arrived
    auto order = [&trackIds](auto lhs, auto rhs) {
      const auto l = find(trackIds.begin(), trackIds.end(), lhs->trackId());
      const auto r = find(trackIds.begin(), trackIds.end(), rhs->trackId());
      return distance(r, l) < 0;
    };
    stable_sort(stream.begin(), stream.end(), order);
    // keep first state (best due to previous sorts) per track id
    stream.erase(
        unique(stream.begin(), stream.end(), [](State* lhs, State* rhs) { return lhs->trackId() == rhs->trackId(); }),
        stream.end());
  }

  // updates state
  void KalmanFilter::update(State*& state) {
    if (state->isSkip()) {
      if (layer_ == setup_->numLayers() - 1)
        throw cms::Exception("logic");
      else if (state->trackPattern()[layer_ + 1])
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
    // helix inv2R wrt input helix
    double x0 = state->x0();
    // helix phi at radius ChosenRofPhi wrt input helix
    double x1 = state->x1();
    // helix cot(Theta) wrt input helix
    double x2 = state->x2();
    // helix z at radius chosenRofZ wrt input helix
    double x3 = state->x3();
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofPhi
    const double H00 = state->H00();
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofZ
    const double H12 = state->H12();
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
    const double R00C = S01_->digi(v0 + S01);
    const double R00 = R00_->digi(R00C + H00 * S00);
    const double R11C = S13_->digi(v1 + S13);
    const double R11 = R11_->digi(R11C + H12 * S12);
    // improved dynamic cancelling
    const int msb0 = max(0, (int)ceil(log2(R00 / R00_->base())));
    const int msb1 = max(0, (int)ceil(log2(R11 / R11_->base())));
    const double R00Rough = R00Rough_->digi(R00 * pow(2., 16 - msb0));
    const double invR00Approx = invR00Approx_->digi(1. / R00Rough);
    const double invR00Cor = invR00Cor_->digi(2. - invR00Approx * R00Rough);
    const double invR00 = invR00_->digi(invR00Approx * invR00Cor * pow(2., 16 - msb0));
    const double R11Rough = R11Rough_->digi(R11 * pow(2., 16 - msb1));
    const double invR11Approx = invR11Approx_->digi(1. / R11Rough);
    const double invR11Cor = invR11Cor_->digi(2. - invR11Approx * R11Rough);
    const double invR11 = invR11_->digi(invR11Approx * invR11Cor * pow(2., 16 - msb1));
    // Kalman gain matrix K = S*R(inv)
    const double K00 = K00_->digi(S00 * invR00);
    const double K10 = K10_->digi(S01 * invR00);
    const double K21 = K21_->digi(S12 * invR11);
    const double K31 = K31_->digi(S13 * invR11);
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
    m0_->updateRangeActual(m0);
    m1_->updateRangeActual(m1);
    v0_->updateRangeActual(v0);
    v1_->updateRangeActual(v1);
    H00_->updateRangeActual(H00);
    H12_->updateRangeActual(H12);
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
    states_.emplace_back(State(state, {x0, x1, x2, x3, C00, C11, C22, C33, C01, C23}));
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
