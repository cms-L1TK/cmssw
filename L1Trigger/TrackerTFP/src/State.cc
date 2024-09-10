#include "L1Trigger/TrackerTFP/interface/State.h"

using namespace std;
using namespace tt;

namespace trackerTFP {

  // proto state constructor
  State::State(KalmanFilterFormats* formats,
               TrackCTB* track,
               const vector<vector<StubCTB*>>& stubs,
               const TTBV& maybePattern,
               int trackId)
      : formats_(formats),
        setup_(formats->setup()),
        track_(track),
        stubs_(stubs),
        maybePattern_(maybePattern),
        trackId_(trackId),
        parent_(nullptr),
        stub_(nullptr),
        hitPattern_(0, setup_->numLayers()),
        trackPattern_(0, setup_->numLayers()) {
    // first stub from first layer on input track with stubs
    for (int layer = setup_->numLayers() - 1; layer >= 0; layer--) {
      const vector<StubCTB*>& stubs = stubs_[layer];
      if (stubs.empty())
        continue;
      trackPattern_.set(layer);
      stub_ = stubs.front();
      layer_ = layer;
    }
    DataFormatKF& dfH00 = formats_->format(VariableKF::H00);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    // stub parameter
    H12_ = dfH00.digi(stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ()));
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

  // combinatoric state constructor
  State::State(State* state, StubCTB* stub, int layer) : State(state) {
    DataFormatKF& dfH00 = formats_->format(VariableKF::H00);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    parent_ = state->parent();
    stub_ = stub;
    layer_ = layer;
    H12_ = dfH00.digi(stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ()));
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

  // updated state constructor
  State::State(State* state, const vector<double>& doubles) : State(state) {
    DataFormatKF& dfH00 = formats_->format(VariableKF::H00);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    parent_ = state;
    // updated track parameter and uncertainties
    x0_ = doubles[0];
    x1_ = doubles[1];
    x2_ = doubles[2];
    x3_ = doubles[3];
    chi20_ = doubles[4];
    chi21_ = doubles[5];
    C00_ = doubles[6];
    C11_ = doubles[7];
    C22_ = doubles[8];
    C33_ = doubles[9];
    C01_ = doubles[10];
    C23_ = doubles[11];
    // update maps
    hitPattern_.set(layer_);
    // pick next stub (first stub in next layer with stub)
    if (hitPattern_.count() >= setup_->kfMinLayers()) {
      layer_ = 0;
      return;
    }
    stub_ = nullptr;
    for (int nextLayer = layer_ + 1; nextLayer < setup_->numLayers(); nextLayer++) {
      if (trackPattern_[nextLayer]) {
        stub_ = stubs_[nextLayer].front();
        layer_ = nextLayer;
        break;
      }
    }
    if (!stub_)
      return;
    H12_ = dfH00.digi(stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ()));
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

  // seed building state constructor
  State::State(State* state, int layer) : State(state) {
    DataFormatKF& dfH00 = formats_->format(VariableKF::H00);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    parent_ = state;
    hitPattern_.set(layer);
    for (int nextLayer = layer + 1; nextLayer < setup_->numLayers(); nextLayer++) {
      if (trackPattern_[nextLayer]) {
        stub_ = stubs_[nextLayer].front();
        layer_ = nextLayer;
        break;
      }
    }
    H12_ = dfH00.digi(stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ()));
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

  //
  State* State::update(deque<State>& states, int layer) {
    if (layer_ != layer || hitPattern_.count() == setup_->kfNumSeedStubs())
      return this;
    states.emplace_back(this, layer);
    return &states.back();
  }

  //
  State* State::combSeed(deque<State>& states, int layer) {
    // handle trivial state
    if (layer_ != layer || hitPattern_.count() == setup_->kfNumSeedStubs())
      return nullptr;
    // pick next stub on layer
    const vector<StubCTB*>& stubs = stubs_[layer];
    const int pos = distance(stubs.begin(), find(stubs.begin(), stubs.end(), stub_)) + 1;
    if (pos < (int)stubs.size()) {
      states.emplace_back(this, stubs[pos], layer);
      return &states.back();
    }
    // skip this layer
    if (trackPattern_[layer + 1] &&
        (hitPattern_.count() + trackPattern_.count(layer + 1, setup_->kfMaxSeedingLayer(), '1') >=
         setup_->kfNumSeedStubs()) &&
        trackPattern_.count() > setup_->kfMinLayers()) {
      states.emplace_back(this, stubs_[layer + 1].front(), layer + 1);
      return &states.back();
    }
    return nullptr;
  }

  //
  State* State::comb(deque<State>& states, int layer) {
    // handle skipping
    if (layer_ > layer)
      return nullptr;
    // handle max reached
    if (hitPattern_.count() == setup_->kfMaxLayers())
      return nullptr;
    // handle end reached
    if (!stub_)
      return nullptr;
    // handle min reached
    const vector<StubCTB*>& stubs = stubs_[layer];
    if (layer_ == 0) {
      if (trackPattern_[layer] && gapCheck(layer)) {
        states.emplace_back(this, stubs.front(), layer);
        return &states.back();
      }
      return nullptr;
    }
    // handle multiple stubs on layer
    const int pos = distance(stubs.begin(), find(stubs.begin(), stubs.end(), stub_)) + 1;
    if (pos < (int)stubs.size()) {
      states.emplace_back(this, stubs[pos], layer);
      return &states.back();
    }
    // handle skip
    int nextLayer = layer + 1;
    for (; nextLayer < setup_->numLayers(); nextLayer++)
      if (trackPattern_[nextLayer])
        break;
    if (gapCheck(nextLayer)) {
      states.emplace_back(this, stubs_[nextLayer].front(), nextLayer);
      return &states.back();
    }
    return nullptr;
  }

  //
  bool State::gapCheck(int layer) const {
    if (layer >= setup_->numLayers())
      return false;
    bool gap(false);
    bool doubleGap(false);
    int hits(0);
    int gaps(0);
    int available(0);
    for (int k = 0; k < layer; k++) {
      if (hitPattern_[k]) {
        hits++;
        gap = false;
      } else if (!maybePattern_[k]) {
        gaps++;
        if (gap)
          doubleGap = true;
        gap = true;
      }
    }
    for (int k = setup_->numLayers() - 1; k >= layer; k--)
      if (trackPattern_[k])
        available++;
    const int needed = setup_->kfMinLayers() - hits;
    if (doubleGap)
      return false;
    if (gaps > setup_->kfMaxGaps())
      return false;
    if (available < needed)
      return false;
    return true;
  }

  // copy constructor
  State::State(State* state)
      : formats_(state->formats_),
        setup_(state->setup_),
        track_(state->track_),
        stubs_(state->stubs_),
        maybePattern_(state->maybePattern_),
        trackId_(state->trackId_),
        parent_(state->parent_),
        stub_(state->stub_),
        layer_(state->layer_),
        hitPattern_(state->hitPattern_),
        trackPattern_(state->trackPattern_),
        x0_(state->x0_),
        x1_(state->x1_),
        x2_(state->x2_),
        x3_(state->x3_),
        chi20_(state->chi20_),
        chi21_(state->chi21_),
        C00_(state->C00_),
        C01_(state->C01_),
        C11_(state->C11_),
        C22_(state->C22_),
        C23_(state->C23_),
        C33_(state->C33_),
        H12_(state->H12_),
        v0_(state->v0_),
        v1_(state->v1_) {}

}  // namespace trackerTFP
