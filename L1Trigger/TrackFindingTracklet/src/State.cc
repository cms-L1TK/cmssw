#include "L1Trigger/TrackFindingTracklet/interface/State.h"

#include <cmath>
#include <vector>
#include <deque>
#include <algorithm>
#include <iterator>

using namespace std;
using namespace tt;
using namespace trackerTFP;

namespace trklet {

  // proto state constructor
  State::State(KalmanFilterFormats* formats,
               TrackCTB* track,
               const vector<StubCTB*>& stubs,
               const TTBV& maybePattern,
               int trackId)
      : formats_(formats),
        setup_(formats->setup()),
        track_(track),
        maybePattern_(maybePattern),
        trackId_(trackId),
        parent_(nullptr),
        stub_(nullptr),
        layer_(-1),
        hitPattern_(0, setup_->numLayers()),
        trackPattern_(0, setup_->numLayers()),
        x0_(0.),
        x1_(0.),
        x2_(0.),
        x3_(0.),
        x4_(0.),
        C00_(9.e9),
        C01_(0.),
        C11_(9.e9),
        C22_(9.e9),
        C23_(0.),
        C33_(9.e9),
        C44_(pow(setup_->maxD0(), 2)),
        C40_(0.),
        C41_(0.) {
    stubs_ = vector<vector<StubCTB*>>(stubs.size());
    for (int layer = 0; layer < setup_->numLayers(); layer++) {
      StubCTB* stub = stubs[layer];
      if (stub) {
        stubs_[layer].push_back(stubs[layer]);
        trackPattern_.set(layer);
        if (!stub_) {
          stub_ = stub;
          layer_ = layer;
        }
      }
    }
    DataFormatKF& dfH00 = formats_->format(VariableKF::H00);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    // stub parameter
    H12_ = dfH00.digi(stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ()));
    H04_ = stub_->r() + setup_->chosenRofPhi();
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
    H12_ = stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
    H04_ = stub_->r() + setup_->chosenRofPhi();
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
    x4_ = doubles[4];
    chi20_ = doubles[5];
    chi21_ = doubles[6];
    C00_ = doubles[7];
    C11_ = doubles[8];
    C22_ = doubles[9];
    C33_ = doubles[10];
    C01_ = doubles[11];
    C23_ = doubles[12];
    C44_ = doubles[13];
    C40_ = doubles[14];
    C41_ = doubles[15];
    // update maps
    hitPattern_.set(layer_);
    // pick next stub (first stub in next layer with stub)
    if (hitPattern_.count() >= setup_->kfMinLayers()) {
      layer_ = -1;
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
    H04_ = stub_->r() + setup_->chosenRofPhi();
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
    H04_ = stub_->r() + setup_->chosenRofPhi();
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
    // skip layers
    for (int nextLayer = layer + 1; nextLayer < setup_->kfMaxSeedingLayer(); nextLayer++) {
      if (!trackPattern_[nextLayer])
        continue;
      const int maxSeedStubs = hitPattern_.count() + trackPattern_.count(nextLayer, setup_->kfMaxSeedingLayer(), '1');
      if (maxSeedStubs < setup_->kfNumSeedStubs())
        continue;
      const int maxStubs = maxSeedStubs + trackPattern_.count(setup_->kfMaxSeedingLayer(), setup_->numLayers(), '1');
      if (maxStubs < setup_->kfMinLayers())
        continue;
      states.emplace_back(this, stubs_[nextLayer].front(), nextLayer);
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
    if (layer_ == -1) {
      if (trackPattern_[layer]) {
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
    // not enough layer left
    if (hitPattern_.count() + trackPattern_.count(nextLayer, setup_->numLayers()) < setup_->kfMinLayers())
      return nullptr;
    if (nextLayer < setup_->numLayers() && trackPattern_[nextLayer]) {
      states.emplace_back(this, stubs_[nextLayer].front(), nextLayer);
      return &states.back();
    }
    return nullptr;
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
        x4_(state->x4_),
        chi20_(state->chi20_),
        chi21_(state->chi21_),
        C00_(state->C00_),
        C01_(state->C01_),
        C11_(state->C11_),
        C22_(state->C22_),
        C23_(state->C23_),
        C33_(state->C33_),
        C44_(state->C44_),
        C40_(state->C40_),
        C41_(state->C41_),
        H12_(state->H12_),
        H04_(state->H04_),
        v0_(state->v0_),
        v1_(state->v1_) {}

}  // namespace trklet
