#include "L1Trigger/TrackFindingTracklet/interface/State.h"

using namespace std;
using namespace tt;
using namespace trackerTFP;

namespace trklet {

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
        layer_(-1),
        hitPattern_(0, setup_->numLayers()),
        trackPattern_(0, setup_->numLayers()),
        type_(combSkip) {
    DataFormatKF& dfX0 = formats_->format(VariableKF::x0);
    DataFormatKF& dfX1 = formats_->format(VariableKF::x1);
    DataFormatKF& dfX2 = formats_->format(VariableKF::x2);
    DataFormatKF& dfX3 = formats_->format(VariableKF::x3);
    DataFormatKF& dfC00 = formats_->format(VariableKF::C00);
    DataFormatKF& dfC11 = formats_->format(VariableKF::C11);
    DataFormatKF& dfC22 = formats_->format(VariableKF::C22);
    DataFormatKF& dfC33 = formats_->format(VariableKF::C33);
    DataFormatKF& dfC01 = formats_->format(VariableKF::C01);
    DataFormatKF& dfC23 = formats_->format(VariableKF::C23);
    DataFormatKF& dfH00 = formats_->format(VariableKF::H00);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    // initial track parameter residuals w.r.t. found track
    x0_ = dfX0.digi(0.);
    x1_ = dfX1.digi(0.);
    x2_ = dfX2.digi(0.);
    x3_ = dfX3.digi(0.);
    // initial uncertainties
    C00_ = dfC00.range() - .5 * dfC00.base();
    C11_ = dfC11.range() - .5 * dfC11.base();
    C22_ = dfC22.range() - .5 * dfC22.base();
    C33_ = dfC33.range() - .5 * dfC33.base();
    C01_ = dfC01.digi(0.);
    C23_ = dfC23.digi(0.);
    // first stub from first layer on input track with stubs
    for (int layer = setup_->numLayers() - 1; layer >= 0; layer--) {
      const vector<StubCTB*>& stubs = stubs_[layer];
      if (stubs.empty())
        continue;
      trackPattern_.set(layer);
      stub_ = stubs.front();
      layer_ = layer;
    }
    H12_ = stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
    if (layer_ > 0)
      type_ = skip;
  }

  // combinatoric state constructor
  State::State(State* state, StubCTB* stub, int layer, Type type) : State(state) {
    DataFormatKF& dfH00 = formats_->format(VariableKF::H00);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    parent_ = state->parent();
    stub_ = stub;
    layer_ = layer;
    type_ = type;
    H12_ = stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
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
    C00_ = doubles[4];
    C11_ = doubles[5];
    C22_ = doubles[6];
    C33_ = doubles[7];
    C01_ = doubles[8];
    C23_ = doubles[9];
    // update maps
    hitPattern_.set(layer_);
    if (hitPattern_.count() >= setup_->kfMinLayers())
      type_ = done;
    // pick next stub (first stub in next layer with stub)
    stub_ = nullptr;
    if (type_ == done)
      return;
    for (int nextLayer = layer_ + 1; nextLayer < setup_->numLayers(); nextLayer++) {
      if (trackPattern_[nextLayer]) {
        stub_ = stubs_[nextLayer].front();
        if (nextLayer > layer_ + 1)
          type_ = skip;
        layer_ = nextLayer;
        break;
      }
    }
    if (!stub_)
      return;
    H12_ = stub_->r() + dfH00.digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

  //
  State* State::comb(deque<State>& states, int layer) {
    // handle trivial state
    if (type_ == skip)
      return nullptr;
    // prepare others state handling
    //bool pre(true);
    //bool post(false);
    //bool doubleGap(false);
    int hits(0);
    int gaps(0);
    int available(0);
    for (int k = 0; k < layer; k++) {
      if (hitPattern_[k]) {
        hits++;
        //pre = true;
      } else if (!maybePattern_[k]) {
        gaps++;
        //if (layer > 0 && !pre)
          //doubleGap = true;
        //pre = false;
      }
    }
    for (int k = setup_->numLayers() - 1; k > layer; k--) {
      if (trackPattern_[k]) {
        available++;
        //post = true;
      } //else if (!maybePattern_[k])
        //post = false;
    }
    const int needed = setup_->kfMinLayers() - hits;
    const vector<StubCTB*>& stubs = stubs_[layer];
    const int pos = distance(stubs.begin(), find(stubs.begin(), stubs.end(), stub_)) + 1;
    // non trivial state handling
    if (type_ == done) {
      // pick first stub on layer if available and gap criteria by adding this still met
      //if (trackPattern_[layer] && gaps <= setup_->kfMaxGaps() && !doubleGap) {
      if (trackPattern_[layer] && gaps <= setup_->kfMaxGaps()) {
        states.emplace_back(this, stubs.front(), layer, combDone);
        return &states.back();
      }
    } else if (type_ == combDone) {
      // pick next stub on layer if available
      if (pos < (int)stubs.size()) {
        states.emplace_back(this, stubs[pos], layer, combDone);
        return &states.back();
      }
    } else if (type_ == combSkip) {
      // pick next stub on layer if available
      if (pos < (int)stubs.size()) {
        states.emplace_back(this, stubs[pos], layer, combSkip);
        return &states.back();
      //} else if (available > 0 && available >= needed && gaps < setup_->kfMaxGaps() && pre && post) {
      } else if (available > 0 && available >= needed && gaps < setup_->kfMaxGaps()) {
        // pick first stub on next layer with stubs
        StubCTB* stub = nullptr;
        int nextLayer = layer + 1;
        for (; nextLayer < setup_->numLayers(); nextLayer++) {
          if (trackPattern_[nextLayer]) {
            stub = stubs_[nextLayer].front();
            break;
          }
        }
        states.emplace_back(this, stub, nextLayer, skip);
        return &states.back();
      }
    }
    return nullptr;
  }

  //
  State* State::unskip(deque<State>& states, int layer) {
    states.emplace_back(this, stubs_[layer].front(), layer, combSkip);
    return &states.back();
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
        C00_(state->C00_),
        C01_(state->C01_),
        C11_(state->C11_),
        C22_(state->C22_),
        C23_(state->C23_),
        C33_(state->C33_),
        H12_(state->H12_),
        v0_(state->v0_),
        v1_(state->v1_),
        type_(state->type_) {}

}  // namespace trklet
