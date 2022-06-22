#include "L1Trigger/TrackerTFP/interface/State.h"

using namespace std;
using namespace tt;

namespace trackerTFP {

  // default constructor
  State::State(State* state)
      : dataFormats_(state->dataFormats_),
        setup_(state->setup_),
        track_(state->track_),
        trackId_(state->trackId_),
        parent_(state->parent_),
        stub_(state->stub_),
        layerMap_(state->layerMap_),
        hitPattern_(state->hitPattern_),
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
        numSkippedLayers_(state->numSkippedLayers_),
        numConsistentLayers_(state->numConsistentLayers_) {}

  // proto state constructor
  State::State(const DataFormats* dataFormats, TrackKFin* track, int trackId)
      : dataFormats_(dataFormats),
        setup_(dataFormats->setup()),
        track_(track),
        trackId_(trackId),
        parent_(nullptr),
        stub_(nullptr),
        layerMap_(setup_->numLayers()),
        hitPattern_(0, setup_->numLayers()),
        numSkippedLayers_(0),
        numConsistentLayers_(0) {
    // initial track parameter residuals w.r.t. found track
    x0_ = 0.;
    x1_ = 0.;
    x2_ = 0.;
    x3_ = 0.;
    // initial uncertainties
    C00_ = pow(dataFormats_->base(Variable::inv2R, Process::kfin), 2) - 1.e-12;
    C11_ = pow(dataFormats_->base(Variable::phiT, Process::kfin), 2) - 1.e-12;
    C22_ = pow(dataFormats_->base(Variable::cot, Process::kfin), 2) - 1.e-12;
    C33_ = pow(dataFormats_->base(Variable::zT, Process::kfin), 2) - 1.e-12;
    C01_ = 0.;
    C23_ = 0.;
    // first stub from first layer on input track with stubs
    stub_ = track->layerStub(track->hitPattern().plEncode());
  }

  // combinatoric state constructor
  State::State(State* state, StubKFin* stub) : State(state) {
    parent_ = state->parent();
    stub_ = stub;
  }

  // updated state constructor
  State::State(State* state, const std::vector<double>& doubles) : State(state) {
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
    const int layer = stub_->layer();
    hitPattern_.set(layer);
    const vector<StubKFin*>& stubs = track_->layerStubs(layer);
    layerMap_[layer] = distance(stubs.begin(), find(stubs.begin(), stubs.end(), stub_));
    // pick next stub (first stub in next layer with stub)
    stub_ = nullptr;
    if (hitPattern_.count() == setup_->kfMaxLayers())
      return;
    for (int nextLayer = layer + 1; nextLayer < setup_->numLayers(); nextLayer++) {
      if (track_->hitPattern(nextLayer)) {
        stub_ = track_->layerStub(nextLayer);
        break;
      }
    }
  }

}  // namespace trackerTFP
