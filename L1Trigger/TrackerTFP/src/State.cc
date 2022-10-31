#include "L1Trigger/TrackerTFP/interface/State.h"

using namespace std;
using namespace tt;

namespace trackerTFP {

  // default constructor
  State::State(State* state)
      : formats_(state->formats_),
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
        H12_(state->H12_),
        v0_(state->v0_),
        v1_(state->v1_) {}

  // proto state constructor
  State::State(KalmanFilterFormats* formats, TrackKFin* track, int trackId)
      : formats_(formats),
        setup_(formats->setup()),
        track_(track),
        trackId_(trackId),
        parent_(nullptr),
        stub_(nullptr),
        layerMap_(setup_->numLayers()),
        hitPattern_(0, setup_->numLayers()) {
    const DataFormats* dataFormats = formats_->dataFormats();
    const DataFormat& dfInv2R = dataFormats->format(Variable::inv2R, Process::ht);
    const DataFormat& dfPhiT = dataFormats->format(Variable::phiT, Process::ht);
    const DataFormat& dfCot = dataFormats->format(Variable::cot, Process::gp);
    const DataFormat& dfZT = dataFormats->format(Variable::zT, Process::gp);
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
    DataFormatKF& dfH12 = formats_->format(VariableKF::H12);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    // initial track parameter residuals w.r.t. found track
    x0_ = dfX0.digi(0.);
    x1_ = dfX1.digi(0.);
    x2_ = dfX2.digi(0.);
    x3_ = dfX3.digi(0.);
    // initial uncertainties
    C00_ = dfC00.digi(pow(dfInv2R.base(), 2) * pow(2, setup_->kfShiftInitialC00()) - 1.e-12);
    C11_ = dfC11.digi(pow(dfPhiT.base(), 2) * pow(2, setup_->kfShiftInitialC11()) - 1.e-12);
    C22_ = dfC22.digi(pow(dfCot.base(), 2) * pow(2, setup_->kfShiftInitialC22()) - 1.e-12);
    C33_ = dfC33.digi(pow(dfZT.base(), 2) * pow(2, setup_->kfShiftInitialC33()) - 1.e-12);
    C01_ = dfC01.digi(0.);
    C23_ = dfC23.digi(0.);
    // first stub from first layer on input track with stubs
    stub_ = track->layerStub(track->hitPattern().plEncode());
    H12_ = dfH12.digi(stub_->r() + setup_->chosenRofPhi() - setup_->chosenRofZ());
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

  // combinatoric state constructor
  State::State(State* state, StubKFin* stub) : State(state) {
    DataFormatKF& dfH12 = formats_->format(VariableKF::H12);
    DataFormatKF& dfv0 = formats_->format(VariableKF::v0);
    DataFormatKF& dfv1 = formats_->format(VariableKF::v1);
    parent_ = state->parent();
    stub_ = stub;
    if (!stub)
      return;
    H12_ = dfH12.digi(stub_->r() + setup_->chosenRofPhi() - setup_->chosenRofZ());
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

  // updated state constructor
  State::State(State* state, const std::vector<double>& doubles) : State(state) {
    DataFormatKF& dfH12 = formats_->format(VariableKF::H12);
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
    if (!stub_)
      return;
    H12_ = dfH12.digi(stub_->r() + setup_->chosenRofPhi() - setup_->chosenRofZ());
    v0_ = dfv0.digi(pow(stub_->dPhi(), 2));
    v1_ = dfv1.digi(pow(stub_->dZ(), 2));
  }

}  // namespace trackerTFP
