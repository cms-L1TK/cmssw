#ifndef L1Trigger_TrackerTFP_State_h
#define L1Trigger_TrackerTFP_State_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/KalmanFilterFormats.h"

#include <vector>
#include <numeric>

namespace trackerTFP {

  // Class to represent a Kalman Filter State
  class State {
  public:
    // default constructor
    State(State* state);
    // proto state constructor
    State(KalmanFilterFormats* formats, TrackKFin* track, int trackId);
    // combinatoric state constructor
    State(State* state, StubKFin* stub);
    // updated state constructor
    State(State* state, const std::vector<double>& doubles);
    ~State() {}
    // input track
    TrackKFin* track() const { return track_; }
    // parent state (nullpointer if no parent available)
    State* parent() const { return parent_; }
    // stub to add to state
    StubKFin* stub() const { return stub_; }
    // hitPattern of so far added stubs
    const TTBV& hitPattern() const { return hitPattern_; }
    // track id of input track
    int trackId() const { return trackId_; }
    // pattern of maybe layers for input track
    const TTBV& maybePattern() const { return track_->maybePattern(); }
    // layer id of the current stub to add
    int layer() const { return (stub_ ? stub_->layer() : -1); }
    // helix inv2R wrt input helix
    double x0() const { return x0_; }
    // helix phi at radius ChosenRofPhi wrt input helix
    double x1() const { return x1_; }
    // helix cot(Theta) wrt input helix
    double x2() const { return x2_; }
    // helix z at radius chosenRofZ wrt input helix
    double x3() const { return x3_; }
    // cov. matrix element
    double C00() const { return C00_; }
    // cov. matrix element
    double C01() const { return C01_; }
    // cov. matrix element
    double C11() const { return C11_; }
    // cov. matrix element
    double C22() const { return C22_; }
    // cov. matrix element
    double C23() const { return C23_; }
    // cov. matrix element
    double C33() const { return C33_; }
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofPhi
    double H00() const { return stub_->r(); }
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofZ
    double H12() const { return H12_; }
    // stub phi residual wrt input helix
    double m0() const { return stub_->phi(); }
    // stub z residual wrt input helix
    double m1() const { return stub_->z(); }
    // stub projected phi uncertainty
    double dPhi() const { return stub_->dPhi(); }
    // stub projected z uncertainty
    double dZ() const { return stub_->dZ(); }
    // squared stub projected phi uncertainty instead of wheight (wrong but simpler)
    double v0() const { return v0_; }
    // squared stub projected z uncertainty instead of wheight (wrong but simpler)
    double v1() const { return v1_; }
  private:
    // provides data fomats
    KalmanFilterFormats* formats_;
    // provides run-time constants
    const tt::Setup* setup_;
    // input track
    TrackKFin* track_;
    // track id
    int trackId_;
    // previous state, nullptr for first states
    State* parent_;
    // stub to add
    StubKFin* stub_;
    // shows which stub on each layer has been added so far
    std::vector<int> layerMap_;
    // shows which layer has been added so far
    TTBV hitPattern_;
    // helix inv2R wrt input helix
    double x0_;
    // helix phi at radius ChosenRofPhi wrt input helix
    double x1_;
    // helix cot(Theta) wrt input helix
    double x2_;
    // helix z at radius chosenRofZ wrt input helix
    double x3_;
    // cov. matrix
    double C00_;
    double C01_;
    double C11_;
    double C22_;
    double C23_;
    double C33_;
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofZ
    double H12_;
    // squared stub projected phi uncertainty instead of wheight (wrong but simpler)
    double v0_;
    // squared stub projected z uncertainty instead of wheight (wrong but simpler)
    double v1_;
  };

}  // namespace trackerTFP

#endif