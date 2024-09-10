#ifndef L1Trigger_TrackFindingTracklet_State_h
#define L1Trigger_TrackFindingTracklet_State_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/KalmanFilterFormats.h"

#include <vector>
#include <numeric>

namespace trklet {

  // Class to represent a Kalman Filter State
  class State {
  public:
    // copy constructor
    State(State* state);
    // proto state constructor
    State(trackerTFP::KalmanFilterFormats* formats,
          trackerTFP::TrackCTB* track,
          const std::vector<trackerTFP::StubCTB*>& stubs,
          const TTBV& maybe,
          int trackId);
    // updated state constructor
    State(State* state, const std::vector<double>& doubles);
    // combinatoric state constructor
    State(State* state, trackerTFP::StubCTB* stub, int layer);
    // seed building state constructor
    State(State* state, int layer);
    ~State() {}
    //
    State* comb(std::deque<State>& states, int layer);
    //
    State* combSeed(std::deque<State>& states, int layer);
    //
    State* update(std::deque<State>& states, int layer);
    // input track
    trackerTFP::TrackCTB* track() const { return track_; }
    // parent state (nullpointer if no parent available)
    State* parent() const { return parent_; }
    // stub to add to state
    trackerTFP::StubCTB* stub() const { return stub_; }
    // hitPattern of so far added stubs
    const TTBV& hitPattern() const { return hitPattern_; }
    // shows which layer the found track has stubs on
    const TTBV& trackPattern() const { return trackPattern_; }
    // track id of input track
    int trackId() const { return trackId_; }
    // pattern of maybe layers for input track
    const TTBV& maybePattern() const { return maybePattern_; }
    // layer id of the current stub to add
    int layer() const { return layer_; }
    // helix inv2R wrt input helix
    double x0() const { return x0_; }
    // helix phi at radius ChosenRofPhi wrt input helix
    double x1() const { return x1_; }
    // helix cot(Theta) wrt input helix
    double x2() const { return x2_; }
    // helix z at radius chosenRofZ wrt input helix
    double x3() const { return x3_; }
    //
    double x4() const { return x4_; }
    //
    double chi20() const { return chi20_; }
    //
    double chi21() const { return chi21_; }
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
    double C44() const { return C44_; }
    double C40() const { return C40_; }
    double C41() const { return C41_; }
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofPhi
    double H00() const { return stub_->r(); }
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofZ
    double H12() const { return H12_; }
    //
    double H04() const { return H04_; }
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
    //const std::vector<trackerTFP::StubCTB*>& stubs() const { return stubs_; }

  private:
    // provides data fomats
    trackerTFP::KalmanFilterFormats* formats_;
    // provides run-time constants
    const tt::Setup* setup_;
    // input track
    trackerTFP::TrackCTB* track_;
    // input track stubs
    std::vector<std::vector<trackerTFP::StubCTB*>> stubs_;
    // pattern of maybe layers for input track
    TTBV maybePattern_;
    // track id
    int trackId_;
    // previous state, nullptr for first states
    State* parent_;
    // stub to add
    trackerTFP::StubCTB* stub_;
    // layer id of the current stub to add
    int layer_;
    // shows which layer has been added so far
    TTBV hitPattern_;
    // shows which layer the found track has stubs on
    TTBV trackPattern_;
    // helix inv2R wrt input helix
    double x0_;
    // helix phi at radius ChosenRofPhi wrt input helix
    double x1_;
    // helix cot(Theta) wrt input helix
    double x2_;
    // helix z at radius chosenRofZ wrt input helix
    double x3_;
    // impact parameter in 1/cm
    double x4_;
    //
    double chi20_;
    //
    double chi21_;
    // cov. matrix
    double C00_;
    double C01_;
    double C11_;
    double C22_;
    double C23_;
    double C33_;
    double C44_;
    double C40_;
    double C41_;
    // Derivative of predicted stub coords wrt helix params: stub radius minus chosenRofZ
    double H12_;
    double H04_;
    // squared stub projected phi uncertainty instead of wheight (wrong but simpler)
    double v0_;
    // squared stub projected z uncertainty instead of wheight (wrong but simpler)
    double v1_;
  };

}  // namespace trklet

#endif