#include "L1Trigger/TrackFindingTMTT/interface/KalmanState.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include <TMatrixD.h>

using namespace std;

namespace tmtt {

  KalmanState::KalmanState()
      : settings_(nullptr),
        kLayer_(0),
        last_state_(nullptr),
        stub_(nullptr),
        chi2rphi_(0),
        chi2rz_(0),
        barrel_(true),
        nSkipped_(0),
        hitPattern_(0) {}

  KalmanState::KalmanState(const Settings* settings,
			   const L1track3D &candidate,
                           unsigned nSkipped,
                           int kLayer,
                           const KalmanState *last_state,
                           const TVectorD &vecX,
                           const TMatrixD &matC,
                           const TMatrixD &matK,
                           const TMatrixD &matV,
                           const Stub *stub,
                           double chi2rphi,
                           double chi2rz)
      : settings_(settings),
        kLayer_(kLayer),
        last_state_(last_state),
        vecX_(vecX),
        stub_(stub),
        chi2rphi_(chi2rphi),
        chi2rz_(chi2rz),
        nSkipped_(nSkipped),
        l1track3D_(candidate) {
    matC_.Clear();
    matC_.ResizeTo(matC.GetNrows(), matC.GetNcols());
    matC_ = matC;
    matK_.ResizeTo(matK.GetNrows(), matK.GetNcols());
    matK_ = matK;
    matV_.ResizeTo(matV.GetNrows(), matV.GetNcols());
    matV_ = matV;
    kalmanChi2RphiScale_ = settings_->kalmanChi2RphiScale();

    hitPattern_ = 0;
    if (last_state != nullptr)
      hitPattern_ = last_state->hitPattern();  // Bit encoded list of hit layers
    if (stub != nullptr && kLayer_ >= 0)
      hitPattern_ |= (1 << (kLayer_));

    r_ = 0.1;
    z_ = 0;
    barrel_ = true;
    endcapRing_ = 0;

    if (stub != nullptr) {
      r_ = stub->r();
      z_ = stub->z();
      barrel_ = stub->barrel();
      endcapRing_ = stub->endcapRing();
    }

    n_stubs_ = 1 + kLayer_ - nSkipped_;
  }

  KalmanState::KalmanState(const KalmanState &p)
      : settings_(p.settings()),
        kLayer_(p.layer()),
        endcapRing_(p.endcapRing()),
        r_(p.r()),
        z_(p.z()),
        last_state_(p.last_state()),
        vecX_(p.vectorX()),
        matC_(p.matrixC()),
        matK_(p.matrixK()),
        matV_(p.matrixV()),
        stub_(p.stub()),
        chi2rphi_(p.chi2rphi()),
        chi2rz_(p.chi2rz()),
        n_stubs_(p.nStubLayers()),
        barrel_(p.barrel()),
        nSkipped_(p.nSkippedLayers()),
        l1track3D_(p.candidate()) {}

  KalmanState &KalmanState::operator=(const KalmanState &other) {
    if (&other == this)
      return *this;

    settings_ = other.settings();
    kLayer_ = other.layer();
    endcapRing_ = other.endcapRing();
    r_ = other.r();
    z_ = other.z();
    last_state_ = other.last_state();
    vecX_ = other.vectorX();
    matC_ = other.matrixC();
    matK_ = other.matrixK();
    matV_ = other.matrixV();
    stub_ = other.stub();
    chi2rphi_ = other.chi2rphi();
    chi2rz_ = other.chi2rz();
    n_stubs_ = other.nStubLayers();
    barrel_ = other.barrel();
    nSkipped_ = other.nSkippedLayers();
    l1track3D_ = other.candidate();
    return *this;
  }

  bool KalmanState::good(const TP *tp) const {
    const KalmanState *state = this;
    while (state) {
      const Stub *stub = state->stub();
      if (stub != nullptr) {
        const set<const TP *> &tps = stub->assocTPs();
        if (tps.find(tp) == tps.end())
          return false;
      }
      state = state->last_state();
    }
    return true;
  }

  double KalmanState::reducedChi2() const {
    if (2 * n_stubs_ - vecX_.GetNrows() > 0)
      return (this->chi2()) / (2 * n_stubs_ - vecX_.GetNrows());
    else
      return 0;
  }

  const KalmanState *KalmanState::last_update_state() const {
    const KalmanState *state = this;
    while (state) {
      if (state->stub() != nullptr)
        return state;
      state = state->last_state();
    }
    return 0;
  }

  std::vector<const Stub *> KalmanState::stubs() const {
    std::vector<const Stub *> all_stubs;

    const KalmanState *state = this;
    while (state) {
      const Stub *stub = state->stub();
      if (stub != nullptr) 
        all_stubs.push_back(stub);
      state = state->last_state();
    }
    std::reverse(all_stubs.begin(), all_stubs.end());  // Put innermost stub first.
    return all_stubs;
  }

}  // namespace tmtt
