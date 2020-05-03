#ifndef L1Trigger_TrackFindingTMTT_KalmanState_h
#define L1Trigger_TrackFindingTMTT_KalmanState_h

#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/KFbase.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <map>

///=== Represents helix state & last associated stub.
///=== All variable names & equations come from Fruhwirth KF paper
///=== http://dx.doi.org/10.1016/0168-9002%2887%2990887-4

namespace tmtt {

  class KFbase;
  class KalmanState;
  class Stub;
  class Settings;

  class KalmanState {
  public:
    KalmanState();
    KalmanState(const Settings *settings,
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
                double chi2rz);

    KalmanState(const KalmanState &p);
    ~KalmanState() {}

    KalmanState &operator=(const KalmanState &other);

    const Settings *settings() const { return settings_; }
    // Layer where next stub to extend this state should be sought.
    unsigned nextLayer() const { return (1 + kLayer_); }
    // Layer of last added stub. (-1 if no stubs yet).
    int layer() const { return kLayer_; }
    unsigned endcapRing() const { return endcapRing_; }
    bool barrel() const { return barrel_; }
    unsigned nSkippedLayers() const { return nSkipped_; }
    // Hit coordinates.
    double r() const { return r_; }
    double z() const { return z_; }
    const KalmanState *last_state() const { return last_state_; }
    // Helix parameters (1/2R, phi relative to sector, z0, tanLambda)
    const TVectorD &vectorX() const { return vecX_; }
    // Covariance matrix on helix params.
    const TMatrixD &matrixC() const { return matC_; }
    // Kalman Gain matrix
    const TMatrixD &matrixK() const { return matK_; }
    // Hit position covariance matrix.
    const TMatrixD &matrixV() const { return matV_; }
    // Last added stub
    const Stub *stub() const { return stub_; }
    // Track used to seed KF.
    const L1track3D &candidate() const { return l1track3D_; }

    double chi2() const { return chi2rphi_ + chi2rz_; }
    double chi2scaled() const { return chi2rphi_ / kalmanChi2RphiScale_ + chi2rz_; }  // Improves electron performance.
    double chi2rphi() const { return chi2rphi_; }
    double chi2rz() const { return chi2rz_; }
    unsigned nStubLayers() const { return n_stubs_; }
    unsigned int hitPattern() const { return hitPattern_; }  // Bit-encoded KF layers the fitted track has stubs in.

    bool good(const TP *tp) const;
    double reducedChi2() const;
    const KalmanState *last_update_state() const;
    std::vector<const Stub *> stubs() const;

    void setChi2(double chi2rphi, double chi2rz) {
      chi2rphi_ = chi2rphi;
      chi2rz_ = chi2rz;
    }

    // If using HLS, note/get additional output produced by HLS core.
    //void setHLSselect(unsigned int mBinHelix, unsigned int cBinHelix, bool consistent) { mBinHelixHLS_ = mBinHelix; cBinHelixHLS_ = cBinHelix; consistentHLS_ = consistent;}
    //void getHLSselect(unsigned int& mBinHelix, unsigned int& cBinHelix, bool& consistent) const { mBinHelix = mBinHelixHLS_; cBinHelix = cBinHelixHLS_; consistent = consistentHLS_;}

  private:
    const Settings *settings_;
    int kLayer_;
    unsigned endcapRing_;
    double r_;
    double z_;
    const KalmanState *last_state_;
    TVectorD vecX_;
    TMatrixD matC_;
    TMatrixD matK_;
    TMatrixD matV_;
    const Stub *stub_;
    double chi2rphi_;
    double chi2rz_;
    unsigned int kalmanChi2RphiScale_;
    unsigned n_stubs_;
    bool barrel_;
    unsigned nSkipped_;
    L1track3D l1track3D_;
    unsigned int hitPattern_;

    // Additional output from HLS if using it.
    unsigned int mBinHelixHLS_;
    unsigned int cBinHelixHLS_;
    bool consistentHLS_;
  };

}  // namespace tmtt

#endif
