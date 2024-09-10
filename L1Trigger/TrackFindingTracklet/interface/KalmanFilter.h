#ifndef L1Trigger_TrackFindingTracklet_KalmanFilter_h
#define L1Trigger_TrackFindingTracklet_KalmanFilter_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "L1Trigger/TrackerTFP/interface/KalmanFilterFormats.h"
#include "L1Trigger/TrackFindingTracklet/interface/State.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <deque>

namespace trklet {

  /*! \class  trklet::KalmanFilter
   *  \brief  Class to do helix fit to all tracks in a region.
   *          All variable names & equations come from Fruhwirth KF paper
   *          http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
   *          Summary of variables:
   *          m = hit position (phi,z)
   *          V = hit position 2x2 covariance matrix in (phi,z).
   *          x = helix params
   *          C = helix params 4x4 covariance matrix
   *          r = residuals
   *          H = 2x4 derivative matrix (expected stub position w.r.t. helix params)
   *          K = KF gain 2x2 matrix
   *          x' & C': Updated values of x & C after KF iteration
   *          Boring: F = unit matrix; pxcov = C
   *          Summary of equations:
   *          S = H*C (2x4 matrix); St = Transpose S
   *          R = V + H*C*Ht (KF paper) = V + H*St (used here at simpler): 2x2 matrix
   *          Rinv = Inverse R
   *          K = St * Rinv : 2x2 Kalman gain matrix * det(R)
   *          r = m - H*x
   *          x' = x + K*r
   *          C' = C - K*H*C (KF paper) = C - K*S (used here as simpler)
   *  \author Thomas Schuh
   *  \date   2024, Sep
   */
  class KalmanFilter {
  public:
    KalmanFilter(const edm::ParameterSet& iConfig,
                 const tt::Setup* setup,
                 const trackerTFP::DataFormats* dataFormats,
                 const trackerTFP::LayerEncoding* layerEncoding,
                 trackerTFP::KalmanFilterFormats* kalmanFilterFormats,
                 int region,
                 tt::TTTracks& ttTracks);
    ~KalmanFilter() {}
    // read in and organize input tracks and stubs
    void consume(const tt::StreamsTrack& streamsTrack, const tt::StreamsStub& streamsStub);
    // fill output products
    void produce(tt::StreamsStub& streamsStub,
                 tt::StreamsTrack& streamsTrack,
                 int& numAcceptedStates,
                 int& numLostStates);

  private:
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;
    //
    double digi(double val, double base) const { return (floor(val / base + 1.e-11) + .5) * base; }

    // calculates the helix params & their cov. matrix from a pair of stubs
    void calcSeeds();
    // apply final cuts
    void finalize();
    // Transform States into output products
    void conv(tt::StreamsStub& streamsStub, tt::StreamsTrack& streamsTrack);
    // adds a layer to states
    void addLayer();
    // adds a layer to states to build seeds
    void addSeedLayer();
    // Assign next combinatoric (i.e. not first in layer) stub to state
    void comb(State*& state);
    // best state selection
    void accumulator();
    // updates state
    void update(State*& state) { use5ParameterFit_ ? update5(state) : update4(state); }
    // updates state using 4 paramter fit
    void update4(State*& state);
    // updates state using 5 parameter fit
    void update5(State*& state);

    // true if truncation is enbaled
    bool enableTruncation_;
    //
    bool use5ParameterFit_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const trackerTFP::DataFormats* dataFormats_;
    // provides layer Encoding
    const trackerTFP::LayerEncoding* layerEncoding_;
    // provides dataformats of Kalman filter internals
    trackerTFP::KalmanFilterFormats* kalmanFilterFormats_;
    // processing region
    int region_;
    //
    tt::TTTracks& ttTracks_;
    // container of tracks
    std::vector<trackerTFP::TrackCTB> tracks_;
    // container of stubs
    std::vector<trackerTFP::StubCTB> stubs_;
    // container of all Kalman Filter states
    std::deque<State> states_;
    // processing stream
    std::deque<State*> stream_;
    // current layer used during state propagation
    int layer_;

    // dataformats of Kalman filter internals

    trackerTFP::DataFormatKF* x0_;
    trackerTFP::DataFormatKF* x1_;
    trackerTFP::DataFormatKF* x2_;
    trackerTFP::DataFormatKF* x3_;
    trackerTFP::DataFormatKF* H00_;
    trackerTFP::DataFormatKF* H12_;
    trackerTFP::DataFormatKF* m0_;
    trackerTFP::DataFormatKF* m1_;
    trackerTFP::DataFormatKF* v0_;
    trackerTFP::DataFormatKF* v1_;
    trackerTFP::DataFormatKF* r0_;
    trackerTFP::DataFormatKF* r1_;
    trackerTFP::DataFormatKF* S00_;
    trackerTFP::DataFormatKF* S01_;
    trackerTFP::DataFormatKF* S12_;
    trackerTFP::DataFormatKF* S13_;
    trackerTFP::DataFormatKF* K00_;
    trackerTFP::DataFormatKF* K10_;
    trackerTFP::DataFormatKF* K21_;
    trackerTFP::DataFormatKF* K31_;
    trackerTFP::DataFormatKF* R00_;
    trackerTFP::DataFormatKF* R11_;
    trackerTFP::DataFormatKF* R00Rough_;
    trackerTFP::DataFormatKF* R11Rough_;
    trackerTFP::DataFormatKF* invR00Approx_;
    trackerTFP::DataFormatKF* invR11Approx_;
    trackerTFP::DataFormatKF* invR00Cor_;
    trackerTFP::DataFormatKF* invR11Cor_;
    trackerTFP::DataFormatKF* invR00_;
    trackerTFP::DataFormatKF* invR11_;
    trackerTFP::DataFormatKF* C00_;
    trackerTFP::DataFormatKF* C01_;
    trackerTFP::DataFormatKF* C11_;
    trackerTFP::DataFormatKF* C22_;
    trackerTFP::DataFormatKF* C23_;
    trackerTFP::DataFormatKF* C33_;
  };

}  // namespace trklet

#endif