#ifndef L1Trigger_TrackerTFP_KalmanFilter_h
#define L1Trigger_TrackerTFP_KalmanFilter_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "L1Trigger/TrackerTFP/interface/KalmanFilterFormats.h"
#include "L1Trigger/TrackerTFP/interface/State.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <deque>
#include <utility>

namespace trackerTFP {

  // Class to do helix fit to all tracks in a region.
  class KalmanFilter {
  public:
    KalmanFilter(const edm::ParameterSet& iConfig,
                 const tt::Setup* setup,
                 const DataFormats* dataFormats,
                 const LayerEncoding* layerEncoding,
                 KalmanFilterFormats* kalmanFilterFormats,
                 std::vector<TrackKF>& tracks,
                 std::vector<StubKF>& stubs);
    ~KalmanFilter() {}

    // fill output products
    void produce(const std::vector<std::vector<TrackCTB*>>& tracksIn,
                 const std::vector<std::vector<StubCTB*>>& stubsIn,
                 std::vector<std::vector<TrackKF*>>& tracksOut,
                 std::vector<std::vector<std::vector<StubKF*>>>& stubsOut,
                 int& numAcceptedStates,
                 int& numLostStates,
                 std::deque<std::pair<double, double>>& chi2s);

  private:
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;

    // apply final cuts
    void finalize(std::deque<State*>& stream);
    // Transform States into Tracks
    void conv(const std::deque<State*>& states,
              std::vector<TrackKF*>& tracks,
              std::vector<std::vector<StubKF*>>& stubs);
    // adds a layer to states
    void addLayer(std::deque<State*>& stream);
    // Assign next combinatoric (i.e. not first in layer) stub to state
    void comb(State*& state);
    // best state selection
    void accumulator(std::deque<State*>& stream);
    // updates state
    void update(State*& state);

    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // provides layer Encoding
    const LayerEncoding* layerEncoding_;
    // provides dataformats of Kalman filter internals
    KalmanFilterFormats* kalmanFilterFormats_;
    // container of output tracks
    std::vector<TrackKF>& tracks_;
    // container of output stubs
    std::vector<StubKF>& stubs_;
    // container of all Kalman Filter states
    std::deque<State> states_;
    // current layer used during state propagation
    int layer_;

    // dataformats of Kalman filter internals

    DataFormatKF* x0_;
    DataFormatKF* x1_;
    DataFormatKF* x2_;
    DataFormatKF* x3_;
    DataFormatKF* H00_;
    DataFormatKF* H12_;
    DataFormatKF* m0_;
    DataFormatKF* m1_;
    DataFormatKF* v0_;
    DataFormatKF* v1_;
    DataFormatKF* r0_;
    DataFormatKF* r1_;
    DataFormatKF* S00_;
    DataFormatKF* S01_;
    DataFormatKF* S12_;
    DataFormatKF* S13_;
    DataFormatKF* K00_;
    DataFormatKF* K10_;
    DataFormatKF* K21_;
    DataFormatKF* K31_;
    DataFormatKF* R00_;
    DataFormatKF* R11_;
    DataFormatKF* R00Rough_;
    DataFormatKF* R11Rough_;
    DataFormatKF* invR00Approx_;
    DataFormatKF* invR11Approx_;
    DataFormatKF* invR00Cor_;
    DataFormatKF* invR11Cor_;
    DataFormatKF* invR00_;
    DataFormatKF* invR11_;
    DataFormatKF* C00_;
    DataFormatKF* C01_;
    DataFormatKF* C11_;
    DataFormatKF* C22_;
    DataFormatKF* C23_;
    DataFormatKF* C33_;
    DataFormatKF* r02_;
    DataFormatKF* r12_;
    DataFormatKF* chi20_;
    DataFormatKF* chi21_;
  };

}  // namespace trackerTFP

#endif