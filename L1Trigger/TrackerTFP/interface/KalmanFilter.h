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

    std::vector<DataFormatKF*> x0_;
    std::vector<DataFormatKF*> x1_;
    std::vector<DataFormatKF*> x2_;
    std::vector<DataFormatKF*> x3_;
    std::vector<DataFormatKF*> H00_;
    std::vector<DataFormatKF*> H12_;
    std::vector<DataFormatKF*> m0_;
    std::vector<DataFormatKF*> m1_;
    std::vector<DataFormatKF*> v0_;
    std::vector<DataFormatKF*> v1_;
    std::vector<DataFormatKF*> r0_;
    std::vector<DataFormatKF*> r1_;
    std::vector<DataFormatKF*> S00_;
    std::vector<DataFormatKF*> S01_;
    std::vector<DataFormatKF*> S12_;
    std::vector<DataFormatKF*> S13_;
    std::vector<DataFormatKF*> K00_;
    std::vector<DataFormatKF*> K10_;
    std::vector<DataFormatKF*> K21_;
    std::vector<DataFormatKF*> K31_;
    std::vector<DataFormatKF*> R00_;
    std::vector<DataFormatKF*> R11_;
    std::vector<DataFormatKF*> R00Rough_;
    std::vector<DataFormatKF*> R11Rough_;
    std::vector<DataFormatKF*> invR00Approx_;
    std::vector<DataFormatKF*> invR11Approx_;
    std::vector<DataFormatKF*> invR00Cor_;
    std::vector<DataFormatKF*> invR11Cor_;
    std::vector<DataFormatKF*> invR00_;
    std::vector<DataFormatKF*> invR11_;
    std::vector<DataFormatKF*> C00_;
    std::vector<DataFormatKF*> C01_;
    std::vector<DataFormatKF*> C11_;
    std::vector<DataFormatKF*> C22_;
    std::vector<DataFormatKF*> C23_;
    std::vector<DataFormatKF*> C33_;
    std::vector<DataFormatKF*> r02_;
    std::vector<DataFormatKF*> r12_;
    std::vector<DataFormatKF*> chi20_;
    std::vector<DataFormatKF*> chi21_;
  };

}  // namespace trackerTFP

#endif