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

  // Class to do helix fit to all tracks in a region.
  class KalmanFilter {
  public:
    KalmanFilter(const edm::ParameterSet& iConfig,
                 const tt::Setup* setup,
                 const trackerTFP::DataFormats* dataFormats,
                 const trackerTFP::LayerEncoding* layerEncoding,
                 trackerTFP::KalmanFilterFormats* kalmanFilterFormats,
                 std::vector<trackerTFP::TrackKF>& tracks,
                 std::vector<trackerTFP::StubKF>& stubs);
    ~KalmanFilter() {}

    // fill output products
    void produce(const std::vector<std::vector<trackerTFP::TrackCTB*>>& tracksIn,
                 const std::vector<std::vector<trackerTFP::StubCTB*>>& stubsIn,
                 std::vector<std::vector<trackerTFP::TrackKF*>>& tracksOut,
                 std::vector<std::vector<std::vector<trackerTFP::StubKF*>>>& stubsOut,
                 int& numAcceptedStates,
                 int& numLostStates);

  private:
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;

    // apply final cuts
    void finalize(std::deque<State*>& stream);
    // Transform States into Tracks
    void conv(const std::deque<State*>& states,
              std::vector<trackerTFP::TrackKF*>& tracks,
              std::vector<std::vector<trackerTFP::StubKF*>>& stubs);
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
    const trackerTFP::DataFormats* dataFormats_;
    // provides layer Encoding
    const trackerTFP::LayerEncoding* layerEncoding_;
    // provides dataformats of Kalman filter internals
    trackerTFP::KalmanFilterFormats* kalmanFilterFormats_;
    // container of output tracks
    std::vector<trackerTFP::TrackKF>& tracks_;
    // container of output stubs
    std::vector<trackerTFP::StubKF>& stubs_;
    // container of all Kalman Filter states
    std::deque<State> states_;
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