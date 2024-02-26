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
                 int region);
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
    // apply final cuts
    void finalize();
    // Transform States into output products
    void conv(tt::StreamsStub& streamsStub, tt::StreamsTrack& streamsTrack) const;
    // adds a layer to states
    void addLayer();
    // Assign next combinatoric (i.e. not first in layer) stub to state
    void comb(State*& state);
    // best state selection
    void accumulator();
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
    // processing region
    int region_;
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

    std::vector<trackerTFP::DataFormatKF*> x0_;
    std::vector<trackerTFP::DataFormatKF*> x1_;
    std::vector<trackerTFP::DataFormatKF*> x2_;
    std::vector<trackerTFP::DataFormatKF*> x3_;
    std::vector<trackerTFP::DataFormatKF*> H00_;
    std::vector<trackerTFP::DataFormatKF*> H12_;
    std::vector<trackerTFP::DataFormatKF*> m0_;
    std::vector<trackerTFP::DataFormatKF*> m1_;
    std::vector<trackerTFP::DataFormatKF*> v0_;
    std::vector<trackerTFP::DataFormatKF*> v1_;
    std::vector<trackerTFP::DataFormatKF*> r0_;
    std::vector<trackerTFP::DataFormatKF*> r1_;
    std::vector<trackerTFP::DataFormatKF*> S00_;
    std::vector<trackerTFP::DataFormatKF*> S01_;
    std::vector<trackerTFP::DataFormatKF*> S12_;
    std::vector<trackerTFP::DataFormatKF*> S13_;
    std::vector<trackerTFP::DataFormatKF*> K00_;
    std::vector<trackerTFP::DataFormatKF*> K10_;
    std::vector<trackerTFP::DataFormatKF*> K21_;
    std::vector<trackerTFP::DataFormatKF*> K31_;
    std::vector<trackerTFP::DataFormatKF*> R00_;
    std::vector<trackerTFP::DataFormatKF*> R11_;
    std::vector<trackerTFP::DataFormatKF*> R00Rough_;
    std::vector<trackerTFP::DataFormatKF*> R11Rough_;
    std::vector<trackerTFP::DataFormatKF*> invR00Approx_;
    std::vector<trackerTFP::DataFormatKF*> invR11Approx_;
    std::vector<trackerTFP::DataFormatKF*> invR00Cor_;
    std::vector<trackerTFP::DataFormatKF*> invR11Cor_;
    std::vector<trackerTFP::DataFormatKF*> invR00_;
    std::vector<trackerTFP::DataFormatKF*> invR11_;
    std::vector<trackerTFP::DataFormatKF*> C00_;
    std::vector<trackerTFP::DataFormatKF*> C01_;
    std::vector<trackerTFP::DataFormatKF*> C11_;
    std::vector<trackerTFP::DataFormatKF*> C22_;
    std::vector<trackerTFP::DataFormatKF*> C23_;
    std::vector<trackerTFP::DataFormatKF*> C33_;
    std::vector<trackerTFP::DataFormatKF*> r02_;
    std::vector<trackerTFP::DataFormatKF*> r12_;
    std::vector<trackerTFP::DataFormatKF*> chi20_;
    std::vector<trackerTFP::DataFormatKF*> chi21_;
  };

}  // namespace trklet

#endif