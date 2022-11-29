#ifndef L1Trigger_TrackerTFP_KFin_h
#define L1Trigger_TrackerTFP_KFin_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <deque>

namespace trackerTFP {

  // Class to interconnect Process::zht with Process::kf
  class KFin {
  public:
    KFin(const edm::ParameterSet& iConfig, const tt::Setup* setup, const DataFormats* dataFormats, const LayerEncoding* layerEncoding, int region);
    ~KFin() {}
    // read in and organize input product
    void consume(const tt::StreamsStub& streams);
    // read in and remove TTTrackRefs from argument
    void consume(std::vector<TTTrackRef>& ttTrackRefs);
    // fill output products
    void produce(tt::StreamsTrack& acceptedTracks, tt::StreamsStub& acceptedStubs, tt::StreamsTrack& lostTracks, tt::StreamsStub& lostStubs);
  private:
    //
    void prepare(std::vector<std::deque<TrackKFin*>>& streamsTracks, std::vector<std::deque<StubKFin*>>& streamsStubs);
    //
    void route(std::vector<std::deque<TrackKFin*>>& inputs, std::deque<TrackKFin*>& output) const;
    //
    void route(std::vector<std::deque<StubKFin*>>& input, std::vector<std::deque<StubKFin*>>& outputs) const;
    // remove and return last element of vector, returns nullRef if empty
    TTTrackRef pop_back(std::vector<TTTrackRef>& ttTrackRefs) const;
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;
    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // provides layer Encoding
    const LayerEncoding* layerEncoding_;
    // processing region (0 - 8)
    int region_;
    // container of input stubs
    std::vector<StubZHT> stubsZHT_;
    // container of output stubs
    std::vector<StubKFin> stubsKFin_;
    // container of output TTtracks
    std::vector<TTTrackRef> ttTrackRefs_;
    // container of output tracks
    std::vector<TrackKFin> tracksKFin_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<StubZHT*>> input_;
  };

}  // namespace trackerTFP

#endif