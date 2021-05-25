#ifndef L1Trigger_TrackerTFP_SeedFilter_h
#define L1Trigger_TrackerTFP_SeedFilter_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <deque>

namespace trackerTFP {

  // Class to clean up MHT track candidates with rough tracking in r-z
  class SeedFilter {
  public:
    SeedFilter(const edm::ParameterSet& iConfig, const tt::Setup* setup, const DataFormats* dataFormats, const LayerEncoding* layerEncoding, int region);
    ~SeedFilter(){}

    // read in and organize input product
    void consume(const tt::StreamsStub& streams);
    // fill output products
    void produce(tt::StreamsStub& accepted, tt::StreamsStub& lost);

  private:
    // stub layerId (barrel: 1-6, endcap: 11-15)
    int layerId(StubMHT* stub) const;
    // remove and return first element of deque, returns nullptr if empty
    template<class T>
    T* pop_front(std::deque<T*>& ts) const;

    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // provides layer encoding
    const LayerEncoding* layerEncoding_;
    // processing region (0 - 8)
    int region_;
    // dataformat of cot
    DataFormat cot_;
    // dataformat of zT
    DataFormat zT_;
    // dataformat of z
    DataFormat z_;
    // container of input stubs
    std::vector<StubMHT> stubsMHT_;
    // container of output stubs
    std::deque<StubSF> stubsSF_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<StubMHT*>> input_;
  };

}

#endif