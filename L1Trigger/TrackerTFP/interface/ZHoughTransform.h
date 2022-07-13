#ifndef L1Trigger_TrackerTFP_ZHoughTransform_h
#define L1Trigger_TrackerTFP_ZHoughTransform_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <deque>

namespace trackerTFP {

  // Class to refine MHT track candidates in r-z
  class ZHoughTransform {
  public:
    ZHoughTransform(const edm::ParameterSet& iConfig,
                    const tt::Setup* setup,
                    const DataFormats* dataFormats,
                    int region);
    ~ZHoughTransform() {}

    // read in and organize input product (fill vector input_)
    void consume(const tt::StreamsStub& streams);
    // fill output products
    void produce(tt::StreamsStub& accepted, tt::StreamsStub& lost) const;

  private:
    // run single track through r-z hough transform ans store result
    void produce(const std::vector<StubMHT*>& track, tt::StreamStub& accepted) const;

    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // processing region (0 - 8)
    int region_;
    // container of input stubs
    std::vector<StubMHT> stubs_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<StubMHT*>> input_;
  };

}  // namespace trackerTFP

#endif