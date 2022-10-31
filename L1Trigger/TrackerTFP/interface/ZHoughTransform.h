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
    void produce(tt::StreamsStub& accepted, tt::StreamsStub& lost);
  private:
    // run single track through r-z hough transform and store result
    void array(const std::vector<StubMHT*>& track, std::deque<StubZHT*>& stream);
    // run single track through final r-z cell and store result
    void cell(const std::vector<StubZHT*>& track, std::deque<StubZHT*>& stream) const;
    //
    double digi(double val, double base) const { return (std::floor(val / base + 1.e-12) + .5) * base; }
    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // processing region (0 - 8)
    int region_;
    // container of input stubs
    std::vector<StubMHT> stubsMHT_;
    // container of output stubs
    std::vector<StubZHT> stubsZHT_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<StubMHT*>> input_;
    //
    double baseZ_;
    //
    double baseCot_;
    //
    double baseZT_;
    //
    std::vector<TTBV> validZ0_;
  };

}  // namespace trackerTFP

#endif