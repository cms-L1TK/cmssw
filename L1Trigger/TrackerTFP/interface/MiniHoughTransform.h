#ifndef L1Trigger_TrackerTFP_MiniHoughTransform_h
#define L1Trigger_TrackerTFP_MiniHoughTransform_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <deque>

namespace trackerTFP {

  // Class to refine HT track candidates in r-phi, by subdividing each HT cell into a finer granularity array
  class MiniHoughTransform {
  public:
    MiniHoughTransform(const edm::ParameterSet& iConfig,
                       const tt::Setup* setup,
                       const DataFormats* dataFormats,
                       int region);
    ~MiniHoughTransform() {}

    // read in and organize input product (fill vector input_)
    void consume(const tt::StreamsStub& streams);
    // fill output products
    void produce(tt::StreamsStub& accepted, tt::StreamsStub& lost);

  private:
    // run single track through r-phi hough transform and store result
    void array(const std::vector<StubHT*>& track, std::deque<StubMHT*>& stream);
    // run single track through final r-phi cell and store result
    void cell(const std::vector<StubMHT*>& track, std::deque<StubMHT*>& stream) const;
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
    std::vector<StubHT> stubsHT_;
    // container of output stubs
    std::vector<StubMHT> stubsMHT_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<StubHT*>> input_;
    //
    double basePhi_;
    //
    double baseInv2R_;
    //
    double basePhiT_;
  };

}  // namespace trackerTFP

#endif