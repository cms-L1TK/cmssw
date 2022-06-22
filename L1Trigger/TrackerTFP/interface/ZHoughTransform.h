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
    struct State {
      State(const StubMHT& stub, const tt::Setup* setup) : z_(stub.z()), dZ_(stub.dZ()), layer_(stub.layer()), cot_(0), zT_(0), trackId_(stub.trackId()), stub_(&stub) {
        r_ = stub.r() + setup->chosenRofPhi() - setup->chosenRofZ();
      }
      void update(int cot, int zT, const DataFormats* df) {
        cot_ += cot;
        zT_ += zT;
        z_ -= df->base(Variable::zT, Process::zht) * zT + df->base(Variable::cot, Process::zht) * cot * r_;
      }
      double r_;
      double z_;
      double dZ_;
      int layer_;
      int cot_;
      int zT_;
      int trackId_;
      const StubMHT* stub_;
    };
    // perform finer pattern recognition per track
    void stage(int iter, std::vector<State*>& stream);

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
    // container of intermediate stubs
    std::vector<State> states_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<State*>> input_;
  };

}  // namespace trackerTFP

#endif