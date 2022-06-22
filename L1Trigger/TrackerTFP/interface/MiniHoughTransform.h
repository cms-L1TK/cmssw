#ifndef L1Trigger_TrackerTFP_MiniHoughTransform_h
#define L1Trigger_TrackerTFP_MiniHoughTransform_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <set>
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
    struct State {
      State(const StubHT& stub, const tt::Setup* setup) : r_(stub.r()), phi_(stub.phi()), dPhi_(stub.dPhi()), layer_(stub.layer()), trackId_(stub.trackId()), stubHT_(&stub) {
        const double p = pow(2., setup->mhtNumStages());
        inv2R_ = floor((stub.inv2R() + .5) * setup->mhtNumBinsInv2R() * p);
        phiT_ = floor((stub.phiT() + .5) * setup->mhtNumBinsPhiT() * p);
      }
      void update(int inv2R, int phiT, const DataFormats* df) {
        inv2R_ += inv2R;
        phiT_ += phiT;
        phi_ -= df->base(Variable::phiT, Process::mht) * phiT + df->base(Variable::inv2R, Process::mht) * inv2R * r_;
      }
      double r_;
      double phi_;
      double dPhi_;
      int layer_;
      int inv2R_;
      int phiT_;
      int trackId_;
      const StubHT* stubHT_;
    };
    // perform finer pattern recognition per track
    void stage(int iter, std::vector<State*>& stream);

    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // dataformat of inv2R
    DataFormat inv2R_;
    // dataformat of phiT
    DataFormat phiT_;
    // processing region (0 - 8)
    int region_;
    // number of inv2R bins used in HT
    int numBinsInv2R_;
    // number of chained mhts
    int numStages_;
    // number of cells used in MHT
    int numCells_;
    // container of input stubs
    std::vector<StubHT> stubs_;
    // container of intermediate stubs
    std::vector<State> states_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<State*>> input_;
  };

}  // namespace trackerTFP

#endif