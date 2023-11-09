#ifndef L1Trigger_TrackerTFP_HoughTransform_h
#define L1Trigger_TrackerTFP_HoughTransform_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <vector>
#include <deque>

namespace trackerTFP {

  // Class to find initial rough candidates in r-phi in a region
  class HoughTransform {
  public:
    HoughTransform(const edm::ParameterSet& iConfig,
                   const tt::Setup* setup,
                   const DataFormats* dataFormats,
                   std::deque<StubHT>& stubs);
    ~HoughTransform() {}
    // fill output products
    void produce(const std::vector<std::vector<StubGP*>>& streamsIn,
                 std::vector<std::deque<StubHT*>>& streamsOut,
                 std::vector<std::deque<StubHT*>>& streamsTrunc);

  private:
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;
    // associate stubs with phiT bins in this inv2R column
    void fillIn(int inv2R,
                int sector,
                const std::vector<StubGP*>& input,
                std::deque<StubHT*>& output,
                std::deque<StubHT*>& tuncated);
    // identify tracks
    void readOut(const std::deque<StubHT*>& input, std::deque<StubHT*>& output) const;
    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // data format of inv2R
    const DataFormat* inv2R_;
    // data format of phiT
    const DataFormat* phiT_;
    // data format of phi
    const DataFormat* phi_;
    // container of output stubs
    std::deque<StubHT>& stubsHT_;
  };

}  // namespace trackerTFP

#endif