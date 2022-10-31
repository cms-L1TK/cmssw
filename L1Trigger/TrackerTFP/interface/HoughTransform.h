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
    HoughTransform(const edm::ParameterSet& iConfig, const tt::Setup* setup, const DataFormats* dataFormats, int region);
    ~HoughTransform() {}
    // read in and organize input product
    void consume(const tt::StreamsStub& streams);
    // fill output products
    void produce(tt::StreamsStub& accepted, tt::StreamsStub& lost);
  private:
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;
    // remove and return last element of vector, returns nullptr if empty
    template <class T>
    T* pop_back(std::vector<T*>& ts) const;
    // associate stubs with phiT bins in this inv2R column
    void fillIn(int inv2R, std::vector<StubGP*>& input, std::deque<StubHT*>& accepted, std::deque<StubHT*>& lost);
    // identify tracks
    void readOut(const std::deque<StubHT*>& input, std::deque<StubHT*>& output ) const;
    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // data format of inv2R
    DataFormat inv2R_;
    // data format of phiT
    DataFormat phiT_;
    // processing region (0 - 8)
    int region_;
    // container of input stubs
    std::vector<StubGP> stubsGP_;
    // container of output stubs
    std::vector<StubHT> stubsHT_;
    // h/w liked organized pointer to input stubs
    std::vector<std::vector<std::vector<StubGP*>>> input_;
  };

}  // namespace trackerTFP

#endif