#ifndef L1Trigger_TrackFindingTracklet_AllStubProcessor_h
#define L1Trigger_TrackFindingTracklet_AllStubProcessor_h

#include "L1Trigger/TrackFindingTracklet/interface/Setup.h"
#include "L1Trigger/TrackFindingTracklet/interface/DataFormats.h"

#include <vector>
#include <deque>

namespace trklet {

  /*! \class  trklet::AllStubProcessor
   *  \brief  Class to process All Stubs
   *  \author Thomas Schuh
   *  \date   2026, Nov
   */
  class AllStubProcessor {
  public:
    AllStubProcessor(const Setup*, const DataFormats*, int);
    ~AllStubProcessor() = default;
    // read in and organize input
    void consume(const tt::StreamsStub&);
    // fill output
    void produce(tt::StreamsStub&);

  private:
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;
    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const Setup* setup_;
    // provides data formats
    const DataFormats* dataFormats_;
    // processing region (0 - 8) aka processing phi nonant
    const int region_;
    // input data container
    std::vector<StubASP> stubs_;
    // stub streams
    std::vector<std::deque<StubASP*>> streams_;
  };

}  // namespace trklet

#endif
