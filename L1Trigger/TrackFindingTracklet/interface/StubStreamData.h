#ifndef L1Trigger_TrackFindingTracklet_interface_StubStreamData_h
#define L1Trigger_TrackFindingTracklet_interface_StubStreamData_h

#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"

#include <string>

namespace trklet {

  class L1TStub;

  class StubStreamData {
  public:
    StubStreamData() {}
    
    StubStreamData(int iSeed, const L1TStub& stub, const std::string& residuals):
       iSeed_(iSeed), stub_(stub), residuals_(residuals) {}

    ~StubStreamData() = default;

    int iSeed() const {return iSeed_;}
    const L1TStub& stub() const {return stub_;}
    const std::string& residuals() const {return residuals_;}

  private:
    int iSeed_{-1};
    L1TStub stub_;
    std::string residuals_{""};
   
  };
};  // namespace trklet
#endif
