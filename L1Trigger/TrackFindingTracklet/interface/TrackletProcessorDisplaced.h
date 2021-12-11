#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletProcessorDisplaced_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletProcessorDisplaced_h

#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletLUT.h"
#include "L1Trigger/TrackFindingTracklet/interface/CircularBuffer.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletEngineUnit.h"

#include <vector>
#include <tuple>
#include <map>

namespace trklet {

  class Settings;
  class Globals;
  class MemoryBase;
  class AllStubsMemory;
  class AllInnerStubsMemory;
  class VMStubsTEMemory;
  class StubPairsMemory;

  class TrackletProcessorDisplaced : public TrackletCalculatorDisplaced {
  public:
    TrackletProcessorDisplaced(std::string name, Settings const& settings, Globals* globals);

    ~TrackletProcessorDisplaced() override = default;

    void addOutputProjection(TrackletProjectionsMemory*& outputProj, MemoryBase* memory);

    void addOutput(MemoryBase* memory, std::string output) override;

    void addInput(MemoryBase* memory, std::string input) override;

    void execute(unsigned int iSector, double phimin, double phimax);

  private:
    int iTC_;
    int iAllStub_;
    unsigned int maxStep_;
    int count_;
    unsigned int layerdisk_;
    /* VMStubsTEMemory* outervmstubs_; */

    int TCIndex_;
    int layer_;
    int disk_;

    int layer1_;
    int layer2_;
    int layer3_;
    int disk1_;
    int disk2_;
    int disk3_;

    int firstphibits_;
    int secondphibits_;
    int thirdphibits_;

    int iSeed_;

    double rproj_[N_LAYER - 2];
    int lproj_[N_LAYER - 2];
    double zproj_[N_DISK - 2];
    int dproj_[N_DISK - 2];
    double rzmeanInv_[N_DISK - 2];

    unsigned int iSector_;
    double phimin_, phimax_;

    int nbitszfinebintable_;
    int nbitsrfinebintable_;

    /*                                 istub          imem          start imem    end imem */
    /* std::tuple<CircularBuffer<TEData>, unsigned int, unsigned int, unsigned int, unsigned int> tebuffer_; */

    /* std::vector<TrackletEngineUnit> teunits_; */

    TrackletLUT innerTable_;         //projection to next layer/disk
    TrackletLUT innerOverlapTable_;  //projection to disk from layer
    TrackletLUT innerThirdTable_;    //projection to disk1 for extended - iseed=10

    std::vector<double> toR_;
    std::vector<double> toZ_;

    std::vector<AllStubsMemory*> innerallstubs_;
    std::vector<AllStubsMemory*> middleallstubs_;
    std::vector<AllStubsMemory*> outerallstubs_;
    std::vector<StubPairsMemory*> stubpairs_;
    /* std::vector<StubTripletsMemory*> stubtriplets_; */
    std::vector<VMStubsTEMemory*> innervmstubs_;
    /* std::vector<VMStubsTEMemory*> outervmstubs_; */
    VMStubsTEMemory* outervmstubs_;

    TrackletParametersMemory* trackletpars_;
    StubTripletsMemory* stubtriplets_;
                                                         
    std::vector<std::vector<TrackletProjectionsMemory*> > trackletprojlayers_;
    std::vector<std::vector<TrackletProjectionsMemory*> > trackletprojdisks_;

    std::map<std::string, std::vector<std::vector<std::string> > > tmpSPTable_;
    std::map<std::string, std::vector<std::map<std::string, unsigned> > > spTable_;
    std::vector<bool> table_;

  };

};  // namespace trklet
#endif
