// TrackletProcessor: this class is an evolved version, performing the tasks of the TrackletEngine+TrackletCalculator.
// It will combine TEs that feed into a TC to a single module.
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletProcessor_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
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
  class TrackletParametersMemory;

  class TrackletProcessor : public ProcessBase {
  public:
    TrackletProcessor(std::string name, Settings const& settings, Globals* globals);

    ~TrackletProcessor() override = default;

    void addOutput(MemoryBase* memory, std::string output) override;

    void addInput(MemoryBase* memory, std::string input) override;

    void execute(unsigned int iSector, double phimin, double phimax);

    void init(int iSeed);

    void exacttracklet(double r1,
                       double z1,
                       double phi1,
                       double r2,
                       double z2,
                       double phi2,
                       double,
                       double& rinv,
                       double& phi0,
                       double& t,
                       double& z0);

    void calcPars(unsigned int idr,
                  int iphi1,
                  int ir1,
                  int iz1,
                  int iphi2,
                  int ir2,
                  int iz2,
                  int& irinv_new,
                  int& iphi0_new,
                  int& iz0_new,
                  int& it_new,
                  bool print = false);

    bool goodTrackPars(bool goodrinv, bool goodz0);

    bool inSector(int iphi0, int irinv, double phi0approx, double rinvapprox);

    bool processStubPair(const Stub* innerFPGAStub,
                         const L1TStub* innerStub,
                         const Stub* outerFPGAStub,
                         const L1TStub* outerStub,
                         bool print);

  private:
    unsigned int iSeed_;
    unsigned int layerdisk1_;
    unsigned int layerdisk2_;
    unsigned int iSector_;

    int TCIndex_;

    double phimin_, phimax_;
    double phiHG_;

    TrackletParametersMemory* trackletpars_;

    int iTC_;
    int iAllStub_;

    unsigned int maxStep_;

    VMStubsTEMemory* outervmstubs_;

    //                                 istub          imem          start imem    end imem
    std::tuple<CircularBuffer<TEData>, unsigned int, unsigned int, unsigned int, unsigned int> tebuffer_;

    std::vector<TrackletEngineUnit> teunits_;

    std::vector<AllInnerStubsMemory*> innerallstubs_;
    std::vector<AllStubsMemory*> outerallstubs_;

    TrackletLUT pttableinner_;
    TrackletLUT pttableouter_;
    TrackletLUT useregiontable_;

    int nbitsfinephi_;
    int nbitsfinephidiff_;

    int innerphibits_;
    int outerphibits_;

    unsigned int nbitszfinebintable_;
    unsigned int nbitsrfinebintable_;

    unsigned int nbitsrzbin_;

    TrackletLUT innerTable_;         //projection to next layer/disk
    TrackletLUT innerOverlapTable_;  //projection to disk from layer

    //Constants for coordinates and track parameter definitions
    int n_phi_;
    int n_r_;
    int n_z_;
    int n_z0_;
    int n_phi0_;
    int n_rinv_;
    int n_t_;

    //Constants used for tracklet parameter calculations
    int n_Deltar_;
    int n_delta0_;
    int n_deltaz_;
    int n_delta1_;
    int n_delta2_;
    int n_delta12_;
    int n_a_;
    int n_r6_;
    int n_ifact_;
    int n_delta02_;
    int n_x6_;
    int n_it1_;
    int n_HG_;

    std::vector<int> LUT_idrinv_;
  };

};  // namespace trklet
#endif
