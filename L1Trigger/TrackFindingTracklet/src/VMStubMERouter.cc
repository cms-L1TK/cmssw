#include "L1Trigger/TrackFindingTracklet/interface/VMStubMERouter.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace std;
using namespace trklet;

VMStubMERouter::VMStubMERouter(string name, Settings const& settings, Globals* global)
    : ProcessBase(name, settings, global) {
  layerdisk_ = initLayerDisk(7);

  unsigned int region = name[12] - 'A';
  assert(region < settings_.nallstubs(layerdisk_));

  vmstubsMEPHI_.resize(1, nullptr);

}

void VMStubMERouter::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }

  if (output.substr(0, 10) == "allstubout") {
    AllStubsMemory* tmp = dynamic_cast<AllStubsMemory*>(memory);
    allstubs_.push_back(tmp);
    return;
  }

  if (output == "vmstubout") {
    VMStubsMEMemory* tmp = dynamic_cast<VMStubsMEMemory*>(memory);
    assert(tmp != nullptr);
    tmp->resize(16 * settings_.nvmme(layerdisk_));
    assert(vmstubsMEPHI_[0] == nullptr);
    vmstubsMEPHI_[0] = tmp;
    return;
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find output : " << output;
}

void VMStubMERouter::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }
  if (input == "allstubin") {
    AllStubsMemory* tmp1 = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp1 != nullptr);
    if (tmp1 != nullptr) {
      stubinputs_.push_back(tmp1);
    }
    return;
  }
  if (input == "vmstubin") {
    VMStubsMEMemory* tmp1 = dynamic_cast<VMStubsMEMemory*>(memory);
    assert(tmp1 != nullptr);
    if (tmp1 != nullptr) {
      vmstubsinput_.push_back(tmp1);
    }
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find input : " << input;
}

void VMStubMERouter::execute(unsigned int) {
  unsigned int allStubCounter = 0;

  //bool print = getName() == "VMR_D1PHIB" && iSector == 3;
  //print = false;

  std::cout << "Will process allstubs:" << getName() << std::endl;

  //Loop over the input stubs
  for (auto& stubinput : stubinputs_) {
    for (unsigned int i = 0; i < stubinput->nStubs(); i++) {
      if (allStubCounter > settings_.maxStep("VMR"))
        continue;
      if (allStubCounter >= (1 << N_BITSMEMADDRESS))
        continue;

      const Stub* stub = stubinput->getStub(i);

      allStubCounter++;

      for (auto& allstub : allstubs_) {
        allstub->addStub(stub);
      }
    }
  }

  const unsigned int nbin = vmstubsinput_[0]->size();

  std::cout << "Will process vmstubs:" << getName() << std::endl;

  std::cout << vmstubsinput_[0]->getName() << " " << nbin <<std::endl;

  for (unsigned int ivmbin = 0; ivmbin < nbin; ivmbin++) {
    const unsigned int nvmstub = vmstubsinput_[0]->nStubsBin(ivmbin);
    for (unsigned int istub = 0; istub < nvmstub; istub++) {
      vmstubsMEPHI_[0]->addStub(vmstubsinput_[0]->getVMStubMEBin(ivmbin, istub), ivmbin);
    }
  }
}
