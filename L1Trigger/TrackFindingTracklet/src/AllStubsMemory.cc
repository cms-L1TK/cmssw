#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "DataFormats/L1TrackTrigger/interface/TTBV.h"
#include <iomanip>

using namespace std;
using namespace trklet;

AllStubsMemory::AllStubsMemory(string name, Settings const& settings) : MemoryBase(name, settings) {}

void AllStubsMemory::writeStubs(bool first, unsigned int iSector) {
  iSector_ = iSector;
  const string dirS = settings_.memPath() + "Stubs/";
  openFile(first, dirS, "AllStubs_");

  for (unsigned int j = 0; j < stubs_.size(); j++) {
    string stub = stubs_[j]->str();
    out_ << hexstr(j) << " " << stub << " " << hexFormat(stub) << endl;
  }
  out_.close();
}

void AllStubsMemory::writeStubs(std::vector<StubStreamData>& streamsAllStubRaw) {
  streamsAllStubRaw.reserve(stubs_.size());
  for (const Stub* stub : stubs_) {
    if (!stub) {
      streamsAllStubRaw.emplace_back();
      continue;
    }
    static constexpr int widthDisk2Sidentifier = 3;
    static constexpr int widthDiskPSradius = 12;
    string r = stub->r().str();
    if (!stub->isBarrel()) {
      if (stub->isPSmodule())
        r = FPGAWord(stub->rvalue(), widthDiskPSradius).str();
      else
        r = string(widthDisk2Sidentifier, '0') + r;
    }
    std::string s;
    s += r;
    s += stub->z().str();
    s += stub->phi().str();
    if (stub->alpha().nbits() > 0)
      s += stub->alpha().str();
    s += stub->bend().str();
    streamsAllStubRaw.emplace_back(0, *(stub->l1tstub()), s);
  }
}
