#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iomanip>
#include <filesystem>

using namespace std;
using namespace trklet;

FullMatchMemory::FullMatchMemory(string name, Settings const& settings) : MemoryBase(name, settings) {
  size_t pos = find_nth(name, 0, "_", 1);
  assert(pos != string::npos);
  initLayerDisk(pos + 1, layer_, disk_);
}

void FullMatchMemory::addMatch(Tracklet* tracklet, const Stub* stub) {
  if (!settings_.doKF() || !settings_.doMultipleMatches()) {
    // When allowing only one stub per track per layer (or no KF implying same).
    for (auto& match : matches_) {
      if (match.first == tracklet) {  //Better match: replace existing one
        std::cout << "Replacing stub=" << trklet::hexFormat(match.second->str()) << " with " << trklet::hexFormat(stub->str()) << std::endl;
        std::pair<Tracklet*, const Stub*> tmp(tracklet, stub);
        if(barrel_test<trklet::N_LAYER && layer_==(int)barrel_test+1)
        std::cout << "FullMatch=" << trklet::hexFormat(tmp.first->fullmatchstr(layer_)) << std::endl;
        if(disk_test<trklet::N_LAYER+N_DISK && disk_==(int)disk_test)
        std::cout << "FullMatch=" << trklet::hexFormat(tmp.first->fullmatchdiskstr(disk_)) << std::endl;
        match.second = stub;
        return;
      }
    }
  }
  std::pair<Tracklet*, const Stub*> tmp(tracklet, stub);
  //Check that we have the right TCID order
  if (!matches_.empty()) {
    if ((!settings_.doKF() && matches_[matches_.size() - 1].first->TCID() >= tracklet->TCID()) ||
        (settings_.doKF() && matches_[matches_.size() - 1].first->TCID() > tracklet->TCID())) {
      edm::LogPrint("Tracklet") << "Wrong TCID ordering in " << getName() << " : "
                                << matches_[matches_.size() - 1].first->TCID() << " " << tracklet->TCID() << " "
                                << matches_[matches_.size() - 1].first->trackletIndex() << " "
                                << tracklet->trackletIndex();
    }
  }
  matches_.push_back(tmp);
  if(barrel_test<trklet::N_LAYER && layer_==(int)barrel_test+1)
  std::cout << "FullMatch=" << trklet::hexFormat(tmp.first->fullmatchstr(layer_)) << std::endl;
  if(disk_test<trklet::N_LAYER+N_DISK && disk_==(int)disk_test)
  std::cout << "FullMatch=" << trklet::hexFormat(tmp.first->fullmatchdiskstr(disk_)) << std::endl;
}

void FullMatchMemory::writeMC(bool first, unsigned int iSector) {
  iSector_ = iSector;
  const string dirM = settings_.memPath() + "Matches/";

  std::ostringstream oss;
  oss << dirM << "FullMatches_" << getName() << "_" << std::setfill('0') << std::setw(2) << (iSector_ + 1) << ".dat";
  auto const& fname = oss.str();

  openfile(out_, first, dirM, fname, __FILE__, __LINE__);

  out_ << "BX = " << (bitset<3>)bx_ << " Event : " << event_ << endl;

  for (unsigned int j = 0; j < matches_.size(); j++) {
    string match = (layer_ > 0) ? matches_[j].first->fullmatchstr(layer_) : matches_[j].first->fullmatchdiskstr(disk_);
    out_ << hexstr(j) << " " << match << " " << trklet::hexFormat(match) << endl;
    std::cout << "Wrote FullMatch=" << trklet::hexFormat(match) << std::endl;
  }
  out_.close();

  bx_++;
  event_++;
  if (bx_ > 7)
    bx_ = 0;
}
