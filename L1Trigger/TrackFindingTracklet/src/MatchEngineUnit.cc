#include "L1Trigger/TrackFindingTracklet/interface/MatchEngineUnit.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletLUT.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include <bitset>

using namespace std;
using namespace trklet;

MatchEngineUnit::MatchEngineUnit(const Settings& settings,
                                 bool barrel,
                                 unsigned int layerdisk,
                                 const TrackletLUT& luttable)
    : settings_(settings), luttable_(luttable), candmatches_(3) {
  idle_ = true;
  print_ = false;
  imeu_ = -1;
  barrel_ = barrel;
  layerdisk_ = layerdisk;
  good__ = false;
  good__t = false;
  good___ = false;
}

void MatchEngineUnit::setAlmostFull() { almostfullsave_ = candmatches_.nearfull(); }

void MatchEngineUnit::init(VMStubsMEMemory* vmstubsmemory,
                           unsigned int istep,
                           unsigned int nrzbins,
                           unsigned int rzbin,
                           unsigned int phibin,
                           int shift,
                           int projrinv,
                           int projfinerz,
                           int projfinephi,
                           bool usefirstMinus,
                           bool usefirstPlus,
                           bool usesecondMinus,
                           bool usesecondPlus,
                           bool isPSseed,
                           Tracklet* proj) {
  vmstubsmemory_ = vmstubsmemory;
  istep_ = istep;
  idle_ = false;
  nrzbins_ = nrzbins;
  rzbin_ = rzbin; // Remove rzfirst != rzlast flag
  phibin_ = phibin;
  shift_ = shift;
  istub_ = 0;
  iuse_ = 0;
  projrinv_ = projrinv;
  projfinerz_ = projfinerz;
  projfinephi_ = projfinephi;
  use_.clear();
  int sign = (proj->t() > 0.0) ? 1 : -1;
  int disk = sign * (layerdisk_ - N_LAYER + 1);
  if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_)) {
  if (barrel_test<N_LAYER && layerdisk_ < N_LAYER)
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj->trackletprojstr(layerdisk_+1)) << std::endl;
  else if (disk_test < N_DISK && layerdisk_ >= N_LAYER && abs(disk) == disk_test)
  std::cout << std::hex << "proj=" << trklet::hexFormat(proj->trackletprojstrdisk(disk)) << std::endl;
  }
  std::cout << "Initializing rzbin=" << rzbin_ << " for ";
  if (layerdisk_ < N_LAYER)
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj->trackletprojstr(layerdisk_+1)) << std::endl;
  else
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj->trackletprojstrdisk(disk)) << std::endl;
  if (usefirstMinus) {
    use_.emplace_back(0, 0);
    std::cout << "HERE usefirstMinus" << std::endl;
  }
  if (usesecondMinus) {
    use_.emplace_back(1, 0);
    std::cout << "HERE usesecondMinus" << std::endl;
  }
  if (usefirstPlus) {
    use_.emplace_back(0, 1);
    std::cout << "HERE usefirstPlus" << std::endl;
  }
  if (usesecondPlus) {
    use_.emplace_back(1, 1);
    std::cout << "HERE usesecondPlus" << std::endl;
  }
  assert(!use_.empty());
  std::cout << " nuse=" << use_.size() << std::endl;
  isPSseed_ = isPSseed;
  proj_ = proj;

  good__ = false;
  if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_test)) {
  std::cout << std::hex << "Initializing MEU " << " with ";
  if (barrel_test<N_LAYER && layerdisk_ < N_LAYER)
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstr(layerdisk_+1)) << std::endl;
  else if (disk_test < N_DISK && layerdisk_ >= N_LAYER && abs(disk) == disk_test)
  std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstrdisk(disk)) << std::endl;
    std::cout << "received" << " zbin=" << rzbin_ << std::endl;
        std::cout << "usefirstMinus=" << usefirstMinus << "\t" <<
                     "usesecondMinus=" << usesecondMinus << "\t" <<
                     "usefirstPlus=" << usefirstPlus << "\t" <<
                     "usesecondPlus=" << usesecondPlus << std::endl;
  }
}

void MatchEngineUnit::step(unsigned int istep) {
  good__ = !idle() && !almostfullsave_;

  if (!good__)
    return;

  unsigned int slot = (phibin_ + use_[iuse_].second) * nrzbins_ + rzbin_ + use_[iuse_].first;

  int sign = (proj_->t() > 0.0) ? 1 : -1;
  int disk = sign * (layerdisk_ - N_LAYER + 1);
  if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_test)) {
  if (barrel_test<N_LAYER && layerdisk_ < N_LAYER)
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstr(layerdisk_+1));
  else if (disk_test < N_DISK && layerdisk_ >= N_LAYER && abs(disk) == disk_test)
  std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstrdisk(disk));
  std::cout << " before projfinerz__=" << projfinerz_ << std::endl;
  }
  projfinerz__ = projfinerz_ - (1 << NFINERZBITS) * use_[iuse_].first;
  if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_test)) {
  std::cout << "after projfinerz__=" << projfinerz__ << " use=(" << use_[iuse_].first << "," << use_[iuse_].second << ")" << std::endl;
  }
  projfinephi__ = projfinephi_;
  if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_test)) {
  if (barrel_test<N_LAYER && layerdisk_ < N_LAYER)
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstr(layerdisk_+1));
  else if (disk_test < N_DISK && layerdisk_ >= N_LAYER && abs(disk) == disk_test)
  std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstrdisk(disk));
  std::cout << " before projfinephi__=" << projfinephi__ << std::endl;
  }
  if (use_[iuse_].second == 0) {
    if (shift_ == -1) {
      projfinephi__ -= (1 << NFINEPHIBITS);
    }
  } else {
    //When we get here shift_ is either 1 or -1
    if (shift_ == 1) {
      projfinephi__ += (1 << NFINEPHIBITS);
    }
  }
  if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_test)) {
  std::cout << "after projfinephi__=" << projfinephi__ << " use=(" << use_[iuse_].first << "," << use_[iuse_].second << ")" << std::endl;
  }

  vmstub__ = vmstubsmemory_->getVMStubMEBin(slot, istub_);
  rzbin__ = rzbin_;
  if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_test)) {
  if (barrel_test<N_LAYER && layerdisk_ < N_LAYER)
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstr(layerdisk_+1)) << std::endl;
  else if (disk_test < N_DISK && layerdisk_ >= N_LAYER && abs(disk) == disk_test)
  std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstrdisk(disk)) << std::endl;
  std::cout << "zbin=" << rzbin_ << std::endl;
  std::cout << "use_=" << use_[iuse_].first << ", " << use_[iuse_].second << std::endl;
  std::cout << "stubid=" << trklet::hexFormat(vmstub__.stub()->stubindex().str()) << " (" << trklet::hexFormat(vmstub__.stub()->str()) << ")" << std::endl;
  std::cout << "first=" << use_[iuse_].first << "\tsecond=" << use_[iuse_].second << std::endl;
  std::cout << "building slot " << (phibin_ + use_[iuse_].second) << "*" << nrzbins_ << "+" << rzbin_ << "+"
            << use_[iuse_].first << std::endl;
  std::cout << "building slot " << "(" << phibin_ << "+" << use_[iuse_].second << ")" << "*" << nrzbins_ << "+" << rzbin_ << "+"
            << use_[iuse_].first << " with nstubs=" << vmstubsmemory_->nStubsBin(slot) << " and iuse=" << iuse_ << std::endl;
  std::cout << "istub_=" << istub_ << std::endl;
  std::cout << "stubadd=" << (slot << 4) + istub_ << std::endl;
  }
  string stub = vmstub__.stub()->stubindex().str();
  stub += "|";
  if (!barrel_ && vmstub__.isPSmodule())
    stub += "0";
  stub += vmstub__.stub()->bend().str();

  FPGAWord finephipos = vmstub__.finephi();
  stub += "|" + finephipos.str();
  FPGAWord finepos = vmstub__.finerz();
  stub += "|" + finepos.str();
  if (abs(disk) == disk_test || layerdisk_ == barrel_test) {
  if (barrel_test<N_LAYER && layerdisk_ < N_LAYER) {
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstr(layerdisk_+1)) << " use_=(" << use_[iuse_].first << ", " << use_[iuse_].second << ") nuse=" << use_.size();
  }
  else if (disk_test < N_DISK && layerdisk_ >= N_LAYER && abs(disk) == disk_test)
    std::cout << std::hex << "proj=" << trklet::hexFormat(proj_->trackletprojstrdisk(disk)) << " use_[" << iuse_ << "]=(" << use_[iuse_].first << ", " << use_[iuse_].second << ") nuse=" << use_.size();
  std::cout << std::hex << " " << (layerdisk_ < N_LAYER ? "Barrel" : "Disk") << " MEU " << imeu_ << " pipelining vmstub=" << trklet::hexFormat(stub) << " isPS=" << vmstub__.isPSmodule() << " bend=" << vmstub__.stub()->bend().str() << " stub=" << trklet::hexFormat(vmstub__.stub()->str()) << " with stubaddr=(" << slot << "," << istub_ << "/" << vmstubsmemory_->nStubsBin(slot) << ")" << " stubid=" << trklet::hexFormat(vmstub__.stub()->stubindex().str()) << "\tistep=" << istep << std::endl;
  std::cout << trklet::hexFormat(stub) << "\t" << stub << std::endl;
  }

  isPSseed__ = isPSseed_;
  projrinv__ = projrinv_;
  proj__ = proj_;

  istub_++;
  if (istub_ >= vmstubsmemory_->nStubsBin(slot)) {
    iuse_++;
    if (iuse_ < use_.size()) {
      istub_ = 0;
    } else {
      idle_ = true;
    }
  }
}

void MatchEngineUnit::processPipeline() {
  if (good___) {
    //bool isPSmodule = vmstub___.isPSmodule();
    int stubfinerz = vmstub___.finerz().value();
    int stubfinephi = vmstub___.finephi().value();
    int stubbend = vmstub___.bend().value();
    bool isPSmodule = false;

    if (barrel_) {
      isPSmodule = layerdisk_ < N_PSLAYER;
    } else {
      const int absz = (1 << settings_.MEBinsBits()) - 1;
      std::cout << "absz=" << absz << std::endl;
      if (layerdisk_ < N_LAYER + 2) {
        isPSmodule = ((rzbin___ & absz) < 3) || ((rzbin___ & absz) == 3 && stubfinerz <= 3);
      } else {
        isPSmodule = ((rzbin___ & absz) < 3) || ((rzbin___ & absz) == 3 && stubfinerz <= 2);
      }
    }
    isPSmodule = barrel_ ? isPSmodule : stubbend < (1<<N_BENDBITS_2S)-1;

    int sign = (proj___->t() > 0.0) ? 1 : -1;
    int disk = sign * (layerdisk_ - N_LAYER + 1);
    if(!barrel_)
    std::cout << "rzbin=" << rzbin___ << " " << std::bitset<4>(rzbin___) << " r=" << std::dec << vmstub___.stub()->r().value() << std::hex << " stubfinerz=" << stubfinerz << " " << std::bitset<trklet::NFINERZBITS>(stubfinerz) << " stubfinephi=" << stubfinephi << " " << std::bitset<trklet::NFINEPHIBITS>(stubfinephi) << " stubindex=" << vmstub___.stubindex().value() << " bend=" << stubbend << " " << std::bitset<N_BENDBITS_2S>(stubbend) << " isPS=" << isPSmodule << " true isPS=" << vmstub___.isPSmodule() << " matches=" << (isPSmodule == vmstub___.isPSmodule()) << " DISK=" << layerdisk_ - N_LAYER + 1 << " for ";
    if (layerdisk_ < N_LAYER)
      std::cout << std::hex << "proj=" << trklet::hexFormat(proj___->trackletprojstr(layerdisk_+1)) << std::endl;
    else
      std::cout << std::hex << "proj=" << trklet::hexFormat(proj___->trackletprojstrdisk(disk)) << std::endl;
    //isPSmodule = vmstub___.isPSmodule();
    assert(isPSmodule == vmstub___.isPSmodule());
    /*
    */

    int deltaphi = stubfinephi - projfinephi___;

    constexpr int idphicut = 3;

    bool dphicut = (abs(deltaphi) < idphicut);
    std::cout << "dphicut=" << dphicut << " stubfinephi=" << stubfinephi <<  " =" << " projfinephi___=" << projfinephi___ << " idphicut=" << idphicut << std::endl;

    int nbits = isPSmodule ? N_BENDBITS_PS : N_BENDBITS_2S;

    int diskps = (!barrel_) && isPSmodule;

    //here we always use the larger number of bits for the bend
    unsigned int index = (diskps << (nbits + NRINVBITS)) + (projrinv___ << nbits) + stubbend;//vmstub___.bend().value();

    //Check if stub z position consistent
    int idrz = stubfinerz - projfinerz___;
    bool pass;

    if (barrel_) {
      if (isPSseed___) {
        constexpr int drzcut = 1;
        pass = std::abs(idrz) <= drzcut;
      } else {
        constexpr int drzcut = 5;
        pass = std::abs(idrz) <= drzcut;
      }
    } else {
      if (isPSmodule) {
        constexpr int drzcut = 1;
        pass = std::abs(idrz) <= drzcut;
      } else {
        constexpr int drzcut = 3;
        pass = std::abs(idrz) <= drzcut;
      }
    }

    bool goodpair = (pass && dphicut) && luttable_.lookup(index);
    std::cout << proj___ << std::endl;
    if (disk_test< N_DISK && (abs(disk) == disk_test || layerdisk_ == barrel_test)) {
      int drzcut = 0;
      if (barrel_) {
        if (isPSseed___) {
          drzcut = 1;
        }
        else {
          drzcut = 5;
        }
      }
      else {
        if (isPSmodule) {
          drzcut = 1;
        }
        else {
          drzcut = 3;
        }
      }
      if (layerdisk_ < N_LAYER) {
        std:: cout << "allpass=" << goodpair << "\tproj=" << trklet::hexFormat(proj___->trackletprojstr(layerdisk_+1)) << "\tgood___=" << good___ << "\tpassphi=" << dphicut << "(projfinephi=" << projfinephi___ << " - stubfinephi=" << stubfinephi << ") < " << idphicut << "\tpass=" << pass << " (";
        std::cout << "projfinezadj___=" << projfinerz___ << " - stubfinez=" << stubfinerz << ") < " << drzcut;
        std::cout << "\ttable[" << index << "]=" << luttable_.lookup(index) << " (diskps=" << diskps << " projrinv___=" << projrinv___ << " stubbend=" << stubbend << ")";
        //std::cout << "\ttable[" << index << "]=" << luttable_.lookup(index) << " (diskps=" << diskps << " projrinv___=" << projrinv___ << " stubbend=" << vmstub___.bend().value() << ")";
      }
      else {
        std:: cout << "allpass=" << goodpair << "\tproj=" << trklet::hexFormat(proj___->trackletprojstrdisk(disk)) << "\tgood___=" << good___ << "\tpassphi=" << dphicut << "(projfinephi=" << projfinephi___ << " - stubfinephi=" << stubfinephi << ") < " << idphicut << "\tpass=" << pass << " (";
        std::cout << "projfinezadj___=" << projfinerz___ << " - stubfinez=" << stubfinerz << ") < " << drzcut;
        std::cout << "\ttable[" << index << "]=" << luttable_.lookup(index) << " (diskps=" << diskps << " projrinv___=" << projrinv___ << " stubbend=" << stubbend << ")";
        //std::cout << "\ttable[" << index << "]=" << luttable_.lookup(index) << " (diskps=" << diskps << " projrinv___=" << projrinv___ << " stubbend=" << vmstub___.bend().value() << ")";
      }
      std::cout << "\tfor stub=" << trklet::hexFormat(vmstub___.stub()->str()) << " isPS=" << isPSmodule << " diskps=" << diskps << std::endl;
    }

    std::pair<Tracklet*, const Stub*> tmppair(proj___, vmstub___.stub());

    if (goodpair) {
      candmatches_.store(tmppair);
    }
  }

  proj___ = proj__t;
  projfinephi___ = projfinephi__t;
  projfinerz___ = projfinerz__t;
  projrinv___ = projrinv__t;
  isPSseed___ = isPSseed__t;
  good___ = good__t;
  vmstub___ = vmstub__t;
  rzbin___ = rzbin__t;

  proj__t = proj__;
  projfinephi__t = projfinephi__;
  projfinerz__t = projfinerz__;
  projrinv__t = projrinv__;
  isPSseed__t = isPSseed__;
  good__t = good__;
  vmstub__t = vmstub__;
  rzbin__t = rzbin__;
}

void MatchEngineUnit::print() {
  std::cout << "Content of MEU" << std::endl;
  for(size_t i = 0; i < candmatches_.size(); i++) {
    auto *proj = candmatches_[i].first;
    auto *stub = candmatches_[i].second;
    int sign = (proj->t() > 0.0) ? 1 : -1;
    int disk = sign * (layerdisk_ - N_LAYER + 1);
    if (layerdisk_ < N_LAYER)
      std::cout << i << ": proj=" << trklet::hexFormat(proj->trackletprojstr(layerdisk_+1)) << "\tproj=" << trklet::hexFormat(stub->str()) << "\tistep=" << istep_ << std::endl;
    else
      std::cout << i << ": proj=" << trklet::hexFormat(proj->trackletprojstrdisk(disk)) << "\tproj=" << trklet::hexFormat(stub->str()) << "\tistep=" << istep_ << std::endl;
  }
}

void MatchEngineUnit::reset() {
  candmatches_.reset();
  idle_ = true;
  istub_ = 0;
  good__ = false;
  good__t = false;
  good___ = false;
}

int MatchEngineUnit::TCID() const {
  if (!empty()) {
    return peek().first->TCID();
  }

  if (good___) {
    return proj___->TCID();
  }

  if (good__t) {
    return proj__t->TCID();
  }

  if (good__) {
    return proj__->TCID();
  }

  if (idle_) {
    return (1 << (settings_.nbitstrackletindex() + settings_.nbitstcindex())) - 1;
  }

  return proj_->TCID();
}
