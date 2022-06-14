#include "L1Trigger/TrackFindingTracklet/interface/MatchCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"
#include "L1Trigger/TrackFindingTracklet/interface/CandidateMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/HistBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <filesystem>
#include <algorithm>
#include <bitset>

using namespace std;
using namespace trklet;

MatchCalculator::MatchCalculator(string name, Settings const& settings, Globals* global)
    : ProcessBase(name, settings, global),
      phimatchcuttable_(settings),
      zmatchcuttable_(settings),
      rphicutPStable_(settings),
      rphicut2Stable_(settings),
      rcutPStable_(settings),
      rcut2Stable_(settings),
      alphainner_(settings),
      alphaouter_(settings),
      rSSinner_(settings),
      rSSouter_(settings) {
  phiregion_ = name[8] - 'A';
  layerdisk_ = initLayerDisk(3);

  fullMatches_.resize(12, nullptr);

  //TODO - need to sort out constants here
  icorrshift_ = 7;

  if (layerdisk_ < N_PSLAYER) {
    icorzshift_ = -1 - settings_.PS_zderL_shift();
  } else {
    icorzshift_ = -1 - settings_.SS_zderL_shift();
  }
  phi0shift_ = 3;
  fact_ = 1;
  if (layerdisk_ >= N_PSLAYER && layerdisk_ < N_LAYER) {
    fact_ = (1 << (settings_.nzbitsstub(0) - settings_.nzbitsstub(5)));
    icorrshift_ -= (10 - settings_.nrbitsstub(layerdisk_));
    icorzshift_ += (settings_.nzbitsstub(0) - settings_.nzbitsstub(5) + settings_.nrbitsstub(layerdisk_) -
                    settings_.nrbitsstub(0));
    phi0shift_ = 0;
  }

  unsigned int region = getName()[8] - 'A';
  assert(region < settings_.nallstubs(layerdisk_));

  if (layerdisk_ < N_LAYER) {
    phimatchcuttable_.initmatchcut(layerdisk_, TrackletLUT::MatchType::barrelphi, region);
    zmatchcuttable_.initmatchcut(layerdisk_, TrackletLUT::MatchType::barrelz, region);
  } else {
    rphicutPStable_.initmatchcut(layerdisk_, TrackletLUT::MatchType::diskPSphi, region);
    rphicut2Stable_.initmatchcut(layerdisk_, TrackletLUT::MatchType::disk2Sphi, region);
    rcutPStable_.initmatchcut(layerdisk_, TrackletLUT::MatchType::diskPSr, region);
    rcut2Stable_.initmatchcut(layerdisk_, TrackletLUT::MatchType::disk2Sr, region);
    alphainner_.initmatchcut(layerdisk_, TrackletLUT::MatchType::alphainner, region);
    alphaouter_.initmatchcut(layerdisk_, TrackletLUT::MatchType::alphaouter, region);
    rSSinner_.initmatchcut(layerdisk_, TrackletLUT::MatchType::rSSinner, region);
    rSSouter_.initmatchcut(layerdisk_, TrackletLUT::MatchType::rSSouter, region);
  }

  for (unsigned int i = 0; i < N_DSS_MOD * 2; i++) {
    ialphafactinner_[i] = (1 << settings_.alphashift()) * settings_.krprojshiftdisk() * settings_.half2SmoduleWidth() /
                          (1 << (settings_.nbitsalpha() - 1)) / (settings_.rDSSinner(i) * settings_.rDSSinner(i)) /
                          settings_.kphi();
    ialphafactouter_[i] = (1 << settings_.alphashift()) * settings_.krprojshiftdisk() * settings_.half2SmoduleWidth() /
                          (1 << (settings_.nbitsalpha() - 1)) / (settings_.rDSSouter(i) * settings_.rDSSouter(i)) /
                          settings_.kphi();
  }
}

void MatchCalculator::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  if (output.substr(0, 8) == "matchout") {
    auto* tmp = dynamic_cast<FullMatchMemory*>(memory);
    assert(tmp != nullptr);
    unsigned int iSeed = getISeed(memory->getName());
    fullMatches_[iSeed] = tmp;
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find output " << output;
}

void MatchCalculator::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }
  if (input == "allstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    allstubs_ = tmp;
    return;
  }
  if (input == "allprojin") {
    auto* tmp = dynamic_cast<AllProjectionsMemory*>(memory);
    assert(tmp != nullptr);
    allprojs_ = tmp;
    return;
  }
  if (input.substr(0, 5) == "match" && input.substr(input.size() - 2, 2) == "in") {
    auto* tmp = dynamic_cast<CandidateMatchMemory*>(memory);
    assert(tmp != nullptr);
    matches_.push_back(tmp);
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find input " << input;
}

void MatchCalculator::execute(unsigned int iSector, double phioffset) {
  unsigned int countall = 0;
  unsigned int countsel = 0;

  //bool print = getName() == "MC_L4PHIC" && iSector == 3;

  Tracklet* oldTracklet = nullptr;

  std::vector<std::pair<std::pair<Tracklet*, int>, const Stub*> > mergedMatches = mergeMatches(matches_);

  // Number of clock cycles the pipeline in HLS takes to process the projection merging to 
  // produce the first projectio
  unsigned int mergedepth = 3;

  unsigned int maxProc =  std::min(settings_.maxStep("MC") - mergedepth, (unsigned int)mergedMatches.size());

  int best_ideltaphi_barrel = 0xFFFF;
  int best_ideltaz_barrel = 0xFFFF;
  int best_ideltaphi_disk = 0xFFFF;
  int best_ideltar_disk = 0xFFFF;
  unsigned int curr_projid = -1;
  unsigned int next_projid = -1;

  for (unsigned int j = 0; j < maxProc; j++) {
    if (settings_.debugTracklet() && j == 0) {
      edm::LogVerbatim("Tracklet") << getName() << " has " << mergedMatches.size() << " candidate matches";
    }

    countall++;

    const Stub* fpgastub = mergedMatches[j].second;
    Tracklet* tracklet = mergedMatches[j].first.first;
    const L1TStub* stub = fpgastub->l1tstub();

    //check that the matches are orderd correctly
    //allow equal here since we can have more than one cadidate match per tracklet projection
    if (oldTracklet != nullptr) {
      assert(oldTracklet->TCID() <= tracklet->TCID());
    }
    oldTracklet = tracklet;

    if (layerdisk_ < N_LAYER) {
      //Integer calculation

      const Projection& proj = tracklet->proj(layerdisk_);

      int ir = fpgastub->r().value();
      int iphi = proj.fpgaphiproj().value();
      int icorr = (ir * proj.fpgaphiprojder().value()) >> icorrshift_;
      iphi += icorr;

      int iz = proj.fpgarzproj().value();
      int izcor = (ir * proj.fpgarzprojder().value() + (1 << (icorzshift_ - 1))) >> icorzshift_;
      iz += izcor;

      int ideltaz = fpgastub->z().value() - iz;
      int ideltaphi = (fpgastub->phi().value() << phi0shift_) - (iphi << (settings_.phi0bitshift() - 1 + phi0shift_));

      //Floating point calculations

      double phi = stub->phi() - phioffset;
      double r = stub->r();
      double z = stub->z();

      if (settings_.useapprox()) {
        double dphi = reco::reduceRange(phi - fpgastub->phiapprox(0.0, 0.0));
        assert(std::abs(dphi) < 0.001);
        phi = fpgastub->phiapprox(0.0, 0.0);
        z = fpgastub->zapprox();
        r = fpgastub->rapprox();
      }

      if (phi < 0)
        phi += 2 * M_PI;

      double dr = r - settings_.rmean(layerdisk_);
      assert(std::abs(dr) < settings_.drmax());

      double dphi = reco::reduceRange(phi - (proj.phiproj() + dr * proj.phiprojder()));

      double dz = z - (proj.rzproj() + dr * proj.rzprojder());

      double dphiapprox = reco::reduceRange(phi - (proj.phiprojapprox() + dr * proj.phiprojderapprox()));

      double dzapprox = z - (proj.rzprojapprox() + dr * proj.rzprojderapprox());

      int seedindex = tracklet->getISeed();
      unsigned int projindex = mergedMatches[j].first.second;               // Allproj index
      curr_projid = next_projid;
      next_projid = projindex;
      bool newtracklet = (j == 0 || projindex != curr_projid);
      if (j == 0)  best_ideltar_disk = (1 << (fpgastub->r().nbits() - 1));
      if (newtracklet) {
        best_ideltaphi_barrel = (int)phimatchcuttable_.lookup(seedindex);
        best_ideltaz_barrel = (int)zmatchcuttable_.lookup(seedindex);
      }
      /*
      */

      assert(phimatchcuttable_.lookup(seedindex) > 0);
      assert(zmatchcuttable_.lookup(seedindex) > 0);

      if (settings_.bookHistos()) {
        bool truthmatch = tracklet->stubtruthmatch(stub);

        HistBase* hists = globals_->histograms();
        hists->FillLayerResidual(layerdisk_ + 1,
                                 seedindex,
                                 dphiapprox * settings_.rmean(layerdisk_),
                                 ideltaphi * settings_.kphi1() * settings_.rmean(layerdisk_),
                                 ideltaz * fact_ * settings_.kz(),
                                 dz,
                                 truthmatch);
      }

      //This would catch significant consistency problems in the configuration - helps to debug if there are problems.
      if (std::abs(dphi) > 0.5 * settings_.dphisectorHG() || std::abs(dphiapprox) > 0.5 * settings_.dphisectorHG()) {
        throw cms::Exception("LogicError")
            << "WARNING dphi and/or dphiapprox too large : " << dphi << " " << dphiapprox << endl;
      }

      if (settings_.writeMonitorData("Residuals")) {
        double pt = 0.01 * settings_.c() * settings_.bfield() / std::abs(tracklet->rinv());

        globals_->ofstream("layerresiduals.txt")
            << layerdisk_ + 1 << " " << seedindex << " " << pt << " "
            << ideltaphi * settings_.kphi1() * settings_.rmean(layerdisk_) << " "
            << dphiapprox * settings_.rmean(layerdisk_) << " "
            << phimatchcuttable_.lookup(seedindex) * settings_.kphi1() * settings_.rmean(layerdisk_) << "   "
            << ideltaz * fact_ * settings_.kz() << " " << dz << " "
            << zmatchcuttable_.lookup(seedindex) * settings_.kz() << endl;
      }

      /*
      bool imatch = (std::abs(ideltaphi) <= (int)phimatchcuttable_.lookup(seedindex)) &&
                    (ideltaz * fact_ < (int)zmatchcuttable_.lookup(seedindex)) &&
                    (ideltaz * fact_ >= - (int)zmatchcuttable_.lookup(seedindex));
      */
      bool imatch = (std::abs(ideltaphi) <= best_ideltaphi_barrel) &&
                    (ideltaz * fact_ < best_ideltaz_barrel) &&
                    (ideltaz * fact_ >= - best_ideltaz_barrel);
      if (imatch) {
        best_ideltaphi_barrel = std::abs(ideltaphi);
        best_ideltaz_barrel = std::abs(ideltaz * fact_);
      }

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << getName() << " imatch = " << imatch << " ideltaphi cut " << ideltaphi << " "
                                     << phimatchcuttable_.lookup(seedindex) << " ideltaz*fact cut " << ideltaz * fact_
                                     << " " << zmatchcuttable_.lookup(seedindex);
      }

      if (imatch) {
        countsel++;

        tracklet->addMatch(layerdisk_,
                           ideltaphi,
                           ideltaz,
                           dphi,
                           dz,
                           dphiapprox,
                           dzapprox,
                           (phiregion_ << 7) + fpgastub->stubindex().value(),
                           mergedMatches[j].second);

        if (settings_.debugTracklet()) {
          edm::LogVerbatim("Tracklet") << "Accepted full match in layer " << getName() << " " << tracklet;
        }

        fullMatches_[seedindex]->addMatch(tracklet, mergedMatches[j].second);
      }
    } else {  //disk matches
      {
        int projindex = mergedMatches[j].first.second;               // Allproj index
        string stubid = mergedMatches[j].second->stubindex().str();  // stub ID
        int stubindex = 0;
        for(unsigned int i = 0; i < stubid.length();  i++) {
          if(stubid[i] == '1')
            stubindex += pow(2, stubid.length()-1-i);
        }
        FPGAWord tmp;
        if (projindex >= (1 << 7)) {
          projindex = (1 << 7) - 1;
        }
        tmp.set(projindex, 7, true, __LINE__, __FILE__);
        std::cout << "found a cmatch=" << std::bitset<7>(projindex) << "|" << stub->isPSmodule() << "|" << std::bitset<7>(stubindex) << "\t" << trklet::hexFormat(tmp.str() + stubid) << std::endl;
      }

      //check that stubs and projections in same half of detector
      assert(stub->z() * tracklet->t() > 0.0);

      int sign = (tracklet->t() > 0.0) ? 1 : -1;
      int disk = sign * (layerdisk_ - (N_LAYER - 1));
      assert(disk != 0);
        std::cout << std::hex << "stub=" << trklet::hexFormat(fpgastub->str()) << std::endl;
        std::cout << std::hex << "proj=" << trklet::hexFormat(tracklet->trackletprojstrD(disk)) << std::endl;

      //Perform integer calculations here

      const Projection& proj = tracklet->proj(layerdisk_);
      string disks[] = { "", "D1", "D2", "D3", "D4", "D5" };

      int iz = fpgastub->z().value();
      std::cout << "iz=" << iz << "\t" << std::bitset<7>(iz) << std::endl;
      int iphi = proj.fpgaphiproj().value();
      std::cout << "iphi=" << iphi << "\t" << std::bitset<14>(iphi) << std::endl;

      //TODO - need to express interms of constants
      int shifttmp = 6;
      int iphicorr = (iz * proj.fpgaphiprojder().value()) >> shifttmp;
      std::cout << "iphicorr=" << iphicorr << std::endl;

      iphi += iphicorr;
      std::cout << "iphi=" << iphi << "\t" << std::bitset<14>(iphi) << std::endl;

      int ir = proj.fpgarzproj().value();
      std::cout << "ir=" << ir << "\t" << std::bitset<12>(ir) << std::endl;

      //TODO - need to express interms of constants
      int shifttmp2 = 7;
      int ircorr = (iz * proj.fpgarzprojder().value()) >> shifttmp2;
      std::cout << "ircorr=" << ircorr << "\t" << std::bitset<7>(ircorr) << std::endl;

      ir += ircorr;
      std::cout << "ir=" << ir << "\t" << std::bitset<14>(ir) << std::endl;

      std::cout << "stubphi=" << fpgastub->phi().value() << std::endl;
      int ideltaphi = fpgastub->phi().value() * settings_.kphi() / settings_.kphi() - iphi;

      int irstub = fpgastub->r().value();
      std::cout << "irstub=" << irstub << std::endl;
      int ialphafact = 0;
      if (!stub->isPSmodule()) {
        assert(irstub < (int)N_DSS_MOD * 2);
        if (abs(disk) <= 2) {
          ialphafact = ialphafactinner_[irstub];
          irstub = settings_.rDSSinner(irstub) / settings_.kr();
        } else {
          ialphafact = ialphafactouter_[irstub];
          irstub = settings_.rDSSouter(irstub) / settings_.kr();
        }
      }
      std::cout << "irstub (lookup)=" << irstub << std::endl;

      //TODO stub and projection r should not use different # bits...
      int ideltar = (irstub >> 1) - ir;
      std::cout << "ideltar=" << ideltar << "\t" << std::bitset<7>(ideltar) << std::endl;
      std::cout << "ideltaphi=" << ideltaphi << "\t" << std::bitset<7>(ideltaphi) << std::endl;

      if (!stub->isPSmodule()) {
        int ialpha = fpgastub->alpha().value();
        std::cout << "ialpha=" << ialpha << "\t" << std::bitset<4>(ialpha) << std::endl;
        int iphialphacor = ((ideltar * ialpha * ialphafact) >> settings_.alphashift());
        ideltaphi += iphialphacor;
        std::cout << "ialphafact=" << ialphafact << "\t" << std::bitset<10>(ialphafact) << std::endl;
        std::cout << "iphialphacor=" << iphialphacor << std::endl;
      }
      std::cout << "ideltaphi=" << ideltaphi << "\t" << std::bitset<7>(ideltaphi) << std::endl;

      //Perform floating point calculations here

      double phi = stub->phi() - phioffset;
      double z = stub->z();
      double r = stub->r();

      if (settings_.useapprox()) {
        double dphi = reco::reduceRange(phi - fpgastub->phiapprox(0.0, 0.0));
        assert(std::abs(dphi) < 0.001);
        phi = fpgastub->phiapprox(0.0, 0.0);
        z = fpgastub->zapprox();
        r = fpgastub->rapprox();
      }

      if (phi < 0)
        phi += 2 * M_PI;

      double dz = z - sign * settings_.zmean(layerdisk_ - N_LAYER);

      if (std::abs(dz) > settings_.dzmax()) {
        throw cms::Exception("LogicError") << __FILE__ << " " << __LINE__ << " " << name_ << " " << tracklet->getISeed()
                                           << "\n stub " << stub->z() << " disk " << disk << " " << dz;
      }

      double phiproj = proj.phiproj() + dz * proj.phiprojder();

      double rproj = proj.rzproj() + dz * proj.rzprojder();

      double deltar = r - rproj;

      double dr = stub->r() - rproj;

      double dphi = reco::reduceRange(phi - phiproj);

      double dphiapprox = reco::reduceRange(phi - (proj.phiprojapprox() + dz * proj.phiprojderapprox()));

      double drapprox = stub->r() - (proj.rzprojapprox() + dz * proj.rzprojderapprox());

      double drphi = dphi * stub->r();
      double drphiapprox = dphiapprox * stub->r();

      if (!stub->isPSmodule()) {
        double alphanorm = stub->alphanorm();
        dphi += dr * alphanorm * settings_.half2SmoduleWidth() / stub->r2();
        dphiapprox += drapprox * alphanorm * settings_.half2SmoduleWidth() / stub->r2();

        drphi += dr * alphanorm * settings_.half2SmoduleWidth() / stub->r();
        drphiapprox += dr * alphanorm * settings_.half2SmoduleWidth() / stub->r();
      }

      int seedindex = tracklet->getISeed();

      int idrphicut = rphicutPStable_.lookup(seedindex);
      int idrcut = rcutPStable_.lookup(seedindex);
      if (!stub->isPSmodule()) {
        idrphicut = rphicut2Stable_.lookup(seedindex);
        idrcut = rcut2Stable_.lookup(seedindex);
      }

      unsigned int projindex = mergedMatches[j].first.second;               // Allproj index
      curr_projid = next_projid;
      next_projid = projindex;
      bool newtracklet = (j == 0 || projindex != curr_projid);
      if (j == 0)  best_ideltar_disk = (1 << (fpgastub->r().nbits() - 1));
      if (newtracklet) {
        best_ideltaphi_disk = idrphicut;
        best_ideltar_disk = idrcut;
      }

      std::cout << "BEFORE idrphicut=" << idrphicut << std::endl;
      std::cout << "BEFORE idrcut=" << idrcut << std::endl;
      idrphicut = newtracklet ? idrphicut : best_ideltaphi_disk;
      idrcut = newtracklet ? idrcut : best_ideltar_disk;
      std::cout << "idrphicut=" << idrphicut << std::endl;
      std::cout << "idrcut=" << idrcut << std::endl;

      double drphicut = idrphicut * settings_.kphi() * settings_.kr();
      double drcut = idrcut * settings_.krprojshiftdisk();

      bool match, imatch;
      if (std::abs(dphi) < third * settings_.dphisectorHG() &&
          std::abs(dphiapprox) < third * settings_.dphisectorHG()) {  //1/3 of sector size to catch errors
        if (settings_.writeMonitorData("Residuals")) {
          double pt = 0.01 * settings_.c() * settings_.bfield() / std::abs(tracklet->rinv());

          globals_->ofstream("diskresiduals.txt")
              << disk << " " << stub->isPSmodule() << " " << tracklet->layer() << " " << abs(tracklet->disk()) << " "
              << pt << " " << ideltaphi * settings_.kphi() * stub->r() << " " << drphiapprox << " " << drphicut << " "
              << ideltar * settings_.krprojshiftdisk() << " " << deltar << " " << drcut << " " << endl;
        }

        match = (std::abs(drphi) < drphicut) && (std::abs(deltar) < drcut);

        imatch = (std::abs(ideltaphi) * irstub < best_ideltaphi_disk) && (std::abs(ideltar) < best_ideltar_disk);
        std::cout << "std::abs(ideltaphi)=" << std::abs(ideltaphi) << "\tstd::abs(ideltar)=" << std::abs(ideltar) << std::endl;
        std::cout << "std::abs(ideltaphi * irstub) < best_ideltaphi_disk=" << (std::abs(ideltaphi * irstub) <= best_ideltaphi_disk) << "\tstd::abs(ideltar) < best_ideltar_disk=" << (std::abs(ideltar) <= best_ideltar_disk) << std::endl;
        std::cout << "std::abs(ideltaphi) irstub " << std::abs(ideltaphi) << "\t" << irstub << std::endl;
        std::cout << "std::abs(ideltaphi) * irstub=" << std::abs(ideltaphi)*irstub << std::endl;
        std::cout << "best_ideltaphi_disk=" << best_ideltaphi_disk << std::endl;
        std::cout << "best_ideltar_disk=" << best_ideltar_disk << std::endl;
        std::cout << (std::abs(ideltaphi * irstub) < idrphicut) << "\t" << (std::abs(ideltar) < idrcut) << std::endl;
        std::cout << "imatch=" << imatch << std::endl;
        if (imatch) {
          best_ideltaphi_disk = std::abs(ideltaphi) * irstub;
          best_ideltar_disk = std::abs(ideltar);
        }
      } else {
        edm::LogProblem("Tracklet") << "WARNING dphi and/or dphiapprox too large : " << dphi << " " << dphiapprox
                                    << "dphi " << dphi << " Seed / ISeed " << tracklet->getISeed() << endl;
        match = false;
        imatch = false;
      }

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "imatch match disk: " << imatch << " " << match << " " << std::abs(ideltaphi)
                                     << " " << drphicut / (settings_.kphi() * stub->r()) << " " << std::abs(ideltar)
                                     << " " << drcut / settings_.krprojshiftdisk() << " r = " << stub->r();
      }

      if (imatch) {
        countsel++;
        std::cout << "imatch match disk: " << imatch << " " << match << " " << std::abs(ideltaphi)
                  << " " << drphicut / (settings_.kphi() * stub->r()) << " " << std::abs(ideltar)
                  << " " << drcut / settings_.krprojshiftdisk() << " r = " << stub->r() << std::endl;

        if (settings_.debugTracklet()) {
          edm::LogVerbatim("Tracklet") << "MatchCalculator found match in disk " << getName();
        }

        tracklet->addMatch(layerdisk_,
                           ideltaphi,
                           ideltar,
                           drphi / stub->r(),
                           dr,
                           drphiapprox / stub->r(),
                           drapprox,
                           (phiregion_ << 7) + fpgastub->stubindex().value(),
                           fpgastub);

        if (settings_.debugTracklet()) {
          edm::LogVerbatim("Tracklet") << "Accepted full match in disk " << getName() << " " << tracklet;
        }
        string stubid = mergedMatches[j].second->stubindex().str();  // stub ID
        int stubindex = 0;
        for(unsigned int i = 0; i < stubid.length(); i++) {
          if(stubid[i] == '1')
            stubindex += pow(2, i);
        }
        if(name_[3] == 'D' && name_.substr(name_.length()-4, name_.length()) == "PHIC") {
        std::cout << std::hex << "stub=" << trklet::hexFormat(fpgastub->str()) << std::endl;
        std::cout << std::hex << "proj=" << trklet::hexFormat(tracklet->trackletprojstrD(disk)) << std::endl;
	std::cout << name_ << std::endl;
	std::cout << "layerdisk_=" << layerdisk_ << std::endl;
        std::cout << std::dec << "disk=" << abs(disk) << "\t" << disks[abs(disk)] << std::endl;
        std::cout << std::dec << "stubid=" << stubindex << "\t" << std::bitset<15>(stubindex) << std::endl;
        std::cout << std::dec << "projindex=" << projindex << "\t" << std::bitset<15>(projindex) << std::endl;
        std::cout << "cmatch=" << std::bitset<7>(projindex) << "|" << stub->isPSmodule() << "|" << std::bitset<7>(stubindex) << "\t" << trklet::hexFormat(projindex + (stub->isPSmodule() ? "1" : "0") + stubid) << std::endl;
        std::cout << "j=" << j << std::endl;
        std::cout << "iz=" << iz << "\t" << std::bitset<7>(iz) << std::endl;
        std::cout << "iz bits=" << fpgastub->z().nbits() << std::endl;
        std::cout << "iphi=" << iphi - iphicorr << std::endl;
        std::cout << "iphicorr=" << iphicorr << std::endl;
        std::cout << "iphi=" << iphi << std::endl;
        std::cout << "ir=" << ir << std::endl;
        std::cout << "isPSStub=" << stub->isPSmodule() << std::endl;
        std::cout << "ircorr=" << ircorr << std::endl;
        std::cout << "ideltaphi=" << ideltaphi << std::endl;
        std::cout << "irstub=" << irstub << std::endl;
        std::cout << "ideltar=" << ideltar << std::endl;
        std::cout << std::hex << "idrphicut=" << idrphicut << std::endl;
        std::cout << std::hex << "idrcut=" << idrcut << std::endl;
        std::cout << std::dec << "imatch match disk: " << imatch << " " << match << " " << std::abs(ideltaphi)
                  << " " << drphicut / (settings_.kphi() * stub->r()) << " " << std::abs(ideltar)
                  << " " << drcut / settings_.krprojshiftdisk() << " r = " << stub->r() << std::endl;
        std::cout << "phiproj=" << phiproj << std::endl;
        std::cout << "rproj=" << rproj << std::endl;
        std::cout << "deltar=" << deltar << std::endl;
        std::cout << "dr=" << dr << std::endl;
        std::cout << "dphi=" << dphi << std::endl;
        std::cout << "dphiapprox=" << dphiapprox << std::endl;
        std::cout << "drapprox=" << drapprox << std::endl;
        std::cout << "drphi=" << drphi << std::endl;
        std::cout << "drphiapprox=" << drphiapprox << std::endl;
        std::cout << "drphicut=" << drphicut << std::endl;
        std::cout << "drcut=" << drcut << std::endl;
        }

        fullMatches_[seedindex]->addMatch(tracklet, mergedMatches[j].second);
      }
    }
    if (countall >= settings_.maxStep("MC"))
      break;
  }

  if (settings_.writeMonitorData("MC")) {
    globals_->ofstream("matchcalculator.txt") << getName() << " " << countall << " " << countsel << endl;
  }
}

std::vector<std::pair<std::pair<Tracklet*, int>, const Stub*> > MatchCalculator::mergeMatches(
    vector<CandidateMatchMemory*>& candmatch) {
  std::vector<std::pair<std::pair<Tracklet*, int>, const Stub*> > tmp;

  std::vector<unsigned int> indexArray;
  indexArray.reserve(candmatch.size());
  for (unsigned int i = 0; i < candmatch.size(); i++) {
    indexArray.push_back(0);
  }

  int bestIndex = -1;
  do {
    int bestSector = 100;
    int bestTCID = -1;
    bestIndex = -1;
    for (unsigned int i = 0; i < candmatch.size(); i++) {
      if (indexArray[i] >= candmatch[i]->nMatches()) {
        // skip as we were at the end
        continue;
      }
      int TCID = candmatch[i]->getMatch(indexArray[i]).first.first->TCID();
      int dSector = 0;
      if (dSector > 2)
        dSector -= N_SECTOR;
      if (dSector < -2)
        dSector += N_SECTOR;
      assert(abs(dSector) < 2);
      if (dSector == -1)
        dSector = 2;
      if (dSector < bestSector) {
        bestSector = dSector;
        bestTCID = TCID;
        bestIndex = i;
      }
      if (dSector == bestSector) {
        if (TCID < bestTCID || bestTCID < 0) {
          bestTCID = TCID;
          bestIndex = i;
        }
      }
    }
    if (bestIndex != -1) {
      tmp.push_back(candmatch[bestIndex]->getMatch(indexArray[bestIndex]));
      /*
      int stub(0);
      string stubidx = candmatch[bestIndex]->getMatch(indexArray[bestIndex]).second->stubindex().str();
      for(int icmatch = 0; icmatch < 7; icmatch++)
        if(stubidx[icmatch] == '1') stub += pow(2, icmatch);
      int projindex = candmatch[bestIndex]->getMatch(indexArray[bestIndex]).first.second;
      FPGAWord tmp;
      if (projindex >= (1 << 7)) {
        projindex = (1 << 7) - 1;
      }
      tmp.set(projindex, 7, true, __LINE__, __FILE__);
      std::cout << "Found cmatch " << 
                   tmp.str() << "|" <<
                   candmatch[bestIndex]->getMatch(indexArray[bestIndex]).second->stubindex().str() << "\t" <<
                   "projid=" << candmatch[bestIndex]->getMatch(indexArray[bestIndex]).first.second << "\t" <<
                   "stubid=" << stub << std::endl;
      */
      indexArray[bestIndex]++;
    }
  } while (bestIndex != -1);

  if (layerdisk_ < N_LAYER) {
    int lastTCID = -1;
    bool error = false;

    //Allow equal TCIDs since we can have multiple candidate matches
    for (unsigned int i = 1; i < tmp.size(); i++) {
      if (lastTCID > tmp[i].first.first->TCID()) {
        edm::LogProblem("Tracklet") << "Wrong TCID ordering for projections in " << getName() << " last " << lastTCID
                                    << " " << tmp[i].first.first->TCID();
        error = true;
      } else {
        lastTCID = tmp[i].first.first->TCID();
      }
    }

    if (error) {
      for (unsigned int i = 1; i < tmp.size(); i++) {
        edm::LogProblem("Tracklet") << "Wrong order for in " << getName() << " " << i << " " << tmp[i].first.first
                                    << " " << tmp[i].first.first->TCID();
      }
    }
  }

  for (unsigned int i = 0; i < tmp.size(); i++) {
    if (i > 0) {
      //This allows for equal TCIDs. This means that we can e.g. have a track seeded
      //in L1L2 that projects to both L3 and D4. The algorithm will pick up the first hit and
      //drop the second

      assert(tmp[i - 1].first.first->TCID() <= tmp[i].first.first->TCID());
    }
  }

  return tmp;
}
