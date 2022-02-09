#include "L1Trigger/TrackFindingTracklet/interface/MatchEngineUnit.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletLUT.h"

using namespace std;
using namespace trklet;

MatchEngineUnit::MatchEngineUnit(bool barrel, unsigned int layerdisk, const TrackletLUT& luttable)
    : luttable_(luttable), candmatches_(3) {
  idle_ = true;
  barrel_ = barrel;
  layerdisk_ = layerdisk;
  goodpair_ = false;
  goodpair__ = false;
  havepair_ = false;
  good__ = false;
  good___ = false;
}

void MatchEngineUnit::setAlmostFull() {
  almostfullsave_ = candmatches_.nearfull();
}


void MatchEngineUnit::init(VMStubsMEMemory* vmstubsmemory,
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
                           Tracklet* proj,
                           bool,
			   int imeu) {
  vmstubsmemory_ = vmstubsmemory;
  idle_ = false;
  imeu_ = imeu;
  nrzbins_ = nrzbins;
  rzbin_ = rzbin;
  phibin_ = phibin;
  shift_ = shift;
  istub_ = 0;
  iuse_ = 0;
  projrinv_ = projrinv;
  projfinerz_ = projfinerz;
  projfinephi_ = projfinephi;
  use_.clear();
  if (usefirstMinus) {
    use_.emplace_back(0, 0);
  }
  if (usesecondMinus) {
    use_.emplace_back(1, 0);
  }
  if (usefirstPlus) {
    use_.emplace_back(0, 1);
  }
  if (usesecondPlus) {
    use_.emplace_back(1, 1);
  }
  assert(!use_.empty());
  isPSseed_ = isPSseed;
  proj_ = proj;

  //Even when you init a new projection you need to process the pipeline
  //This should be fixed to be done more cleanly - but require synchronizaton
  //with the HLS code

  goodpair__ = goodpair_;
  tmppair__ =  tmppair_;

  havepair_ = false;
  goodpair_ = false;

  good__ = false;

}


void MatchEngineUnit::step(bool print) {

  goodpair__ = goodpair_;
  tmppair__ =  tmppair_;


  havepair_ = false;
  goodpair_ = false;

  good__ = !idle() && !almostfullsave_;

  if (idle() || almostfullsave_)
    return;

  unsigned int slot = (phibin_ + use_[iuse_].second) * nrzbins_ + rzbin_ + use_[iuse_].first;

  projfinerz__ = projfinerz_ - (1 << NFINERZBITS) * use_[iuse_].first;
  projfinephi__ = projfinephi_;
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

  vmstub__ = vmstubsmemory_->getVMStubMEBin(slot, istub_);

  isPSseed__ = isPSseed_;
  projrinv__ = projrinv_;
  proj__ = proj_;

   
  //
  bool isPSmodule = vmstub__.isPSmodule();
  int stubfinerz = vmstub__.finerz().value();
  int stubfinephi = vmstub__.finephi().value();

  int deltaphi = stubfinephi - projfinephi__;

  bool dphicut = (abs(deltaphi) < 3);

  int nbits = isPSmodule ? 3 : 4;

  int diskps = (!barrel_) && isPSmodule;

  unsigned int index = (diskps << (4 + 5)) + (projrinv_ << nbits) + vmstub__.bend().value();

  //cout << "["<<imeu_<<"]OLD index: "<<index<<endl;


  //Check if stub z position consistent
  int idrz = stubfinerz - projfinerz__;
  bool pass;

  if (barrel_) {
    if (isPSseed__) {
      pass = idrz >= -1 && idrz <= 1;
    } else {
      pass = idrz >= -5 && idrz <= 5;
    }
  } else {
    if (isPSmodule) {
      pass = idrz >= -1 && idrz <= 1;
    } else {
      pass = idrz >= -3 && idrz <= 3;
    }
  }

  // Detailed printout for comparison with HLS code
  if (print)
    edm::LogVerbatim("Tracklet") << "MEU TrkId stubindex : " << 128 * proj_->TCIndex() + proj_->trackletIndex() << " "
                                 << vmstub__.stubindex().value() << "   "
                                 << ((pass && dphicut) && luttable_.lookup(index)) << " index=" << index
                                 << " projrinv bend : " << projrinv_ << " " << vmstub__.bend().value()
                                 << "  shift_ isPSseed_ :" << shift_ << " " << isPSseed_ << " slot=" << slot;

  //Check if stub bend and proj rinv consistent

  goodpair_ = (pass && dphicut) && luttable_.lookup(index);
  havepair_ = true;

  if (havepair_) {
    std::pair<Tracklet*, const Stub*> tmppair(proj_, vmstub__.stub());
    tmppair_ = tmppair;
  }

  //
  

  if (good__) {
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

}


void MatchEngineUnit::processPipeline() {

  if (good___) {
  
    bool isPSmodule = vmstub___.isPSmodule();
    int stubfinerz = vmstub___.finerz().value();
    int stubfinephi = vmstub___.finephi().value();
    
    int deltaphi = stubfinephi - projfinephi___;

    bool dphicut = (abs(deltaphi) < 3);

    int nbits = isPSmodule ? 3 : 4;

    int diskps = (!barrel_) && isPSmodule;

    unsigned int index = (diskps << (4 + 5)) + (projrinv___ << nbits) + vmstub___.bend().value();

    //Check if stub z position consistent
    int idrz = stubfinerz - projfinerz___;
    bool pass;

    if (barrel_) {
      if (isPSseed___) {
	pass = idrz >= -1 && idrz <= 1;
      } else {
	pass = idrz >= -5 && idrz <= 5;
      }
    } else {
      if (isPSmodule) {
	pass = idrz >= -1 && idrz <= 1;
      } else {
	pass = idrz >= -3 && idrz <= 3;
      }
    }

    // Detailed printout for comparison with HLS code
    //if (print)
    //  edm::LogVerbatim("Tracklet") << "MEU TrkId stubindex : " << 128 * proj_->TCIndex() + proj_->trackletIndex() << " "
    //                               << vmstub__.stubindex().value() << "   "
    //                               << ((pass && dphicut) && luttable_.lookup(index)) << " index=" << index
    //                               << " projrinv bend : " << projrinv_ << " " << vmstub__.bend().value()
    //                               << "  shift_ isPSseed_ :" << shift_ << " " << isPSseed_ << " slot=" << slot;
    
    //Check if stub bend and proj rinv consistent

    //cout << "["<<imeu_<<"]good___ index size : "<<good___<<" "<<index<<" "<<luttable_.size()<<endl;
    
    bool goodpair = (pass && dphicut) && luttable_.lookup(index);
    //cout << "goodpair : "<<goodpair<<endl;

    std::pair<Tracklet*, const Stub*> tmppair(proj___, vmstub___.stub());
    
    if (goodpair) {
      candmatches_.store(tmppair);
    }

 
  }

  proj___ = proj__;
  projfinephi___ = projfinephi__;
  projfinerz___ = projfinerz__;
  projrinv___ = projrinv__;
  isPSseed___ = isPSseed__;
  good___ =  good__;
  vmstub___ = vmstub__;


}

void MatchEngineUnit::reset() {
  candmatches_.reset();
  idle_ = true;
  istub_ = 0;
  goodpair_ = false;
  havepair_ = false;
  good__ = false;
  good___ = false;
}

int MatchEngineUnit::TCID() const {
  if (!empty()) {
    return peek().first->TCID();
  }

  if (good___) {
    return proj___->TCID();
  }

  if (good__) {
    return proj__->TCID();
  }

  //if (idle_ && !havepair_) {
  //  return 16383;
  //}

  //if (havepair_) {
  //  return tmppair_.first->TCID();
  //}


  if (idle_) {
    return 16383;
  }


  return proj_->TCID();
}
