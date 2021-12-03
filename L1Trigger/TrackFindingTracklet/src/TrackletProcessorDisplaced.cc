#include "L1Trigger/TrackFindingTracklet/interface/TrackletProcessorDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllInnerStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPairsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"
#include "L1Trigger/TrackFindingTracklet/interface/IMATH_TrackletCalculator.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <utility>
#include <tuple>

using namespace std;
using namespace trklet;

TrackletProcessorDisplaced::TrackletProcessorDisplaced(string name, Settings const& settings, Globals* globals)
    : TrackletCalculatorDisplaced(name, settings, globals)
      // tebuffer_(CircularBuffer<TEData>(3), 0, 0, 0, 0),
 {
  // iAllStub_ = -1;

  for (unsigned int ilayer = 0; ilayer < N_LAYER; ilayer++) {
    vector<TrackletProjectionsMemory*> tmp(settings_.nallstubs(ilayer), nullptr);
    trackletprojlayers_.push_back(tmp);
  }

  for (unsigned int idisk = 0; idisk < N_DISK; idisk++) {
    vector<TrackletProjectionsMemory*> tmp(settings_.nallstubs(idisk + N_LAYER), nullptr);
    trackletprojdisks_.push_back(tmp);
  }

  // outervmstubs_ = nullptr; 
  // initLayerDisksandISeed(layerdisk1_, layerdisk2_, iSeed_);

  layer_ = 0;
  disk_ = 0;
  layer1_ = 0;
  layer2_ = 0;
  layer3_ = 0;
  disk1_ = 0;
  disk2_ = 0;
  disk3_ = 0;

  string name1 = name.substr(1);  //this is to correct for "TPD" having one more letter then "TP"
  if (name1[3] == 'L')
    layer_ = name1[4] - '0';
  if (name1[3] == 'D')
    disk_ = name1[4] - '0';

  if (name1[3] == 'L')
    layer1_ = name1[4] - '0';
  if (name_[3] == 'D')
    disk1_ = name1[4] - '0';
  if (name_[5] == 'L')
    layer2_ = name1[6] - '0';
  if (name_[5] == 'D')
    disk2_ = name1[6] - '0';

  std::cout<< "Name: " << name1 << " |Layer: " << layer_ << " |Disk: " << disk_ << std::endl;
 
  // set TC index
  iSeed_ = 0;

  int iTC = name1[9] - 'A';

  if (name1.substr(3, 6) == "L3L4L2"){
    iSeed_ = 8;
    layer3_ = 2;
  }
  else if (name1.substr(3, 6) == "L5L6L4"){
    iSeed_ = 9;
    layer3_ = 4;
  }
  else if (name1.substr(3, 6) == "L2L3D1"){
    iSeed_ = 10;
    disk3_ = 1;
  }
  else if (name1.substr(3, 6) == "D1D2L2"){
    iSeed_ = 11;
    layer3_ = 2;
  }
  assert(iSeed_ != 0);

  firstphibits_ = settings_.nfinephi(0, iSeed_);

  if ((layer2_ == 4 && layer3_ == 2) || (layer2_ == 6 && layer3_ == 4)) {
    secondphibits_ = settings_.nfinephi(1, iSeed_);
    thirdphibits_ = settings_.nfinephi(2, iSeed_);
  }
  if ((layer2_ == 3 && disk3_ == 1) || (disk2_ == 2 && layer3_ == 2)) {
    secondphibits_ = settings_.nfinephi(1, iSeed_);
    thirdphibits_ = settings_.nfinephi(2, iSeed_);
  }


  TCIndex_ = (iSeed_ << 4) + iTC;
  assert(TCIndex_ >= 128 && TCIndex_ < 191);

  assert((layer_ != 0) || (disk_ != 0));

  toR_.clear();
  toZ_.clear();

  if (iSeed_ == 8 || iSeed_ == 9) {
    if (layer_ == 3) {
      rzmeanInv_[0] = 1.0 / settings_.rmean(2 - 1);
      rzmeanInv_[1] = 1.0 / settings_.rmean(3 - 1);
      rzmeanInv_[2] = 1.0 / settings_.rmean(4 - 1);

      rproj_[0] = settings_.rmean(0);
      rproj_[1] = settings_.rmean(4);
      rproj_[2] = settings_.rmean(5);
      lproj_[0] = 1;
      lproj_[1] = 5;
      lproj_[2] = 6;

      dproj_[0] = 1;
      dproj_[1] = 2;
      dproj_[2] = 0;
      toZ_.push_back(settings_.zmean(0));
      toZ_.push_back(settings_.zmean(1));
    }
    if (layer_ == 5) {
      rzmeanInv_[0] = 1.0 / settings_.rmean(4 - 1);
      rzmeanInv_[1] = 1.0 / settings_.rmean(5 - 1);
      rzmeanInv_[2] = 1.0 / settings_.rmean(6 - 1);

      rproj_[0] = settings_.rmean(0);
      rproj_[1] = settings_.rmean(1);
      rproj_[2] = settings_.rmean(2);
      lproj_[0] = 1;
      lproj_[1] = 2;
      lproj_[2] = 3;

      dproj_[0] = 0;
      dproj_[1] = 0;
      dproj_[2] = 0;
    }
    for (unsigned int i = 0; i < N_LAYER - 3; ++i)
      toR_.push_back(rproj_[i]);
  }

  if (iSeed_ == 10 || iSeed_ == 11) {
    if (layer_ == 2) {
      rzmeanInv_[0] = 1.0 / settings_.rmean(2 - 1);
      rzmeanInv_[1] = 1.0 / settings_.rmean(3 - 1);
      rzmeanInv_[2] = 1.0 / settings_.zmean(1 - 1);

      rproj_[0] = settings_.rmean(0);
      lproj_[0] = 1;
      lproj_[1] = -1;
      lproj_[2] = -1;

      zproj_[0] = settings_.zmean(1);
      zproj_[1] = settings_.zmean(2);
      zproj_[2] = settings_.zmean(3);
      dproj_[0] = 2;
      dproj_[1] = 3;
      dproj_[2] = 4;
    }
    if (disk_ == 1) {
      rzmeanInv_[0] = 1.0 / settings_.rmean(2 - 1);
      rzmeanInv_[1] = 1.0 / settings_.zmean(1 - 1);
      rzmeanInv_[2] = 1.0 / settings_.zmean(2 - 1);

      rproj_[0] = settings_.rmean(0);
      lproj_[0] = 1;
      lproj_[1] = -1;
      lproj_[2] = -1;

      zproj_[0] = settings_.zmean(2);
      zproj_[1] = settings_.zmean(3);
      zproj_[2] = settings_.zmean(4);
      dproj_[0] = 3;
      dproj_[1] = 4;
      dproj_[2] = 5;
    }
    toR_.push_back(settings_.rmean(0));
    for (unsigned int i = 0; i < N_DISK - 2; ++i)
      toZ_.push_back(zproj_[i]);
  }
 }

void TrackletProcessorDisplaced::addOutputProjection(TrackletProjectionsMemory*& outputProj, MemoryBase* memory) {
  outputProj = dynamic_cast<TrackletProjectionsMemory*>(memory);
  assert(outputProj != nullptr);
}

void TrackletProcessorDisplaced::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output " << output;
  }

  if (output == "trackpar") {
    auto* tmp = dynamic_cast<TrackletParametersMemory*>(memory);
    assert(tmp != nullptr);
    trackletpars_ = tmp;
    return;
  }

  if (output.substr(0, 7) == "projout") {
    //output is on the form 'projoutL2PHIC' or 'projoutD3PHIB'
    auto* tmp = dynamic_cast<TrackletProjectionsMemory*>(memory);
    assert(tmp != nullptr);

    unsigned int layerdisk = output[8] - '1';   //layer or disk counting from 0
    unsigned int phiregion = output[12] - 'A';  //phiregion counting from 0

    if (output[7] == 'L') {
      assert(layerdisk < N_LAYER);
      assert(phiregion < trackletprojlayers_[layerdisk].size());
      //check that phiregion not already initialized
      assert(trackletprojlayers_[layerdisk][phiregion] == nullptr);
      trackletprojlayers_[layerdisk][phiregion] = tmp;
      return;
    }

    if (output[7] == 'D') {
      assert(layerdisk < N_DISK);
      assert(phiregion < trackletprojdisks_[layerdisk].size());
      //check that phiregion not already initialized
      assert(trackletprojdisks_[layerdisk][phiregion] == nullptr);
      trackletprojdisks_[layerdisk][phiregion] = tmp;
      return;
    }
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find output : " << output;

}



void TrackletProcessorDisplaced::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input " << input;
  }

  if (input == "thirdallstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    innerallstubs_.push_back(tmp);
    return;
  }
  if (input == "firstallstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    middleallstubs_.push_back(tmp);
    return;
  }
  if (input == "secondallstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    outerallstubs_.push_back(tmp);
    return;
  }
  if (input.substr(0, 8) == "stubpair") {
    auto *tmp = dynamic_cast<StubPairsMemory*>(memory);
    assert(tmp != nullptr);
    stubpairs_.push_back(tmp);
    return;
  }

  // if (input.find("stubtriplet") == 0) {
  //   auto* tmp = dynamic_cast<StubTripletsMemory*>(memory);
  //   assert(tmp != nullptr);
  //   stubtriplets_.push_back(tmp);
  //   return;
  // }

  if (input == "thirdvmstubin") {
    auto* tmp = dynamic_cast<VMStubsTEMemory*>(memory);
    assert(tmp != nullptr);
    innervmstubs_.push_back(tmp);
    return;
  }
  if (input == "secondvmstubin") {
    auto* tmp = dynamic_cast<VMStubsTEMemory*>(memory);
    assert(tmp != nullptr);
    outervmstubs_.push_back(tmp);
    return;
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find input : " << input;

}


void TrackletProcessorDisplaced::execute(unsigned int iSector, double phimin, double phimax) {

  unsigned int countall = 0;
  unsigned int countsel = 0;
  unsigned int countpass = 0;
  unsigned int nThirdStubs = 0;
  unsigned int nInnerStubs = 0;
  count_ = 0;
  
  phimin_ = phimin;
  phimax_ = phimax;
  iSector_ = iSector;

  for (unsigned int iInnerMem = 0; iInnerMem < outervmstubs_.size();
       nInnerStubs += outervmstubs_.at(iInnerMem)->nVMStubs(), iInnerMem++
       );

  for (unsigned int iThirdMem = 0; iThirdMem < innervmstubs_.size();
       nThirdStubs += innervmstubs_.at(iThirdMem)->nVMStubs(), iThirdMem++
       );

  assert(!innerallstubs_.empty());
  assert(!middleallstubs_.empty());
  assert(!outerallstubs_.empty());
  assert(!innervmstubs_.empty());
  assert(!outervmstubs_.empty());
  assert(stubpairs_.empty());


  std::cout << "CHECKING EXECUTE TPD MODULE" << std::endl;

  for (auto& iInnerMem : outervmstubs_) {
    assert(iInnerMem->nVMStubs() == iInnerMem->nVMStubs());
    for (unsigned int i = 0; i < iInnerMem->nVMStubs(); i++) {
      const VMStubTE& firstvmstub = iInnerMem->getVMStubTE(i);
      if (settings_.debugTracklet()) {
	// edm::LogVerbatim("Tracklet") 
	std::cout
	  << "In " << getName() << " have first stub";
      }

      if ((layer1_ == 3 && layer2_ == 4) || (layer1_ == 5 && layer2_ == 6)) {
        int lookupbits = firstvmstub.vmbits().value() & 1023;
        int zdiffmax = (lookupbits >> 7);
        int newbin = (lookupbits & 127);
        int bin = newbin / 8;

        int zbinfirst = newbin & 7;

        int start = (bin >> 1);
        int last = start + (bin & 1);

        assert(last < 8);

        if (settings_.debugTracklet()) {
	  // edm::LogVerbatim("Tracklet")
	  std::cout
	    << "Will look in zbins " << start << " to " << last;
        }

	// for (auto& iOuterMem : middleallstubs_) {
	//   assert(iOuterMem->nStubs() == iOuterMem->nStubs());
	//   for (unsigned int j = 0; j < iOuterMem->nStubs(); j++) {

	//     if (settings_.debugTracklet()) {
	//       edm::LogVerbatim("Tracklet") << "In " << getName() << " have second stub(1) " << " " << j;
        //     }

	//     if (countall >= settings_.maxStep("TE"))
        //       break;
        //     countall++;
        //     const Stub& secondallstub = middleallstubs_->getStub(j);
	//   }

        for (int ibin = start; ibin <= last; ibin++) {
	  
	  for (auto& iOuterMem : middleallstubs_){
	    assert(iOuterMem->nStubs() == iOuterMem->nStubs());
	    for (unsigned int j = 0; j < iOuterMem->nStubs(); j++) {
	    
	      if (settings_.debugTracklet()) {
		// edm::LogVerbatim("Tracklet")
		std::cout
		  << "In " << getName() << " have second stub(1) " << ibin << " " << j;
	      }

	      if (countall >= settings_.maxStep("TE"))
		break;
	      countall++;
	      const Stub* secondallstub = iOuterMem->getStub(j);

	      int zbin = (secondallstub->z().value() & 7);
	      
	      if (start != ibin)
		zbin += 8;
	      if (zbin < zbinfirst || zbin - zbinfirst > zdiffmax) {
		if (settings_.debugTracklet()) {
		  // edm::LogVerbatim("Tracklet")
		  std::cout
		    << "Stubpair rejected because of wrong zbin";
		}
		continue;
	      }

	      assert(firstphibits_ != -1);
	      assert(secondphibits_ != -1);

	      // FPGAWord iphifirstbin = firstvmstub.finephi();
	      // FPGAWord iphisecondbin = secondallstub.finephi();

	      // unsigned int index = (iphifirstbin.value() << secondphibits_) + iphisecondbin.value();
	      
	      FPGAWord firstbend = firstvmstub.bend();
	      FPGAWord secondbend = secondallstub->bend();

	      unsigned int index = (index << firstbend.nbits()) + firstbend.value();
	      index = (index << secondbend.nbits()) + secondbend.value();

            if ((settings_.enableTripletTables() && !settings_.writeTripletTables())
		// && (index >= table_.size() || table_.at(index).empty())) 
		){
              if (settings_.debugTracklet()) {
	    	// edm::LogVerbatim("Tracklet")
		std::cout
	    	  << "Stub pair rejected because of stub pt cut bends : "
	    	  << settings_.benddecode(firstvmstub.bend().value(), layer1_ - 1, firstvmstub.isPSmodule()) << " "
	    	  << settings_.benddecode(secondallstub->bend().value(), layer2_ - 1, secondallstub->isPSmodule());
              }

              // FIXME temporarily commented out until stub bend table fixed
              // if (!settings_.writeTripletTables())
              //  continue;
            }

            // if (settings_.debugTracklet())
	    //   edm::LogVerbatim("Tracklet") << "Adding layer-layer pair in " << getName();
            // for (unsigned int isp = 0; isp < stubpairs_.size(); ++isp) {
            //   if (!settings_.enableTripletTables() || settings_.writeTripletTables() || table_.at(index).count(isp)) {
            //     if (settings_.writeMonitorData("Seeds")) {
            //       ofstream fout("seeds.txt", ofstream::app);
            //       fout << __FILE__ << ":" << __LINE__ << " " << name_ << " " << iSeed_ << endl;
            //       fout.close();
            //     }
            //     stubpairs_.at(isp)->addStubPair(firstvmstub, secondvmstub, index, getName());
            //   }
            // }

            countpass++;
	    }

	  }

	}

      } else if (layer1_ == 2 && layer2_ == 3) {
        int lookupbits = firstvmstub.vmbits().value() & 1023;
        int zdiffmax = (lookupbits >> 7);
        int newbin = (lookupbits & 127);
        int bin = newbin / 8;

        int zbinfirst = newbin & 7;

        int start = (bin >> 1);
        int last = start + (bin & 1);

        assert(last < 8);

	if (settings_.debugTracklet()) {
	  edm::LogVerbatim("Tracklet") << "Will look in zbins " << start << " to " << last;
        }
       
	for (int ibin = start; ibin <= last; ibin++) {
	  for (auto& iOuterMem : middleallstubs_){
            assert(iOuterMem->nStubs() == iOuterMem->nStubs());
            for (unsigned int j = 0; j < iOuterMem->nStubs(); j++) {

	      if (settings_.debugTracklet()) {
		edm::LogVerbatim("Tracklet") << "In " << getName() << " have second stub(2) ";
	      }

	      if (countall >= settings_.maxStep("TE"))
		break;
	      countall++;

	      const Stub* secondallstub = iOuterMem->getStub(j);

              int zbin = (secondallstub->z().value() & 7);

              if (start != ibin)
                zbin += 8;

	      if (zbin < zbinfirst || zbin - zbinfirst > zdiffmax) {
		if (settings_.debugTracklet()) {
		  edm::LogVerbatim("Tracklet") << "Stubpair rejected because of wrong zbin";
		}
		continue;
	      }

	      assert(firstphibits_ != -1);
	      assert(secondphibits_ != -1);

	      // FPGAWord iphifirstbin = firstvmstub.finephi();
	      // FPGAWord iphisecondbin = secondallstub.finephi();

	      // unsigned int index = (iphifirstbin.value() << secondphibits_) + iphisecondbin.value();
	      
	      FPGAWord firstbend = firstvmstub.bend();
	      FPGAWord secondbend = secondallstub->bend();

	      unsigned int index = (index << firstbend.nbits()) + firstbend.value();
	      index = (index << secondbend.nbits()) + secondbend.value();

	      if ((settings_.enableTripletTables() && !settings_.writeTripletTables())
		  // && (index >= table_.size() || table_.at(index).empty())) 
		  ){
		if (settings_.debugTracklet()) {
		  // edm::LogVerbatim("Tracklet")
		  std::cout
		    << "Stub pair rejected because of stub pt cut bends : "
		    << settings_.benddecode(firstvmstub.bend().value(), layer1_ - 1, firstvmstub.isPSmodule()) << " "
		    << settings_.benddecode(secondallstub->bend().value(), layer2_ - 1, secondallstub->isPSmodule());
		}

		continue;
	      }
		// if (settings_.debugTracklet())
		//   edm::LogVerbatim("Tracklet") << "Adding layer-layer pair in " << getName();
		// for (unsigned int isp = 0; isp < stubpairs_.size(); ++isp) {
		//   if (!settings_.enableTripletTables() || settings_.writeTripletTables() || table_.at(index).count(isp)) {
		//     if (settings_.writeMonitorData("Seeds")) {
		//       ofstream fout("seeds.txt", ofstream::app);
		//       fout << __FILE__ << ":" << __LINE__ << " " << name_ << " " << iSeed_ << endl;
		//       fout.close();
		//     }
		//     stubpairs_.at(isp)->addStubPair(firstvmstub, secondvmstub, index, getName());
		//   }
		// }

		countpass++;
	    }

	  }

	}

      }




    }

  }



}

 //  bool print = true;
//   int hacksum = 0;

//   if (print) {
//     edm::LogVerbatim("Tracklet") << "In TrackletProcessorDisplaced::execute : " << getName() << " " << nThirdStubs << ":";
//     std::cout << "In TrackletProcessorDisplaced::execute : " << getName() << " " << nThirdStubs << ":";
//     for (unsigned int i = 0; i < innervmstubs_.size(); ++i) {
//       edm::LogVerbatim("Tracklet") << innervmstubs_.at(i)->getName() << " " << innervmstubs_.at(i)->nVMStubs();
//       std::cout << innervmstubs_.at(i)->getName() << " " << innervmstubs_.at(i)->nVMStubs();
//     }

//     int s = 0;
//     std::string oss = "";
//     for (unsigned int i = 0; i < stubpairs_.size(); ++i) {
//       oss += std::to_string(stubpairs_.at(i)->nStubPairs());
//       oss += " ";
//       s += stubpairs_.at(i)->nStubPairs();
//     }
//     hacksum += nThirdStubs * s;
//     edm::LogVerbatim("Tracklet") << oss;
//     for (unsigned int i = 0; i < stubpairs_.size(); ++i) {
//       edm::LogVerbatim("Tracklet") << "                                          " << stubpairs_.at(i)->getName();
//       std::cout<< "                                          " << stubpairs_.at(i)->getName();
//     }
//   }

//   tmpSPTable_.clear();

//   for (unsigned int i = 0; i < stubpairs_.size(); ++i) {
//     for (unsigned int j = 0; j < stubpairs_.at(i)->nStubPairs(); ++j) {
//       if (print)
// 	edm::LogVerbatim("Tracklet") << "     *****    " << stubpairs_.at(i)->getName() << " "
//                                      << stubpairs_.at(i)->nStubPairs();

//       auto firstvmstub = stubpairs_.at(i)->getVMStub1(j);
//       auto secondvmstub = stubpairs_.at(i)->getVMStub2(j);

//       if ((layer2_ == 4 && layer3_ == 2) || (layer2_ == 6 && layer3_ == 4)) {
//         constexpr unsigned int vmbitshift = 10;
//         int lookupbits = (int)((firstvmstub.vmbits().value() >> vmbitshift) & 1023);  //1023=2^vmbitshift-1
//         int newbin = (lookupbits & 127);
//         int bin = newbin / 8;

//         int start = (bin >> 1);
//         int last = start + (bin & 1);

// 	for (int ibin = start; ibin <= last; ibin++) {
//           for (unsigned int k = 0; k < innervmstubs_.size(); k++) {
//             string vmsteSuffix = innervmstubs_.at(k)->getLastPartOfName();
//             vmsteSuffix = vmsteSuffix.substr(0, vmsteSuffix.find_last_of('n'));
//             if (stubpairs_.at(i)->getLastPartOfName() != vmsteSuffix)
//               continue;
//             for (unsigned int l = 0; l < innervmstubs_.at(k)->nVMStubsBinned(ibin); l++) {
//               if (settings_.debugTracklet()) {
// 		edm::LogVerbatim("Tracklet") << "In " << getName() << " have third stub";
//               }

//               if (countall >= settings_.maxStep("TRE"))
//                 break;
//               countall++;

//               const VMStubTE &thirdvmstub = innervmstubs_.at(k)->getVMStubTEBinned(ibin, l);

//               assert(secondphibits_ != -1);
//               assert(thirdphibits_ != -1);

//               unsigned int nvmsecond = settings_.nallstubs(layer2_ - 1) * settings_.nvmte(1, iSeed_);
//               unsigned int nvmbitssecond = nbits(nvmsecond);

//               FPGAWord iphisecondbin = secondvmstub.stub()->iphivmFineBins(nvmbitssecond, secondphibits_);

//               //currently not using same number of bits as in the TED
//               //assert(iphisecondbin==(int)secondvmstub.finephi());
//               FPGAWord iphithirdbin = thirdvmstub.finephi();

//               unsigned int index = (iphisecondbin.value() << thirdphibits_) + iphithirdbin.value();

//               FPGAWord secondbend = secondvmstub.bend();
//               FPGAWord thirdbend = thirdvmstub.bend();

//               index = (index << secondbend.nbits()) + secondbend.value();
//               index = (index << thirdbend.nbits()) + thirdbend.value();


// 	      if ((settings_.enableTripletTables() && !settings_.writeTripletTables()) &&
//                   (index >= table_.size() || !table_[index])) {
//                 if (settings_.debugTracklet()) {
// 		  edm::LogVerbatim("Tracklet")
// 		    << "Stub pair rejected because of stub pt cut bends : "
// 		    << settings_.benddecode(secondvmstub.bend().value(), layer2_ - 1, secondvmstub.isPSmodule())
// 		    << " " << settings_.benddecode(thirdvmstub.bend().value(), layer3_ - 1, thirdvmstub.isPSmodule());
//                 }

// 	      }
// 	      if (settings_.writeTripletTables()) {
//                 if (index >= table_.size())
//                   table_.resize(index + 1, false);
//                 table_[index] = true;

//                 const unsigned spIndex = stubpairs_.at(i)->getIndex(j);
//                 const string &tedName = stubpairs_.at(i)->getTEDName(j);
//                 if (!tmpSPTable_.count(tedName))
//                   tmpSPTable_[tedName];
//                 if (spIndex >= tmpSPTable_.at(tedName).size())
//                   tmpSPTable_.at(tedName).resize(spIndex + 1);
//                 tmpSPTable_.at(tedName).at(spIndex).push_back(stubpairs_.at(i)->getName());
//               }

//               if (settings_.debugTracklet())
// 		edm::LogVerbatim("Tracklet") << "Adding layer-layer pair in " << getName();
//               if (settings_.writeMonitorData("Seeds")) {
//                 ofstream fout("seeds.txt", ofstream::app);
//                 fout << __FILE__ << ":" << __LINE__ << " " << name_ << " " << iSeed_ << endl;
//                 fout.close();
//               }
//               stubtriplets_->addStubs(thirdvmstub.stub(),
//                                       (stubpairs_.at(i))->getVMStub1(j).stub(),
//                                       (stubpairs_.at(i))->getVMStub2(j).stub());

//               countpass++;
//             }
//           }
//         }

//       }


//       else if (disk2_ == 2 && layer3_ == 2) {
//         int lookupbits = (int)((firstvmstub.vmbits().value() >> 10) & 1023);
//         int newbin = (lookupbits & 127);
//         int bin = newbin / 8;

//         int start = (bin >> 1);
//         int last = start + (bin & 1);

//         if (firstvmstub.stub()->disk().value() < 0) {  //TODO - negative disk should come from memory
//           start = settings_.NLONGVMBINS() - last - 1;
//           last = settings_.NLONGVMBINS() - start - 1;
//         }

//         for (int ibin = start; ibin <= last; ibin++) {
//           for (unsigned int k = 0; k < innervmstubs_.size(); k++) {
//             string vmsteSuffix = innervmstubs_.at(k)->getLastPartOfName();
//             vmsteSuffix = vmsteSuffix.substr(0, vmsteSuffix.find_last_of('n'));
//             if (stubpairs_.at(i)->getLastPartOfName() != vmsteSuffix)
//               continue;

//             for (unsigned int l = 0; l < innervmstubs_.at(k)->nVMStubsBinned(ibin); l++) {
//               if (countall >= settings_.maxStep("TRE"))
//                 break;
//               countall++;

//               const VMStubTE &thirdvmstub = innervmstubs_.at(k)->getVMStubTEBinned(ibin, l);

//               assert(secondphibits_ != -1);
//               assert(thirdphibits_ != -1);

//               FPGAWord iphisecondbin = secondvmstub.finephi();
//               FPGAWord iphithirdbin = thirdvmstub.finephi();

//               unsigned int index = (iphisecondbin.value() << thirdphibits_) + iphithirdbin.value();

//               FPGAWord secondbend = secondvmstub.bend();
//               FPGAWord thirdbend = thirdvmstub.bend();

//               index = (index << secondbend.nbits()) + secondbend.value();
//               index = (index << thirdbend.nbits()) + thirdbend.value();


// 	      if ((settings_.enableTripletTables() && !settings_.writeTripletTables()) &&
//                   (index >= table_.size() || !table_[index])) {
//                 if (settings_.debugTracklet()) {
// 		  edm::LogVerbatim("Tracklet")
// 		    << "Stub triplet rejected because of stub pt cut bends : "
// 		    << settings_.benddecode(secondvmstub.bend().value(), disk2_ + 5, secondvmstub.isPSmodule()) << " "
// 		    << settings_.benddecode(thirdvmstub.bend().value(), layer3_ - 1, thirdvmstub.isPSmodule());
//                 }
//                 continue;
//               }
//               if (settings_.writeTripletTables()) {
//                 if (index >= table_.size())
//                   table_.resize(index + 1, false);
//                 table_[index] = true;

//                 const unsigned spIndex = stubpairs_.at(i)->getIndex(j);
//                 const string &tedName = stubpairs_.at(i)->getTEDName(j);
//                 if (!tmpSPTable_.count(tedName))
//                   tmpSPTable_[tedName];
//                 if (spIndex >= tmpSPTable_.at(tedName).size())
//                   tmpSPTable_.at(tedName).resize(spIndex + 1);
//                 tmpSPTable_.at(tedName).at(spIndex).push_back(stubpairs_.at(i)->getName());
//               }

//               if (settings_.debugTracklet())
// 		edm::LogVerbatim("Tracklet") << "Adding layer-disk pair in " << getName();
//               if (settings_.writeMonitorData("Seeds")) {
//                 ofstream fout("seeds.txt", ofstream::app);
//                 fout << __FILE__ << ":" << __LINE__ << " " << name_ << " " << iSeed_ << endl;
//                 fout.close();
//               }
//               stubtriplets_->addStubs(thirdvmstub.stub(),
//                                       (stubpairs_.at(i))->getVMStub1(j).stub(),
//                                       (stubpairs_.at(i))->getVMStub2(j).stub());
//               countpass++;
//             }
//           }
//         }
//       }

//       else if (layer2_ == 3 && disk3_ == 1) {
//         int lookupbits = (int)((firstvmstub.vmbits().value() >> 10) & 1023);

//         int newbin = (lookupbits & 127);
//         int bin = newbin / 8;

//         int start = (bin >> 1);
//         int last = start + (bin & 1);

//         for (int ibin = start; ibin <= last; ibin++) {
//           for (unsigned int k = 0; k < innervmstubs_.size(); k++) {
//             string vmsteSuffix = innervmstubs_.at(k)->getLastPartOfName();
//             vmsteSuffix = vmsteSuffix.substr(0, vmsteSuffix.find_last_of('n'));
//             if (stubpairs_.at(i)->getLastPartOfName() != vmsteSuffix)
//               continue;
//             for (unsigned int l = 0; l < innervmstubs_.at(k)->nVMStubsBinned(ibin); l++) {
//               if (countall >= settings_.maxStep("TRE"))
//                 break;
//               countall++;

//               const VMStubTE &thirdvmstub = innervmstubs_.at(k)->getVMStubTEBinned(ibin, l);

//               assert(secondphibits_ != -1);
//               assert(thirdphibits_ != -1);

//               unsigned int nvmsecond;

//               nvmsecond = settings_.nallstubs(layer2_ - 1) * settings_.nvmte(1, iSeed_);
//               unsigned int nvmbitssecond = nbits(nvmsecond);

//               FPGAWord iphisecondbin = secondvmstub.stub()->iphivmFineBins(nvmbitssecond, secondphibits_);

//               //currentlty not using same number of bits as in the TED
//               //assert(iphisecondbin==(int)secondvmstub.finephi());
//               FPGAWord iphithirdbin = thirdvmstub.finephi();

//               unsigned int index = (iphisecondbin.value() << thirdphibits_) + iphithirdbin.value();

//               FPGAWord secondbend = secondvmstub.bend();
//               FPGAWord thirdbend = thirdvmstub.bend();

//               index = (index << secondbend.nbits()) + secondbend.value();
//               index = (index << thirdbend.nbits()) + thirdbend.value();

// 	      if ((settings_.enableTripletTables() && !settings_.writeTripletTables()) &&
//                   (index >= table_.size() || !table_[index])) {
//                 if (settings_.debugTracklet()) {
// 		  edm::LogVerbatim("Tracklet")
// 		    << "Stub pair rejected because of stub pt cut bends : "
// 		    << settings_.benddecode(secondvmstub.bend().value(), layer2_ - 1, secondvmstub.isPSmodule())
// 		    << " " << settings_.benddecode(thirdvmstub.bend().value(), disk3_ + 5, thirdvmstub.isPSmodule());
//                 }
//                 continue;
//               }
//               if (settings_.writeTripletTables()) {
//                 if (index >= table_.size())
//                   table_.resize(index + 1, false);
//                 table_[index] = true;

//                 const unsigned spIndex = stubpairs_.at(i)->getIndex(j);
//                 const string &tedName = stubpairs_.at(i)->getTEDName(j);
//                 if (!tmpSPTable_.count(tedName))
//                   tmpSPTable_[tedName];
//                 if (spIndex >= tmpSPTable_.at(tedName).size())
//                   tmpSPTable_.at(tedName).resize(spIndex + 1);
//                 tmpSPTable_.at(tedName).at(spIndex).push_back(stubpairs_.at(i)->getName());
//               }

//               if (settings_.debugTracklet())
// 		edm::LogVerbatim("Tracklet") << "Adding layer-disk pair in " << getName();
//               if (settings_.writeMonitorData("Seeds")) {
//                 ofstream fout("seeds.txt", ofstream::app);
//                 fout << __FILE__ << ":" << __LINE__ << " " << name_ << " " << iSeed_ << endl;
//                 fout.close();
//               }
//               stubtriplets_->addStubs(thirdvmstub.stub(),
//                                       (stubpairs_.at(i))->getVMStub1(j).stub(),
//                                       (stubpairs_.at(i))->getVMStub2(j).stub());
//               countpass++;
//             }
//           }
//         }
//       }

//     }
//   }

//   for (const auto &tedName : tmpSPTable_) {
//     for (unsigned spIndex = 0; spIndex < tedName.second.size(); spIndex++) {
//       if (tedName.second.at(spIndex).empty())
//         continue;
//       vector<string> entry(tedName.second.at(spIndex));
//       sort(entry.begin(), entry.end());
//       entry.erase(unique(entry.begin(), entry.end()), entry.end());
//       const string &spName = entry.at(0);

//       if (!spTable_.count(tedName.first))
//         spTable_[tedName.first];
//       if (spIndex >= spTable_.at(tedName.first).size())
//         spTable_.at(tedName.first).resize(spIndex + 1);
//       if (!spTable_.at(tedName.first).at(spIndex).count(spName))
//         spTable_.at(tedName.first).at(spIndex)[spName] = 0;
//       spTable_.at(tedName.first).at(spIndex)[spName]++;
//     }
//   }

//   if (settings_.writeMonitorData("TPD")) {
//     globals_->ofstream("trackletprocessordisplaced.txt") << getName() << " " << countall << " " << countpass << endl;
//   }


// }

