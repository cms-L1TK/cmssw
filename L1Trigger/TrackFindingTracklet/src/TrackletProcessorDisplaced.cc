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
  : TrackletCalculatorDisplaced(name, settings, globals),
    // tebuffer_(CircularBuffer<TEData>(3), 0, 0, 0, 0),
    innerTable_(settings),
    innerOverlapTable_(settings),
    innerThirdTable_(settings) 
 {
   middleallstubs_.clear();
   outervmstubs_ = nullptr;

  // iAllStub_ = -1;
   layerdisk_ = initLayerDisk(4);

   unsigned int region = name.substr(1)[9] - 'A';
   std::cout << "region: " << region <<std::endl;
   // assert(region < settings_.nallstubs(layerdisk_));

   if (layerdisk_ == LayerDisk::L1 || layerdisk_ == LayerDisk::L2 || layerdisk_ == LayerDisk::L3 ||
       layerdisk_ == LayerDisk::L5 || layerdisk_ == LayerDisk::D1 || layerdisk_ == LayerDisk::D3) {
     innerTable_.initVMRTable(layerdisk_, TrackletLUT::VMRTableType::inner, region);  //projection to next layer/disk
   }

   if (layerdisk_ == LayerDisk::L1 || layerdisk_ == LayerDisk::L2) {
     innerOverlapTable_.initVMRTable(layerdisk_, TrackletLUT::VMRTableType::inneroverlap, region);  //projection to disk from layer
   }

   if (layerdisk_ == LayerDisk::L2 || layerdisk_ == LayerDisk::L3 || layerdisk_ == LayerDisk::L5 ||
       layerdisk_ == LayerDisk::D1) {
     innerThirdTable_.initVMRTable(layerdisk_, TrackletLUT::VMRTableType::innerthird, region);  //projection to third layer/disk
   }

   nbitszfinebintable_ = settings_.vmrlutzbits(layerdisk_);
   nbitsrfinebintable_ = settings_.vmrlutrbits(layerdisk_);

  for (unsigned int ilayer = 0; ilayer < N_LAYER; ilayer++) {
    vector<TrackletProjectionsMemory*> tmp(settings_.nallstubs(ilayer), nullptr);
    trackletprojlayers_.push_back(tmp);
  }

  for (unsigned int idisk = 0; idisk < N_DISK; idisk++) {
    vector<TrackletProjectionsMemory*> tmp(settings_.nallstubs(idisk + N_LAYER), nullptr);
    trackletprojdisks_.push_back(tmp);
  }

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
    outervmstubs_ = tmp;
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
  unsigned int nOuterStubs = 0;
  count_ = 0;
  
  phimin_ = phimin;
  phimax_ = phimax;
  iSector_ = iSector;

  for (unsigned int iInnerMem = 0; iInnerMem < middleallstubs_.size();
       nInnerStubs += middleallstubs_.at(iInnerMem)->nStubs(), iInnerMem++
       );

  // for (unsigned int iThirdMem = 0; iThirdMem < innervmstubs_.size();
  //      nThirdStubs += innervmstubs_.at(iThirdMem)->nVMStubs(), iThirdMem++
  //      );

  assert(!innerallstubs_.empty());
  assert(!middleallstubs_.empty());
  assert(!outerallstubs_.empty());
  assert(!innervmstubs_.empty());
  assert(outervmstubs_ != nullptr);
  assert(stubpairs_.empty());


  std::cout << "\nRUNNING ON TPD::EXECUTE\n" << std::endl;

  std::cout << "layer 1: " << layer1_ << " layer 2: " << layer2_ << " layer 3: " <<  layer3_ << std::endl;
  std::cout << "disk 1: " << disk1_ << " disk 2: " << disk2_ << " disk 3: " <<  disk3_ << std::endl;

  for (auto& iInnerMem : middleallstubs_){
    assert(iInnerMem->nStubs() == iInnerMem->nStubs());
    for (unsigned int j = 0; j < iInnerMem->nStubs(); j++) {
      const Stub* firstallstub = iInnerMem->getStub(j);

      if (settings_.debugTracklet()) {
        // edm::LogVerbatim("Tracklet")                                                         
	std::cout
          << "In " << getName() << " have first stub\n";
      }

      int inner = 0;
      bool negdisk = (firstallstub->disk().value() < 0);
      int indexz = (((1 << (firstallstub->z().nbits() - 1)) + firstallstub->z().value()) >> (firstallstub->z().nbits() - nbitszfinebintable_));
      int indexr = -1;
      if (layerdisk_ > (N_LAYER - 1)) {
        if (negdisk) {
          indexz = (1 << nbitszfinebintable_) - indexz;
        }
        indexr = firstallstub->r().value();
        if (firstallstub->isPSmodule()) {
          indexr = firstallstub->r().value() >> (firstallstub->r().nbits() - nbitsrfinebintable_);
        }
      } else {
        //Take the top nbitsfinebintable_ bits of the z coordinate. The & is to handle the negative z values.
        indexr = (((1 << (firstallstub->r().nbits() - 1)) + firstallstub->r().value()) >> (firstallstub->r().nbits() - nbitsrfinebintable_));
      }

      assert(indexz >= 0);
      assert(indexr >= 0);
      assert(indexz < (1 << nbitszfinebintable_));
      assert(indexr < (1 << nbitsrfinebintable_));

      // int melut = meTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
      // assert(melut >= 0);

      unsigned int lutwidth = settings_.lutwidthtab(inner, iSeed_);
      if (settings_.extended()) {
	lutwidth = settings_.lutwidthtabextended(inner, iSeed_);
      }

      int lutval = -999;

      if (iSeed_ < Seed::L1D1 || iSeed_ > Seed::L2D1) {
	lutval = innerTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
      } else {
	lutval = innerOverlapTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
      }

      if (lutval == -1)
	continue;
      if (settings_.extended() &&
	  (iSeed_ == Seed::L3L4 || iSeed_ == Seed::L5L6 || iSeed_ == Seed::D1D2 || iSeed_ == Seed::L2L3D1)) {
	int lutval2 = innerThirdTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
	if (lutval2 != -1)
	  lutval += (lutval2 << 10);
      }

      assert(lutval >= 0);
      // assert(lutwidth > 0);

      std::cout << "lutval: " << lutval << " lutwidth: " << lutwidth <<std::endl;

      FPGAWord binlookup(lutval, lutwidth, true, __LINE__, __FILE__);
      FPGAWord finephi = firstallstub->iphivmFineBins(settings_.nphireg(inner, iSeed_), settings_.nfinephi(inner, iSeed_));

      if ((layer1_ == 3 && layer2_ == 4) || (layer1_ == 5 && layer2_ == 6)) {
	
        if (settings_.debugTracklet())
	  // edm::LogVerbatim("Tracklet")
	  std::cout
	    << getName() << " Layer-layer pair";

	int lookupbits = binlookup.value() & 1023;
        int zdiffmax = (lookupbits >> 7);
        int newbin = (lookupbits & 127);
        int bin = newbin / 8;

        int zbinfirst = newbin & 7;

        int start = (bin >> 1);
        int last = start + (bin & 1);

	// std::cout << "start" << start << "last" << last << std::endl;

        assert(last < 8);

	if (settings_.debugTracklet()) {
          // edm::LogVerbatim("Tracklet")                                                       
	  std::cout
            << "Will look in zbins " << start << " to " << last;
        }

	for (int ibin = start; ibin <= last; ibin++) {
	  for (unsigned int j = 0; j < outervmstubs_->nVMStubsBinned(ibin); j++) {

	    if (settings_.debugTracklet()) {
	      // edm::LogVerbatim("Tracklet")
	      std::cout
		<< "In " << getName() << " have second stub(1) " << ibin << " " << j;
	    }

	    if (countall >= settings_.maxStep("TE"))
	      break;
	    countall++;
	    const VMStubTE& secondvmstub = outervmstubs_->getVMStubTEBinned(ibin,j);

	    int zbin = (secondvmstub.vmbits().value() & 7);
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

	    FPGAWord iphifirstbin = finephi;
	    FPGAWord iphisecondbin = secondvmstub.finephi();

	    unsigned int index = (iphifirstbin.value() << secondphibits_) + iphisecondbin.value();
	      
	    FPGAWord firstbend = firstallstub->bend();
	    FPGAWord secondbend = secondvmstub.bend();

	    index = (index << firstbend.nbits()) + firstbend.value();
	    index = (index << secondbend.nbits()) + secondbend.value();

	    if ((settings_.enableTripletTables() && !settings_.writeTripletTables())){  
	      // && (index >= table_.size() || table_.at(index).empty())) {
	      if (settings_.debugTracklet()) {
		// edm::LogVerbatim("Tracklet")
		std::cout
		  << "Stub pair rejected because of stub pt cut bends : "
		  << settings_.benddecode(firstallstub->bend().value(), layer1_ - 1, firstallstub->isPSmodule()) << " "
		  << settings_.benddecode(secondvmstub.bend().value(), layer2_ - 1, secondvmstub.isPSmodule());
	      }

	      //FIXME temporarily commented out until stub bend table fixed
	      //if (!settings_.writeTripletTables())
	      //  continue;
	    }

	    countpass++;
		
	  }

	}

      } else if (layer1_ == 2 && layer2_ == 3) {

        if (settings_.debugTracklet())
	  // edm::LogVerbatim("Tracklet")
	  std::cout
	    << getName() << " Layer-layer pair";

        int lookupbits = binlookup.value() & 1023;
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

	for (int ibin = start; ibin <= last; ibin++) {
	  for (unsigned int j = 0; j < outervmstubs_->nVMStubsBinned(ibin); j++) {

	    if (settings_.debugTracklet()) {
	      // edm::LogVerbatim("Tracklet")
	      std::cout
		<< "In " << getName() << " have second stub(1) " << ibin << " " << j;
	    }

	    if (countall >= settings_.maxStep("TE"))
	      break;
	    countall++;
	    const VMStubTE& secondvmstub = outervmstubs_->getVMStubTEBinned(ibin,j);

	    int zbin = (secondvmstub.vmbits().value() & 7);
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

	    FPGAWord iphifirstbin = finephi;
	    FPGAWord iphisecondbin = secondvmstub.finephi();

	    unsigned int index = (iphifirstbin.value() << secondphibits_) + iphisecondbin.value();
	    
	    FPGAWord firstbend = firstallstub->bend();
	    FPGAWord secondbend = secondvmstub.bend();

	    index = (index << firstbend.nbits()) + firstbend.value();
	    index = (index << secondbend.nbits()) + secondbend.value();

	    if ((settings_.enableTripletTables() && !settings_.writeTripletTables())) { 
	      // && (index >= table_.size() || table_.at(index).empty())) {
	      if (settings_.debugTracklet()) {
		// edm::LogVerbatim("Tracklet")
		std::cout
		  << "Stub pair rejected because of stub pt cut bends : "
		  << settings_.benddecode(firstallstub->bend().value(), layer1_ - 1, firstallstub->isPSmodule()) << " "
		  << settings_.benddecode(secondvmstub.bend().value(), layer2_ - 1, secondvmstub.isPSmodule());
		
	      }
	
	      continue;
	  
	    }

	  }
      
	}

      } else if (disk1_ == 1 && disk2_ == 2) {
        if (settings_.debugTracklet())
	  // edm::LogVerbatim("Tracklet")
	  std::cout
	    << getName() << " Disk-disk pair";


	int lookupbits = binlookup.value() & 511;
        bool negdisk = firstallstub->disk().value() < 0;
        int rdiffmax = (lookupbits >> 6);
        int newbin = (lookupbits & 63);
        int bin = newbin / 8;

        int rbinfirst = newbin & 7;

        int start = (bin >> 1);
        if (negdisk)
          start += 4;
        int last = start + (bin & 1);
        assert(last < 8);

	for (int ibin = start; ibin <= last; ibin++) {
          if (settings_.debugTracklet()) {
	    std::cout
	    // edm::LogVerbatim("Tracklet") 
	      << getName() << " looking for matching stub in " << outervmstubs_->getName()
	      << " in bin = " << ibin << " with " << outervmstubs_->nVMStubsBinned(ibin)
	      << " stubs";
          }
          for (unsigned int j = 0; j < outervmstubs_->nVMStubsBinned(ibin); j++) {
            if (countall >= settings_.maxStep("TE"))
              break;
            countall++;

            const VMStubTE& secondvmstub = outervmstubs_->getVMStubTEBinned(ibin, j);

            int rbin = (secondvmstub.vmbits().value() & 7);
            if (start != ibin)
              rbin += 8;
            if (rbin < rbinfirst)
              continue;
            if (rbin - rbinfirst > rdiffmax)
              continue;

            unsigned int irsecondbin = secondvmstub.vmbits().value() >> 2;

            FPGAWord iphifirstbin = finephi;
            FPGAWord iphisecondbin = secondvmstub.finephi();

            unsigned int index = (irsecondbin << (secondphibits_ + firstphibits_)) +
	      (iphifirstbin.value() << secondphibits_) + iphisecondbin.value();
	    
            FPGAWord firstbend = firstallstub->bend();
            FPGAWord secondbend = secondvmstub.bend();

            index = (index << firstbend.nbits()) + firstbend.value();
            index = (index << secondbend.nbits()) + secondbend.value();

            if ((settings_.enableTripletTables() && !settings_.writeTripletTables())) { 
		// && (index >= table_.size() || table_.at(index).empty())) {
              if (settings_.debugTracklet()) {
		// edm::LogVerbatim("Tracklet")
		std::cout
		  << "Stub pair rejected because of stub pt cut bends : "
		  << settings_.benddecode(firstallstub->bend().value(), disk1_ + 5, firstallstub->isPSmodule()) << " "
		  << settings_.benddecode(secondvmstub.bend().value(), disk2_ + 5, secondvmstub.isPSmodule());
              }
              continue;
            }
	    
	  }

	}
	
      }
	  
    }

  }

}


