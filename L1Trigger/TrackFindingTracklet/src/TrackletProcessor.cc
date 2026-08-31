#include "L1Trigger/TrackFindingTracklet/interface/TrackletProcessor.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllInnerStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/HistBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "L1Trigger/L1TCommon/interface/BitShift.h"

#include <utility>
#include <tuple>
#include <sstream>

using namespace std;
using namespace trklet;

TrackletProcessor::TrackletProcessor(string name, Settings const& settings, Globals* globals)
  : ProcessBase(name, settings, globals),
    tebuffer_(CircularBuffer<TEData>(3), 0, 0, 0, 0),
    pttableinner_(settings),
    pttableouter_(settings),
    useregiontable_(settings),
    innerTable_(settings),
    innerOverlapTable_(settings) {
  iAllStub_ = -1;

  outervmstubs_ = nullptr;

  initLayerDisksandISeed(layerdisk1_, layerdisk2_, iSeed_);

  double rmin = -1.0;
  double rmax = -1.0;

  if (iSeed_ == Seed::L1L2 || iSeed_ == Seed::L2L3 || iSeed_ == Seed::L3L4 || iSeed_ == Seed::L5L6) {
    rmin = settings_.rmean(layerdisk1_);
    rmax = settings_.rmean(layerdisk2_);
  } else {
    if (iSeed_ == Seed::L1D1) {
      rmax = settings_.rmaxdiskl1overlapvm();
      rmin = settings_.rmean(layerdisk1_);
    } else if (iSeed_ == Seed::L2D1) {
      rmax = settings_.rmaxdiskvm();
      rmin = settings_.rmean(layerdisk1_);
    } else {
      rmax = settings_.rmaxdiskvm();
      rmin = rmax * settings_.zmean(layerdisk2_ - N_LAYER - 1) / settings_.zmean(layerdisk2_ - N_LAYER);
    }
  }

  init(iSeed_);

  double dphimax = asin(0.5 * settings_.maxrinv() * rmax) - asin(0.5 * settings_.maxrinv() * rmin);

  //number of fine phi bins in sector
  int nfinephibins =
      settings_.nallstubs(layerdisk2_) * settings_.nvmte(1, iSeed_) * (1 << settings_.nfinephi(1, iSeed_));
  double dfinephi = settings_.dphisectorHG() / nfinephibins;

  nbitsfinephi_ = settings_.nbitsallstubs(layerdisk2_) + settings_.nbitsvmte(1, iSeed_) + settings_.nfinephi(1, iSeed_);

  int nbins = 2.0 * (dphimax / dfinephi + 1.0);

  nbitsfinephidiff_ = log(nbins) / log(2.0) + 1;

  nbitszfinebintable_ = settings_.vmrlutzbits(layerdisk1_);
  nbitsrfinebintable_ = settings_.vmrlutrbits(layerdisk1_);

  nbitsrzbin_ = N_RZBITS;
  if (iSeed_ == Seed::D1D2 || iSeed_ == Seed::D3D4)
    nbitsrzbin_ = 2;

  innerphibits_ = settings_.nfinephi(0, iSeed_);
  outerphibits_ = settings_.nfinephi(1, iSeed_);

  if (layerdisk1_ == LayerDisk::L1 || layerdisk1_ == LayerDisk::L2 || layerdisk1_ == LayerDisk::L3 ||
      layerdisk1_ == LayerDisk::L5 || layerdisk1_ == LayerDisk::D1 || layerdisk1_ == LayerDisk::D3) {
    innerTable_.initVMRTable(layerdisk1_, TrackletLUT::VMRTableType::inner);  //projection to next layer/disk
  }

  if (layerdisk1_ == LayerDisk::L1 || layerdisk1_ == LayerDisk::L2) {
    innerOverlapTable_.initVMRTable(layerdisk1_,
                                    TrackletLUT::VMRTableType::inneroverlap);  //projection to disk from layer
  }

  // set TC index
  iTC_ = name_[7] - 'A';
  assert(iTC_ >= 0 && iTC_ < 14);

  TCIndex_ = (iSeed_ << 4) + iTC_;
  assert(TCIndex_ >= 0 && TCIndex_ <= (int)settings_.ntrackletmax());

  maxStep_ = settings_.maxStep("TP");
}

void TrackletProcessor::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  if (output == "trackpar") {
    auto* tmp = dynamic_cast<TrackletParametersMemory*>(memory);
    assert(tmp != nullptr);
    trackletpars_ = tmp;
    return;
  }

  //if (output.substr(0,4) == "proj") {
  //Hack to keep proj output in config - but ignore in application
  //  return;
  //}

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find output : " << output;
}

void TrackletProcessor::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }

  if (input == "outervmstubin") {
    auto* tmp = dynamic_cast<VMStubsTEMemory*>(memory);
    assert(tmp != nullptr);
    outervmstubs_ = tmp;
    iAllStub_ = tmp->getName()[11] - 'A';
    if (iSeed_ == Seed::L2L3)
      iAllStub_ = tmp->getName()[11] - 'I';
    if (iSeed_ == Seed::L1D1 || iSeed_ == Seed::L2D1) {
      if (tmp->getName()[11] == 'X')
        iAllStub_ = 0;
      if (tmp->getName()[11] == 'Y')
        iAllStub_ = 1;
      if (tmp->getName()[11] == 'Z')
        iAllStub_ = 2;
      if (tmp->getName()[11] == 'W')
        iAllStub_ = 3;
    }

    unsigned int iTP = getName()[7] - 'A';

    pttableinner_.initTPlut(true, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephidiff_, iTP);
    pttableouter_.initTPlut(false, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephidiff_, iTP);

    //need iAllStub_ set before building the table

    useregiontable_.initTPregionlut(
        iSeed_, layerdisk1_, layerdisk2_, iAllStub_, nbitsfinephidiff_, nbitsfinephi_, pttableinner_, iTP);

    TrackletEngineUnit teunit(&settings_,
                              nbitsfinephi_,
                              layerdisk1_,
                              layerdisk2_,
                              iSeed_,
                              nbitsfinephidiff_,
                              iAllStub_,
                              &pttableinner_,
                              &pttableouter_,
                              outervmstubs_);

    teunits_.resize(settings_.teunits(iSeed_), teunit);

    return;
  }

  if (input == "innerallstubin") {
    auto* tmp = dynamic_cast<AllInnerStubsMemory*>(memory);
    assert(tmp != nullptr);
    if (innerallstubs_.size() == 2) {  //FIXME this should be done with better logic with reading the input stubs
      innerallstubs_.insert(innerallstubs_.begin(), tmp);
    } else {
      innerallstubs_.push_back(tmp);
    }

    //FIXME should be done once after all inputs are added
    tebuffer_ = tuple<CircularBuffer<TEData>, unsigned int, unsigned int, unsigned int, unsigned int>(
        CircularBuffer<TEData>(3), 0, 0, 0, innerallstubs_.size());

    return;
  }
  if (input == "outerallstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    outerallstubs_.push_back(tmp);
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find input : " << input;
}

void TrackletProcessor::execute(unsigned int iSector, double phimin, double phimax) {
  //bool print = (iSector == 3) && (getName() == "TP_D1D2C");
  bool print = false;

  phimin_ = phimin;
  phimax_ = phimax;
  iSector_ = iSector;

  if (!settings_.useSeed(iSeed_))
    return;

  //Not most elegant solution; but works - firmware will run for all clocks, but to save
  //CPU we stop processig when there is nothing more to do
  int donecount = 0;

  //Consistency checks
  assert(iAllStub_ >= 0);
  assert(iAllStub_ < (int)settings_.nallstubs(layerdisk2_));
  assert(outervmstubs_ != nullptr);

  //used to collect performance data
  unsigned int countsel = 0;

  unsigned int countteall = 0;
  unsigned int stubpairs = 0;

  unsigned int ninnerstubs = 0;

  //Actual implemenation starts here

  //Reset the tebuffer
  std::get<0>(tebuffer_).reset();
  std::get<1>(tebuffer_) = 0;
  std::get<2>(tebuffer_) = std::get<3>(tebuffer_);

  //Reset the teunits
  for (auto& teunit : teunits_) {
    teunit.reset();
  }

  TEData tedata;
  TEData tedata__;
  TEData tedata___;
  bool goodtedata = false;
  bool goodtedata__ = false;
  bool goodtedata___ = false;

  bool tebuffernearfull;

  for (unsigned int istep = 0; istep < maxStep_; istep++) {
    // These print statements are not on by default but can be enabled for the
    // comparison with HLS code to track differences.
    if (print) {
      CircularBuffer<TEData>& tedatabuffer = std::get<0>(tebuffer_);
      unsigned int& istub = std::get<1>(tebuffer_);
      unsigned int& imem = std::get<2>(tebuffer_);
      std::stringstream mess;
      mess << "istep=" << istep << " TEBuffer: " << istub << " " << imem << " " << tedatabuffer.rptr() << " "
           << tedatabuffer.wptr();
      int k = -1;
      for (auto& teunit : teunits_) {
        k++;
        mess << " [" << k << " " << teunit.rptr() << " " << teunit.wptr() << " " << teunit.idle() << "]";
      }
      edm::LogVerbatim("Tracklet") << mess.str();
    }

    CircularBuffer<TEData>& tedatabuffer = std::get<0>(tebuffer_);
    tebuffernearfull = tedatabuffer.nearfull();

    //
    // First block here checks if there is a teunit with data that should should be used
    // to calculate the tracklet parameters
    //

    TrackletEngineUnit* teunitptr = nullptr;

    for (auto& teunit : teunits_) {
      teunit.setNearFull();
      if (!teunit.empty()) {
        teunitptr = &teunit;
      }
    }

    if (teunitptr != nullptr) {
      auto stubpair = teunitptr->read();
      stubpairs++;

      if (trackletpars_->nTracklets() >= settings_.ntrackletmax()) {
        edm::LogVerbatim("Tracklet") << "Will break on too many tracklets in " << getName();
        break;
      }
      const Stub* innerFPGAStub = stubpair.first;
      const L1TStub* innerStub = innerFPGAStub->l1tstub();

      const Stub* outerFPGAStub = stubpair.second;
      const L1TStub* outerStub = outerFPGAStub->l1tstub();

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "TrackletProcessor execute " << getName() << "[" << iSector_ << "]";
      }

      bool accept = processStubPair(innerFPGAStub, innerStub, outerFPGAStub, outerStub, print);

      /*      
      if (iSeed_ == Seed::L1L2 || iSeed_ == Seed::L2L3 || iSeed_ == Seed::L3L4 || iSeed_ == Seed::L5L6) {
        accept = barrelSeeding(innerFPGAStub, innerStub, outerFPGAStub, outerStub, print);
      } else if (iSeed_ == Seed::D1D2 || iSeed_ == Seed::D3D4) {
        accept = diskSeeding(innerFPGAStub, innerStub, outerFPGAStub, outerStub, print);
      } else {
        accept = overlapSeeding(innerFPGAStub, innerStub, outerFPGAStub, outerStub);
      }
      */

      if (accept)
        countsel++;

      if (trackletpars_->nTracklets() >= settings_.ntrackletmax()) {
        edm::LogVerbatim("Tracklet") << "Will break on number of tracklets in " << getName();
        assert(0);
        break;
      }

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "TrackletProcessor execute done";
      }
    }

    //
    // The second block fills the teunit if data in buffer and process TEUnit step
    //
    //

    bool notemptytebuffer = !tedatabuffer.empty();

    int ite = -1;
    for (auto& teunit : teunits_) {
      ite++;
      if (teunit.idle()) {
        if (notemptytebuffer) {
          teunit.init(std::get<0>(tebuffer_).read());
          notemptytebuffer = false;  //prevent initialzing another TE unit
        }
      }
      teunit.step(print, istep, ite);
    }

    //
    // The third block here checks if we have input stubs to process
    //
    //

    if (goodtedata___)
      tedatabuffer.store(tedata___);

    goodtedata = false;

    unsigned int& istub = std::get<1>(tebuffer_);
    unsigned int& imem = std::get<2>(tebuffer_);
    unsigned int imemend = std::get<4>(tebuffer_);

    if ((!tebuffernearfull) && imem < imemend && istub < innerallstubs_[imem]->nStubs()) {
      ninnerstubs++;

      const Stub* stub = innerallstubs_[imem]->getStub(istub);

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << getName() << " Have stub in " << innerallstubs_[imem]->getName();
      }

      bool negdisk = (stub->disk().value() < 0);  //FIXME stub needs to contain bit for +/- z disk

      FPGAWord phicorr = stub->phicorr();
      int innerfinephi = phicorr.bits(phicorr.nbits() - nbitsfinephi_, nbitsfinephi_);
      FPGAWord innerbend = stub->bend();

      //Take the top nbitszfinebintable_ bits of the z coordinate
      int indexz = (stub->z().value() >> (stub->z().nbits() - nbitszfinebintable_)) & ((1 << nbitszfinebintable_) - 1);
      int indexr = -1;
      if (layerdisk1_ > (N_LAYER - 1)) {
        if (negdisk) {
          indexz = ((1 << nbitszfinebintable_) - 1) - indexz;
        }
        indexr =
            stub->rvalue() >>
            (stub->r().nbits() + 1 -
             nbitsrfinebintable_);  // + 1 required to offset artificial decrease in # of diskps r bits from 12 -> 11 to make space for negDisk bit
      } else {  //Take the top nbitsfinebintable_ bits of the z coordinate
        indexr = (stub->rvalue() >> (stub->r().nbits() - nbitsrfinebintable_)) & ((1 << nbitsrfinebintable_) - 1);
      }

      int lutval = -1;
      if (iSeed_ == Seed::L1D1 ||
          iSeed_ == Seed::L2D1) {  //FIXME should only be one table - but will need coordination with HLS code.
        lutval = innerOverlapTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
      } else {
        lutval = innerTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
      }

      if (lutval != -1) {
        unsigned int lutwidth = settings_.lutwidthtab(0, iSeed_);
        FPGAWord lookupbits(lutval, lutwidth, true, __LINE__, __FILE__);

        int rzfinebinfirst = lookupbits.bits(0, NFINERZBITS);       //finerz
        int next = lookupbits.bits(NFINERZBITS, 1);                 //use next r/z bin
        int start = lookupbits.bits(NFINERZBITS + 1, nbitsrzbin_);  //rz bin
        int rzdiffmax = lookupbits.bits(NFINERZBITS + 1 + nbitsrzbin_, NFINERZBITS);

        if ((iSeed_ == Seed::D1D2 || iSeed_ == Seed::D3D4) && negdisk) {  //TODO - need to store negative disk
          start += (1 << nbitsrzbin_);
        }
        int last = start + next;

        int nbins = (1 << N_RZBITS);

        unsigned int useregindex = (innerfinephi << innerbend.nbits()) + innerbend.value();
        if (iSeed_ == Seed::D1D2 || iSeed_ == Seed::D3D4 || iSeed_ == Seed::L1D1 || iSeed_ == Seed::L2D1) {
          //FIXME If the lookupbits were rationally organized this would be much simpler
          unsigned int nrbits = 3;
          int ir = ((start & ((1 << (nrbits - 1)) - 1)) << 1) + (rzfinebinfirst >> (NFINERZBITS - 1));
          useregindex = (useregindex << nrbits) + ir;
        }

        unsigned int usereg = useregiontable_.lookup(useregindex);

        tedata.regions_.clear();
        tedata.stub_ = stub;
        tedata.rzbinfirst_ = rzfinebinfirst;
        tedata.start_ = start;
        tedata.innerfinephi_ = innerfinephi;
        tedata.rzdiffmax_ = rzdiffmax;
        tedata.innerbend_ = innerbend;

        std::string mask = "";

        for (int ibin = start; ibin <= last; ibin++) {
          for (unsigned int ireg = 0; ireg < settings_.nvmte(1, iSeed_); ireg++) {
            if (!(usereg & (1 << ireg))) {
              mask = "0" + mask;
              continue;
            }

            if (settings_.debugTracklet()) {
              edm::LogVerbatim("Tracklet") << getName() << " looking for matching stub in bin " << ibin << " with "
                                           << outervmstubs_->nVMStubsBinned(ireg * nbins + ibin) << " stubs";
            }
            assert(ireg * nbins + ibin < outervmstubs_->nBin());
            int nstubs = outervmstubs_->nVMStubsBinned(ireg * nbins + ibin);

            if (nstubs > 0) {
              mask = "1" + mask;
              tedata.regions_.emplace_back(tuple<int, int, int>(ibin - start, ireg, nstubs));
              countteall += nstubs;
            } else {
              mask = "0" + mask;
            }
          }
        }

        if (!tedata.regions_.empty()) {
          goodtedata = true;
        }
      }
      istub++;
      if (istub >= innerallstubs_[imem]->nStubs()) {
        istub = 0;
        imem++;
      }
    } else if ((!tebuffernearfull) && imem < imemend && istub == 0) {
      imem++;
    }

    goodtedata___ = goodtedata__;
    goodtedata__ = goodtedata;

    tedata___ = tedata__;
    tedata__ = tedata;

    //
    // stop looping over istep if done
    //

    bool done = true;

    if (imem < imemend || (!tedatabuffer.empty())) {
      done = false;
    }

    for (auto& teunit : teunits_) {
      if (!(teunit.idle() && teunit.empty()))
        done = false;
    }

    if (done) {
      donecount++;
    }

    //FIXME This should be done cleaner... Not too hard, but need to check fully the TEBuffer and TEUnit buffer.
    if (donecount > 4) {
      break;
    }
  }

  //
  // Done with processing - collect performance statistics
  //

  if (settings_.writeMonitorData("TP")) {
    globals_->ofstream("trackletprocessor.txt") << getName() << " " << ninnerstubs   //# inner stubs
                                                << " " << outervmstubs_->nVMStubs()  //# outer stubs
                                                << " " << countteall                 //# pairs tried in TE
                                                << " " << stubpairs                  //# stubs pairs
                                                << " " << countsel                   //# tracklets found
                                                << endl;
  }
}


void TrackletProcessor::init(int iSeed) {
  phiHG_ = settings_.dphisectorHG();

  //Constants used for tracklet parameter calculations
  //These controls the number of bits that are kept in the intermediate
  //calculations. They are in some sense arbitrary, but should be
  //tuned to get an accurate result without using too much
  //resources
  std::map<int, int> na = {{Seed::L1L2, 15},
                           {Seed::L2L3, 15},
                           {Seed::L3L4, 15},
                           {Seed::L5L6, 15},
                           {Seed::D1D2, 15},
                           {Seed::D3D4, 15},
                           {Seed::L1D1, 14},
                           {Seed::L2D1, 14}};
  std::map<int, int> ndelta0 = {{Seed::L1L2, 12},
                                {Seed::L2L3, 12},
                                {Seed::L3L4, 12},
                                {Seed::L5L6, 12},
                                {Seed::D1D2, 12},
                                {Seed::D3D4, 12},
                                {Seed::L1D1, 12},
                                {Seed::L2D1, 12}};
  std::map<int, int> ndelta02 = {{Seed::L1L2, 15},
                                 {Seed::L2L3, 15},
                                 {Seed::L3L4, 15},
                                 {Seed::L5L6, 15},
                                 {Seed::D1D2, 15},
                                 {Seed::D3D4, 15},
                                 {Seed::L1D1, 15},
                                 {Seed::L2D1, 15}};
  std::map<int, int> ndelta1 = {{Seed::L1L2, 11},
                                {Seed::L2L3, 11},
                                {Seed::L3L4, 11},
                                {Seed::L5L6, 11},
                                {Seed::D1D2, 11},
                                {Seed::D3D4, 11},
                                {Seed::L1D1, 11},
                                {Seed::L2D1, 11}};
  std::map<int, int> ndelta2 = {{Seed::L1L2, 15},
                                {Seed::L2L3, 15},
                                {Seed::L3L4, 15},
                                {Seed::L5L6, 15},
                                {Seed::D1D2, 15},
                                {Seed::D3D4, 15},
                                {Seed::L1D1, 15},
                                {Seed::L2D1, 15}};
  std::map<int, int> ndelta12 = {{Seed::L1L2, 15},
                                 {Seed::L2L3, 15},
                                 {Seed::L3L4, 15},
                                 {Seed::L5L6, 15},
                                 {Seed::D1D2, 15},
                                 {Seed::D3D4, 15},
                                 {Seed::L1D1, 15},
                                 {Seed::L2D1, 15}};
  std::map<int, int> ndeltaz = {{Seed::L1L2, 9},
                                {Seed::L2L3, 9},
                                {Seed::L3L4, 9},
                                {Seed::L5L6, 9},
                                {Seed::D1D2, 9},
                                {Seed::D3D4, 9},
                                {Seed::L1D1, 7},
                                {Seed::L2D1, 7}};
  std::map<int, int> nr6 = {{Seed::L1L2, 0},
                            {Seed::L2L3, 0},
                            {Seed::L3L4, 0},
                            {Seed::L5L6, 0},
                            {Seed::D1D2, 0},
                            {Seed::D3D4, 0},
                            {Seed::L1D1, 0},
                            {Seed::L2D1, 0}};
  std::map<int, int> nifact = {{Seed::L1L2, 11},
                               {Seed::L2L3, 12},
                               {Seed::L3L4, 12},
                               {Seed::L5L6, 13},
                               {Seed::D1D2, 11},
                               {Seed::D3D4, 12},
                               {Seed::L1D1, 9},
                               {Seed::L2D1, 9}};
  std::map<int, int> nx6 = {{Seed::L1L2, 16},
                            {Seed::L2L3, 16},
                            {Seed::L3L4, 16},
                            {Seed::L5L6, 16},
                            {Seed::D1D2, 16},
                            {Seed::D3D4, 16},
                            {Seed::L1D1, 16},
                            {Seed::L2D1, 16}};
  std::map<int, int> nit1 = {{Seed::L1L2, 5},
                             {Seed::L2L3, 5},
                             {Seed::L3L4, 5},
                             {Seed::L5L6, 5},
                             {Seed::D1D2, 5},
                             {Seed::D3D4, 5},
                             {Seed::L1D1, 5},
                             {Seed::L2D1, 5}};
  std::map<int, int> nHG = {{Seed::L1L2, 10},
                            {Seed::L2L3, 11},
                            {Seed::L3L4, 11},
                            {Seed::L5L6, 11},
                            {Seed::D1D2, 11},
                            {Seed::D3D4, 11},
                            {Seed::L1D1, 10},
                            {Seed::L2D1, 11}};
  std::map<int, int> ndeltar = {{Seed::L1L2, 25},
                                {Seed::L2L3, 25},
                                {Seed::L3L4, 25},
                                {Seed::L5L6, 25},
                                {Seed::D1D2, 25},
                                {Seed::D3D4, 25},
                                {Seed::L1D1, 25},
                                {Seed::L2D1, 25}};

  //The offests that are commented out here are useful for adjusting the number of bits
  n_delta0_ = ndelta0[iSeed];    // + stoi(getenv("N_DELTA0"));
  n_deltaz_ = ndeltaz[iSeed];    // + stoi(getenv("N_DELTAZ"));
  n_delta1_ = ndelta1[iSeed];    // + stoi(getenv("N_DELTA1"));
  n_delta2_ = ndelta2[iSeed];    // + stoi(getenv("N_DELTA2"));
  n_delta12_ = ndelta12[iSeed];  // + stoi(getenv("N_DELTA12"));
  n_a_ = na[iSeed];              // + stoi(getenv("N_A"));
  n_r6_ = nr6[iSeed];            // + stoi(getenv("N_R6"));
  n_ifact_ = nifact[iSeed];      // + stoi(getenv("N_IFACT"));
  n_delta02_ = ndelta02[iSeed];  // + stoi(getenv("N_DELTA02"));
  n_x6_ = nx6[iSeed];            // + stoi(getenv("N_X6"));
  n_it1_ = nit1[iSeed];          // + stoi(getenv("N_IT1"));
  n_HG_ = nHG[iSeed];            // + stoi(getenv("N_HG"));
  n_Deltar_ = ndeltar[iSeed];    // + stoi(getenv("N_DELTAR"));

  //Constants for coordinates and track parameter definitions
  n_phi_ = settings_.nphibitsstub(N_LAYER - 1);
  n_r_ = settings_.nrbitsstub(N_LAYER) + 1;
  n_z_ = settings_.nzbitsstub(0) - 1;
  n_z0_ = settings_.nbitsz0() + 1;
  n_phi0_ = settings_.nbitsphi0() - 2;
  n_rinv_ = settings_.nbitsrinv() - 1;
  n_t_ = settings_.nbitst() - 5;

  if (settings_.barrelSeed(iSeed)) {
    LUT_idrinv_.resize(512);
    for (int idr = -256; idr < 256; idr++) {
      int uidr = idr;
      if (uidr < 0)
        uidr += 512;
      int idrabs = idr + settings_.irmean(layerdisk2_) - settings_.irmean(layerdisk1_);
      LUT_idrinv_[uidr] = 0.5 + float(1 << n_Deltar_) / idrabs;
    }
  }

  if (settings_.diskSeed(iSeed)) {
    LUT_idrinv_.resize(512);
    for (unsigned int idr = 1; idr < 512; idr++) {
      LUT_idrinv_[idr] = 0.5 + float(1 << n_Deltar_) / idr;
    }
  }

  if (settings_.overlapSeed(iSeed)) {
    LUT_idrinv_.resize(1024);
    for (unsigned int idr = 1; idr < 1024; idr++) {
      LUT_idrinv_[idr] = 0.5 + float(1 << n_Deltar_) / idr;
    }
  }
}

void TrackletProcessor::exacttracklet(double r1,
				      double z1,
				      double phi1,
				      double r2,
				      double z2,
				      double phi2,
				      double,
				      double& rinv,
				      double& phi0,
				      double& t,
				      double& z0) {
  double deltaphi = reco::reducePhiRange(phi1 - phi2);

  double dist = sqrt(r2 * r2 + r1 * r1 - 2 * r1 * r2 * cos(deltaphi));

  rinv = 2 * sin(deltaphi) / dist;

  if (r1 > r2) {
    //std::cout << "Flipping rinv" << std::endl;
    rinv = -rinv;
  }

  double phi1tmp = phi1 - phimin_;

  phi0 = reco::reducePhiRange(phi1tmp + asin(0.5 * r1 * rinv));

  double rhopsi1 = 2 * asin(0.5 * r1 * rinv) / rinv;
  double rhopsi2 = 2 * asin(0.5 * r2 * rinv) / rinv;

  t = (z1 - z2) / (rhopsi1 - rhopsi2);

  z0 = z1 - t * rhopsi1;
}

void TrackletProcessor::calcPars(unsigned int idr,
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
				 bool print) {
  long int idz = iz2 - iz1;

  assert(idr < LUT_idrinv_.size());
  long int invdr = LUT_idrinv_[idr];

  long int idelta0 = ((iphi2 - iphi1) * invdr) >> n_delta0_;
  long int ideltaz = (idz * invdr) >> n_deltaz_;

  long int idelta1 = (ir1 * idelta0) >> n_delta1_;
  long int idelta2 = (ir2 * idelta0) >> n_delta2_;

  long int idelta12 = (idelta1 * idelta2) >> n_delta12_;

  long int iHG = phiHG_ * phiHG_ * (1 << n_HG_);

  long int ia = ((1 << n_a_) - ((idelta12 * iHG) >> (2 * n_Deltar_ + 2 * n_phi_ + n_HG_ - 2 * n_delta0_ - n_delta1_ -
                                                     n_delta2_ - n_delta12_ + 1 - n_a_)));

  long int ifact = (1 << n_ifact_) * phiHG_ * phiHG_ / 6.0;

  long ir6 = ((ir1 + ir2) * ifact) >> (n_ifact_ - n_r6_);

  long int idelta02 = (idelta0 * idelta2) >> n_delta02_;

  long int ix6 =
      (-(1 << n_x6_) +
       ((ir6 * idelta02) >> (n_r6_ + 2 * n_Deltar_ + 2 * n_phi_ - n_x6_ - n_delta2_ - n_delta02_ - 2 * n_delta0_)));

  //Temporary hack here
  long int it1 = (ir1 * ideltaz) >> (n_Deltar_ - n_deltaz_ - n_it1_);

  irinv_new = (((-idelta0 * ia) >> (n_phi_ + n_a_ - n_rinv_ + n_Deltar_ - n_delta0_ - n_r_ - 1 - 1)) + 1) >> 1;

  iphi0_new =
      (((iphi1 + ((idelta1 * ix6) >> (n_Deltar_ + n_x6_ - n_delta0_ - n_delta1_))) >> (n_phi_ - n_phi0_ - 1)) + 1) >> 1;

  long int shift_tmp = n_Deltar_ + n_a_ + n_z_ - n_t_ - n_deltaz_ - n_r_;
  it_new = (((ideltaz * ia) >> (shift_tmp - 1)) + 1) >> 1;

  iz0_new = (iz1 << 1) + ((((it1 * ix6) >> (n_x6_ + n_it1_ - 2)) + 1) >> 1);

  if (print) {
    std::cout << "=================================" << std::endl;
    std::cout << "ir1 iz1: " << ir1 << " " << iz1 << std::endl;
    std::cout << "ir2 iz2: " << ir2 << " " << iz2 << std::endl;
    std::cout << "ir6 idelta02 : " << ir6 << " " << idelta02 << std::endl;
    std::cout << "it1 ix6 : " << it1 << " " << ix6 << std::endl;
    std::cout << "it iz0 : " << it_new << " " << iz0_new << std::endl;
  }
}

bool TrackletProcessor::goodTrackPars(bool goodrinv, bool goodz0) {
  bool success = true;
  if (!goodrinv) {
    if (settings_.debugTracklet()) {
      edm::LogVerbatim("Tracklet") << getName() << " TrackletCalculatorBase irinv too large";
    }
    success = false;
  }
  if (!goodz0) {
    if (settings_.debugTracklet()) {
      edm::LogVerbatim("Tracklet") << getName() << " TrackletCalculatorBase z0 cut to large";
    }
    success = false;
  }
  return success;
}

bool TrackletProcessor::inSector(int iphi0, int irinv, double phi0approx, double rinvapprox) {
  double phicritapprox = phi0approx - asin(0.5 * settings_.rcrit() * rinvapprox);

  int ifactor = 0.5 * settings_.rcrit() * settings_.krinvpars() / settings_.kphi0pars() * (1 << 8);
  int iphicrit = iphi0 - (irinv >> 8) * ifactor;

  int iphicritmincut = settings_.phicritminmc() / settings_.kphi0pars();
  int iphicritmaxcut = settings_.phicritmaxmc() / settings_.kphi0pars();

  bool keepapprox = (phi0approx >= 0) && (phicritapprox > settings_.phicritminmc()) &&
                    (phicritapprox < settings_.phicritmaxmc()),
       keep = (iphi0 >= 0) && (iphicrit > iphicritmincut) && (iphicrit < iphicritmaxcut);
  if (settings_.debugTracklet())
    if (keepapprox && !keep)
      edm::LogVerbatim("Tracklet") << getName()
                                   << " Tracklet kept with exact phicrit cut but not approximate, phicritapprox: "
                                   << phicritapprox;
  if (settings_.usephicritapprox()) {
    return keepapprox;
  } else {
    return keep;
  }

  return true;
}


bool TrackletProcessor::processStubPair(const Stub* innerFPGAStub,
					const L1TStub* innerStub,
					const Stub* outerFPGAStub,
					const L1TStub* outerStub,
					bool print) {
  if (settings_.debugTracklet()) {
    edm::LogVerbatim("Tracklet") << "TrackletCalculatorBase " << getName()
                                 << " trying stub pair in layer (inner outer): " << innerFPGAStub->layer().value()
                                 << " " << outerFPGAStub->layer().value();
  }

  double r1 = innerStub->r();
  double z1 = innerStub->z();
  double phi1 = innerStub->phi();

  double r2 = outerStub->r();
  double z2 = outerStub->z();
  double phi2 = outerStub->phi();

  double rinv, phi0, t, z0;

  exacttracklet(r1, z1, phi1, r2, z2, phi2, outerStub->sigmaz(), rinv, phi0, t, z0);

  //now binary

  int ir1 = innerFPGAStub->rvalue();
  int iphi1 = innerFPGAStub->phi().value();
  int iz1 = innerFPGAStub->z().value();

  int ir2 = outerFPGAStub->rvalue();
  int iphi2 = outerFPGAStub->phi().value();
  int iz2 = outerFPGAStub->z().value();

  int idr, ir1abs, ir2abs, iz1abs, iz2abs;
  
  //barrel

  if (settings_.barrelSeed(iSeed_)) {
  
    iphi1 <<= (settings_.nphibitsstub(5) - settings_.nphibitsstub(layerdisk1_));
    iphi2 <<= (settings_.nphibitsstub(5) - settings_.nphibitsstub(layerdisk2_));
    ir1 <<= (8 - settings_.nrbitsstub(layerdisk1_));
    ir2 <<= (8 - settings_.nrbitsstub(layerdisk2_));

    iz1abs = iz1 << (settings_.nzbitsstub(0) - settings_.nzbitsstub(layerdisk1_));
    iz2abs = iz2 << (settings_.nzbitsstub(0) - settings_.nzbitsstub(layerdisk2_));

    //Each of ir1 and ir2 are signed 8 bit integers. idr is signed 9 bit integer
    idr = ir2 - ir1;

    if (idr < 0)
      idr += 512;


    unsigned int ir1mean = settings_.irmean(layerdisk1_);
    unsigned int ir2mean = settings_.irmean(layerdisk2_);

    ir1abs = ir1 + ir1mean;
    ir2abs = ir2 + ir2mean;
  } else  if (settings_.diskSeed(iSeed_)) {
    
    //disk

    //To get same precission as for layers.
    iphi1 <<= (settings_.nphibitsstub(5) - settings_.nphibitsstub(0));
    iphi2 <<= (settings_.nphibitsstub(5) - settings_.nphibitsstub(0));
    
    //Each of ir1 and ir2 are signed 8 bit integers. idr is signed 9 bit integer
    idr = ir2 - ir1;

    int sign = (innerFPGAStub->disk().value() < 0) ? -1 : 1;
    
    unsigned int iz1mean = sign * settings_.izmean(layerdisk1_ - N_LAYER);
    unsigned int iz2mean = sign * settings_.izmean(layerdisk2_ - N_LAYER);
    
    iz1abs = iz1 + iz1mean;
    iz2abs = iz2 + iz2mean;
    
    ir1abs = ir1;
    ir2abs = ir2;
  } else {
    
    //overlap

    //Protection for wrong radii. Could be handled cleaner to avoid problem with floating point calculation and with overflows in the integer calculation.
    if (r2 < r1 + 1.5) {
      return false;
    }

    
    //To get global precission
    int ll = innerFPGAStub->layer().value() + 1;
    ir1 = l1t::bitShift(ir1, (8 - settings_.nrbitsstub(ll - 1)));
    iphi1 <<= (settings_.nphibitsstub(5) - settings_.nphibitsstub(0));
    iphi2 <<= (settings_.nphibitsstub(5) - settings_.nphibitsstub(0));

    unsigned int ir1mean = settings_.irmean(layerdisk1_);

    ir1abs = ir1 + ir1mean;
    ir2abs = ir2;
    
    idr = ir2 - ir1abs;

    if (idr >= (int)LUT_idrinv_.size()) {
      return false;
    }

    int iz2mean = settings_.izmean(layerdisk2_ - N_LAYER);

    if (iz1 < 0) {
      iz2mean = -iz2mean;
    }

    iz1abs = iz1;
    iz2abs = iz2 + iz2mean;
  }
  
  int irinv_new, iphi0_new, iz0_new, it_new;

  calcPars(idr, iphi1, ir1abs, iz1abs, iphi2, ir2abs, iz2abs, irinv_new, iphi0_new, iz0_new, it_new, print);

  bool rinvcut = abs(irinv_new) < settings_.rinvcut() * (settings_.zlength() * (1 << n_rinv_)) / phiHG_;
  bool z0cut = abs(iz0_new) < settings_.z0cut() * (1 << n_z0_) / settings_.zlength();
  if (settings_.barrelSeed(iSeed_) && (iSeed_ != 0)) {
    z0cut = abs(iz0_new) < 1.5 * settings_.z0cut() * (1 << n_z0_) / settings_.zlength();
  }

  if (!goodTrackPars(rinvcut, z0cut)) {
    if (settings_.debugTracklet()) {
      edm::LogVerbatim("Tracklet") << getName() << " Failed rinv or z0 cut";
    }
    return false;
  }

  if (!inSector(iphi0_new, irinv_new, phi0, rinv)) {
    if (settings_.debugTracklet()) {
      edm::LogVerbatim("Tracklet") << getName() << " Failed in sector check";
    }
    return false;
  }

  if (settings_.writeMonitorData("TPars")) {
    globals_->ofstream("trackletpars.txt")
        << "Trackpars " << layerdisk1_ + 1 << "   "
        << " " << rinv << " " << irinv_new << "   " << phi0 << " " << iphi0_new << "   " << t << " " << it_new << "   "
        << " " << iz0_new << endl;
  }

  Tracklet* tracklet = new Tracklet(settings_,
                                    globals_,
                                    iSeed_,
                                    innerFPGAStub,
                                    nullptr,
                                    outerFPGAStub,
                                    rinv,
                                    phi0,
                                    0.0,
                                    z0,
                                    t,
                                    rinv,
                                    phi0,
                                    0.0,
                                    z0,
                                    t,
                                    irinv_new,
                                    iphi0_new,
                                    0,
                                    iz0_new,
                                    it_new,
                                    settings_.diskSeed(iSeed_),
				    settings_.overlapSeed(iSeed_));

  if (settings_.debugTracklet()) {
    edm::LogVerbatim("Tracklet") << "TrackletCalculator " << getName() << " Found tracklet for seed = " << iSeed_ << " "
                                 << iSector_ << " phi0 = " << phi0;
  }

  tracklet->setTrackletIndex(trackletpars_->nTracklets());
  tracklet->setTCIndex(TCIndex_);

  if (settings_.writeMonitorData("Seeds")) {
    ofstream fout("seeds.txt", ofstream::app);
    fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector_ << " " << tracklet->getISeed() << endl;
    fout.close();
  }
  trackletpars_->addTracklet(tracklet);

  if (settings_.bookHistos()) {
    HistBase* hists = globals_->histograms();
    int tp = tracklet->tpseed();
    hists->fillTrackletParams(settings_,
                              globals_,
                              iSeed_,
                              iSector_,
                              rinv,
                              irinv_new * settings_.krinvpars(),
                              phi0,
                              iphi0_new * settings_.kphi0pars(),
                              asinh(t),
                              asinh(it_new * settings_.ktpars()),
                              z0,
                              iz0_new * settings_.kz0pars(),
                              tp);
  }

  return true;
}

