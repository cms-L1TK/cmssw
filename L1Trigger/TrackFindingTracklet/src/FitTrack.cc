#include "L1Trigger/TrackFindingTracklet/interface/FitTrack.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/HybridFit.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubStreamData.h"
#include "DataFormats/L1TrackTrigger/interface/TTBV.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <sstream>

using namespace std;
using namespace trklet;

FitTrack::FitTrack(string name, Settings const& settings, Globals* global)
    : ProcessBase(name, settings, global), trackfit_(nullptr) {
  fullmatch_.resize(N_LAYER + N_DISK);
}

void FitTrack::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  if (output == "trackout") {
    TrackFitMemory* tmp = dynamic_cast<TrackFitMemory*>(memory);
    assert(tmp != nullptr);
    trackfit_ = tmp;
    return;
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " addOutput, output = " << output << " not known";
}

void FitTrack::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }
  if (input.substr(0, 4) == "tpar") {
    auto* tmp = dynamic_cast<TrackletParametersMemory*>(memory);
    assert(tmp != nullptr);
    seedtracklet_.push_back(tmp);
    return;
  }

  for (unsigned int i = 0; i < N_LAYER + N_DISK; i++) {
    std::ostringstream oss;
    oss << "fullmatch" << i << "in";
    ;
    auto const& str = oss.str();
    if (input.substr(0, 11) == str.substr(0, 11)) {
      auto* tmp = dynamic_cast<FullMatchMemory*>(memory);
      assert(tmp != nullptr);
      fullmatch_[i].push_back(tmp);
      return;
    }
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " input = " << input << " not found";
}

#ifdef USEHYBRID
void FitTrack::trackFitKF(Tracklet* tracklet,
                          std::vector<const Stub*>& trackstublist,
                          std::vector<std::pair<int, int>>& stubidslist) {
  if (settings_.doKF()) {
    // From full match lists, collect all the stubs associated with the tracklet seed

    // Get seed stubs first
    trackstublist.emplace_back(tracklet->innerFPGAStub());
    if (tracklet->getISeed() >= (int)N_TRKLSEED + 1)
      trackstublist.emplace_back(tracklet->middleFPGAStub());
    trackstublist.emplace_back(tracklet->outerFPGAStub());

    // Now get ALL matches (can have multiple per layer)
    for (unsigned int k = 0; k < fullmatch_.size(); k++) {
      for (const auto& i : fullmatch_[k]) {
        for (unsigned int j = 0; j < i->nMatches(); j++) {
          if (i->getTracklet(j)->TCID() == tracklet->TCID()) {
            trackstublist.push_back(i->getMatch(j).second);
          }
        }
      }
    }

    // For merge removal, loop through the resulting list of stubs to calculate their stubids
    if (settings_.removalType() == "merge") {
      for (const auto& it : trackstublist) {
        int layer = it->layer().value() + 1;  // Assume layer (1-6) stub first
        if (it->layer().value() < 0) {        // if disk stub, though...
          layer = it->disk().value() + 10 * it->disk().value() / abs(it->disk().value());  //disk = +/- 11-15
        }
        stubidslist.push_back(std::make_pair(layer, it->phiregionaddress()));
      }

      // And that's all we need! The rest is just for track fit (done in PurgeDuplicate)

    } else {
      // Track fit only called here if not running duplicate removal
      // before fit. (e.g. If skipping duplicate removal).
      HybridFit hybridFitter(iSector_, settings_, globals_);
      hybridFitter.Fit(tracklet, trackstublist);
    }
  }
}

#else  // Code for pure Tracklet algo -- now deleted...


#endif

void FitTrack::trackFitFake(Tracklet* tracklet, std::vector<const Stub*>&, std::vector<std::pair<int, int>>&) {
  tracklet->setFitPars(tracklet->rinvapprox(),
                       tracklet->phi0approx(),
                       tracklet->d0approx(),
                       tracklet->tapprox(),
                       tracklet->z0approx(),
                       0.0,
                       0.0,
                       tracklet->rinv(),
                       tracklet->phi0(),
                       tracklet->d0(),
                       tracklet->t(),
                       tracklet->z0(),
                       0.0,
                       0.0,
                       tracklet->fpgarinv().value(),
                       tracklet->fpgaphi0().value(),
                       tracklet->fpgad0().value(),
                       tracklet->fpgat().value(),
                       tracklet->fpgaz0().value(),
                       0,
                       0,
                       0);
  return;
}

std::vector<Tracklet*> FitTrack::orderedMatches(vector<FullMatchMemory*>& fullmatch) {
  std::vector<Tracklet*> tmp;

  std::vector<unsigned int> indexArray;
  for (auto& imatch : fullmatch) {
    //check that we have correct order
    if (imatch->nMatches() > 1) {
      for (unsigned int j = 0; j < imatch->nMatches() - 1; j++) {
        assert(imatch->getTracklet(j)->TCID() <= imatch->getTracklet(j + 1)->TCID());
      }
    }

    if (settings_.debugTracklet() && imatch->nMatches() != 0) {
      edm::LogVerbatim("Tracklet") << "orderedMatches: " << imatch->getName() << " " << imatch->nMatches();
    }

    indexArray.push_back(0);
  }

  int bestIndex = -1;
  do {
    int bestTCID = -1;
    bestIndex = -1;
    for (unsigned int i = 0; i < fullmatch.size(); i++) {
      if (indexArray[i] >= fullmatch[i]->nMatches()) {
        //skip as we were at the end
        continue;
      }
      int TCID = fullmatch[i]->getTracklet(indexArray[i])->TCID();
      if (TCID < bestTCID || bestTCID < 0) {
        bestTCID = TCID;
        bestIndex = i;
      }
    }
    if (bestIndex != -1) {
      tmp.push_back(fullmatch[bestIndex]->getTracklet(indexArray[bestIndex]));
      indexArray[bestIndex]++;
    }
  } while (bestIndex != -1);

  for (unsigned int i = 0; i < tmp.size(); i++) {
    if (i > 0) {
      //This allows for equal TCIDs. This means that we can e.g. have a track seeded in L1L2 that projects to both L3 and D4.
      //The algorithm will pick up the first hit and drop the second.
      if (tmp[i - 1]->TCID() > tmp[i]->TCID()) {
        edm::LogVerbatim("Tracklet") << "Wrong TCID ordering in " << getName() << " : " << tmp[i - 1]->TCID() << " "
                                     << tmp[i]->TCID();
      }
    }
  }

  return tmp;
}

// Adds the fitted track to the output memories to be used by pure Tracklet algo.
// (Also used by Hybrid algo with non-exact Old KF emulation)
// Also create output streams, that bypass these memories, (so can include gaps in time),
// to be used by Hybrid case with exact New KF emulation.

void FitTrack::execute(deque<string>& streamTrackRaw,
                       vector<deque<StubStreamData>>& streamsStubRaw,
                       unsigned int iSector) {
  // merge
  std::vector<std::vector<Tracklet*>> matches;

  matches.reserve(fullmatch_.size());
  for (unsigned int i = 0; i < fullmatch_.size(); i++) {
    matches.push_back(orderedMatches(fullmatch_[i]));
  }

  bool print = getName() == "FT_D1D2" && iSector == 3;
  print = false;

  iSector_ = iSector;

  unsigned int indexArray[N_LAYER + N_DISK];
  for (unsigned int i = 0; i < N_LAYER + N_DISK; i++) {
    indexArray[i] = 0;
  }

  unsigned int count = 0;
  unsigned int countFit = 0;
  unsigned int countAll = 0;

  int istep = -1;

  Tracklet* bestTracklet = nullptr;
  do {
    istep++;
    count++;
    bestTracklet = nullptr;

    for (unsigned int i = 0; i < matches.size(); i++) {
      if (indexArray[i] < matches[i].size()) {
        if (bestTracklet == nullptr) {
          bestTracklet = matches[i][indexArray[i]];
        } else {
          if (matches[i][indexArray[i]]->TCID() < bestTracklet->TCID()) {
            bestTracklet = matches[i][indexArray[i]];
          }
        }
      }
    }

    if (bestTracklet == nullptr)
      break;

    countAll++;

    //Counts total number of matched hits
    int nMatches = 0;

    //Counts unique hits in each layer
    int nMatchesUniq = 0;

    std::stringstream mess;
    if (print)
      mess << "istep = " << istep;

    for (unsigned int i = 0; i < matches.size(); i++) {
      bool match = false;
      while (indexArray[i] < matches[i].size() && matches[i][indexArray[i]] == bestTracklet) {
        if (print)
          mess << " match" << i;
        indexArray[i]++;
        nMatches++;
        match = true;
      }

      if (match)
        nMatchesUniq++;
    }

    if (print) {
      mess << " nMatchesUniq = " << nMatchesUniq;
      edm::LogVerbatim("Tracklet") << mess.str();
    }

    if (settings_.debugTracklet()) {
      edm::LogVerbatim("Tracklet") << getName() << " : nMatches = " << nMatches << " nMatchesUniq = " << nMatchesUniq;
    }

    std::vector<const Stub*> trackstublist;
    std::vector<std::pair<int, int>> stubidslist;
    // Track Builder cut of >= 4 layers with stubs.
    if ((bestTracklet->getISeed() >= (int)N_SEED_PROMPT && nMatchesUniq >= 1) ||
        nMatchesUniq >= 2) {  //For seeds index >=8 (triplet seeds), there are three stubs associated from start.
      countFit++;

#ifdef USEHYBRID
      if (settings_.fakefit()) {
        trackFitKF(bestTracklet, trackstublist, stubidslist);
      } else {
        trackFitKF(bestTracklet, trackstublist, stubidslist);
      }
#else
      if (settings_.fakefit()) {
        trackFitFake(bestTracklet, trackstublist, stubidslist);
      } else {
        trackFitChisq(bestTracklet, trackstublist, stubidslist);
      }
#endif

      if (settings_.removalType() == "merge") {
        trackfit_->addStubList(trackstublist);
        trackfit_->addStubidsList(stubidslist);
        bestTracklet->setTrackIndex(trackfit_->nTracks());
        trackfit_->addTrack(bestTracklet);
      } else if (bestTracklet->fit()) {
        assert(trackfit_ != nullptr);
        if (settings_.writeMonitorData("Seeds")) {
          ofstream fout("seeds.txt", ofstream::app);
          fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_"
               << " " << bestTracklet->getISeed() << endl;
          fout.close();
        }
        bestTracklet->setTrackIndex(trackfit_->nTracks());
        trackfit_->addTrack(bestTracklet);
      }
    }

    // store bit and clock accurate TB output
    if (settings_.storeTrackBuilderOutput() && bestTracklet) {
      // add gap if TrackBuilder rejected track (due to too few stub layers).
      if (!bestTracklet->fit()) {
        static const string invalid = "0";
        streamTrackRaw.emplace_back(invalid);
        for (auto& stream : streamsStubRaw)
          stream.emplace_back();
        continue;
      }
      // convert Track word
      const string rinv = bestTracklet->fpgarinv().str();
      const string phi0 = bestTracklet->fpgaphi0().str();
      const string z0 = bestTracklet->fpgaz0().str();
      const string t = bestTracklet->fpgat().str();
      const int seedType = bestTracklet->getISeed();
      const string seed = TTBV(seedType, settings_.nbitsseed()).str();
      const string valid("1");
      streamTrackRaw.emplace_back(valid + seed + rinv + phi0 + z0 + t);

      // convert projected stubs
      TTBV hitPattern(0, N_LAYER + N_DISK);
      for (unsigned int ilayer = 0; ilayer < N_LAYER + N_DISK; ilayer++) {
        if (!bestTracklet->match(ilayer))
          continue;
        hitPattern.set(ilayer);
        const Residual& resid = bestTracklet->resid(ilayer);
        // create bit accurate 64 bit word
        // Need to extract the corrected r value
        string r = resid.stubptr()->r().str();
        const string& phi = resid.fpgaphiresid().str();
        const string& rz = resid.fpgarzresid().str();
        const L1TStub* stub = resid.stubptr()->l1tstub();
        static constexpr int widthDisk2Sidentifier = 8;
        bool disk2S = (stub->disk() != 0) && (stub->isPSmodule() == 0);
        if (disk2S) {
          r = string(widthDisk2Sidentifier, '0') + r;
          // edm::LogVerbatim("Tracklet") << "2s stub : " << r;
        }
        bool diskPS = (stub->disk() != 0) && (stub->isPSmodule() != 0);
        if (diskPS) {
          //  string oldr = r;
          FPGAWord tmp(resid.stubptr()->rvalue(), 12);
          r = tmp.str();
          // edm::LogVerbatim("Tracklet") << "old r: "<< oldr << "new r: " << r;
        }
        const string& stubId = resid.fpgastubid().str();
        // store seed, L1TStub, and bit accurate 64 bit word in clock accurate output
        streamsStubRaw[ilayer].emplace_back(seedType, *stub, valid + stubId + r + phi + rz);
      }
      // convert seed stubs
      const string& stubId0 = bestTracklet->innerFPGAStub()->phiregionaddressstr();
      const L1TStub* stub0 = bestTracklet->innerFPGAStub()->l1tstub();
      const int layer0 = bestTracklet->innerFPGAStub()->layerdisk();
      hitPattern.set(layer0);
      streamsStubRaw[layer0].emplace_back(seedType, *stub0, valid + stubId0);
      const string& stubId1 = bestTracklet->outerFPGAStub()->phiregionaddressstr();
      const L1TStub* stub1 = bestTracklet->outerFPGAStub()->l1tstub();
      const int layer1 = bestTracklet->outerFPGAStub()->layerdisk();
      hitPattern.set(layer1);
      streamsStubRaw[layer1].emplace_back(seedType, *stub1, valid + stubId1);
      // fill all layers that have no stubs with gaps
      for (int ilayer : hitPattern.ids(false))
        streamsStubRaw[ilayer].emplace_back();
    }

  } while (bestTracklet != nullptr && count < settings_.maxStep("TB"));

  if (settings_.writeMonitorData("FT")) {
    globals_->ofstream("fittrack.txt") << getName() << " " << countAll << " " << countFit << endl;
  }
}
