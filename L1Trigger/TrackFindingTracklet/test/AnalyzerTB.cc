#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/TrackFindingTracklet/interface/DataFormats.h"
#include "SimDataFormats/Associations/interface/TTTypes.h"

#include <TProfile.h>
#include <TH1F.h>

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <cmath>
#include <numeric>
#include <sstream>

namespace trklet {

  /*! \class  trklet::AnalyzerTB
   *  \brief  Class to analyze emulated All Stubs
   *  \author Thomas Schuh
   *  \date   2025, Nov
   */
  class AnalyzerTB : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
  public:
    AnalyzerTB(const edm::ParameterSet& iConfig);
    void beginJob() override {}
    void beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override;
    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    void endRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override {}
    void endJob() override {}

  private:
    // ED input token of all Stubs
    edm::EDGetTokenT<tt::StreamsStub> edGetToken_;
    // Setup token
    edm::ESGetToken<Setup, trackerDTC::SetupRcd> esGetToken_;
    // helper class to store configurations
    const Setup* setup_;
    //
    std::vector<TH1F*> hisIn_;
    std::vector<TH1F*> hisMid_;
    std::vector<TH1F*> hisOut_;
    TProfile* profIn_;
    TProfile* profMid_;
    TProfile* profOut_;
  };

  AnalyzerTB::AnalyzerTB(const edm::ParameterSet& iConfig) {
    usesResource("TFileService");
    const std::string& label = iConfig.getParameter<std::string>("InputLabelTM");
    const std::string& branch = iConfig.getParameter<std::string>("BranchAllStubs");
    // book  ED products
    edGetToken_ = consumes(edm::InputTag(label, branch));
    // book ES products
    esGetToken_ = esConsumes<edm::Transition::BeginRun>();
  }

  void AnalyzerTB::beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetToken_);
    // book histograms
    edm::Service<TFileService> fs;
    TFileDirectory dir;
    dir = fs->mkdir("in");
    constexpr int numChanIn = 48;
    constexpr int maxOccIn = 110;
    hisIn_ = std::vector<TH1F*>(numChanIn);
    for (int i = 0; i < numChanIn; i++)
      hisIn_[i] = dir.make<TH1F>(
          ("His Channel " + std::to_string(i) + " Occupancy").c_str(), ";", maxOccIn, -.5, maxOccIn - .5);
    profIn_ = dir.make<TProfile>("Prof Channel Occupancy", ";", numChanIn, -.5, numChanIn - .5);
    dir = fs->mkdir("mid");
    constexpr int numChanMid = 12;
    constexpr int maxOccMid = 490;
    hisMid_ = std::vector<TH1F*>(numChanMid);
    for (int i = 0; i < numChanMid; i++)
      hisMid_[i] = dir.make<TH1F>(
          ("His Channel " + std::to_string(i) + " Occupancy").c_str(), ";", maxOccMid, -.5, maxOccMid - .5);
    profMid_ = dir.make<TProfile>("Prof Channel Occupancy", ";", numChanMid, -.5, numChanMid - .5);
    dir = fs->mkdir("out");
    const int numChanOut = setup_->asMux().size();
    constexpr int maxOccOut = 170;
    hisOut_ = std::vector<TH1F*>(numChanOut);
    for (int i = 0; i < numChanOut; i++)
      hisOut_[i] = dir.make<TH1F>(
          ("His Channel " + std::to_string(i) + " Occupancy").c_str(), ";", maxOccOut, -.5, maxOccOut - .5);
    profOut_ = dir.make<TProfile>("Prof Channel Occupancy", ";", numChanOut, -.5, numChanOut - .5);
  }

  void AnalyzerTB::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // read in All Stubs
    const tt::StreamsStub& streams = iEvent.get(edGetToken_);
    // count stubs
    for (int iRegion = 0; iRegion < setup_->sysNumRegion(); iRegion++) {
      const int offset = iRegion * 48;
      std::vector<int> nums(48, 0);
      for (int iChan = 0; iChan < 48; iChan++) {
        int& num = nums[iChan];
        for (const tt::FrameStub& frame : streams[offset + iChan])
          if (frame.first.isNonnull())
            num++;
        hisIn_[iChan]->Fill(num);
        profIn_->Fill(iChan, num);
      }
      // mid
      for (int i = 0; i < 12; i++) {
        int num(0);
        for (int j = 0; j < 4; j++)
          num += nums[4 * i + j];
        hisMid_[i]->Fill(num);
        profMid_->Fill(i, num);
      }
      // out
      int i(0);
      for (const std::vector<int>& mux : setup_->asMux()) {
        int num(0);
        for (int j : mux)
          num += nums[j];
        if (num > 161)
          num = 161;
        hisOut_[i]->Fill(num);
        profOut_->Fill(i++, num);
      }
    }
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::AnalyzerTB);
