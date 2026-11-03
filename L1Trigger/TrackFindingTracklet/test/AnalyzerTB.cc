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
    std::vector<TH1F*> his_;
    TProfile* prof_;
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
    // stub z postion from tracklet residuals plus tracklet seed parameter vs stub z psotion from TTStubs projected to stub radius from tracklet using tracklet seed parameter
    dir = fs->mkdir("");
    constexpr int numChan = 48;
    constexpr int maxOcc = 180;
    his_ = std::vector<TH1F*>(numChan);
    for (int i = 0; i < numChan; i++)
      his_[i] = dir.make<TH1F>(("His Channel Occupancy " + std::to_string(i)).c_str(), ";", maxOcc, -.5, maxOcc - .5);
    prof_ = dir.make<TProfile>("Prof Channel Occupancy", ";", numChan, -.5, numChan - .5);
  }

  void AnalyzerTB::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // read in All Stubs
    const tt::StreamsStub& streams = iEvent.get(edGetToken_);
    // count stubs
    for (int iRegion = 0; iRegion < setup_->sysNumRegion(); iRegion++) {
      const int offset = iRegion * 48;
      for (int iChan = 0; iChan < 48; iChan++) {
        int num(0);
        for (const tt::FrameStub& frame : streams[offset + iChan])
          if (frame.first.isNonnull())
            num++;
        his_[iChan]->Fill(num);
        prof_->Fill(iChan, num);
      }
    }
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::AnalyzerTB);
