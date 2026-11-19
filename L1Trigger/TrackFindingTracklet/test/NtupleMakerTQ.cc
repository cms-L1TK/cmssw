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
#include "SimDataFormats/Associations/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/StubAssociation.h"
#include "L1Trigger/TrackTrigger/interface/Associator.h"
#include "L1Trigger/TrackFindingTracklet/interface/DataFormats.h"

#include <TTree.h>

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <cmath>
#include <numeric>
#include <sstream>

namespace trklet {

  /*! \class  trklet::NtupleMakerTQ
   *  \brief  Class to write data for Track Quality training
   *  \author Thomas Schuh
   *  \date   2026, Oct
   */
  class NtupleMakerTQ : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
  public:
    NtupleMakerTQ(const edm::ParameterSet& iConfig);
    void beginJob() override {}
    void beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override;
    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    void endRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override {}
    void endJob() override {}

  private:
    // ED input token of stubs
    edm::EDGetTokenT<tt::StreamsStub> edGetTokenStubs_;
    // ED input token of tracks
    edm::EDGetTokenT<tt::StreamsTrack> edGetTokenTracks_;
    // ED output token for stub association for fake rate
    edm::EDGetTokenT<tt::StubAssociation> edGetToken_;
    // Associator token
    edm::ESGetToken<tt::Associator, trackerDTC::SetupRcd> esGetTokenAssociator_;
    // Setup token
    edm::ESGetToken<Setup, trackerDTC::SetupRcd> esGetTokenSetup_;
    // DataFormats token
    edm::ESGetToken<DataFormats, trackerDTC::SetupRcd> esGetTokenDataFormats_;
    // TTree
    TTree* tTree_;
    // training data
    bool fake_;
    int hits_;
    int cot_;
    int z0_;
    int chi20_;
    int chi21_;
  };

  NtupleMakerTQ::NtupleMakerTQ(const edm::ParameterSet& iConfig) {
    usesResource("TFileService");
    // book in- and output ED products
    const std::string& labelKF = iConfig.getParameter<std::string>("OutputLabelKF");
    const std::string& labelTQ = iConfig.getParameter<std::string>("OutputLabelTQ");
    const std::string& branchStubs = iConfig.getParameter<std::string>("BranchStubs");
    const std::string& branchTracks = iConfig.getParameter<std::string>("BranchTracks");
    edGetTokenStubs_ = consumes(edm::InputTag(labelKF, branchStubs));
    edGetTokenTracks_ = consumes(edm::InputTag(labelTQ, branchTracks));
    const std::string& labelMC = iConfig.getParameter<std::string>("LabelMC");
    const std::string& branchMC = iConfig.getParameter<std::string>("BranchFake");
    edGetToken_ = consumes(edm::InputTag(labelMC, branchMC));
    // book ES products
    esGetTokenAssociator_ = esConsumes();
    esGetTokenSetup_ = esConsumes();
    esGetTokenDataFormats_ = esConsumes();
  }

  void NtupleMakerTQ::beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) {
    // book histograms
    edm::Service<TFileService> fs;
    TFileDirectory dir;
    dir = fs->mkdir("");
    tTree_ = dir.make<TTree>("events", "event features for ML");
    tTree_->Branch("fake", &fake_, "Fake/B");
    tTree_->Branch("hits", &hits_, "Hits/I");
    tTree_->Branch("cot", &cot_, "Cot/I");
    tTree_->Branch("z0", &z0_, "Z0/I");
    tTree_->Branch("chi20", &chi20_, "Chi20/I");
    tTree_->Branch("chi21", &chi21_, "Chi21/I");
  }

  void NtupleMakerTQ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // read in tracks and stubs products
    const tt::StreamsStub& streamsStub = iEvent.get(edGetTokenStubs_);
    const tt::StreamsTrack& streamsTrack = iEvent.get(edGetTokenTracks_);
    const Setup* setup = &iSetup.getData(esGetTokenSetup_);
    const DataFormats* df = &iSetup.getData(esGetTokenDataFormats_);
    // read in MCTruth
    tt::Associator assoc = iSetup.getData(esGetTokenAssociator_);
    assoc.consume(iEvent.get(edGetToken_));
    // analyze and associate tracks with TrackingParticles
    for (int region = 0; region < setup->sysNumRegion(); region++) {
      const int offset = region * setup->kfNumLayers();
      const tt::StreamTrack& streamTrackTQ = streamsTrack[2 * region + 1];
      const tt::StreamTrack& streamTrackKF = streamsTrack[2 * region + 0];
      const int numFrames = streamTrackTQ.size();
      for (int frame = 0; frame < numFrames; frame++) {
        if (streamTrackTQ[frame].first.isNull())
          continue;
        const TrackKF trackKF(streamTrackKF[frame], df);
        const TrackTQ trackTQ(streamTrackTQ[frame], df);
        std::vector<TTStubRef> ttStubRefs;
        ttStubRefs.reserve(setup->kfNumLayers());
        for (int layer = 0; layer < setup->kfNumLayers(); layer++) {
          const TTStubRef& ttStubRef = streamsStub[offset + layer][frame].first;
          if (ttStubRef.isNonnull())
            ttStubRefs.push_back(ttStubRef);
        }
        const std::vector<TPPtr>& tpPtrs = assoc.associateFinal(ttStubRefs);
        // fill TTree
        fake_ = tpPtrs.empty();
        hits_ = trackTQ.hitPattern().val();
        cot_ = df->format(Variable::cot, Process::kf).integer(trackKF.cot());
        z0_ = df->format(Variable::z0, Process::kf).integer(trackKF.z0());
        chi20_ = df->format(Variable::chi20, Process::tq).integer(trackTQ.chi20());
        chi21_ = df->format(Variable::chi21, Process::tq).integer(trackTQ.chi21());
        tTree_->Fill();
      }
    }
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::NtupleMakerTQ);
