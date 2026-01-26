#include <memory>
#include <string>
#include <vector>
#include <cmath>

// CMSSW core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Data formats
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Common/interface/RefVector.h"

// TFileService + ROOT
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

class L1TrackZ0NtupleProducer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  using TTTrackType = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection = std::vector<TTTrackType>;
  using TTTrackRefCollection = edm::RefVector<TTTrackCollection>;

  explicit L1TrackZ0NtupleProducer(const edm::ParameterSet& iConfig)
      : verticesTag_(iConfig.getParameter<edm::InputTag>("vertices")),
        associatedTracksTag_(iConfig.getParameter<edm::InputTag>("associatedTracks")),
        verticesToken_(consumes<l1t::VertexCollection>(verticesTag_)),
        associatedTracksToken_(consumes<TTTrackRefCollection>(associatedTracksTag_)),
        debug_(iConfig.getUntrackedParameter<bool>("debug", false)),
        storeDz_(iConfig.getUntrackedParameter<bool>("storeDz", true)) {
    usesResource("TFileService");

    edm::Service<TFileService> fs;

    tree_ = fs->make<TTree>("z0Tree", "L1 PV and associated track z0/mva");

    tree_->Branch("run",   &run_,   "run/i");
    tree_->Branch("lumi",  &lumi_,  "lumi/i");
    tree_->Branch("event", &event_, "event/l");

    tree_->Branch("leading_vertex_z0", &leading_vertex_z0_);
    tree_->Branch("trk_z0s",           &trk_z0s_);
    tree_->Branch("trk_mvas",          &trk_mvas_);

    // Handy extras
    tree_->Branch("trk_pts",  &trk_pts_);
    tree_->Branch("trk_etas", &trk_etas_);

    if (storeDz_) {
      tree_->Branch("trk_dz", &trk_dz_);
    }
  }

  void analyze(const edm::Event& iEvent, const edm::EventSetup&) override {
    // Clear per-event buffers
    leading_vertex_z0_.clear();
    trk_z0s_.clear();
    trk_mvas_.clear();
    trk_pts_.clear();
    trk_etas_.clear();
    trk_dz_.clear();

    run_   = iEvent.id().run();
    lumi_  = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();

    const auto& vertices = iEvent.get(verticesToken_);
    if (vertices.empty()) {
      if (debug_) {
        edm::LogPrint("L1TrackZ0NtupleProducer") << "No L1 vertices in event " << iEvent.id();
      }
      return;
    }

    // Assume leading PV is front()
    const float pvZ0 = vertices.front().z0();
    leading_vertex_z0_.push_back(pvZ0);

    const auto& assocTracks = iEvent.get(associatedTracksToken_);

    if (debug_) {
      edm::LogPrint("L1TrackZ0NtupleProducer")
          << "Event " << iEvent.id()
          << " pvZ0=" << pvZ0
          << " nAssoc=" << assocTracks.size();
    }

    for (const auto& trkRef : assocTracks) {
      const auto& trk = *trkRef;

      const float z0  = trk.z0();
      const float pt  = trk.momentum().perp();
      const float eta = trk.momentum().eta();
      const float mva = trk.trkMVA1();

      trk_z0s_.push_back(z0);
      trk_mvas_.push_back(mva);
      trk_pts_.push_back(pt);
      trk_etas_.push_back(eta);

      if (storeDz_) {
        trk_dz_.push_back(z0 - pvZ0);
      }
    }

    tree_->Fill();
  }

private:
  // Inputs
  edm::InputTag verticesTag_;
  edm::InputTag associatedTracksTag_;

  // Tokens
  edm::EDGetTokenT<l1t::VertexCollection> verticesToken_;
  edm::EDGetTokenT<TTTrackRefCollection> associatedTracksToken_;

  // Config
  bool debug_;
  bool storeDz_;

  // Output tree
  TTree* tree_{nullptr};

  // Event IDs
  UInt_t run_{0};
  UInt_t lumi_{0};
  ULong64_t event_{0};

  // Branch buffers
  std::vector<float> leading_vertex_z0_; // size = 1 per event
  std::vector<float> trk_z0s_;
  std::vector<float> trk_mvas_;
  std::vector<float> trk_pts_;
  std::vector<float> trk_etas_;
  std::vector<float> trk_dz_;            // optional
};

DEFINE_FWK_MODULE(L1TrackZ0NtupleProducer);
