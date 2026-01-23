#include <memory>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1Trigger/interface/Vertex.h"

// Optional histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"

class L1VertexAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1VertexAnalyzer(const edm::ParameterSet& iConfig)
      : verticesTag_(iConfig.getParameter<edm::InputTag>("vertices")),
        verticesToken_(consumes<l1t::VertexCollection>(verticesTag_)),
        debug_(iConfig.getUntrackedParameter<bool>("debug", false)) {
    usesResource("TFileService");

    edm::Service<TFileService> fs;
    h_nVtx_   = fs->make<TH1F>("nVtx",   "L1 vertices per event;N(vertices);Events", 101, -0.5, 100.5);
    h_vtxZ0_  = fs->make<TH1F>("vtxZ0",  "L1 vertex z0;z0 [cm];Vertices",          240, -30.0, 30.0);
    h_vtxPt_  = fs->make<TH1F>("vtxPt",  "L1 vertex p_{T};p_{T} [GeV];Vertices",   200, 0.0, 200.0);
    h_ptVsZ0_ = fs->make<TH2F>("ptVsZ0", "L1 vertex p_{T} vs z0;z0 [cm];p_{T} [GeV]", 240, -30.0, 30.0, 200, 0.0, 200.0);
    h_nTrk_   = fs->make<TH1F>("nTrk",   "Tracks per L1 vertex;N(tracks);Vertices", 101, -0.5, 100.5);
  }

  void analyze(const edm::Event& iEvent, const edm::EventSetup&) override {
    const auto& vertices = iEvent.get(verticesToken_);

    h_nVtx_->Fill(vertices.size());

    if (debug_) {
      edm::LogPrint("L1VertexAnalyzer") << "Event " << iEvent.id() << " has " << vertices.size() << " L1 vertices";
    }

    for (const auto& vtx : vertices) {
      const float z0 = vtx.z0();
      const float pt = vtx.pt();
      const unsigned ntrk = vtx.tracks().size();

      h_vtxZ0_->Fill(z0);
      h_vtxPt_->Fill(pt);
      h_ptVsZ0_->Fill(z0, pt);
      h_nTrk_->Fill(ntrk);

      if (debug_) {
        edm::LogPrint("L1VertexAnalyzer")
            << "  vtx: z0=" << z0 << " cm, pt=" << pt << " GeV, nTracks=" << ntrk;
      }
    }
  }

private:
  edm::InputTag verticesTag_;
  edm::EDGetTokenT<l1t::VertexCollection> verticesToken_;
  bool debug_;

  TH1F* h_nVtx_{nullptr};
  TH1F* h_vtxZ0_{nullptr};
  TH1F* h_vtxPt_{nullptr};
  TH2F* h_ptVsZ0_{nullptr};
  TH1F* h_nTrk_{nullptr};
};

DEFINE_FWK_MODULE(L1VertexAnalyzer);
