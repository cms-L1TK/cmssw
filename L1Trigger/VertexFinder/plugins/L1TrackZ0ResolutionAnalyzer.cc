#include <memory>
#include <string>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// L1 vertices (the one you already use successfully)
#include "DataFormats/L1Trigger/interface/Vertex.h"

// L1 tracks (RefVector of TTTracks)
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Common/interface/RefVector.h"

// Optional histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"

class L1TrackZ0ResolutionAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  using TTTrackType = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection = std::vector<TTTrackType>;
  using TTTrackRefCollection = edm::RefVector<TTTrackCollection>;

  explicit L1TrackZ0ResolutionAnalyzer(const edm::ParameterSet& iConfig)
      : verticesTag_(iConfig.getParameter<edm::InputTag>("vertices")),
        selectedTracksTag_(iConfig.getParameter<edm::InputTag>("selectedTracks")),
        associatedTracksTag_(iConfig.getParameter<edm::InputTag>("associatedTracks")),
        verticesToken_(consumes<l1t::VertexCollection>(verticesTag_)),
        selectedTracksToken_(consumes<TTTrackRefCollection>(selectedTracksTag_)),
        associatedTracksToken_(consumes<TTTrackRefCollection>(associatedTracksTag_)),
        debug_(iConfig.getUntrackedParameter<bool>("debug", false)) {
    usesResource("TFileService");

    edm::Service<TFileService> fs;

    h_dz_sel_ = fs->make<TH1F>("dz_selected",
                              "Selected tracks: z0(track) - z0(L1 PV);#Delta z_{0} [cm];Tracks",
                              400, -10.0, 10.0);
    h_dz_assoc_ = fs->make<TH1F>("dz_associated",
                                "Associated tracks: z0(track) - z0(L1 PV);#Delta z_{0} [cm];Tracks",
                                400, -10.0, 10.0);

    h_dz_sel_vs_eta_ = fs->make<TH2F>("dz_selected_vs_abseta",
                                     "Selected: #Delta z_{0} vs |#eta|;|#eta|;#Delta z_{0} [cm]",
                                     60, 0.0, 3.0, 400, -10.0, 10.0);
    h_dz_assoc_vs_eta_ = fs->make<TH2F>("dz_associated_vs_abseta",
                                       "Associated: #Delta z_{0} vs |#eta|;|#eta|;#Delta z_{0} [cm]",
                                       60, 0.0, 3.0, 400, -10.0, 10.0);
    h_dz_assoc_vs_mva_ = fs->make<TH2F>("dz_associated_vs_trk_mva1",
                                    "Associated: #Delta z_{0} vs MVA;MVA;#Delta z_{0} [cm]",
                                    50, 0.0, 1.0, 400, -10.0, 10.0);

    h_dz_sel_vs_pt_ = fs->make<TH2F>("dz_selected_vs_pt",
                                    "Selected: #Delta z_{0} vs p_{T};p_{T} [GeV];#Delta z_{0} [cm]",
                                    60, 0.0, 120.0, 400, -10.0, 10.0);
    h_dz_assoc_vs_pt_ = fs->make<TH2F>("dz_associated_vs_pt",
                                      "Associated: #Delta z_{0} vs p_{T};p_{T} [GeV];#Delta z_{0} [cm]",
                                      60, 0.0, 120.0, 400, -10.0, 10.0);

    h_nSel_ = fs->make<TH1F>("nSelected", "N selected tracks;N;Events", 301, -0.5, 300.5);
    h_nAssoc_ = fs->make<TH1F>("nAssociated", "N associated tracks;N;Events", 301, -0.5, 300.5);
    h_nVertices_ = fs->make<TH1F>("nVertices", "N associated tracks;N;Events", 301, -0.5, 300.5);
  }

  void analyze(const edm::Event& iEvent, const edm::EventSetup&) override {
    const auto& vertices = iEvent.get(verticesToken_);
    if (vertices.empty()) {
      if (debug_) edm::LogPrint("L1TrackZ0ResolutionAnalyzer") << "No L1 vertices in event " << iEvent.id();
      return;
    }

    const float pvZ0 = vertices.front().z0(); // assumes collection is sorted so [0] is leading PV

    const auto& selTracks = iEvent.get(selectedTracksToken_);
    const auto& assocTracks = iEvent.get(associatedTracksToken_);

    h_nSel_->Fill(selTracks.size());
    h_nAssoc_->Fill(assocTracks.size());
    h_nVertices_->Fill(vertices.front().tracks().size());

    if (debug_) {
      edm::LogPrint("L1TrackZ0ResolutionAnalyzer")
          << "Event " << iEvent.id()
          << " PV z0=" << pvZ0
          << " nSel=" << selTracks.size()
          << " nAssoc=" << assocTracks.size();
    }

    // Selected tracks residuals
    for (const auto& trkRef : selTracks) {
      const auto& trk = *trkRef;
      const float dz = trk.z0() - pvZ0;
      const float abseta = std::abs(trk.momentum().eta());
      const float pt = trk.momentum().perp();

      h_dz_sel_->Fill(dz);
      h_dz_sel_vs_eta_->Fill(abseta, dz);
      h_dz_sel_vs_pt_->Fill(pt, dz);
    }

    // Associated tracks residuals
    for (const auto& trkRef : assocTracks) {
      const auto& trk = *trkRef;
      // apply mva cut
      // if (trk.trkMVA1() < 0.5) continue;
      const float dz = trk.z0() - pvZ0;
      const float abseta = std::abs(trk.momentum().eta());
      const float pt = trk.momentum().perp();

      h_dz_assoc_vs_mva_->Fill(trk.trkMVA1(), dz);
      h_dz_assoc_->Fill(dz);
      h_dz_assoc_vs_eta_->Fill(abseta, dz);
      h_dz_assoc_vs_pt_->Fill(pt, dz);
    }
  }

private:
    edm::InputTag verticesTag_;
    edm::InputTag selectedTracksTag_;
    edm::InputTag associatedTracksTag_;

    edm::EDGetTokenT<l1t::VertexCollection> verticesToken_;
    edm::EDGetTokenT<TTTrackRefCollection> selectedTracksToken_;
    edm::EDGetTokenT<TTTrackRefCollection> associatedTracksToken_;

    bool debug_;

    TH1F* h_dz_sel_{nullptr};
    TH1F* h_dz_assoc_{nullptr};
    TH2F* h_dz_sel_vs_eta_{nullptr};
    TH2F* h_dz_assoc_vs_eta_{nullptr};
    TH2F* h_dz_sel_vs_pt_{nullptr};
    TH2F* h_dz_assoc_vs_pt_{nullptr};
    TH2F* h_dz_assoc_vs_mva_{nullptr};
    TH1F* h_nSel_{nullptr};
    TH1F* h_nAssoc_{nullptr};
    TH1F* h_nVertices_{nullptr};

};

DEFINE_FWK_MODULE(L1TrackZ0ResolutionAnalyzer);
