#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include <TTree.h>
#include <vector>

class VertexWordAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit VertexWordAnalyzer(const edm::ParameterSet&);
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<l1t::VertexWordCollection> vertexToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVertexToken_;

  TTree* tree_;

  // l1t vertex-word branches (one entry per event, vectors per event)
  std::vector<int>*          v_valid_;
  std::vector<double>*       v_z0_;
  std::vector<unsigned int>* v_multiplicity_;
  std::vector<double>*       v_pt_;
  std::vector<unsigned int>* v_quality_;
  std::vector<unsigned int>* v_inverseMultiplicity_;
  std::vector<unsigned int>* v_unassigned_;

  // all sim vertices (vectors per event)
  std::vector<float>* sv_x_;
  std::vector<float>* sv_y_;
  std::vector<float>* sv_z_;
  std::vector<float>* sv_t_;

  // "PV" = first sim vertex (scalars per event)
  float pv_x_;
  float pv_y_;
  float pv_z_;
  float pv_t_;
  int   pv_found_;
};

VertexWordAnalyzer::VertexWordAnalyzer(const edm::ParameterSet& iConfig)
    : vertexToken_(consumes<l1t::VertexWordCollection>(iConfig.getParameter<edm::InputTag>("vertexTag"))),
      simVertexToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertexTag"))) {
  usesResource(TFileService::kSharedResource);
}

void VertexWordAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs.isAvailable())
    return;

  tree_ = fs->make<TTree>("vertexTree", "Vertex word tree");

  v_valid_               = new std::vector<int>;
  v_z0_                  = new std::vector<double>;
  v_multiplicity_        = new std::vector<unsigned int>;
  v_pt_                  = new std::vector<double>;
  v_quality_             = new std::vector<unsigned int>;
  v_inverseMultiplicity_ = new std::vector<unsigned int>;
  v_unassigned_          = new std::vector<unsigned int>;

  sv_x_ = new std::vector<float>;
  sv_y_ = new std::vector<float>;
  sv_z_ = new std::vector<float>;
  sv_t_ = new std::vector<float>;

  tree_->Branch("vtx_valid",               &v_valid_);
  tree_->Branch("vtx_z0",                  &v_z0_);
  tree_->Branch("vtx_multiplicity",        &v_multiplicity_);
  tree_->Branch("vtx_pt",                  &v_pt_);
  tree_->Branch("vtx_quality",             &v_quality_);
  tree_->Branch("vtx_inverseMultiplicity", &v_inverseMultiplicity_);
  tree_->Branch("vtx_unassigned",          &v_unassigned_);

  tree_->Branch("sim_vtx_x", &sv_x_);
  tree_->Branch("sim_vtx_y", &sv_y_);
  tree_->Branch("sim_vtx_z", &sv_z_);
  tree_->Branch("sim_vtx_t", &sv_t_);

  tree_->Branch("sim_pv_x", &pv_x_, "sim_pv_x/F");
  tree_->Branch("sim_pv_y", &pv_y_, "sim_pv_y/F");
  tree_->Branch("sim_pv_z", &pv_z_, "sim_pv_z/F");
  tree_->Branch("sim_pv_t", &pv_t_, "sim_pv_t/F");
  tree_->Branch("sim_pv_found", &pv_found_, "sim_pv_found/I");
}

void VertexWordAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  v_valid_->clear();
  v_z0_->clear();
  v_multiplicity_->clear();
  v_pt_->clear();
  v_quality_->clear();
  v_inverseMultiplicity_->clear();
  v_unassigned_->clear();

  sv_x_->clear();
  sv_y_->clear();
  sv_z_->clear();
  sv_t_->clear();

  pv_x_ = 0.f;
  pv_y_ = 0.f;
  pv_z_ = 0.f;
  pv_t_ = 0.f;
  pv_found_ = 0;

  edm::Handle<l1t::VertexWordCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  edm::Handle<std::vector<SimVertex>> simVertices;
  iEvent.getByToken(simVertexToken_, simVertices);

  // --- l1t vertex words ---
  for (const auto& v : *vertices) {
    v_valid_->push_back(v.valid() ? 1 : 0);
    v_z0_->push_back(v.z0());
    v_multiplicity_->push_back(v.multiplicity());
    v_pt_->push_back(v.pt());
    v_quality_->push_back(v.quality());
    v_inverseMultiplicity_->push_back(v.inverseMultiplicity());
    v_unassigned_->push_back(v.unassigned());
  }

  // --- sim vertices ---
  if (simVertices.isValid()) {
    if (!simVertices->empty()) {
      for (const auto& sv : *simVertices) {
        sv_x_->push_back(sv.position().x());
        sv_y_->push_back(sv.position().y());
        sv_z_->push_back(sv.position().z());
        sv_t_->push_back(sv.position().t());
      }

      // "PV" = first sim vertex
      const SimVertex& simPV = simVertices->front();
      pv_x_ = simPV.position().x();
      pv_y_ = simPV.position().y();
      pv_z_ = simPV.position().z();
      pv_t_ = simPV.position().t();
      pv_found_ = 1;
    } else {
      edm::LogWarning("MissingCollectionEntries") << "Warning: SimVertex collection is empty!";
    }
  } else {
    edm::LogWarning("DataNotFound") << "Warning: SimVertexHandle not found in the event";
  }

  tree_->Fill();
}

DEFINE_FWK_MODULE(VertexWordAnalyzer);