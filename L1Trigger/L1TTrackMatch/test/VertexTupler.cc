#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TTClusterAssociationMap.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <cmath>

using namespace edm;

class SimpleVertexAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimpleVertexAnalyzer(const edm::ParameterSet& iConfig);
  ~SimpleVertexAnalyzer() override {}

private:
  void beginJob() override;
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override {}

  // Input tags
  edm::InputTag SimVertexInputTag;
  edm::InputTag TrackingParticleInputTag;
  edm::InputTag MCTruthClusterInputTag;

  // Tokens
  edm::EDGetTokenT<std::vector<SimVertex>> SimVertexToken_;
  edm::EDGetTokenT<std::vector<TrackingParticle>> TrackingParticleToken_;
  edm::EDGetTokenT<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> ttClusterMCTruthToken_;

  // ES tokens
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

  // ROOT tree and branches
  TTree* eventTree;
  
  // SimVertex branches
  std::vector<float>* m_pv_x;
  std::vector<float>* m_pv_y;
  std::vector<float>* m_pv_z;
  std::vector<float>* m_pv_t;
  
  // TrackingParticle branches
  std::vector<float>* m_tp_pt;
  std::vector<float>* m_tp_eta;
  std::vector<float>* m_tp_phi;
  std::vector<float>* m_tp_lz;
  std::vector<float>* m_tp_lxy;
  std::vector<float>* m_tp_z0;
  std::vector<float>* m_tp_d0;
  std::vector<int>* m_tp_pdgid;
  std::vector<int>* m_tp_eventid;
  std::vector<int>* m_tp_charge;

  // Parameters
  int MyProcess;
  double TP_minPt;
  double TP_maxEta;
  double TP_maxZ0;
  bool DebugMode;
};

SimpleVertexAnalyzer::SimpleVertexAnalyzer(const edm::ParameterSet& iConfig) {
  SimVertexInputTag = iConfig.getParameter<edm::InputTag>("SimVertexInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  MCTruthClusterInputTag = iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  
  // Get parameters
  MyProcess = iConfig.getParameter<int>("MyProcess");
  TP_minPt = iConfig.getParameter<double>("TP_minPt");
  TP_maxEta = iConfig.getParameter<double>("TP_maxEta");
  TP_maxZ0 = iConfig.getParameter<double>("TP_maxZ0");
  DebugMode = iConfig.getParameter<bool>("DebugMode");
  
  SimVertexToken_ = consumes<std::vector<SimVertex>>(SimVertexInputTag);
  TrackingParticleToken_ = consumes<std::vector<TrackingParticle>>(TrackingParticleInputTag);
  ttClusterMCTruthToken_ = consumes<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>>(MCTruthClusterInputTag);
  
  bFieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord>(edm::ESInputTag("", ""));
  
  usesResource(TFileService::kSharedResource);
}

void SimpleVertexAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs.isAvailable()) return;

  // Initialize SimVertex branches
  m_pv_x = new std::vector<float>;
  m_pv_y = new std::vector<float>;
  m_pv_z = new std::vector<float>;
  m_pv_t = new std::vector<float>;
  
  // Initialize TrackingParticle branches
  m_tp_pt = new std::vector<float>;
  m_tp_eta = new std::vector<float>;
  m_tp_phi = new std::vector<float>;
  m_tp_lz = new std::vector<float>;
  m_tp_lxy = new std::vector<float>;
  m_tp_z0 = new std::vector<float>;
  m_tp_d0 = new std::vector<float>;
  m_tp_pdgid = new std::vector<int>;
  m_tp_eventid = new std::vector<int>;
  m_tp_charge = new std::vector<int>;

  eventTree = fs->make<TTree>("vertexTree", "Vertex tree");
  
  // SimVertex branches
  eventTree->Branch("mc_vertex_x", &m_pv_x);
  eventTree->Branch("mc_vertex_y", &m_pv_y);
  eventTree->Branch("mc_vertex_z", &m_pv_z);
  eventTree->Branch("mc_vertex_t", &m_pv_t);
  
  // TrackingParticle branches
  eventTree->Branch("tp_pt", &m_tp_pt);
  eventTree->Branch("tp_eta", &m_tp_eta);
  eventTree->Branch("tp_phi", &m_tp_phi);
  eventTree->Branch("tp_lz", &m_tp_lz);
  eventTree->Branch("tp_lxy", &m_tp_lxy);
  eventTree->Branch("tp_z0", &m_tp_z0);
  eventTree->Branch("tp_d0", &m_tp_d0);
  eventTree->Branch("tp_pdgid", &m_tp_pdgid);
  eventTree->Branch("tp_eventid", &m_tp_eventid);
  eventTree->Branch("tp_charge", &m_tp_charge);
}

void SimpleVertexAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Clear SimVertex vectors
  m_pv_x->clear();
  m_pv_y->clear();
  m_pv_z->clear();
  m_pv_t->clear();
  
  // Clear TrackingParticle vectors
  m_tp_pt->clear();
  m_tp_eta->clear();
  m_tp_phi->clear();
  m_tp_lz->clear();
  m_tp_lxy->clear();
  m_tp_z0->clear();
  m_tp_d0->clear();
  m_tp_pdgid->clear();
  m_tp_eventid->clear();
  m_tp_charge->clear();

  // Get MC truth vertex
  edm::Handle<std::vector<SimVertex>> SimVertexHandle;
  iEvent.getByToken(SimVertexToken_, SimVertexHandle);
  if (SimVertexHandle.isValid()) {
    if (!SimVertexHandle->empty()) {
      const SimVertex simPVh = *(SimVertexHandle->begin());
      m_pv_x->push_back(simPVh.position().x());
      m_pv_y->push_back(simPVh.position().y());
      m_pv_z->push_back(simPVh.position().z());
      m_pv_t->push_back(simPVh.position().t());
    }
  }

  // Get MC truth cluster association map
  edm::Handle<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);

  // Get magnetic field
  const MagneticField& bField = iSetup.getData(bFieldToken_);

  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles with cuts
  // ----------------------------------------------------------------------------------------------
  edm::Handle<std::vector<TrackingParticle>> TrackingParticleHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);
  
  if (TrackingParticleHandle.isValid()) {
    if (DebugMode)
      edm::LogVerbatim("Tracklet") << "\n Loop over tracking particles!";

    int this_tp = 0;
    std::vector<TrackingParticle>::const_iterator iterTP;
    for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
      edm::Ptr<TrackingParticle> tp_ptr(TrackingParticleHandle, this_tp);
      this_tp++;

      int tmp_eventid = iterTP->eventId().event();
      if (MyProcess != 1 && tmp_eventid > 0)
        continue;  //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)

      float tmp_tp_pt = iterTP->pt();
      float tmp_tp_charge = iterTP->charge();
      float tmp_tp_eta = iterTP->eta();

      if (tmp_tp_pt < TP_minPt)  // Save CPU by applying these cuts here.
        continue;
      if (tmp_tp_charge == 0.)
        continue;
      if (std::abs(tmp_tp_eta) > TP_maxEta)
        continue;

      float tmp_tp_phi = iterTP->phi();
      float tmp_tp_t = iterTP->tanl();
      float tmp_tp_vz = iterTP->vz();
      float tmp_tp_vx = iterTP->vx();
      float tmp_tp_vy = iterTP->vy();
      int tmp_tp_pdgid = iterTP->pdgId();

      // ----------------------------------------------------------------------------------------------
      // get d0/z0 accurately propagated back to the IP

      float delx = -tmp_tp_vx;
      float dely = -tmp_tp_vy;

      float b_field = bField.inTesla(GlobalPoint(0, 0, 0)).z();
      float c_converted = CLHEP::c_light / 1.0E5;
      float r2_inv = tmp_tp_charge * c_converted * b_field / tmp_tp_pt / 2.0;

      float tmp_tp_x0p = delx - (1. / (2. * r2_inv) * sin(tmp_tp_phi));
      float tmp_tp_y0p = dely + (1. / (2. * r2_inv) * cos(tmp_tp_phi));
      float tmp_tp_rp = sqrt(tmp_tp_x0p * tmp_tp_x0p + tmp_tp_y0p * tmp_tp_y0p);
      float tmp_tp_d0 = tmp_tp_charge * tmp_tp_rp - (1. / (2. * r2_inv));

      float delphi = reco::deltaPhi(tmp_tp_phi, atan2(-r2_inv * tmp_tp_x0p, r2_inv * tmp_tp_y0p));
      float tmp_tp_z0 = tmp_tp_vz + tmp_tp_t * delphi / (2.0 * r2_inv);

      // ----------------------------------------------------------------------------------------------

      if (MyProcess == 13 && abs(tmp_tp_pdgid) != 13)
        continue;
      if (MyProcess == 11 && abs(tmp_tp_pdgid) != 11)
        continue;
      if ((MyProcess == 6 || MyProcess == 15 || MyProcess == 211) && abs(tmp_tp_pdgid) != 211)
        continue;

      if (std::abs(tmp_tp_z0) > TP_maxZ0)
        continue;

      float tmp_tp_lxy = sqrt(tmp_tp_vx * tmp_tp_vx + tmp_tp_vy * tmp_tp_vy);
      float tmp_tp_lz = tmp_tp_vz;

      if (DebugMode)
        edm::LogVerbatim("Tracklet") << "Tracking particle, pt: " << tmp_tp_pt << " eta: " << tmp_tp_eta
                                     << " phi: " << tmp_tp_phi << " z0: " << tmp_tp_z0 << " d0: " << tmp_tp_d0
                                     << " pdgid: " << tmp_tp_pdgid << " eventID: " << iterTP->eventId().event()
                                     << " ttclusters " << MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size();

      // ----------------------------------------------------------------------------------------------
      // only consider TPs associated with >= 1 cluster

      if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).empty()) {
        if (DebugMode)
          edm::LogVerbatim("Tracklet") << "No matching TTClusters for TP, continuing...";
        continue;
      }

      // Fill TrackingParticle branches
      m_tp_pt->push_back(tmp_tp_pt);
      m_tp_eta->push_back(tmp_tp_eta);
      m_tp_phi->push_back(tmp_tp_phi);
      m_tp_lz->push_back(tmp_tp_lz);
      m_tp_lxy->push_back(tmp_tp_lxy);
      m_tp_z0->push_back(tmp_tp_z0);
      m_tp_d0->push_back(tmp_tp_d0);
      m_tp_pdgid->push_back(tmp_tp_pdgid);
      m_tp_eventid->push_back(tmp_eventid);
      m_tp_charge->push_back(tmp_tp_charge);
    }
  }

  // Fill the tree
  eventTree->Fill();
}

DEFINE_FWK_MODULE(SimpleVertexAnalyzer);