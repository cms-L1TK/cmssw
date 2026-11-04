//===========================================================================
// Author: Ian Tomalin
// Date: 2026
//
// Plots P2 TTStubs (=stubs used for L1 tracking),
// including comparisons with truth,
// e.g. for stub bend resolution plots.
//
// N.B. This code is based on that used for offline P2 clusters,
// RecoLocalTracker/SiPhase2Clusterizer/test/ClustersValidation.cc
//
// Can be run from AnalyzerClusterStub_cfg.py
//===========================================================================

// system includes
#include <map>
#include <vector>
#include <algorithm>

// user includes
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "SimDataFormats/Associations/interface/TTStubAssociationMap.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "FWCore/Utilities/interface/Exception.h"

// ROOT includes
#include <TH2F.h>
#include <TH1F.h>
#include <THStack.h>

struct StubHistos {
  // use TH1D instead of TH1F to avoid stauration at 2^31
  // above this increments with +1 don't work for float, need double

  TH1D* numberStubs[3];
  TH1D* numberStubMatch[3];
  // Bend histos made only for genuine stubs from prompt TP with Pt > 2 GeV.
  TH1F* bend[3];
  TH1F* bendRes[3];
};

class AnalyzerStub : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit AnalyzerStub(const edm::ParameterSet&);
  ~AnalyzerStub() override;
  void beginJob() override;
  void beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override;
  void endRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override {}
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  std::map<unsigned int, StubHistos>::iterator createLayerHistograms(unsigned int);

  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> esTokenGeom_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> esTokenTopo_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> esTokenBfield_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  const edm::EDGetTokenT<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>> tokenStubs_;
  edm::EDGetTokenT<TTStubAssociationMap<Ref_Phase2TrackerDigi_>> tokenStubAssoc_;

  const bool catECasRings_;

  TH2F* trackerLayout_;
  TH2F* trackerLayoutXY_;
  TH2F* trackerLayoutXYBar_;
  TH2F* trackerLayoutXYEC_;

  std::map<unsigned int, StubHistos> histograms_;

  // Cuts for TrackingParticles of most interest.
  const float maxEta_ = 2.5, minPt_ = 2.0, maxVertR_ = 1.0, maxVertZ_ = 15.0;
  const bool useSignal_ = true, useIntime_ = true, useCharged_ = true, useStable_ = true;
  TrackingParticleSelector tpSel_;

  double bField_;
};

AnalyzerStub::AnalyzerStub(const edm::ParameterSet& conf)
    : esTokenGeom_(esConsumes()),
      esTokenTopo_(esConsumes()),
      esTokenBfield_(esConsumes<edm::Transition::BeginRun>()),
      tpToken_(consumes<std::vector<TrackingParticle>>(conf.getParameter<edm::InputTag>("trackingParticles"))),
      tokenStubs_(
          consumes<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>>(conf.getParameter<edm::InputTag>("ttstubs"))),
      tokenStubAssoc_(
          consumes<TTStubAssociationMap<Ref_Phase2TrackerDigi_>>(conf.getParameter<edm::InputTag>("ttstubassoc"))),
      catECasRings_(conf.getParameter<bool>("ECasRings")),
      tpSel_(
          minPt_, 9.9e9, -maxEta_, maxEta_, maxVertR_, maxVertZ_, 0, useSignal_, useIntime_, useCharged_, useStable_),
      bField_(0.) {
  usesResource(TFileService::kSharedResource);
}

AnalyzerStub::~AnalyzerStub() = default;

void AnalyzerStub::beginJob() {
  edm::Service<TFileService> fs;
  fs->file().cd("/");
  TFileDirectory td = fs->mkdir("Common");
  // Create common histograms
  trackerLayout_ = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 1200, 0.0, 120.0);
  trackerLayoutXY_ = td.make<TH2F>("XVsY", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
  trackerLayoutXYBar_ = td.make<TH2F>("XVsYBar", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
  trackerLayoutXYEC_ = td.make<TH2F>("XVsYEC", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
}

void AnalyzerStub::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) {
  /// Get CMS Bfield
  const MagneticField* theMagneticField = &iSetup.getData(esTokenBfield_);
  bField_ = theMagneticField->inTesla(GlobalPoint(0, 0, 0)).z();
}

void AnalyzerStub::analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {
  // Get tracking particles (truth)
  edm::Handle<std::vector<TrackingParticle>> tParts;
  event.getByToken(tpToken_, tParts);

  // Get the stubs
  edm::Handle<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>> stubs;
  event.getByToken(tokenStubs_, stubs);

  // Get association map
  edm::Handle<TTStubAssociationMap<Ref_Phase2TrackerDigi_>> stubsAssoc;
  event.getByToken(tokenStubAssoc_, stubsAssoc);

  // Get the geometry
  const TrackerGeometry* tkGeom = &eventSetup.getData(esTokenGeom_);
  const TrackerTopology* tkTopo = &eventSetup.getData(esTokenTopo_);

  /*
   * Validation   
   */

  // Number of stubs vs layer
  std::map<unsigned int, unsigned int> nStubs[3], nStubMatch[3];

  // Loop over modules
  for (auto DSViter = stubs->begin(); DSViter != stubs->end(); ++DSViter) {
    // Get the detector unit's id
    unsigned int rawid(DSViter->detId());
    DetId detId(rawid);
    bool barrel = tkTopo->side(detId) == 0;
    unsigned int layer = (not barrel) * 1000;  // don't split up endcap sides
    if (!layer) {
      layer += tkTopo->layer(detId);
    } else {
      layer += (catECasRings_ ? tkTopo->endcapRingP2(detId) * 10 : tkTopo->layer(detId));
    }
    TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
    unsigned int det = 0;
    if (mType == TrackerGeometry::ModuleType::Ph2PSP) {
      det = 1;
    } else if (mType == TrackerGeometry::ModuleType::Ph2SS) {
      det = 2;
    } else {
      throw cms::Exception("LogicError") << "Unknown detector type! " << det;
    }

    // initialize the nstub counters if they don't exist for this layer
    if (nStubs[det].find(layer) == nStubs[det].end()) {
      nStubs[det].emplace(layer, 0);
    }
    if (nStubMatch[det].find(layer) == nStubMatch[det].end()) {
      nStubMatch[det].emplace(layer, 0);
    }

    // Create histograms for the layer if they do not yet exist
    std::map<unsigned int, StubHistos>::iterator histogramLayer(histograms_.find(layer));
    if (histogramLayer == histograms_.end())
      histogramLayer = createLayerHistograms(layer);

    // Loop over the stubs in the detector unit
    for (auto stubIt = DSViter->begin(); stubIt != DSViter->end(); ++stubIt) {
      TTStubRef ttStubRef = edmNew::makeRefTo(stubs, stubIt);

      // MC truth association
      bool genuine = stubsAssoc->isGenuine(ttStubRef);

      // stub count
      ++(nStubs[det].at(layer));
      if (genuine) {
        ++(nStubMatch[det].at(layer));
      }

      // Get stub position
      const DetId detId1 = detId + 1;  // seed sensor
      const GeomDetUnit* detUnit1 = tkGeom->idToDetUnit(detId1);
      const PixelTopology* pixTopo =
          dynamic_cast<const PixelTopology*>(&(dynamic_cast<const PixelGeomDetUnit*>(detUnit1)->specificTopology()));
      const Plane& plane = dynamic_cast<const PixelGeomDetUnit*>(detUnit1)->surface();
      const MeasurementPoint& mp = ttStubRef->clusterRef(0)->findAverageLocalCoordinatesCentered();
      GlobalPoint posStub = plane.toGlobal(pixTopo->localPosition(mp));

      // Fill the position histograms
      trackerLayout_->Fill(posStub.z(), posStub.perp());
      trackerLayoutXY_->Fill(posStub.x(), posStub.y());
      if (layer < 1000) {
        trackerLayoutXYBar_->Fill(posStub.x(), posStub.y());
      } else {
        trackerLayoutXYEC_->Fill(posStub.x(), posStub.y());
      }

      if (genuine) {
        const TrackingParticlePtr& tp = stubsAssoc->findTrackingParticlePtr(ttStubRef);
        if (tpSel_(*tp)) {
          // Here we have genuine stub matched to prompt TP with Pt > 2 GeV.

          auto histogramLayer(histograms_.find(layer));
          if (histogramLayer == histograms_.end())
            throw cms::Exception("LogicError") << "No histo for an existing thing! This should not happen!";

          // Plot stub bend & bend resolution
          // First calculate true bend (based on code in SensorModule.cc & StubGL.cc).
          const float cSpeed = 1.0e8 * CLHEP::c_light;  // units of cm/s
          const float invPtToDphi = bField_ * cSpeed / 2.0e13;
          // Info about tracker module & two sensors in it.
          const GlobalPoint pos1 = GlobalPoint(detUnit1->position());
          const GeomDetUnit* detUnit2 = tkGeom->idToDetUnit(tkTopo->partnerDetId(detId1));
          const GlobalPoint pos2 = GlobalPoint(detUnit2->position());
          const float sensorSep = (pos2 - pos1).mag();
          const float pitch = pixTopo->pitch().first;
          const float pitchOverSep = pitch / sensorSep;
          const float deltaR = pos2.perp() - pos1.perp();
          const float deltaZ = pos2.z() - pos1.z();
          bool flipped = pos1.mag() > pos2.mag();
          // Module tilt angle w.r.t. beam axis (0=barrel,-+pi/2=+-ve endcap).
          const float tiltAngle = flipped ? atan2(deltaZ, -deltaR) : atan2(-deltaZ, deltaR);
          float dphi = asin(invPtToDphi * posStub.perp() * tp->charge() / tp->pt());
          float lambda = atan2(posStub.z(), posStub.perp());
          float dphiOverBend = pitchOverSep * (cos(tiltAngle) - sin(tiltAngle) * tan(lambda));
          float tpBend = dphi / dphiOverBend;  // Bend of TrackingParticle at stub position

          // +ve charges bend towards -ve phi, so flip sign.
          tpBend = -tpBend;

          // Some modules dont have increasing strip number in direction of increasing phi.
          const auto& localY = plane.rotation().y();
          const double radial = pos1.x() * localY.x() + pos1.y() * localY.y();
          bool signCol = std::signbit(barrel ? localY.z() : radial);
          // Calculation of variable d in StubGL.cc assumes localY points radially
          // inwards in +ve endcap, (where sinTilt=-1), so flip sign there.
          bool side = pos1.z() >= 0;
          if (side && not barrel)
            signCol = not signCol;
          bool signBend = signCol;

          if (signBend)
            tpBend = -tpBend;

          histogramLayer->second.bend[det]->Fill(ttStubRef->bendFE());
          histogramLayer->second.bendRes[det]->Fill(ttStubRef->bendFE() - tpBend);
        }
      }
    }
  }

  // fill the counter histos per layer
  for (unsigned int det = 1; det < 3; ++det) {
    for (auto it : nStubs[det]) {
      auto histogramLayer(histograms_.find(it.first));
      if (histogramLayer == histograms_.end())
        throw cms::Exception("LogicError") << "No histo for an existing counter! This should not happen!";
      histogramLayer->second.numberStubs[det]->Fill(it.second);
    }
    for (auto it : nStubMatch[det]) {
      auto histogramLayer(histograms_.find(it.first));
      if (histogramLayer == histograms_.end())
        throw cms::Exception("LogicError") << "No histo for an existing counter! This should not happen!";
      histogramLayer->second.numberStubMatch[det]->Fill(it.second);
    }
  }
}

// Create the histograms
std::map<unsigned int, StubHistos>::iterator AnalyzerStub::createLayerHistograms(unsigned int ival) {
  std::ostringstream fname1, fname2;

  edm::Service<TFileService> fs;
  fs->file().cd("/");

  std::string tag;
  unsigned int id;
  if (ival < 1000) {
    id = ival;
    fname1 << "Barrel";
    fname2 << "Layer_" << id;
    tag = "_layer_";
  } else {
    int side = ival / 1000;
    id = ival - side * 1000;
    if (ival > 10) {
      id /= 10;
      //        fname1 << "EndCap_Side_" << side;
      fname1 << "EndCap";
      fname2 << "Ring_" << id;
      tag = "_ring_";
    } else {
      id = ival;
      //        fname1 << "EndCap_Side_" << side;
      fname1 << "EndCap";
      fname2 << "Disk_" << id;
      tag = "_disk_";
    }
  }

  TFileDirectory td1 = fs->mkdir(fname1.str());
  TFileDirectory td = td1.mkdir(fname2.str());

  StubHistos local_histos;

  std::ostringstream histoName;

  /*
     * Number of stubs
     */

  histoName.str("");
  histoName << "Number_Stubs_PS" << tag.c_str() << id;
  local_histos.numberStubs[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 0.);
  histoName.str("");
  histoName << "Number_Stubs_2S" << tag.c_str() << id;
  local_histos.numberStubs[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 0.);

  /*
     * Number of stubs matching truth 
     */

  histoName.str("");
  histoName << "Number_Match_Stubs_PS" << tag.c_str() << id;
  local_histos.numberStubMatch[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 0.);
  histoName.str("");
  histoName << "Number_Match_Stubs_2S" << tag.c_str() << id;
  local_histos.numberStubMatch[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 0.);

  /*
     * Stub bend & bend resolution
     */

  histoName.str("");
  histoName << "Stub_Bend_PS" << tag.c_str() << id;
  local_histos.bend[1] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 200, -10., 10.);
  histoName.str("");
  histoName << "Stub_Bend_2S" << tag.c_str() << id;
  local_histos.bend[2] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 200, -10., 10.);

  histoName.str("");
  histoName << "Stub_Bend_Res_PS" << tag.c_str() << id;
  local_histos.bendRes[1] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 200, -20., 20.);
  histoName.str("");
  histoName << "Stub_Bend_Res_2S" << tag.c_str() << id;
  local_histos.bendRes[2] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 200, -20., 20.);

  /*
     * End
     */

  std::pair<std::map<unsigned int, StubHistos>::iterator, bool> insertedIt(
      histograms_.insert(std::make_pair(ival, local_histos)));
  fs->file().cd("/");

  return insertedIt.first;
}

DEFINE_FWK_MODULE(AnalyzerStub);
