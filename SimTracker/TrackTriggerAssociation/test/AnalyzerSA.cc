#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimTracker/TrackTriggerAssociation/interface/StubAssociation.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTypes.h"

#include <TH1F.h>

using namespace std;
using namespace edm;

namespace tt {

  /*! \class  tt::AnalyzerSA
   *  \brief  Class to analyze Tracking particles used by Stub Associator
   *  \author Thomas Schuh
   *  \date   2024, Sep
   */
  class AnalyzerSA : public one::EDAnalyzer<one::WatchRuns, one::SharedResources> {
  public:
    AnalyzerSA(const ParameterSet& iConfig);
    void beginJob() override {}
    void beginRun(const Run& iEvent, const EventSetup& iSetup) override;
    void analyze(const Event& iEvent, const EventSetup& iSetup) override;
    void endRun(const Run& iEvent, const EventSetup& iSetup) override {}
    void endJob() override {}

  private:
    // ED input token
    EDGetTokenT<TrackingParticleCollection> edGetToken_;
    // Histograms
    TH1F* hisMCZ0_;
    TH1F* hisMCD0_;
    TH1F* hisMCVz_;
    TH1F* hisMCVr_;
  };

  AnalyzerSA::AnalyzerSA(const ParameterSet& iConfig) {
    usesResource("TFileService");
    // book input ED product
    edGetToken_ = consumes<TrackingParticleCollection>(InputTag("CleanTP", "AtLeastOneCluster"));
  }

  void AnalyzerSA::beginRun(const Run& iEvent, const EventSetup& iSetup) {
    // book histograms
    Service<TFileService> fs;
    TFileDirectory dir;
    dir = fs->mkdir("MC");
    hisMCZ0_ = dir.make<TH1F>("His TP z0", ";", 600, -300, 300);
    hisMCD0_ = dir.make<TH1F>("His TP d0", ";", 100, -50, 50);
    hisMCVz_ = dir.make<TH1F>("His TV z", ";", 600, -300, 300);
    hisMCVr_ = dir.make<TH1F>("His TV r", ";", 120, 0, 120);
  }

  void AnalyzerSA::analyze(const Event& iEvent, const EventSetup& iSetup) {
    // get TPs
    Handle<TrackingParticleCollection> handle;
    iEvent.getByToken<TrackingParticleCollection>(edGetToken_, handle);
    // analyze TPs
    for (const TrackingParticle& tp : *handle) {
      hisMCZ0_->Fill(tp.z0());
      hisMCD0_->Fill(tp.d0());
      hisMCVz_->Fill(tp.vz());
      hisMCVr_->Fill(sqrt(tp.vx() * tp.vx() + tp.vy() * tp.vy()));
    }
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::AnalyzerSA);
