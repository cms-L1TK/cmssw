#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTypes.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include <string>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;
using namespace edm;

namespace tt {

  /*! \class  tt::CleanTP
   *  \brief  creates a new TPs (with only empty parent vertex, no childs ot TPs stored) with at least one cluster in outer Tracker associated
   *  \author Thomas Schuh
   *  \date   2022, June
   */
  class CleanTP : public stream::EDProducer<> {
  public:
    explicit CleanTP(const ParameterSet&);
    ~CleanTP() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override {}
    void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}

    // ED input token of TTClusterAssociation
    EDGetTokenT<TTClusterAssMap> edGetTokenAss_;
    EDGetTokenT<TrackingVertexCollection> edGetTokenV_;
    // ED output token for accepted stubs
    EDPutTokenT<TrackingParticleCollection> edPutToken_;
  };

  CleanTP::CleanTP(const ParameterSet& iConfig) {
    // book in- and output ed products
    edGetTokenAss_ = consumes<TTClusterAssMap>(iConfig.getParameter<InputTag>("InputTagTTClusterAssMap"));
    edGetTokenV_ = consumes<TrackingVertexCollection>(iConfig.getParameter<InputTag>("InputTagTVs"));
    edPutToken_ = produces<TrackingParticleCollection>(iConfig.getParameter<string>("Branch"));
  }

  void CleanTP::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty output product
    TrackingParticleCollection tps;
    // find TrackingParticle with at least one associated cluster
    Handle<TTClusterAssMap> handleAss;
    iEvent.getByToken<TTClusterAssMap>(edGetTokenAss_, handleAss);
    Handle<TrackingVertexCollection> handleV;
    iEvent.getByToken<TrackingVertexCollection>(edGetTokenV_, handleV);
    const map<TPPtr, vector<TTClusterRef>>& m = handleAss->getTrackingParticleToTTClustersMap();
    tps.reserve(m.size());
    int i(0);
    for (const auto& p : m) {
      const TPPtr& tpPtr = p.first;
      const Ref<TrackingVertexCollection> tvRef(handleV, i++);
      tps.emplace_back(tpPtr->g4Tracks().front(), tvRef);
      tps.back().setNumberOfHits(tpPtr->numberOfHits());
      tps.back().setNumberOfTrackerHits(tpPtr->numberOfTrackerHits());
      tps.back().setNumberOfTrackerLayers(tpPtr->numberOfTrackerLayers());
    }
    // store products
    iEvent.emplace(edPutToken_, move(tps));
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::CleanTP);