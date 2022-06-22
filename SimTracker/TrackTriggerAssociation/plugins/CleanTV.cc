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

using namespace std;
using namespace edm;

namespace tt {

  /*! \class  tt::CleanTV
   *  \brief  creates a copy of parent TrackingVertex of TPs with at least one cluster in outer Tracker associated
   *  \author Thomas Schuh
   *  \date   2022, June
   */
  class CleanTV : public stream::EDProducer<> {
  public:
    explicit CleanTV(const ParameterSet&);
    ~CleanTV() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override {}
    void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}

    // ED input token of TTClusterAssociation
    EDGetTokenT<TTClusterAssMap> edGetToken_;
    // ED output token for accepted stubs
    EDPutTokenT<TrackingVertexCollection> edPutToken_;
  };

  CleanTV::CleanTV(const ParameterSet& iConfig) {
    // book in- and output ed products
    edGetToken_ = consumes<TTClusterAssMap>(iConfig.getParameter<InputTag>("InputTagTTClusterAssMap"));
    edPutToken_ = produces<TrackingVertexCollection>(iConfig.getParameter<string>("Branch"));
  }

  void CleanTV::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty output product
    TrackingVertexCollection tvs;
    // find parent TrackingVertex with at least one associated cluster
    Handle<TTClusterAssMap> handle;
    iEvent.getByToken<TTClusterAssMap>(edGetToken_, handle);
    const map<TPPtr, vector<TTClusterRef>>& m = handle->getTrackingParticleToTTClustersMap();
    tvs.reserve(m.size());
    for (const auto& p : m) {
      const TrackingVertexRef& tvRef = p.first->parentVertex();
      tvs.emplace_back(tvRef->position(), tvRef->inVolume(), tvRef->eventId());
    }
    // store products
    iEvent.emplace(edPutToken_, move(tvs));
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::CleanTV);