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

#include "SimTracker/TrackTriggerAssociation/interface/TTTypes.h"

#include <string>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <algorithm>

using namespace std;
using namespace edm;

namespace tt {

  /*! \class  tt::CleanAssoc
   *  \brief  creates TTClusterAssociationMap with clean TPs
   *  \author Thomas Schuh
   *  \date   2022, June
   */
  class CleanAssoc : public stream::EDProducer<> {
  public:
    explicit CleanAssoc(const ParameterSet&);
    ~CleanAssoc() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override {}
    void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}

    // ED input token of TTClusterAssociation
    EDGetTokenT<TTClusterAssMap> edGetTokenTTClusterAssMap_;
    // ED input token of clean TPs
    EDGetTokenT<TrackingParticleCollection> edGetTokenTrackingParticleCollection_;
    // ED output token
    EDPutTokenT<TTClusterAssMap> edPutToken_;
  };

  CleanAssoc::CleanAssoc(const ParameterSet& iConfig) {
    // book in- and output ed products
    edGetTokenTTClusterAssMap_ = consumes<TTClusterAssMap>(iConfig.getParameter<InputTag>("InputTagTTClusterAssMap"));
    edGetTokenTrackingParticleCollection_ = consumes<TrackingParticleCollection>(iConfig.getParameter<InputTag>("InputTagTPs"));
    edPutToken_ = produces<TTClusterAssMap>(iConfig.getParameter<string>("Branch"));
  }

  void CleanAssoc::produce(Event& iEvent, const EventSetup& iSetup) {
    // recreate TTClusterAssMap  with clean TPs
    map<TTClusterRef, vector<TPPtr>> mapTTClusterRefTPPtrs;
    map<TPPtr, vector<TTClusterRef>> mapTPPtrClusterRefs;
    Handle<TTClusterAssMap> handleTTClusterAssMap;
    iEvent.getByToken<TTClusterAssMap>(edGetTokenTTClusterAssMap_, handleTTClusterAssMap);
    Handle<TrackingParticleCollection> handleTrackingParticleCollection;
    iEvent.getByToken<TrackingParticleCollection>(edGetTokenTrackingParticleCollection_, handleTrackingParticleCollection);
    int i(0);
    for (const auto& p : handleTTClusterAssMap->getTrackingParticleToTTClustersMap()) {
      const TPPtr tpPtr(handleTrackingParticleCollection, i++);
      mapTPPtrClusterRefs.emplace(tpPtr, p.second);
      for (const auto c : p.second)
        mapTTClusterRefTPPtrs[c].push_back(tpPtr);
    }
    // create and store TTClusterAssMap
    TTClusterAssMap ttClusterAssMap;
    ttClusterAssMap.setTTClusterToTrackingParticlesMap(mapTTClusterRefTPPtrs);
    ttClusterAssMap.setTrackingParticleToTTClustersMap(mapTPPtrClusterRefs);
    iEvent.emplace(edPutToken_, move(ttClusterAssMap));
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::CleanAssoc);