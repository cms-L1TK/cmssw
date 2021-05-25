#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimTracker/TrackTriggerAssociation/interface/TTTypes.h"
#include "SimTracker/TrackTriggerAssociation/interface/StubAssociation.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"

#include <vector>
#include <map>
#include <utility>
#include <set>
#include <algorithm>
#include <iterator>
#include <cmath>

using namespace std;
using namespace edm;

namespace tt {

  /*! \class  ttAssociation::TTTrackingParticleAssociator
   *  \brief  Class to associate reconstrucable TrackingParticles with TTStubs
   *  \author Thomas Schuh
   *  \date   2020, Apr
   */
  class StubAssociator : public stream::EDProducer<> {
  public:
    explicit StubAssociator(const ParameterSet&);
    ~StubAssociator() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override;
    void produce(Event&, const EventSetup&) override;
    void endJob() {}

    //
    void fill(set<TPPtr>& tpPtrs, const TTStubRef& ttStubRef, const Handle<TTClusterAssMap>& handle) const;
    //
    void fill(set<TPPtr>& tpPtrs, const TTStubRef& ttStubRef, const Handle<TTStubAssMap>& handle) const;

    // helper classe to store configurations
    const Setup* setup_;
    // ED input token of TTStubs
    EDGetTokenT<TTStubDetSetVec> getTokenTTStubDetSetVec_;
    // ED input token of TTClusterAssociation
    EDGetTokenT<TTClusterAssMap> getTokenTTClusterAssMap_;
    // ED input token of TTStubAssociation
    EDGetTokenT<TTStubAssMap> getTokenTTStubAssMap_;
    // ED output token for recosntructable stub association
    EDPutTokenT<StubAssociation> putTokenReconstructable_;
    // ED output token for selected stub association
    EDPutTokenT<StubAssociation> putTokenSelection_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetToken_;
    //
    bool useTTStubAssMap_;
  };

  StubAssociator::StubAssociator(const ParameterSet& iConfig) :
    useTTStubAssMap_(iConfig.getParameter<bool>("UseTTStubAssMap"))
  {
    // book in- and output ed products
    getTokenTTStubDetSetVec_ = consumes<TTStubDetSetVec>(iConfig.getParameter<InputTag>("InputTagTTStubDetSetVec"));
    if (useTTStubAssMap_)
      getTokenTTStubAssMap_ = consumes<TTStubAssMap>(iConfig.getParameter<InputTag>("InputTagTTStubAssMap"));
    else
      getTokenTTClusterAssMap_ = consumes<TTClusterAssMap>(iConfig.getParameter<InputTag>("InputTagTTClusterAssMap"));
    putTokenReconstructable_ = produces<StubAssociation>(iConfig.getParameter<string>("BranchReconstructable"));
    putTokenSelection_ = produces<StubAssociation>(iConfig.getParameter<string>("BranchSelection"));
    // book ES product
    esGetToken_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
  }

  void StubAssociator::beginRun(const Run& iRun, const EventSetup& iSetup) {
    setup_ = &iSetup.getData(esGetToken_);
  }

  void StubAssociator::produce(Event& iEvent, const EventSetup& iSetup) {
    // associate TTStubs with TrackingParticles
    Handle<TTStubDetSetVec> handleTTStubDetSetVec;
    iEvent.getByToken<TTStubDetSetVec>(getTokenTTStubDetSetVec_, handleTTStubDetSetVec);
    Handle<TTClusterAssMap> handleTTClusterAssMap;
    Handle<TTStubAssMap> handleTTStubAssMap;
    if (useTTStubAssMap_)
      iEvent.getByToken<TTStubAssMap>(getTokenTTStubAssMap_, handleTTStubAssMap);
    else
      iEvent.getByToken<TTClusterAssMap>(getTokenTTClusterAssMap_, handleTTClusterAssMap);
    map<TPPtr, vector<TTStubRef>> mapTPPtrsTTStubRefs;
    for (TTStubDetSetVec::const_iterator ttModule = handleTTStubDetSetVec->begin();
        ttModule != handleTTStubDetSetVec->end();
        ttModule++) {
      for (TTStubDetSet::const_iterator ttStub = ttModule->begin(); ttStub != ttModule->end(); ttStub++) {
        const TTStubRef ttStubRef = makeRefTo(handleTTStubDetSetVec, ttStub);
        set<TPPtr> tpPtrs;
        if (useTTStubAssMap_)
          fill(tpPtrs, ttStubRef, handleTTStubAssMap);
        else
          fill(tpPtrs, ttStubRef, handleTTClusterAssMap);
        for (const TPPtr& tpPtr : tpPtrs)
          mapTPPtrsTTStubRefs[tpPtr].push_back(ttStubRef);
      }
    }
    // associate reconstructable TrackingParticles with TTStubs
    StubAssociation reconstructable(setup_);
    StubAssociation selection(setup_);
    for (const pair<TPPtr, vector<TTStubRef>>& p : mapTPPtrsTTStubRefs) {
      if (!setup_->useForReconstructable(*p.first) || !setup_->reconstructable(p.second))
        continue;
      reconstructable.insert(p.first, p.second);
      if (setup_->useForAlgEff(*p.first))
        selection.insert(p.first, p.second);
    }
    iEvent.emplace(putTokenReconstructable_, move(reconstructable));
    iEvent.emplace(putTokenSelection_, move(selection));
  }

  //
  void StubAssociator::fill(set<TPPtr>& tpPtrs, const TTStubRef& ttStubRef, const Handle<TTClusterAssMap>& handle) const {
    auto isNonnull = [](const TPPtr& tpPtr){ return tpPtr.isNonnull(); };
    for (unsigned int iClus = 0; iClus < 2; iClus++) {
      const vector<TPPtr>& assocPtrs = handle->findTrackingParticlePtrs(ttStubRef->clusterRef(iClus));
      copy_if(assocPtrs.begin(), assocPtrs.end(), inserter(tpPtrs, tpPtrs.begin()), isNonnull);
    }
  }

  //
  void StubAssociator::fill(set<TPPtr>& tpPtrs, const TTStubRef& ttStubRef, const Handle<TTStubAssMap>& handle) const {
    const TPPtr tpPtr = handle->findTrackingParticlePtr(ttStubRef);
    if (tpPtr.isNonnull())
      tpPtrs.insert(tpPtr);
  }

} // namespace tt

DEFINE_FWK_MODULE(tt::StubAssociator);