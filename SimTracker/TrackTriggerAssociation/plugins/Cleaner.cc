#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/StubAssociation.h"

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <utility>

namespace tt {

  /*! \class  tt::Cleaner
   *  \brief  creates TrackingParticleCollection with at least one associated cluster in the TP chain.
   *          adds missing decay vertices to those TPs.
   *          creates TrackingVertexCollection with all vertices for those TPs.
              creates TTClusterAssociationMap with those TPs
   *  \author Thomas Schuh
   *  \date   2025, Aug
   */
  class Cleaner : public edm::stream::EDProducer<> {
  public:
    explicit Cleaner(const edm::ParameterSet&);
    ~Cleaner() override = default;

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;
    // returns start of TP chain
    TPPtr getPrimaryTP(const TPPtr&) const;
    // fills set with TP chain downstream starting from given TPPtr
    void fillTPChainDown(std::set<TPPtr>&, const TPPtr&) const;
    // fills set with TP chain upstream starting from given TPPtr
    void fillTPChainUp(std::set<TPPtr>&, const TPPtr&) const;
    // ED input token of existing TTClusterAssociation
    edm::EDGetTokenT<TTClusterAssMap> edGetToken_;
    // ED output token of clean TPs
    edm::EDPutTokenT<std::vector<TrackingParticle>> edPutTokenTPs_;
    // ED output token of clean TVs
    edm::EDPutTokenT<std::vector<TrackingVertex>> edPutTokenTVs_;
    // ED output token of clean association
    edm::EDPutTokenT<TTClusterAssMap> edPutTokenAssoc_;
  };

  Cleaner::Cleaner(const edm::ParameterSet& iConfig) {
    const edm::InputTag& inputTag = iConfig.getParameter<edm::InputTag>("InputTag");
    const std::string& branch = iConfig.getParameter<std::string>("Branch");
    // book in- and output ed products
    edGetToken_ = consumes<TTClusterAssMap>(inputTag);
    edPutTokenTPs_ = produces<std::vector<TrackingParticle>>(branch);
    edPutTokenTVs_ = produces<std::vector<TrackingVertex>>(branch);
    edPutTokenAssoc_ = produces<TTClusterAssMap>(branch);
  }

  // returns start of TP chain
  TPPtr Cleaner::getPrimaryTP(const TPPtr& tpPtr) const {
    const TVRef& tvRefParent = tpPtr->parentVertex();
    if (tvRefParent->nSourceTracks() > 0) {
      const TPRef& tpRefParent = *tvRefParent->sourceTracks_begin();
      const TPPtr tpPtrParent = edm::refToPtr(tpRefParent);
      return this->getPrimaryTP(tpPtrParent);
    }
    return tpPtr;
  }

  // fills set with TP chain downstream starting from given TPPtr
  void Cleaner::fillTPChainDown(std::set<TPPtr>& tpPtrs, const TPPtr& tpPtr) const {
    tpPtrs.insert(tpPtr);
    if (!tpPtr->decayVertices().empty()) {
      const TVRef& tvRefChild = *tpPtr->decayVertices_begin();
      const TPRef& tpRefChild = *tvRefChild->daughterTracks_begin();
      const TPPtr tpPtrChild = edm::refToPtr(tpRefChild);
      this->fillTPChainDown(tpPtrs, tpPtrChild);
    }
  }

  // fills set with TP chain upstream starting from given TPPtr
  void Cleaner::fillTPChainUp(std::set<TPPtr>& tpPtrs, const TPPtr& tpPtr) const {
    tpPtrs.insert(tpPtr);
    const TVRef& tvRefParent = tpPtr->parentVertex();
    if (tvRefParent->nSourceTracks() > 0) {
      const TPRef& tpRefParent = *tvRefParent->sourceTracks_begin();
      const TPPtr tpPtrParent = edm::refToPtr(tpRefParent);
      this->fillTPChainUp(tpPtrs, tpPtrParent);
    }
  }

  void Cleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    const TTClusterAssMap& ttClusterAssMap = iEvent.get(edGetToken_);
    const auto& refProdTPs = iEvent.getRefBeforePut(edPutTokenTPs_);
    const auto& refProdTVs = iEvent.getRefBeforePut(edPutTokenTVs_);
    // get primay TPs & childs with at least one associated cluster
    std::map<TPPtr, std::set<TPPtr>> primaries;
    for (const auto& p : ttClusterAssMap.getTrackingParticleToTTClustersMap()) {
      if (p.second.empty())
        continue;
      if (std::all_of(p.second.begin(), p.second.end(), [](const auto& c) { return c.isNull(); }))
        continue;
      const TPPtr primary = this->getPrimaryTP(p.first);
      primaries[primary].insert(p.first);
    }
    // fill childs without associated cluster
    std::vector<std::vector<TPPtr>> tpChains;
    tpChains.reserve(primaries.size());
    for (const auto& p : primaries) {
      std::set<TPPtr> chain;
      this->fillTPChainDown(chain, p.first);
      for (const TPPtr& tpPtr : p.second) {
        this->fillTPChainDown(chain, tpPtr);
        this->fillTPChainUp(chain, tpPtr);
      }
      tpChains.emplace_back(chain.begin(), chain.end());
    }
    // count TPs and TVs
    auto acc = [](int sum, const auto& tpPtrs) { return sum + tpPtrs.size(); };
    const int size = std::accumulate(tpChains.begin(), tpChains.end(), 0, acc);
    // create TPs with TVs and associate TPs with TTCluster
    std::vector<TrackingParticle> tps;
    tps.reserve(size);
    std::vector<TrackingVertex> tvs;
    tvs.reserve(size);
    std::map<TPPtr, std::vector<TTClusterRef>> mapTPPtrClusterRefs;
    std::map<TTClusterRef, std::vector<TPPtr>> mapTTClusterRefTPPtrs;
    int index(-1);
    for (const std::vector<TPPtr>& chain : tpChains) {
      for (const TPPtr& tpPtr : chain) {
        index++;
        // construct TV
        const TVRef& tvRef = tpPtr->parentVertex();
        tvs.emplace_back(tvRef->position(), tvRef->inVolume(), tvRef->eventId());
        TrackingVertex& tv = tvs.back();
        tv.addDaughterTrack(TPRef(refProdTPs, index));
        if (tpPtr != chain.front())
          tv.addParentTrack(TPRef(refProdTPs, index - 1));
        // construct TP
        tps.emplace_back(tpPtr->g4Tracks().front(), TVRef(refProdTVs, index));
        TrackingParticle& tp = tps.back();
        tp.setNumberOfHits(tpPtr->numberOfHits());
        tp.setNumberOfTrackerHits(tpPtr->numberOfTrackerHits());
        tp.setNumberOfTrackerLayers(tpPtr->numberOfTrackerLayers());
        if (tpPtr != chain.back())
          tp.addDecayVertex(TVRef(refProdTVs, index + 1));
        // associate
        const TPPtr clean(edm::refToPtr(*tv.daughterTracks_begin()));
        for (const auto& cluster : ttClusterAssMap.findTTClusterRefs(tpPtr)) {
          if (cluster.isNull())
            continue;
          mapTPPtrClusterRefs[clean].push_back(cluster);
          mapTTClusterRefTPPtrs[cluster].push_back(clean);
        }
      }
    }
    // create TTClusterAssMap
    TTClusterAssMap assoc;
    assoc.setTTClusterToTrackingParticlesMap(mapTTClusterRefTPPtrs);
    assoc.setTrackingParticleToTTClustersMap(mapTPPtrClusterRefs);
    // store products
    iEvent.emplace(edPutTokenTPs_, std::move(tps));
    iEvent.emplace(edPutTokenTVs_, std::move(tvs));
    iEvent.emplace(edPutTokenAssoc_, std::move(assoc));
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::Cleaner);
