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
#include <deque>
#include <set>
#include <map>
#include <functional>
#include <utility>

namespace tt {

  /*! \class  tt::Cleaner
   *  \brief  creates TrackingParticleCollection with at least one direct or indirect
   *          (via parent or child) associated cluster.
   *          creates TrackingVertexCollection with all vertices (with corrected
   *          parent TP connection) for those TPs.
   *          creates TTClusterAssociationMap with those TPs
   *  \author Thomas Schuh
   *  \date   2025, Aug
   */
  class Cleaner : public edm::stream::EDProducer<> {
  public:
    explicit Cleaner(const edm::ParameterSet&);
    ~Cleaner() override = default;

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;
    // ED input token of existing TTClusterAssociation
    edm::EDGetTokenT<TTClusterAssMap> edGetTokenAssoc_;
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
    edGetTokenAssoc_ = consumes<TTClusterAssMap>(inputTag);
    edPutTokenTPs_ = produces<std::vector<TrackingParticle>>(branch);
    edPutTokenTVs_ = produces<std::vector<TrackingVertex>>(branch);
    edPutTokenAssoc_ = produces<TTClusterAssMap>(branch);
  }

  void Cleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    const TTClusterAssMap& ttClusterAssMap = iEvent.get(edGetTokenAssoc_);
    const auto& refProdTPs = iEvent.getRefBeforePut(edPutTokenTPs_);
    const auto& refProdTVs = iEvent.getRefBeforePut(edPutTokenTVs_);
    // get primay TPs with at least one direct or indirect (via child) associated cluster
    std::set<TPPtr> tpPtrPrimaries;
    std::function<void(const TPPtr&)> addPrimaryTPPtr;
    addPrimaryTPPtr = [&addPrimaryTPPtr, &tpPtrPrimaries](const TPPtr& tpPtr) {
      const TVRef& tvRefParent = tpPtr->parentVertex();
      if (tvRefParent->nSourceTracks() > 0) {
        const TPRef& tpRefParent = *tvRefParent->sourceTracks_begin();
        const TPPtr tpPtrParent = edm::refToPtr(tpRefParent);
        addPrimaryTPPtr(tpPtrParent);
      }
      tpPtrPrimaries.insert(tpPtr);
    };
    // loop over association map
    for (const auto& p : ttClusterAssMap.getTrackingParticleToTTClustersMap()) {
      // ignore TPs without associated cluster
      if (p.second.empty())
        continue;
      // ignore TPs with only associated null cluster refs
      if (std::all_of(p.second.begin(), p.second.end(), [](const auto& c) { return c.isNull(); }))
        continue;
      // find and fill primary TP into tpPrimaries
      addPrimaryTPPtr(p.first);
    }
    // get all children of those primaries
    std::set<TPPtr> tpPtrChilds;
    std::function<void(const TPPtr&)> addChildTPPtr;
    addChildTPPtr = [&addChildTPPtr, &tpPtrChilds](const TPPtr& tpPtrParent) {
      if (tpPtrParent->decayVertices().empty())
        return;
      const TVRef& tvRef = *tpPtrParent->decayVertices_begin();
      for (auto it = tvRef->daughterTracks_begin(); it != tvRef->daughterTracks_end(); it++) {
        const TPPtr tpPtrChild = edm::refToPtr(*it);
        tpPtrChilds.insert(tpPtrChild);
        addChildTPPtr(tpPtrChild);
      }
    };
    for (const TPPtr& tpPtr : tpPtrPrimaries)
      addChildTPPtr(tpPtr);
    // identify primaries which are actually a child and remove them from primaries
    std::set<TPPtr> tpPtrFakePrimaries;
    std::set_intersection(tpPtrPrimaries.begin(),
                          tpPtrPrimaries.end(),
                          tpPtrChilds.begin(),
                          tpPtrChilds.end(),
                          std::inserter(tpPtrFakePrimaries, tpPtrFakePrimaries.begin()));
    for (const TPPtr& tpPtr : tpPtrFakePrimaries)
      tpPtrPrimaries.erase(tpPtr);
    // create TPs + TVs and associate TPs with TTCluster
    std::deque<TrackingParticle> tps;
    std::deque<TrackingVertex> tvs;
    std::map<TPPtr, std::vector<TTClusterRef>> mapTPPtrClusterRefs;
    std::map<TTClusterRef, std::vector<TPPtr>> mapTTClusterRefTPPtrs;
    std::function<TPRef(const TPPtr&, const TVRef&)> copy;
    copy =
        [&copy, &tps, &tvs, &mapTPPtrClusterRefs, &mapTTClusterRefTPPtrs, &ttClusterAssMap, &refProdTPs, &refProdTVs](
            const TPPtr& tpPtrOld, const TVRef& tvRefNewParent) {
          // constuct TP
          tps.emplace_back(tpPtrOld->g4Tracks().front(), tvRefNewParent);
          TrackingParticle& tpNew = tps.back();
          TPRef tpRefNew(refProdTPs, static_cast<int>(tps.size()) - 1);
          tpNew.setNumberOfHits(tpPtrOld->numberOfHits());
          tpNew.setNumberOfTrackerHits(tpPtrOld->numberOfTrackerHits());
          tpNew.setNumberOfTrackerLayers(tpPtrOld->numberOfTrackerLayers());
          // associate
          const TPPtr tpPtrNew = edm::refToPtr(tpRefNew);
          for (const auto& cluster : ttClusterAssMap.findTTClusterRefs(tpPtrOld)) {
            // ignore null cluster refs
            if (cluster.isNull())
              continue;
            mapTPPtrClusterRefs[tpPtrNew].push_back(cluster);
            mapTTClusterRefTPPtrs[cluster].push_back(tpPtrNew);
          }
          // stop recurse if no childs are present
          if (tpPtrOld->decayVertices().empty())
            return tpRefNew;
          // construct decay TV
          const TVRef& tvRefOld = *tpPtrOld->decayVertices_begin();
          tvs.emplace_back(tvRefOld->position(), tvRefOld->inVolume(), tvRefOld->eventId());
          TrackingVertex& tvNewDecay = tvs.back();
          TVRef tvRefNewDecay(refProdTVs, static_cast<int>(tvs.size()) - 1);
          tvNewDecay.addParentTrack(tpRefNew);
          tpNew.addDecayVertex(tvRefNewDecay);
          // recurse over childs
          for (auto it = tvRefOld->daughterTracks_begin(); it != tvRefOld->daughterTracks_end(); it++) {
            const TPPtr tpPtrOldChild = edm::refToPtr(*it);
            TPRef tpRefNewChild = copy(tpPtrOldChild, tvRefNewDecay);
            tvNewDecay.addDaughterTrack(tpRefNewChild);
          }
          return tpRefNew;
        };
    for (const TPPtr& tpPtrOld : tpPtrPrimaries) {
      // construct primary TV
      const TVRef& tvRefOld = tpPtrOld->parentVertex();
      tvs.emplace_back(tvRefOld->position(), tvRefOld->inVolume(), tvRefOld->eventId());
      TrackingVertex& tvNew = tvs.back();
      TVRef tvRefNew(refProdTVs, static_cast<int>(tvs.size()) - 1);
      // create new TP copies of primaries with children and associate those with TTCluster
      TPRef tpRefNew = copy(tpPtrOld, tvRefNew);
      tvNew.addDaughterTrack(tpRefNew);
    }
    // create TTClusterAssMap
    TTClusterAssMap assoc;
    assoc.setTTClusterToTrackingParticlesMap(mapTTClusterRefTPPtrs);
    assoc.setTrackingParticleToTTClustersMap(mapTPPtrClusterRefs);
    // store products
    iEvent.emplace(edPutTokenTPs_, tps.begin(), tps.end());
    iEvent.emplace(edPutTokenTVs_, tvs.begin(), tvs.end());
    iEvent.emplace(edPutTokenAssoc_, std::move(assoc));
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::Cleaner);
