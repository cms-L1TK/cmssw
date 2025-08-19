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
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/TTTypes.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackTrigger/interface/Associator.h"
#include "SimDataFormats/Associations/interface/StubAssociation.h"

#include <vector>
#include <deque>
#include <map>
#include <utility>
#include <set>
#include <iterator>

namespace tt {

  /*! \class  tt::StubAssociator
   *  \brief  Class to associate reconstrucable TrackingParticles with TTStubs and vice versa
   *          It may associate multiple TPs with a TTStub and can therefore be used to associate
   *          TTTracks with TrackingParticles. This EDProducer creates two StubAssociation EDProducts,
   *          one using only TP that are "reconstructable" (produce stubs in a min. number of layers)
   *          and one using TP that are also "use for the tracking efficiency measurement".
   *  \author Thomas Schuh
   *  \date   2020, Apr
   */
  class StubAssociator : public edm::stream::EDProducer<> {
  public:
    explicit StubAssociator(const edm::ParameterSet&);
    ~StubAssociator() override = default;

  private:
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void produce(edm::Event&, const edm::EventSetup&) override;
    // helper class to store configurations
    const Setup* setup_;
    // helper class to associate TTStubs and TrackingParticle
    const Associator* associator_;
    // ED input token of TTStubs
    edm::EDGetTokenT<TTStubDetSetVec> getTokenTTStubDetSetVec_;
    // ED input token of TTClusterAssociation
    edm::EDGetTokenT<TTClusterAssMap> getTokenTTClusterAssMap_;
    // ED output token for recosntructable stub association
    edm::EDPutTokenT<StubAssociation> putTokenReconstructable_;
    // ED output token for selected stub association
    edm::EDPutTokenT<StubAssociation> putTokenSelection_;
    // Setup token
    edm::ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // Associator token
    edm::ESGetToken<Associator, SetupRcd> esGetTokenAssociator_;
    // pt cut in GeV
    double minPt_;
    // max eta for TP with z0 = 0
    double maxEta0_;
    // half lumi region size in cm
    double maxZ0_;
    // cut on impact parameter in cm
    double maxD0_;
    // cut on vertex pos r in cm
    double maxVertR_;
    // cut on vertex pos z in cm
    double maxVertZ_;
    // cut on TP zT
    double maxZT_;
    // selector to partly select TPs for efficiency measurements
    TrackingParticleSelector tpSelector_;
  };

  StubAssociator::StubAssociator(const edm::ParameterSet& iConfig)
      : minPt_(iConfig.getParameter<double>("MinPt")),
        maxEta0_(iConfig.getParameter<double>("MaxEta0")),
        maxZ0_(iConfig.getParameter<double>("MaxZ0")),
        maxD0_(iConfig.getParameter<double>("MaxD0")),
        maxVertR_(iConfig.getParameter<double>("MaxVertR")),
        maxVertZ_(iConfig.getParameter<double>("MaxVertZ")) {
    // book in- and output ed products
    const auto& ttStubDetSetVec = iConfig.getParameter<edm::InputTag>("InputTagTTStubDetSetVec");
    const auto& ttClusterAssMap = iConfig.getParameter<edm::InputTag>("InputTagTTClusterAssMap");
    const auto& reconstructable = iConfig.getParameter<std::string>("BranchReconstructable");
    const auto& selection = iConfig.getParameter<std::string>("BranchSelection");
    getTokenTTStubDetSetVec_ = consumes<TTStubDetSetVec>(ttStubDetSetVec);
    getTokenTTClusterAssMap_ = consumes<TTClusterAssMap>(ttClusterAssMap);
    putTokenReconstructable_ = produces<StubAssociation>(reconstructable);
    putTokenSelection_ = produces<StubAssociation>(selection);
    // book ES product
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, edm::Transition::BeginRun>();
    esGetTokenAssociator_ = esConsumes<Associator, SetupRcd, edm::Transition::BeginRun>();
  }

  void StubAssociator::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    setup_ = &iSetup.getData(esGetTokenSetup_);
    associator_ = &iSetup.getData(esGetTokenAssociator_);
    maxZT_ = std::sinh(maxEta0_) * setup_->chosenRofZ();
    // configure TrackingParticleSelector
    constexpr double ptMax = 9.e9;
    constexpr int minHit = 0;
    constexpr bool signalOnly = true;
    constexpr bool intimeOnly = true;
    constexpr bool chargedOnly = true;
    constexpr bool stableOnly = false;
    const double maxEta = std::asinh((maxZT_ + maxZ0_) / setup_->chosenRofZ());
    tpSelector_ = TrackingParticleSelector(
        minPt_, ptMax, -maxEta, maxEta, maxVertR_, maxVertZ_, minHit, signalOnly, intimeOnly, chargedOnly, stableOnly);
  }

  void StubAssociator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // associate TTStubs with TrackingParticles
    edm::Handle<TTStubDetSetVec> handle;
    iEvent.getByToken<TTStubDetSetVec>(getTokenTTStubDetSetVec_, handle);
    const TTClusterAssMap& ttClusterAssMap = iEvent.get(getTokenTTClusterAssMap_);
    std::map<TPPtr, std::set<TTStubRef>> mapTPPtrsTTStubRefs;
    for (TTStubDetSetVec::const_iterator ttModule = handle->begin(); ttModule != handle->end(); ttModule++) {
      for (TTStubDetSet::const_iterator ttStub = ttModule->begin(); ttStub != ttModule->end(); ttStub++) {
        const TTStubRef ttStubRef = makeRefTo(handle, ttStub);
        std::set<TPPtr> tpPtrs;
        for (unsigned int iClus = 0; iClus < 2; iClus++)
          for (const TPPtr& tpPtr : ttClusterAssMap.findTrackingParticlePtrs(ttStubRef->clusterRef(iClus)))
            if (tpPtr.isNonnull())
              tpPtrs.insert(tpPtr);
        for (const TPPtr& tpPtr : tpPtrs)
          mapTPPtrsTTStubRefs[tpPtr].insert(ttStubRef);
      }
    }
    // associate TTStubs with primary TrackingParticles
    std::map<TPPtr, std::set<TTStubRef>> mapPrimaryTPPtrsTTStubRefs;
    for (auto& p : mapTPPtrsTTStubRefs) {
      const TPPtr primary = associator_->getPrimaryTP(p.first);
      std::deque<TPPtr> chain;
      associator_->fillTPChain(chain, primary);
      std::set<TTStubRef>& ttStubRefs = mapPrimaryTPPtrsTTStubRefs[primary];
      ttStubRefs.insert(p.second.begin(), p.second.end());
    }
    // associate reconstructable TrackingParticles with TTStubs
    StubAssociation reconstructable;
    StubAssociation selection;
    // fills matched TPPtr into deque, empty if only accumulated stubs matched and not a single TP
    auto fill = [this, &mapTPPtrsTTStubRefs](std::deque<TPPtr>& tpPtrs, const TPPtr& primary) {
      std::deque<TPPtr> chain;
      associator_->fillTPChain(chain, primary);
      for (const TPPtr& tpPtr : chain) {
        if (!mapTPPtrsTTStubRefs.contains(tpPtr))
          continue;
        const std::set<TTStubRef>& stubs = mapTPPtrsTTStubRefs.at(tpPtr);
        if (associator_->reconstructable({stubs.begin(), stubs.end()}))
          tpPtrs.push_back(tpPtr);
      }
    };
    for (const auto& p : mapPrimaryTPPtrsTTStubRefs) {
      const TPPtr& primary = p.first;
      const std::vector<TTStubRef> ttStubRefs(p.second.begin(), p.second.end());
      std::set<int> layers;
      for (const TTStubRef& ttStubRef : ttStubRefs)
        layers.insert(setup_->layerId(ttStubRef));
      // require min layers
      if (!associator_->reconstructable(ttStubRefs))
        continue;
      reconstructable.insert(primary, ttStubRefs);
      std::deque<TPPtr> matched;
      fill(matched, primary);
      if (matched.empty())
        continue;
      // require parameter space and signal only
      bool valid(false);
      for (const TPPtr& tpPtr : matched) {
        if (!tpSelector_(*tpPtr))
          continue;
        const double zT = tpPtr->z0() + tpPtr->tanl() * setup_->chosenRofZ();
        if ((std::abs(tpPtr->d0()) > maxD0_) || (std::abs(tpPtr->z0()) > maxZ0_) || (std::abs(zT) > maxZT_))
          continue;
        valid = true;
      }
      if (valid)
        selection.insert(primary, ttStubRefs);
    }
    iEvent.emplace(putTokenReconstructable_, std::move(reconstructable));
    iEvent.emplace(putTokenSelection_, std::move(selection));
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::StubAssociator);
