#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/Error.h"

#include "L1Trigger/TrackFindingTracklet/interface/Setup.h"
#include "L1Trigger/TrackTrigger/interface/StubPtConsistency.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <numeric>

namespace trklet {

  /*! \class  trklet::ProducerSim
   *  \brief  simulation of Track Processing for prompt or extended track finding
   *  \author Thomas Schuh
   *  \date   2026, June
   */
  class ProducerSim : public edm::stream::EDProducer<edm::stream::WatchRuns> {
  public:
    explicit ProducerSim(const edm::ParameterSet&);
    ~ProducerSim() override = default;

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;
    void beginRun(const edm::Run&, const edm::EventSetup&) override;

    struct Stub {
      double H_;
      double m0_;
      double m1_;
      double v0_;
      double v1_;
    };
    // ED input token of TTTracks
    edm::EDGetTokenT<tt::TTTracks> edGetTokenTracks_;
    // ED output token of TTTracks
    edm::EDPutTokenT<tt::TTTracks> edPutTokenTracks_;
    // Setup token
    edm::ESGetToken<Setup, trackerDTC::SetupRcd> esGetTokenSetup_;
    // helper class to store configurations
    const Setup* setup_;
    //
    std::vector<int> nPer_;
  };

  ProducerSim::ProducerSim(const edm::ParameterSet& iConfig) {
    const edm::InputTag& inputTag = iConfig.getParameter<edm::InputTag>("InputTagTracklet");
    const std::string& branchTracks = iConfig.getParameter<std::string>("BranchTTTracks");
    // book in- and output ED products
    edGetTokenTracks_ = consumes(inputTag);
    edPutTokenTracks_ = produces(branchTracks);
    // book ES products
    esGetTokenSetup_ = esConsumes<edm::Transition::BeginRun>();
  }

  void ProducerSim::beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    // calc permutations for all found track sizes [4 - 7]
    auto fac = [](int n) {
      int f(1);
      for (int i = 1; i <= n; i++)
        f *= i;
      return f;
    };
    auto bc = [fac](int n, int k) { return fac(n) / fac(k) / fac(n - k); };
    nPer_ = std::vector<int>(setup_->kfNumLayers() - setup_->kfMinLayers() + 1, 0);
    for (int i = setup_->kfMinLayers(); i <= setup_->kfNumLayers(); i++)
      for (int j = setup_->kfMinLayers(); j <= i; j++)
        nPer_[i - setup_->kfMinLayers()] += bc(i, j);
  }

  void ProducerSim::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // read input
    edm::Handle<tt::TTTracks> handle;
    iEvent.getByToken(edGetTokenTracks_, handle);
    std::vector<TTTrackRef> ttTrackRefs;
    ttTrackRefs.reserve(handle->size());
    for (int iTrk = 0; iTrk < static_cast<int>(handle->size()); iTrk++)
      ttTrackRefs.emplace_back(handle, iTrk);
    // perform track multiplexinf
    const std::vector<int>& muxOrder = setup_->tmMuxOrder();
    auto order = [&muxOrder](const TTTrackRef& lhs, TTTrackRef& rhs) {
      const auto l = std::find(muxOrder.begin(), muxOrder.end(), lhs->trackSeedType());
      const auto r = std::find(muxOrder.begin(), muxOrder.end(), rhs->trackSeedType());
      return l < r;
    };
    std::sort(ttTrackRefs.begin(), ttTrackRefs.end(), order);
    // perform duplicate removal
    auto equalEnough = [this](const TTTrackRef& lhs, const TTTrackRef& rhs) {
      std::vector<TTStubRef> l = lhs->getStubRefs();
      std::vector<TTStubRef> r = rhs->getStubRefs();
      std::sort(l.begin(), l.end());
      std::sort(r.begin(), r.end());
      std::vector<TTStubRef> same;
      same.reserve(std::min(l.size(), r.size()));
      std::set_intersection(l.begin(), l.end(), r.begin(), r.end(), std::back_inserter(same));
      return static_cast<int>(same.size()) >= setup_->drMinIdenticalStubs();
    };
    std::vector<TTTrackRef*> ptrs;
    ptrs.reserve(ttTrackRefs.size());
    auto toPtr = [](TTTrackRef& ref) { return &ref; };
    std::transform(ttTrackRefs.begin(), ttTrackRefs.end(), std::back_inserter(ptrs), toPtr);
    for (int i = 0; i < static_cast<int>(ptrs.size()); i++) {
      TTTrackRef* master = ptrs[i];
      if (!master)
        continue;
      for (int j = i + 1; j < static_cast<int>(ptrs.size()); j++) {
        TTTrackRef*& slave = ptrs[j];
        if (!slave)
          continue;
        if (equalEnough(*master, *slave))
          slave = nullptr;
      }
    }
    for (int i = static_cast<int>(ptrs.size()) - 1; i >= 0; i--)
      if (!ptrs[i])
        ttTrackRefs.erase(std::next(ttTrackRefs.begin(), i));
    // perform KF
    tt::TTTracks ttTracks;
    ttTracks.reserve(ttTrackRefs.size());
    for (const TTTrackRef& ttTrackRef : ttTrackRefs) {
      const int iRegion = ttTrackRef->phiSector();
      const double phiR = iRegion * setup_->regRangePhiT();
      const double inv2R = -.5 * ttTrackRef->rInv();
      const double cot = ttTrackRef->tanL();
      const std::vector<TTStubRef>& ttStubRefs = ttTrackRef->getStubRefs();
      const int size = ttStubRefs.size();
      std::vector<std::vector<TTStubRef>> permutations;
      permutations.reserve(nPer_[size - setup_->kfMinLayers()]);
      for (int nStubs = setup_->kfMinLayers(); nStubs <= size; nStubs++) {
        // form all nStubs out of size combinations
        std::string bitmask(nStubs, 1);
        bitmask.resize(size, 0);
        do {
          permutations.emplace_back();
          std::vector<TTStubRef>& comb = permutations.back();
          comb.reserve(nStubs);
          for (int i = 0; i < size; ++i)
            if (bitmask[i])
              comb.push_back(ttStubRefs[i]);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }
      ttTracks.emplace_back(0., 0., 0., 0., 0., 9.e3, 9.e3, 0., 0., 0., 0, setup_->simNPar(), setup_->sysBField());
      TTTrack<Ref_Phase2TrackerDigi_>& ttTrack = ttTracks.back();
      // fit all permutations
      for (const std::vector<TTStubRef>& permutation : permutations) {
        TTBV hitPattern(0, setup_->kfNumLayers());
        std::vector<Stub> stubs;
        stubs.reserve(permutation.size());
        for (const TTStubRef& ttStubRef : permutation) {
          const GlobalPoint gp = setup_->stubPosTT(ttStubRef);
          const trackerDTC::SensorModule* sm = setup_->sensorModule(ttStubRef);
          stubs.emplace_back();
          stubs.back().m0_ = tt::deltaPhi(gp.phi() - phiR);
          stubs.back().m1_ = gp.z();
          stubs.back().v0_ = std::pow(sm->dPhi(gp.perp(), inv2R), 2) / 12.;
          stubs.back().v1_ = std::pow(sm->dZ(cot), 2) / 12.;
          stubs.back().H_ = gp.perp();
          hitPattern.set(sm->layerIdReduced());
        }
        double x0(0.);
        double x1(0.);
        double x2(0.);
        double x3(0.);
        double x4(0.);
        double C00(9.e3);
        double C01(0.);
        double C11(9.e3);
        double C22(9.e3);
        double C23(0.);
        double C33(9.e3);
        double C44(setup_->simNPar() == 5 ? 9.e3 : 0.);
        double C40(0.);
        double C41(0.);
        double chi20(0.);
        double chi21(0.);
        // fit twice, first fit without ho corrections, second with
        for (int cor = 0; cor < 2; cor++) {
          // apply ho corrections
          if (cor == 1) {
            const double R = .5 / x0;
            const double R0 = R + x4;
            for (Stub& stub : stubs) {
              const double lin0 = x0 * stub.H_ + x4 / stub.H_;
              const double lin1 = x2 * stub.H_;
              const double nonLin0 = std::asin((stub.H_ * stub.H_ + R0 * R0 - R * R) / 2. / stub.H_ / R0);
              const double nonLin1 = std::abs(R) * x2 * std::acos((R * R + R0 * R0 - stub.H_ * stub.H_) / 2. / R / R0);
              stub.m0_ += lin0 - nonLin0;
              stub.m1_ += lin1 - nonLin1;
            }
            x0 = 0;
            x1 = 0;
            x2 = 0;
            x3 = 0;
            x4 = 0;
            C00 = 9.e3;
            C01 = 0.;
            C11 = 9.e3;
            C22 = 9.e3;
            C23 = 0.;
            C33 = 9.e3;
            C44 = setup_->simNPar() == 5 ? 9.e3 : 0.;
            C40 = 0.;
            C41 = 0.;
            chi20 = 0.;
            chi21 = 0.;
          }
          // add all stubs using KF update maths
          for (const Stub& stub : stubs) {
            const double r0 = stub.m0_ - x1 - x0 * stub.H_ - x4 / stub.H_;
            const double r1 = stub.m1_ - x3 - x2 * stub.H_;
            const double S00 = C01 + stub.H_ * C00 + C40 / stub.H_;
            const double S01 = C11 + stub.H_ * C01 + C41 / stub.H_;
            const double S12 = C23 + stub.H_ * C22;
            const double S13 = C33 + stub.H_ * C23;
            const double S04 = C41 + stub.H_ * C40 + C44 / stub.H_;
            const double R00 = stub.v0_ + S01 + stub.H_ * S00 + S04 / stub.H_;
            const double R11 = stub.v1_ + S13 + stub.H_ * S12;
            const double K00 = S00 / R00;
            const double K10 = S01 / R00;
            const double K21 = S12 / R11;
            const double K31 = S13 / R11;
            const double K40 = S04 / R00;
            x0 += r0 * K00;
            x1 += r0 * K10;
            x2 += r1 * K21;
            x3 += r1 * K31;
            x4 += r0 * K40;
            C00 -= S00 * K00;
            C01 -= S01 * K00;
            C11 -= S01 * K10;
            C22 -= S12 * K21;
            C23 -= S13 * K21;
            C33 -= S13 * K31;
            C44 -= S04 * K40;
            C40 -= S04 * K00;
            C41 -= S04 * K10;
            chi20 += r0 * r0 / R00;
            chi21 += r1 * r1 / R11;
          }
        }
        math::ErrorF<5>::type covMat;
        const std::array<std::array<double, 5>, 5> css{{{{C00, C01, 0., 0., C40}},
                                                        {{C01, C11, 0., 0., C41}},
                                                        {{0., 0., C22, C23, 0.}},
                                                        {{0., 0., C23, C33, 0.}},
                                                        {{C40, C41, 0., 0., C44}}}};
        for (int i = 0; i < 5; i++)
          for (int j = 0; j < 5; j++)
            covMat[i][j] = css[i][j];
        // TTTrack conversion
        TTTrack<Ref_Phase2TrackerDigi_> comb(-2. * x0,
                                             tt::deltaPhi(x1 + phiR),
                                             x2,
                                             x3,
                                             -x4,
                                             chi20,
                                             chi21,
                                             0.,
                                             0.,
                                             0.,
                                             hitPattern.val(),
                                             setup_->simNPar(),
                                             setup_->sysBField(),
                                             iRegion,
                                             ttTrackRef->etaSector(),
                                             0.,
                                             ttTrackRef->trackSeedType(),
                                             covMat);
        // keep best combination
        if (comb.chi2Red() > ttTrack.chi2Red())
          continue;
        ttTrack = comb;
        ttTrack.setStubRefs(permutation);
      }
      // finish TTTrack
      ttTrack.setChi2BendRed(StubPtConsistency::getConsistency(
          ttTrack, setup_->trackerGeometry(), setup_->trackerTopology(), setup_->sysBField(), setup_->simNPar()));
      ttTrack.setTrackWordBits();
    }
    // store products
    iEvent.emplace(edPutTokenTracks_, std::move(ttTracks));
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerSim);
