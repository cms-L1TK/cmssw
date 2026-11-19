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
#include "DataFormats/Math/interface/Error.h"

#include "L1Trigger/TrackFindingTracklet/interface/Setup.h"
#include "L1Trigger/TrackTrigger/interface/StubPtConsistency.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/L1Track.h"

#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <numeric>
#include "conifer.h"

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
    // recalc tracklet parameter + chi2 and cov
    void consume(const L1TrackCollection&, int);
    // perform track multiplexer
    void tm();
    // perform duplicate removal
    void dr();
    // perform KF
    void kf();
    // perform TQ
    void tq();
    // perform tfp
    void tfp(L1TrackCollection&);

    struct Track;
    struct Stub {
      Stub(const TTStubRef&, const Setup*, Track*);
      int layerId() const { return sm_->layerId(); }
      int layerIdReduced() const { return sm_->layerIdReduced(); }
      TTStubRef ttStubRef_;
      const trackerDTC::SensorModule* sm_;
      double H_;
      double m0_;
      double m1_;
      double v0_;
      double v1_;
    };
    struct Track {
      Track() {}
      Track(const L1Track&, const Setup*);
      double chi2Red() const { return (chi20_ + chi21_) / proj_.size(); }
      // recalc seed parameter, cov & chi2
      void calcSeed();
      // add all stubs using KF update maths
      void add(const std::vector<const Stub*>&);
      const L1Track* l1Track_;
      double x0_ = 0.;
      double x1_ = 0.;
      double x2_ = 0.;
      double x3_ = 0.;
      double x4_ = 0.;
      double C00_ = 0.;
      double C01_ = 0.;
      double C11_ = 0.;
      double C22_ = 0.;
      double C23_ = 0.;
      double C33_ = 0.;
      double C44_ = 0.;
      double C40_ = 0.;
      double C41_ = 0.;
      double chi20_ = 0.;
      double chi21_ = 0.;
      double chi20Seed_ = 0.;
      double chi21Seed_ = 0.;
      double mva_ = 0.;
      double phiR_;
      TTBV hitPattern_;
      std::deque<const Stub*> seed_;
      std::deque<const Stub*> proj_;
    };
    // ED input token of TTTracks
    edm::EDGetTokenT<L1TrackCollection> edGetTokenTracks_;
    // ED output token of TTTracks
    edm::EDPutTokenT<L1TrackCollection> edPutTokenTracks_;
    // Setup token
    edm::ESGetToken<Setup, trackerDTC::SetupRcd> esGetTokenSetup_;
    // helper class to store configurations
    const Setup* setup_;
    // BDT
    conifer::BDT<float, float> bdt_;
    // path to model extended tracking
    std::string pathExtended_;
    // number of permutations for number of projections
    std::vector<int> nPer_;
    // stub container
    std::vector<Stub> stubs_;
    // track container
    std::vector<Track> tracks_;
    // data stream
    std::vector<Track*> stream_;
  };

  ProducerSim::ProducerSim(const edm::ParameterSet& iConfig)
      : bdt_(iConfig.getParameter<edm::FileInPath>("BDT").fullPath()),
        pathExtended_(iConfig.getParameter<edm::FileInPath>("BDT5ParSim").fullPath()) {
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
    // calc permutations for all found projection sizes [2 - 5]
    auto fac = [](int n) {
      int f(1);
      for (int i = 1; i <= n; i++)
        f *= i;
      return f;
    };
    auto bc = [fac](int n, int k) { return fac(n) / fac(k) / fac(n - k); };
    nPer_ = std::vector<int>(setup_->kfNumProj() - setup_->kfMinProj() + 1, 0);
    for (int i = setup_->kfMinProj(); i <= setup_->kfNumProj(); i++)
      for (int j = setup_->kfMinProj(); j <= i; j++)
        nPer_[i - setup_->kfMinProj()] += bc(i, j);
    // Load extended BDT model
    if (setup_->simNPar() == 5)
      bdt_ = conifer::BDT<float, float>(pathExtended_);
  }

  void ProducerSim::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // read input
    const L1TrackCollection& input = iEvent.get(edGetTokenTracks_);
    // prep output
    L1TrackCollection output;
    output.reserve(input.size());
    // process per nonant
    for (int iRegion = 0; iRegion < setup_->sysNumRegion(); iRegion++) {
      // recalc tracklet parameter + chi2 and cov
      consume(input, iRegion);
      // perform track multiplexer
      tm();
      // perform duplicate removal
      dr();
      // perform KF
      kf();
      // perform TQ
      tq();
      // perform tfp
      tfp(output);
    }
    // store products
    iEvent.emplace(edPutTokenTracks_, std::move(output));
  }

  // recalc tracklet parameter + chi2 and cov
  void ProducerSim::consume(const L1TrackCollection& l1TrackCollection, int iRegion) {
    // count tracks and stubs
    stubs_.clear();
    tracks_.clear();
    stream_.clear();
    int numTracks(0);
    int numStubs(0);
    for (const L1Track& l1Track : l1TrackCollection) {
      if (static_cast<int>(l1Track.phiSector()) != iRegion)
        continue;
      numTracks++;
      numStubs += l1Track.getStubRefs().size();
    }
    stubs_.reserve(numStubs);
    tracks_.reserve(numTracks);
    stream_.reserve(numTracks);
    for (const L1Track& l1Track : l1TrackCollection) {
      if (static_cast<int>(l1Track.phiSector()) != iRegion)
        continue;
      const std::vector<int>& seedLayers = setup_->tbSeedLayers(l1Track.trackSeedType());
      tracks_.emplace_back(l1Track, setup_);
      Track* track = &tracks_.back();
      stream_.push_back(track);
      for (const TTStubRef& ttStubRef : l1Track.getStubRefs()) {
        stubs_.emplace_back(ttStubRef, setup_, track);
        const Stub* stub = &stubs_.back();
        const auto it = std::find(seedLayers.begin(), seedLayers.end(), stub->layerId());
        if (it == seedLayers.end())
          track->proj_.push_back(stub);
        else {
          track->hitPattern_.set(stub->layerIdReduced());
          track->seed_.push_back(stub);
        }
      }
      // recalc seed parameter, cov & chi2
      track->calcSeed();
    }
  }

  // perform track multiplexer
  void ProducerSim::tm() {
    const std::vector<int>& muxOrder = setup_->tmMuxOrder();
    auto order = [&muxOrder](const Track* lhs, const Track* rhs) {
      const auto l = std::find(muxOrder.begin(), muxOrder.end(), lhs->l1Track_->trackSeedType());
      const auto r = std::find(muxOrder.begin(), muxOrder.end(), rhs->l1Track_->trackSeedType());
      return l < r;
    };
    std::sort(stream_.begin(), stream_.end(), order);
  }

  // perform duplicate removal
  void ProducerSim::dr() {
    auto equalEnough = [this](const Track* lhs, const Track* rhs) {
      std::vector<TTStubRef> l(lhs->l1Track_->getStubRefs());
      std::vector<TTStubRef> r(rhs->l1Track_->getStubRefs());
      std::sort(l.begin(), l.end());
      std::sort(r.begin(), r.end());
      std::vector<TTStubRef> same;
      same.reserve(std::min(l.size(), r.size()));
      std::set_intersection(l.begin(), l.end(), r.begin(), r.end(), std::back_inserter(same));
      return static_cast<int>(same.size()) >= setup_->drMinIdenticalStubs();
    };
    const int size = stream_.size();
    for (int i = 0; i < size; i++) {
      Track* master = stream_[i];
      if (!master)
        continue;
      for (int j = i + 1; j < size; j++) {
        Track*& slave = stream_[j];
        if (!slave)
          continue;
        if (equalEnough(master, slave))
          slave = nullptr;
      }
    }
  }

  // perform KF
  void ProducerSim::kf() {
    for (Track*& track : stream_) {
      if (!track)
        continue;
      const int size = track->proj_.size();
      // get all permutations
      std::vector<std::vector<const Stub*>> permutations;
      permutations.reserve(nPer_[size - setup_->kfMinProj()]);
      for (int nStubs = setup_->kfMinProj(); nStubs <= size; nStubs++) {
        // form all nStubs out of size combinations
        std::string bitmask(nStubs, 1);
        bitmask.resize(size, 0);
        do {
          permutations.emplace_back();
          std::vector<const Stub*>& comb = permutations.back();
          comb.reserve(nStubs);
          for (int i = 0; i < size; ++i)
            if (bitmask[i])
              comb.push_back(track->proj_[i]);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }
      bool first(true);
      Track best;
      // fit all permutations and keep best
      for (const std::vector<const Stub*>& permutation : permutations) {
        Track comb = *track;
        // add all stubs using KF update maths
        comb.add(permutation);
        // keep best combination
        if (!first && (best.chi2Red() < comb.chi2Red()))
          continue;
        first = false;
        best = comb;
      }
      *track = best;
    }
  }

  // perform TQ
  void ProducerSim::tq() {
    for (Track* track : stream_) {
      if (!track)
        continue;
      std::string hp = track->hitPattern_.str();
      std::reverse(hp.begin(), hp.end());
      // BDT input features
      const std::vector<float> feats = {static_cast<float>(track->x3_),
                                        static_cast<float>(track->x2_),
                                        static_cast<float>(track->chi20_),
                                        static_cast<float>(track->chi21_),
                                        static_cast<float>(TTBV(hp).val())};
      // BDT evaluation
      const float mva = bdt_.decision_function(feats)[0];
      // Apply sigmoid activation
      track->mva_ = 1. / (1. + std::exp(-mva));
    }
  }

  // perform tfp
  void ProducerSim::tfp(L1TrackCollection& l1TrackCollection) {
    for (const Track* track : stream_) {
      if (!track)
        continue;
      // convert to TTTrack
      math::ErrorF<5>::type covMat;
      const std::initializer_list<std::initializer_list<double>> css{{track->C00_, track->C01_, 0., 0., track->C40_},
                                                                     {track->C01_, track->C11_, 0., 0., track->C41_},
                                                                     {0., 0., track->C22_, track->C23_, 0.},
                                                                     {0., 0., track->C23_, track->C33_, 0.},
                                                                     {track->C40_, track->C41_, 0., 0., track->C44_}};
      int i(0);
      for (const std::initializer_list<double>& v : css) {
        int j(0);
        for (double d : v)
          covMat[i][j++] = d;
        i++;
      }
      l1TrackCollection.emplace_back(-2. * track->x0_,
                                     tt::deltaPhi(track->x1_ + track->phiR_),
                                     track->x2_,
                                     track->x3_,
                                     -track->x4_,
                                     track->chi20_,
                                     track->chi21_,
                                     track->mva_,
                                     0.,
                                     0.,
                                     track->hitPattern_.val(),
                                     setup_->simNPar(),
                                     setup_->sysBField(),
                                     track->l1Track_->phiSector(),
                                     track->l1Track_->etaSector(),
                                     0.,
                                     track->l1Track_->trackSeedType(),
                                     covMat);
      L1Track& l1Track = l1TrackCollection.back();
      for (const std::deque<const Stub*>& stubs : {track->seed_, track->proj_})
        for (const Stub* stub : stubs)
          l1Track.addStubRef(stub->ttStubRef_);
      l1Track.setChi2BendRed(StubPtConsistency::getConsistency(
          l1Track, setup_->trackerGeometry(), setup_->trackerTopology(), setup_->sysBField(), setup_->simNPar()));
      l1Track.setTrackWordBits();
    }
  }

  ProducerSim::Track::Track(const L1Track& l1Track, const Setup* setup)
      : l1Track_(&l1Track), phiR_(l1Track.phiSector() * setup->regRangePhiT()), hitPattern_(0, setup->kfNumLayers()) {}

  ProducerSim::Stub::Stub(const TTStubRef& ttStubRef, const Setup* setup, Track* track)
      : ttStubRef_(ttStubRef), sm_(setup->sensorModule(ttStubRef)) {
    const double x0 = -.5 * track->l1Track_->rInv();
    const double x1 = track->l1Track_->localPhi();
    const double x2 = track->l1Track_->tanL();
    const double x3 = track->l1Track_->z0();
    const double x4 = -track->l1Track_->d0();
    const GlobalPoint gp = setup->stubPosTT(ttStubRef);
    H_ = gp.perp();
    m0_ = tt::deltaPhi(gp.phi() - track->phiR_);
    m1_ = gp.z();
    v0_ = std::pow(sm_->dPhi(H_, x0), 2) / 3.;
    v1_ = std::pow(sm_->dZ(x2), 2) / 3.;
    // apply HO corrections
    if (x0 == 0. || x4 > H_)
      return;
    const double H2 = std::pow(H_, 2);
    const double R = .5 / x0;
    const double R0 = R + x4;
    // calc precise residual
    m0_ -= x1 + std::asin((H2 + R0 * R0 - R * R) / 2. / H_ / R0);
    m1_ -= x3 + std::abs(R) * x2 * std::acos((R * R + R0 * R0 - H2) / 2. / R / R0);
    // add linear track proj back on
    m0_ += x1 + x0 * H_ + x4 / H_;
    m1_ += x3 + x2 * H_;
  }

  // recalc seed parameter, cov & chi2
  void ProducerSim::Track::calcSeed() {
    const double dH = seed_[1]->H_ - seed_[0]->H_;
    const double dH2 = std::pow(dH, 2);
    // calc r-z parameter using first 2 seed stubs
    x2_ = (seed_[1]->m1_ - seed_[0]->m1_) / dH;
    x3_ = (seed_[1]->H_ * seed_[0]->m1_ - seed_[0]->H_ * seed_[1]->m1_) / dH;
    // calc r-z cov
    C22_ = (seed_[1]->v1_ + seed_[0]->v1_) / dH2;
    C23_ = -(seed_[1]->H_ * seed_[0]->v1_ + seed_[0]->H_ * seed_[1]->v1_) / dH2;
    C33_ = (std::pow(seed_[1]->H_, 2) * seed_[0]->v1_ + std::pow(seed_[0]->H_, 2) * seed_[1]->v1_) / dH2;
    if (seed_.size() == 2) {  // doublets
      // calc r-phi parameter
      x0_ = (seed_[1]->m0_ - seed_[0]->m0_) / dH;
      x1_ = (seed_[1]->H_ * seed_[0]->m0_ - seed_[0]->H_ * seed_[1]->m0_) / dH;
      // calc r-phi cov
      C00_ = (seed_[1]->v0_ + seed_[0]->v0_) / dH2;
      C01_ = -(seed_[1]->H_ * seed_[0]->v0_ + seed_[0]->H_ * seed_[1]->v0_) / dH2;
      C11_ = (std::pow(seed_[1]->H_, 2) * seed_[0]->v0_ + std::pow(seed_[0]->H_, 2) * seed_[1]->v0_) / dH2;
    } else {  // triplets
      const std::vector<double> denumP = {(seed_[0]->m0_ - seed_[1]->m0_) * (seed_[0]->m0_ - seed_[2]->m0_),
                                          (seed_[1]->m0_ - seed_[0]->m0_) * (seed_[1]->m0_ - seed_[2]->m0_),
                                          (seed_[2]->m0_ - seed_[0]->m0_) * (seed_[2]->m0_ - seed_[1]->m0_)};
      const std::vector<double> denumC = {
          std::pow(seed_[0]->H_ - seed_[1]->H_, 2) * std::pow(seed_[0]->H_ - seed_[2]->H_, 2),
          std::pow(seed_[1]->H_ - seed_[0]->H_, 2) * std::pow(seed_[1]->H_ - seed_[2]->H_, 2),
          std::pow(seed_[2]->H_ - seed_[0]->H_, 2) * std::pow(seed_[2]->H_ - seed_[1]->H_, 2)};
      // calc r-phi parameter
      x0_ += seed_[0]->m0_ * seed_[0]->H_ / denumP[0];
      x0_ += seed_[1]->m0_ * seed_[1]->H_ / denumP[1];
      x0_ += seed_[2]->m0_ * seed_[2]->H_ / denumP[2];
      x1_ -= seed_[0]->m0_ * seed_[0]->H_ * (seed_[1]->H_ + seed_[2]->H_) / denumP[0];
      x1_ -= seed_[1]->m0_ * seed_[1]->H_ * (seed_[0]->H_ + seed_[2]->H_) / denumP[1];
      x1_ -= seed_[2]->m0_ * seed_[2]->H_ * (seed_[0]->H_ + seed_[1]->H_) / denumP[2];
      x4_ += seed_[0]->m0_ * seed_[0]->H_ * seed_[1]->H_ * seed_[2]->H_ / denumP[0];
      x4_ += seed_[1]->m0_ * seed_[0]->H_ * seed_[1]->H_ * seed_[2]->H_ / denumP[1];
      x4_ += seed_[2]->m0_ * seed_[0]->H_ * seed_[1]->H_ * seed_[2]->H_ / denumP[2];
      // calc r-phi cov
      C00_ += seed_[0]->v0_ * std::pow(seed_[0]->H_, 2) / denumC[0];
      C00_ += seed_[1]->v0_ * std::pow(seed_[1]->H_, 2) / denumC[1];
      C00_ += seed_[2]->v0_ * std::pow(seed_[2]->H_, 2) / denumC[2];
      C01_ -= seed_[0]->v0_ * std::pow(seed_[0]->H_, 2) * (seed_[1]->H_ + seed_[2]->H_) / denumC[0];
      C01_ -= seed_[1]->v0_ * std::pow(seed_[1]->H_, 2) * (seed_[0]->H_ + seed_[2]->H_) / denumC[1];
      C01_ -= seed_[2]->v0_ * std::pow(seed_[2]->H_, 2) * (seed_[0]->H_ + seed_[1]->H_) / denumC[2];
      C11_ += seed_[0]->v0_ * std::pow(seed_[0]->H_, 2) * std::pow(seed_[1]->H_ + seed_[2]->H_, 2) / denumC[0];
      C11_ += seed_[1]->v0_ * std::pow(seed_[1]->H_, 2) * std::pow(seed_[0]->H_ + seed_[2]->H_, 2) / denumC[1];
      C11_ += seed_[2]->v0_ * std::pow(seed_[2]->H_, 2) * std::pow(seed_[0]->H_ + seed_[1]->H_, 2) / denumC[2];
      C44_ += seed_[0]->v0_ * std::pow(seed_[0]->H_ * seed_[1]->H_ * seed_[2]->H_, 2) / denumC[0];
      C44_ += seed_[1]->v0_ * std::pow(seed_[0]->H_ * seed_[1]->H_ * seed_[2]->H_, 2) / denumC[1];
      C44_ += seed_[2]->v0_ * std::pow(seed_[0]->H_ * seed_[1]->H_ * seed_[2]->H_, 2) / denumC[2];
      C40_ += seed_[0]->v0_ * std::pow(seed_[0]->H_, 2) * seed_[1]->H_ * seed_[2]->H_ / denumC[0];
      C40_ += seed_[1]->v0_ * std::pow(seed_[1]->H_, 2) * seed_[0]->H_ * seed_[2]->H_ / denumC[1];
      C40_ += seed_[2]->v0_ * std::pow(seed_[2]->H_, 2) * seed_[2]->H_ * seed_[1]->H_ / denumC[2];
      C41_ -= seed_[0]->v0_ * std::pow(seed_[0]->H_, 2) * seed_[1]->H_ * seed_[2]->H_ * (seed_[1]->H_ + seed_[2]->H_) /
              denumC[0];
      C41_ -= seed_[1]->v0_ * std::pow(seed_[1]->H_, 2) * seed_[0]->H_ * seed_[2]->H_ * (seed_[0]->H_ + seed_[2]->H_) /
              denumC[1];
      C41_ -= seed_[2]->v0_ * std::pow(seed_[2]->H_, 2) * seed_[0]->H_ * seed_[1]->H_ * (seed_[0]->H_ + seed_[1]->H_) /
              denumC[2];
      // KF update math for r-z
      const double r1 = seed_[2]->m1_ - x3_ - x2_ * seed_[2]->H_;
      const double S12 = C23_ + seed_[2]->H_ * C22_;
      const double S13 = C33_ + seed_[2]->H_ * C23_;
      const double R11 = seed_[2]->v1_ + S13 + seed_[2]->H_ * S12;
      const double K21 = S12 / R11;
      const double K31 = S13 / R11;
      x2_ += r1 * K21;
      x3_ += r1 * K31;
      C22_ -= S12 * K21;
      C23_ -= S13 * K21;
      C33_ -= S13 * K31;
    }
    // calculating chi2s for found track
    for (const std::deque<const Stub*>& stubs : {seed_, proj_}) {
      for (const Stub* stub : stubs) {
        // calc residuals
        const double r0 = stub->m0_ - (x1_ + x0_ * stub->H_ + x4_ / stub->H_);
        const double r1 = stub->m1_ - (x3_ + x2_ * stub->H_);
        // calc stub + track proj uncertainties
        const double S00 = C01_ + stub->H_ * C00_ + C40_ / stub->H_;
        const double S01 = C11_ + stub->H_ * C01_ + C41_ / stub->H_;
        const double S04 = C41_ + stub->H_ * C40_ + C44_ / stub->H_;
        const double S12 = C23_ + stub->H_ * C22_;
        const double S13 = C33_ + stub->H_ * C23_;
        const double R00 = stub->v0_ + S01 + stub->H_ * S00 + S04 / stub->H_;
        const double R11 = stub->v1_ + S13 + stub->H_ * S12;
        // update chi2sSeed
        chi20Seed_ += std::pow(r0, 2) / R00;
        chi21Seed_ += std::pow(r1, 2) / R11;
      }
    }
  }

  // add all stubs using KF update maths
  void ProducerSim::Track::add(const std::vector<const Stub*>& stubs) {
    proj_.clear();
    for (const Stub* stub : stubs) {
      const double r0 = stub->m0_ - x1_ - x0_ * stub->H_ - x4_ / stub->H_;
      const double r1 = stub->m1_ - x3_ - x2_ * stub->H_;
      const double S00 = C01_ + stub->H_ * C00_ + C40_ / stub->H_;
      const double S01 = C11_ + stub->H_ * C01_ + C41_ / stub->H_;
      const double S04 = C41_ + stub->H_ * C40_ + C44_ / stub->H_;
      const double S12 = C23_ + stub->H_ * C22_;
      const double S13 = C33_ + stub->H_ * C23_;
      const double R00 = stub->v0_ + S01 + stub->H_ * S00 + S04 / stub->H_;
      const double R11 = stub->v1_ + S13 + stub->H_ * S12;
      const double K00 = S00 / R00;
      const double K10 = S01 / R00;
      const double K21 = S12 / R11;
      const double K31 = S13 / R11;
      const double K40 = S04 / R00;
      x0_ += r0 * K00;
      x1_ += r0 * K10;
      x2_ += r1 * K21;
      x3_ += r1 * K31;
      x4_ += r0 * K40;
      C00_ -= S00 * K00;
      C01_ -= S01 * K00;
      C11_ -= S01 * K10;
      C22_ -= S12 * K21;
      C23_ -= S13 * K21;
      C33_ -= S13 * K31;
      C44_ -= S04 * K40;
      C40_ -= S04 * K00;
      C41_ -= S04 * K10;
      hitPattern_.set(stub->layerIdReduced());
      proj_.push_back(stub);
    }
    // calc chi2s
    for (const std::deque<const Stub*>& stubs : {seed_, proj_}) {
      for (const Stub* stub : stubs) {
        // calc residuals
        const double r0 = stub->m0_ - (x1_ + x0_ * stub->H_ + x4_ / stub->H_);
        const double r1 = stub->m1_ - (x3_ + x2_ * stub->H_);
        // calc stub + track proj uncertainties
        const double S00 = C01_ + stub->H_ * C00_ + C40_ / stub->H_;
        const double S01 = C11_ + stub->H_ * C01_ + C41_ / stub->H_;
        const double S04 = C41_ + stub->H_ * C40_ + C44_ / stub->H_;
        const double S12 = C23_ + stub->H_ * C22_;
        const double S13 = C33_ + stub->H_ * C23_;
        const double R00 = stub->v0_ + S01 + stub->H_ * S00 + S04 / stub->H_;
        const double R11 = stub->v1_ + S13 + stub->H_ * S12;
        // update chi2s
        chi20_ += std::pow(r0, 2) / R00;
        chi21_ += std::pow(r1, 2) / R11;
      }
    }
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerSim);
