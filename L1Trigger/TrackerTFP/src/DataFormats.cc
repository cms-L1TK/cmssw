#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackTrigger/interface/StubPtConsistency.h"

#include <vector>
#include <deque>
#include <cmath>
#include <tuple>
#include <iterator>
#include <algorithm>
#include <string>
#include <iostream>
#include <numeric>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  // default constructor, trying to needs space as proper constructed object
  DataFormats::DataFormats()
      : numDataFormats_(0),
        formats_(+Variable::end, std::vector<DataFormat*>(+Process::end, nullptr)),
        numUnusedBitsStubs_(+Process::end, TTBV::S_),
        numUnusedBitsTracks_(+Process::end, TTBV::S_),
        numChannel_(+Process::end, 0) {
    setup_ = nullptr;
    countFormats();
    dataFormats_.reserve(numDataFormats_);
    numStreams_.reserve(+Process::end);
    numStreamsStubs_.reserve(+Process::end);
    numStreamsTracks_.reserve(+Process::end);
  }

  // method to count number of unique data formats
  template <Variable v, Process p>
  void DataFormats::countFormats() {
    if constexpr (config_[+v][+p] == p)
      numDataFormats_++;
    if constexpr (++p != Process::end)
      countFormats<v, ++p>();
    else if constexpr (++v != Variable::end)
      countFormats<++v>();
  }

  // proper constructor
  DataFormats::DataFormats(const Setup* setup) : DataFormats() {
    setup_ = setup;
    fillDataFormats();
    for (const Process p : Processes)
      for (const Variable v : stubs_[+p])
        if (p == Process::mht)
        numUnusedBitsStubs_[+p] -= formats_[+v][+p] ? formats_[+v][+p]->width() : 0;
    for (const Process p : Processes)
      for (const Variable v : tracks_[+p])
        numUnusedBitsTracks_[+p] -= formats_[+v][+p] ? formats_[+v][+p]->width() : 0;
    numChannel_[+Process::dtc] = setup_->numDTCsPerRegion();
    numChannel_[+Process::pp] = setup_->numDTCsPerTFP();
    numChannel_[+Process::gp] = setup_->numSectors();
    numChannel_[+Process::ht] = setup_->htNumBinsInv2R();
    numChannel_[+Process::mht] = setup_->htNumBinsInv2R();
    numChannel_[+Process::zht] = setup_->htNumBinsInv2R();
    numChannel_[+Process::kfin] = setup_->kfNumWorker() * setup_->numLayers();
    numChannel_[+Process::kf] = setup_->kfNumWorker();
    transform(numChannel_.begin(), numChannel_.end(), back_inserter(numStreams_), [this](int channel) {
      return channel * setup_->numRegions();
    });
    numStreamsStubs_ = numStreams_;
    numStreamsStubs_[+Process::kf] = numStreams_[+Process::kfin];
    numStreamsTracks_ = vector<int>(+Process::end, 0);
    numStreamsTracks_[+Process::kfin] = numStreams_[+Process::kf];
    numStreamsTracks_[+Process::kf] = numStreams_[+Process::kf];
    // Print digi data format of all variables of any specified algo step
    // (Look at DataFormat.h::tracks_ to see variable names).
    //for (const Variable v : tracks_[+Process::kf]) {
    //  const DataFormat& f = format(v, Process::kf);
    //  cout <<" KF "<< f.base() << " " << f.range() << " " << f.width() << endl;
    //}
  }

  // constructs data formats of all unique used variables and flavours
  template <Variable v, Process p>
  void DataFormats::fillDataFormats() {
    if constexpr (config_[+v][+p] == p) {
      dataFormats_.emplace_back(Format<v, p>(setup_));
      fillFormats<v, p>();
    }
    if constexpr (++p != Process::end)
      fillDataFormats<v, ++p>();
    else if constexpr (++v != Variable::end)
      fillDataFormats<++v>();
  }

  // helper (loop) data formats of all unique used variables and flavours
  template <Variable v, Process p, Process it>
  void DataFormats::fillFormats() {
    if (config_[+v][+it] == p) {
      formats_[+v][+it] = &dataFormats_.back();
    }
    if constexpr (++it != Process::end)
      fillFormats<v, p, ++it>();
  }

  // converts bits to ntuple of variables
  template <typename... Ts>
  void DataFormats::convertStub(Process p, const Frame& bv, tuple<Ts...>& data) const {
    TTBV ttBV(bv);
    extractStub(p, ttBV, data);
  }

  // helper (loop) to convert bits to ntuple of variables
  template <int it, typename... Ts>
  void DataFormats::extractStub(Process p, TTBV& ttBV, std::tuple<Ts...>& data) const {
    Variable v = *next(stubs_[+p].begin(), sizeof...(Ts) - 1 - it);
    formats_[+v][+p]->extract(ttBV, get<sizeof...(Ts) - 1 - it>(data));
    if constexpr (it + 1 != sizeof...(Ts))
      extractStub<it + 1>(p, ttBV, data);
  }

  // converts ntuple of variables to bits
  template <typename... Ts>
  void DataFormats::convertStub(Process p, const std::tuple<Ts...>& data, Frame& bv) const {
    TTBV ttBV(1, numUnusedBitsStubs_[+p]);
    attachStub(p, data, ttBV);
    bv = ttBV.bs();
  }

  // helper (loop) to convert ntuple of variables to bits
  template <int it, typename... Ts>
  void DataFormats::attachStub(Process p, const tuple<Ts...>& data, TTBV& ttBV) const {
    Variable v = *next(stubs_[+p].begin(), it);
    cout << +p << " " << +v << " " << get<it>(data) << endl;
    formats_[+v][+p]->attach(get<it>(data), ttBV);
    if constexpr (it + 1 != sizeof...(Ts))
      attachStub<it + 1>(p, data, ttBV);
  }

  // converts bits to ntuple of variables
  template <typename... Ts>
  void DataFormats::convertTrack(Process p, const Frame& bv, tuple<Ts...>& data) const {
    TTBV ttBV(bv);
    extractTrack(p, ttBV, data);
  }

  // helper (loop) to convert bits to ntuple of variables
  template <int it, typename... Ts>
  void DataFormats::extractTrack(Process p, TTBV& ttBV, std::tuple<Ts...>& data) const {
    Variable v = *next(tracks_[+p].begin(), sizeof...(Ts) - 1 - it);
    formats_[+v][+p]->extract(ttBV, get<sizeof...(Ts) - 1 - it>(data));
    if constexpr (it + 1 != sizeof...(Ts))
      extractTrack<it + 1>(p, ttBV, data);
  }

  // converts ntuple of variables to bits
  template <typename... Ts>
  void DataFormats::convertTrack(Process p, const std::tuple<Ts...>& data, Frame& bv) const {
    TTBV ttBV(1, numUnusedBitsTracks_[+p]);
    attachTrack(p, data, ttBV);
    bv = ttBV.bs();
  }

  // helper (loop) to convert ntuple of variables to bits
  template <int it, typename... Ts>
  void DataFormats::attachTrack(Process p, const tuple<Ts...>& data, TTBV& ttBV) const {
    Variable v = *next(tracks_[+p].begin(), it);
    //cout << +p << " " << +v << " " << get<it>(data) << endl;
    formats_[+v][+p]->attach(get<it>(data), ttBV);
    if constexpr (it + 1 != sizeof...(Ts))
      attachTrack<it + 1>(p, data, ttBV);
  }

  // construct Stub from Frame
  template <typename... Ts>
  Stub<Ts...>::Stub(const FrameStub& frame, const DataFormats* dataFormats, Process p)
      : dataFormats_(dataFormats), p_(p), frame_(frame), trackId_(0) {
    dataFormats_->convertStub(p, frame.second, data_);
  }

  // construct Stub from other Stub
  template <typename... Ts>
  template <typename... Others>
  Stub<Ts...>::Stub(const Stub<Others...>& stub, Ts... data)
      : dataFormats_(stub.dataFormats()),
        p_(++stub.p()),
        frame_(stub.frame().first, Frame()),
        data_(data...),
        trackId_(0) {}

  // construct Stub from TTStubRef
  template <typename... Ts>
  Stub<Ts...>::Stub(const TTStubRef& ttStubRef, const DataFormats* dataFormats, Process p, Ts... data)
      : dataFormats_(dataFormats), p_(p), frame_(ttStubRef, Frame()), data_(data...), trackId_(0) {}

  // construct StubPP from Frame
  StubPP::StubPP(const FrameStub& frame, const DataFormats* formats) : Stub(frame, formats, Process::pp) {
    for (int sectorEta = sectorEtaMin(); sectorEta <= sectorEtaMax(); sectorEta++)
      for (int sectorPhi = 0; sectorPhi < width(Variable::sectorsPhi); sectorPhi++)
        sectors_[sectorEta * width(Variable::sectorsPhi) + sectorPhi] = sectorsPhi()[sectorPhi];
  }

  // construct StubGP from Frame
  StubGP::StubGP(const FrameStub& frame, const DataFormats* formats, int sectorPhi, int sectorEta)
      : Stub(frame, formats, Process::gp), sectorPhi_(sectorPhi), sectorEta_(sectorEta) {
    inv2RBins_ = TTBV(0, setup()->htNumBinsInv2R());
    for (int inv2R = inv2RMin(); inv2R <= inv2RMax(); inv2R++)
      inv2RBins_.set(inv2R + inv2RBins_.size() / 2);
  }

  // construct StubGP from StubPP
  StubGP::StubGP(const StubPP& stub, int sectorPhi, int sectorEta)
      : Stub(stub, stub.r(), stub.phi(), stub.z(), stub.layer(), stub.inv2RMin(), stub.inv2RMax()),
        sectorPhi_(sectorPhi),
        sectorEta_(sectorEta) {
    get<1>(data_) -= (sectorPhi_ - .5) * setup()->baseSector();
    get<2>(data_) -= (r() + setup()->chosenRofPhi()) * setup()->sectorCot(sectorEta_);
    dataFormats_->convertStub(p_, data_, frame_.second);
  }

  // construct StubHT from Frame
  StubHT::StubHT(const FrameStub& frame, const DataFormats* formats, int inv2R)
      : Stub(frame, formats, Process::ht), inv2R_(inv2R) {
    dPhi_ = setup()->dPhi(ttStubRef(), format(Variable::inv2R).floating(inv2R));
    ps_ = setup()->psModule(ttStubRef());
    fillTrackId();
  }

  // construct StubHT from StubGP and HT cell assignment
  StubHT::StubHT(const StubGP& stub, int phiT, int inv2R)
      : Stub(stub, stub.r(), stub.phi(), stub.z(), stub.layer(), stub.sectorPhi(), stub.sectorEta(), phiT), inv2R_(inv2R) {
    get<1>(data_) -= format(Variable::inv2R).floating(inv2R) * r() + format(Variable::phiT).floating(phiT);
    dataFormats_->convertStub(p_, data_, frame_.second);
    fillTrackId();
  }

  void StubHT::fillTrackId() {
    TTBV ttBV(bv());
    trackId_ = ttBV.extract(width(Variable::sectorPhi) + width(Variable::sectorEta) + width(Variable::phiT));
  }

  // construct StubMHT from Frame
  StubMHT::StubMHT(const FrameStub& frame, const DataFormats* formats) : Stub(frame, formats, Process::mht) {
    TTBV ttBV(bv());
    dZ_ = setup()->dZ(ttStubRef());
    ps_ = setup()->psModule(ttStubRef());
    trackId_ = ttBV.extract(width(Variable::sectorPhi) + width(Variable::sectorEta) + width(Variable::inv2R) +
                            width(Variable::phiT));
  }

  StubMHT::StubMHT(const StubHT& stub, double phi, int inv2R, int phiT)
      : Stub(stub, stub.r(), phi, stub.z(), stub.layer(), stub.sectorPhi(), stub.sectorEta(), 0, 0), ps_(stub.ps()) {
    get<6>(data_) = inv2R;
    get<7>(data_) = stub.phiT() * setup()->mhtNumBinsPhiT() + phiT + setup()->mhtNumBinsPhiT() / 2.;
    dataFormats_->convertStub(p_, data_, frame_.second);
  }

  // construct StubZHT from Frame
  StubZHT::StubZHT(const FrameStub& frame, const DataFormats* formats, int inv2R) : Stub(frame, formats, Process::zht) {
    TTBV ttBV(bv());
    inv2R_ = inv2R * setup()->mhtNumBinsInv2R() + this->inv2R() + setup()->mhtNumBinsInv2R() / 2.;
    trackId_ = ttBV.extract(width(Variable::sectorPhi) + width(Variable::sectorEta) + width(Variable::phiT) +
                            width(Variable::inv2R) + width(Variable::zT) + width(Variable::cot));
  }

  StubZHT::StubZHT(const StubMHT& stub, double z, int cot, int zT)
      : Stub(stub, stub.r(), stub.phi(), z, stub.layer(), stub.sectorPhi(), stub.sectorEta(), stub.inv2R(), stub.phiT(), cot, zT) {
    dataFormats_->convertStub(p_, data_, frame_.second);
  }

  // construct StubKFin from Frame
  StubKFin::StubKFin(const FrameStub& frame, const DataFormats* formats, int layer)
      : Stub(frame, formats, Process::kfin), layer_(layer) {}

  // construct StubKFin from StubZHT
  StubKFin::StubKFin(const StubZHT& stub, int layer)
      : Stub(stub, stub.r(), stub.phi(), stub.z(), 0., 0.), layer_(layer) {
    get<3>(data_) = setup()->dPhi(ttStubRef(), format(Variable::inv2R).floating(stub.inv2R()));
    get<4>(data_) = setup()->dZ(ttStubRef());
    dataFormats_->convertStub(p_, data_, frame_.second);
  }

  // construct StubKFin from TTStubRef
  StubKFin::StubKFin(const TTStubRef& ttStubRef,
                     const DataFormats* dataFormats,
                     double r,
                     double phi,
                     double z,
                     double dPhi,
                     double dZ,
                     int layer)
      : Stub(ttStubRef, dataFormats, Process::kfin, r, phi, z, dPhi, dZ), layer_(layer) {
    dataFormats_->convertStub(p_, data_, frame_.second);
  }

  // construct StubKF from Frame
  StubKF::StubKF(const FrameStub& frame, const DataFormats* formats, int layer)
      : Stub(frame, formats, Process::kf), layer_(layer) {}

  // construct StubKF from StubKFin
  StubKF::StubKF(const StubKFin& stub, double inv2R, double phiT, double cot, double zT)
      : Stub(stub, stub.r(), 0., 0., stub.dPhi(), stub.dZ()), layer_(stub.layer()) {
    get<1>(data_) = format(Variable::phi).digi(stub.phi() - (phiT + this->r() * inv2R));
    const double d = setup()->chosenRofPhi() - setup()->chosenRofZ();
    const double rz = format(Variable::r).digi(this->r() + d);
    get<2>(data_) = format(Variable::z).digi(stub.z() - (zT + rz * cot));
    dataFormats_->convertStub(p_, data_, frame_.second);
  }

  // construct Track from Frame
  template <typename... Ts>
  Track<Ts...>::Track(const FrameTrack& frame, const DataFormats* dataFormats, Process p)
      : dataFormats_(dataFormats), p_(p), frame_(frame) {
    dataFormats_->convertTrack(p_, frame.second, data_);
  }

  // construct Track from other Track
  template <typename... Ts>
  template <typename... Others>
  Track<Ts...>::Track(const Track<Others...>& track, Ts... data)
      : dataFormats_(track.dataFormats()), p_(++track.p()), frame_(track.frame().first, Frame()), data_(data...) {}

  // construct Track from Stub
  template <typename... Ts>
  template <typename... Others>
  Track<Ts...>::Track(const Stub<Others...>& stub, const TTTrackRef& ttTrackRef, Ts... data)
      : dataFormats_(stub.dataFormats()), p_(++stub.p()), frame_(ttTrackRef, Frame()), data_(data...) {}

  // construct Track from TTTrackRef
  template <typename... Ts>
  Track<Ts...>::Track(const TTTrackRef& ttTrackRef, const DataFormats* dataFormats, Process p, Ts... data)
      : dataFormats_(dataFormats), p_(p), frame_(ttTrackRef, Frame()), data_(data...) {}

  // construct TrackKFin from Frame
  TrackKFin::TrackKFin(const FrameTrack& frame, const DataFormats* dataFormats, const vector<StubKFin*>& stubs)
      : Track(frame, dataFormats, Process::kfin), stubs_(setup()->numLayers()), hitPattern_(0, setup()->numLayers()) {
    vector<int> nStubs(stubs_.size(), 0);
    for (StubKFin* stub : stubs)
      nStubs[stub->layer()]++;
    for (int layer = 0; layer < dataFormats->setup()->numLayers(); layer++)
      stubs_[layer].reserve(nStubs[layer]);
    for (StubKFin* stub : stubs) {
      const int layer = stub->layer();
      stubs_[layer].push_back(stub);
      hitPattern_.set(layer);
    }
  }

  // construct TrackKFin from StubZHT
  TrackKFin::TrackKFin(const StubZHT& stub, const TTTrackRef& ttTrackRef, const TTBV& maybePattern)
      : Track(stub, ttTrackRef, maybePattern, stub.sectorPhi(), stub.sectorEta(), 0., 0., 0., 0.),
        stubs_(setup()->numLayers()),
        hitPattern_(0, setup()->numLayers()) {
    get<3>(data_) = format(Variable::inv2R, Process::kfin).floating(stub.inv2RFull());
    get<4>(data_) = format(Variable::phiT, Process::mht).floating(stub.phiT());
    get<5>(data_) = format(Variable::cot, Process::zht).floating(stub.cot());
    get<6>(data_) = format(Variable::zT, Process::zht).floating(stub.zT());
    dataFormats_->convertTrack(p_, data_, frame_.second);
  }

  // construct TrackKFin from TTTrackRef
  TrackKFin::TrackKFin(const TTTrackRef& ttTrackRef,
                       const DataFormats* dataFormats,
                       const TTBV& maybePattern,
                       int sectorPhi,
                       int sectorEta,
                       double inv2R,
                       double phiT,
                       double cot,
                       double zT)
      : Track(ttTrackRef, dataFormats, Process::kfin, maybePattern, sectorPhi, sectorEta, inv2R, phiT, cot, zT),
        stubs_(setup()->numLayers()),
        hitPattern_(0, setup()->numLayers()) {
    dataFormats_->convertTrack(p_, data_, frame_.second);
  }

  vector<TTStubRef> TrackKFin::ttStubRefs(const TTBV& hitPattern, const vector<int>& layerMap) const {
    vector<TTStubRef> stubs;
    stubs.reserve(hitPattern.count());
    for (int layer = 0; layer < setup()->numLayers(); layer++)
      if (hitPattern[layer])
        stubs.push_back(stubs_[layer][layerMap[layer]]->ttStubRef());
    return stubs;
  }

  // construct TrackKF from Frame
  TrackKF::TrackKF(const FrameTrack& frame, const DataFormats* dataFormats) : Track(frame, dataFormats, Process::kf) {}

  // construct TrackKF from TrackKfin
  TrackKF::TrackKF(const TrackKFin& track, const vector<StubKF*> stubs, int trackId, double inv2R, double phiT, double cot, double zT)
      : Track(track, false, track.sectorPhi(), track.sectorEta(), track.inv2R(), track.phiT(), track.cot(), track.zT()),
        stubs_(stubs), trackId_(trackId), hitPattern_(0, setup()->numLayers()), numSkippedLayers_(0), numConsistentLayers_(0) {
    get<0>(data_) = abs(inv2R) < dataFormats_->format(Variable::inv2R, Process::zht).base() / 2. &&
                    abs(phiT) < dataFormats_->format(Variable::phiT, Process::zht).base() / 2.;
    get<3>(data_) += inv2R;
    get<4>(data_) += phiT;
    get<5>(data_) += cot;
    get<6>(data_) += zT;
    dataFormats_->convertTrack(p_, data_, frame_.second);
    for (StubKF* stub : stubs_)
      hitPattern_.set(stub->layer());
    auto consistent = [this](int& sum, const StubKF* stub) {
      auto inConsistentRange = [](float v, float r, float d) { return abs(v) <= (r + d) / 2.; };
      // Check stub consistent with helix, allowing for stub & digi uncertainty
      const bool inRange0 = inConsistentRange(stub->phi(), stub->dPhi(), base(Variable::dPhi));
      const bool inRange1 = inConsistentRange(stub->z(), stub->dZ(), base(Variable::dZ));
      return sum += (inRange0 && inRange1 ? 1 : 0);
    };
    numConsistentLayers_ = accumulate(stubs.begin(), stubs.end(), 0, consistent);
    TTBV pattern = hitPattern_;
    pattern |= track.maybePattern();
    // Skipped layers before final stub on state
    numSkippedLayers_ = pattern.count(0, hitPattern_.pmEncode(), false);
  }

  // conversion to TTTrack
  TTTrack<Ref_Phase2TrackerDigi_> TrackKF::ttTrack() const {
    vector<StubKF> stubs;
    stubs.reserve(stubs_.size());
    transform(stubs_.begin(), stubs_.end(), back_inserter(stubs), [](StubKF* stub){ return *stub; });
    return ttTrack(stubs);
  }

  // conversion to TTTrack with given stubs
  TTTrack<Ref_Phase2TrackerDigi_> TrackKF::ttTrack(const vector<StubKF>& stubs) const {
    const double invR = -this->inv2R() * 2.;
    const double phi0 =
        deltaPhi(this->phiT() - this->inv2R() * setup()->chosenRofPhi() +
                 setup()->baseSector() * (this->sectorPhi() - .5) + setup()->baseRegion() * frame_.first->phiSector());
    const double cot = this->cot() + setup()->sectorCot(this->sectorEta());
    const double z0 = this->zT() - this->cot() * setup()->chosenRofZ();
    double chi2phi(0.);
    double chi2z(0.);
    vector<TTStubRef> ttStubRefs;
    ttStubRefs.reserve(stubs.size());
    for (const StubKF& stub : stubs) {
      chi2phi += pow(stub.phi(), 2) / pow(stub.dPhi(), 2);
      chi2z += pow(stub.z(), 2) / pow(stub.dZ(), 2);
      ttStubRefs.push_back(stub.ttStubRef());
    }
    static constexpr int nPar = 4;
    static constexpr double d0 = 0.;
    static constexpr double trkMVA1 = 0.;
    static constexpr double trkMVA2 = 0.;
    static constexpr double trkMVA3 = 0.;
    const int hitPattern = hitPattern_.val();
    const double bField = setup()->bField();
    TTTrack<Ref_Phase2TrackerDigi_> ttTrack(
        invR, phi0, cot, z0, d0, chi2phi, chi2z, trkMVA1, trkMVA2, trkMVA3, hitPattern, nPar, bField);
    ttTrack.setStubRefs(ttStubRefs);
    ttTrack.setPhiSector(frame_.first->phiSector());
    ttTrack.setEtaSector(this->sectorEta());
    ttTrack.setTrackSeedType(frame_.first->trackSeedType());
    ttTrack.setStubPtConsistency(StubPtConsistency::getConsistency(
        ttTrack, setup()->trackerGeometry(), setup()->trackerTopology(), bField, nPar));
    return ttTrack;
  }

  // construct TrackDR from Frame
  TrackDR::TrackDR(const FrameTrack& frame, const DataFormats* dataFormats) : Track(frame, dataFormats, Process::dr) {}

  // construct TrackDR from TrackKF
  TrackDR::TrackDR(const TrackKF& track) : Track(track, 0., 0., 0., 0.) {
    get<0>(data_) = track.inv2R();
    get<1>(data_) = track.phiT() + track.inv2R() * setup()->chosenRofPhi() +
                    dataFormats_->format(Variable::phi, Process::gp).range() * (track.sectorPhi() - .5);
    get<2>(data_) = track.cot() + setup()->sectorCot(track.sectorEta());
    get<3>(data_) = track.zT() - track.cot() * setup()->chosenRofZ();
    dataFormats_->convertTrack(p_, data_, frame_.second);
  }

  // conversion to TTTrack
  TTTrack<Ref_Phase2TrackerDigi_> TrackDR::ttTrack() const {
    const double inv2R = this->inv2R();
    const double phi0 = this->phi0();
    const double cot = this->cot();
    const double z0 = this->z0();
    static constexpr double d0 = 0.;
    static constexpr double chi2phi = 0.;
    static constexpr double chi2z = 0;
    static constexpr double trkMVA1 = 0.;
    static constexpr double trkMVA2 = 0.;
    static constexpr double trkMVA3 = 0.;
    static constexpr int hitPattern = 0.;
    static constexpr int nPar = 4;
    static constexpr double bField = 0.;
    const int sectorPhi = frame_.first->phiSector();
    const int sectorEta = frame_.first->etaSector();
    TTTrack<Ref_Phase2TrackerDigi_> ttTrack(
        inv2R, phi0, cot, z0, d0, chi2phi, chi2z, trkMVA1, trkMVA2, trkMVA3, hitPattern, nPar, bField);
    ttTrack.setPhiSector(sectorPhi);
    ttTrack.setEtaSector(sectorEta);
    return ttTrack;
  }

  template <>
  Format<Variable::inv2R, Process::ht>::Format(const Setup* setup) : DataFormat(true) {
    range_ = 2. * setup->invPtToDphi() / setup->minPt();
    base_ = range_ / (double)setup->htNumBinsInv2R();
    width_ = ceil(log2(setup->htNumBinsInv2R()));
  }
  template <>
  Format<Variable::phiT, Process::ht>::Format(const Setup* setup) : DataFormat(true) {
    range_ = 2. * M_PI / (double)(setup->numRegions() * setup->numSectorsPhi());
    base_ = range_ / (double)setup->htNumBinsPhiT();
    width_ = ceil(log2(setup->htNumBinsPhiT()));
  }
  template <>
  Format<Variable::phi, Process::ht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phiT, Process::ht> phiT(setup);
    range_ = 2. * phiT.base();
    const Format<Variable::phi, Process::dtc> phi(setup);
    base_ = phi.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::dPhi, Process::ht>::Format(const Setup* setup) : DataFormat(false) {
    const Format<Variable::phi, Process::dtc> phi(setup);
    range_ = setup->maxdPhi();
    base_ = phi.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::dZ, Process::ht>::Format(const Setup* setup) : DataFormat(false) {
    const Format<Variable::z, Process::dtc> z(setup);
    range_ = setup->maxdZ();
    base_ = z.base();
    width_ = ceil(log2(range_ / base_));
  }

  template <>
  Format<Variable::phi, Process::gp>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi, Process::dtc> phi(setup);
    base_ = phi.base();
    range_ = phi.range() / setup->numSectorsPhi();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::z, Process::gp>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::dtc> z(setup);
    base_ = z.base();
    range_ = setup->neededRangeChiZ();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::sectorPhi, Process::gp>::Format(const Setup* setup) : DataFormat(false) {
    range_ = setup->numSectorsPhi();
    width_ = ceil(log2(range_));
  }
  template <>
  Format<Variable::sectorEta, Process::gp>::Format(const Setup* setup) : DataFormat(false) {
    range_ = setup->numSectorsEta();
    width_ = ceil(log2(range_));
  }

  template <>
  Format<Variable::r, Process::dtc>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phiT, Process::ht> phiT(setup);
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    const double offset = setup->chosenRofPhi();
    width_ = setup->tmttWidthR();
    range_ = 2. * max(abs(setup->outerRadius() - offset), abs(setup->innerRadius() - offset));
    base_ = phiT.base() / inv2R.base();
    const int shift = ceil(log2(range_ / base_ / pow(2., width_)));
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::phi, Process::dtc>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phiT, Process::ht> phiT(setup);
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    const Format<Variable::r, Process::dtc> r(setup);
    width_ = setup->tmttWidthPhi();
    range_ = 2. * M_PI / (double)setup->numRegions() + inv2R.range() * r.range() / 2.;
    const int shift = ceil(log2(range_ / phiT.base() / pow(2., width_)));
    base_ = phiT.base() * pow(2., shift);
  }
  template <>
  Format<Variable::z, Process::dtc>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::r, Process::dtc> r(setup);
    width_ = setup->tmttWidthZ();
    range_ = 2. * setup->halfLength();
    const int shift = ceil(log2(range_ / r.base())) - width_;
    base_ = r.base() * pow(2., shift);
  }
  template <>
  Format<Variable::layer, Process::dtc>::Format(const Setup* setup) : DataFormat(false) {
    range_ = setup->numLayers();
    width_ = ceil(log2(range_));
  }
  template <>
  Format<Variable::sectorsPhi, Process::dtc>::Format(const Setup* setup) : DataFormat(false) {
    range_ = setup->numSectorsPhi();
    width_ = setup->numSectorsPhi();
  }

  template <>
  Format<Variable::inv2R, Process::mht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::inv2R, Process::ht> ht(setup);
    range_ = ht.base();
    base_ = ht.base() / setup->mhtNumBinsInv2R() / 2.;
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::phiT, Process::mht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phiT, Process::ht> ht(setup);
    range_ = ht.range();
    base_ = ht.base() / setup->mhtNumBinsPhiT() / 2.;
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::phi, Process::mht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::r, Process::dtc> r(setup);
    const Format<Variable::phi, Process::ht> phi(setup);
    const Format<Variable::dPhi, Process::ht> dPhi(setup);
    const Format<Variable::inv2R, Process::mht> inv2R(setup);
    const Format<Variable::phiT, Process::mht> phiT(setup);
    base_ = phi.base();
    range_ = phiT.base() * 4. + r.range() / 2. * inv2R.base() * 4. + dPhi.range();
    width_ = ceil(log2(range_ / base_));
  }

  template <>
  Format<Variable::cot, Process::zht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::zT, Process::zht> zT(setup);
    width_ = ceil(log2(setup->zhtNumBinsCot())) + 1;
    range_ = (zT.range() + 2. * setup->beamWindowZ()) / setup->chosenRofZ();
    const int shift = ceil(log2(range_)) - width_;
    base_ = pow(2., shift);
  }
  template <>
  Format<Variable::zT, Process::zht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::dtc> z(setup);
    width_ = ceil(log2(setup->zhtNumBinsZT())) + 1;
    const vector<int>& quantiles = setup->gpQuantiles();
    const int maxQ = *max_element(quantiles.begin(), quantiles.end());
    range_ = maxQ * setup->baseQuantile();
    const int shift = ceil(log2(range_ / z.base() / pow(2., width_)));
    base_ = z.base() * pow(2., shift);
  }
  template <>
  Format<Variable::z, Process::zht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::gp> z(setup);
    const Format<Variable::dZ, Process::ht> dZ(setup);
    const Format<Variable::cot, Process::zht> cot(setup);
    const Format<Variable::zT, Process::zht> zT(setup);
    const double offset = setup->chosenRofZ();
    const double rangeR = 2. * max(abs(setup->outerRadius() - offset), abs(setup->innerRadius() - offset));
    base_ = z.base();
    range_ = zT.base() * 4. + cot.base() * 4. * rangeR / 2. + dZ.range();
    width_ = ceil(log2(range_ / base_));
  }

  template <>
  Format<Variable::phi, Process::kfin>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi, Process::mht> phi(setup);
    range_ = phi.range() * pow(2, setup->kfinShiftRangePhi());
    base_ = phi.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::z, Process::kfin>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::zht> z(setup);
    range_ = z.range() * pow(2, setup->kfinShiftRangeZ());
    base_ = z.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::hitPattern, Process::kfin>::Format(const Setup* setup) : DataFormat(false) {
    width_ = setup->numLayers();
  }

  template <>
  Format<Variable::inv2R, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::inv2R, Process::dr> dr(setup);
    const Format<Variable::inv2R, Process::mht> mht(setup);
    range_ = mht.range() * setup->kfRangeFactor();
    base_ = dr.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::phiT, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi0, Process::dr> dr(setup);
    const Format<Variable::phiT, Process::mht> mht(setup);
    range_ = mht.range() * setup->kfRangeFactor();
    base_ = dr.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::cot, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::cot, Process::dr> dr(setup);
    const Format<Variable::cot, Process::zht> zht(setup);
    range_ = zht.range() * setup->kfRangeFactor();
    base_ = dr.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::zT, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z0, Process::dr> dr(setup);
    const Format<Variable::zT, Process::zht> zht(setup);
    range_ = zht.range() * setup->kfRangeFactor();
    base_ = dr.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::phi, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi, Process::kfin> phi(setup);
    const Format<Variable::r, Process::dtc> r(setup);
    const Format<Variable::inv2R, Process::kf> inv2R(setup);
    const Format<Variable::phiT, Process::kf> phiT(setup);
    const Format<Variable::dPhi, Process::ht> dPhi(setup);
    range_ = phiT.range() + inv2R.range() * r.range() / 2. + dPhi.range();
    base_ = phi.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::z, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::kfin> z(setup);
    const Format<Variable::dZ, Process::ht> dZ(setup);
    const Format<Variable::zT, Process::kf> zT(setup);
    const Format<Variable::cot, Process::kf> cot(setup);
    const double offset = setup->chosenRofZ();
    const double rangeR = 2. * max(abs(setup->outerRadius() - offset), abs(setup->innerRadius() - offset));
    range_ = zT.range() + cot.range() * rangeR / 2. + dZ.range();
    base_ = z.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::match, Process::kf>::Format(const Setup* setup) : DataFormat(false) {
    width_ = 1;
    range_ = 1.;
  }

  template <>
  Format<Variable::inv2R, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    width_ = setup->tfpWidthInv2R();
    range_ = inv2R.range();
    base_ = inv2R.base();
    const int shift = ceil(log2(range_ / base_ / pow(2., width_)));
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::phi0, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    const Format<Variable::phiT, Process::ht> phiT(setup);
    width_ = setup->tfpWidthPhi0();
    range_ = 2. * M_PI / (double)setup->numRegions() + inv2R.range() * setup->chosenRofPhi();
    base_ = phiT.base();
    const int shift = ceil(log2(range_ / base_ / pow(2., width_)));
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::cot, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::cot, Process::zht> cot(setup);
    width_ = setup->tfpWidthCot();
    range_ = 2. * setup->maxCot();
    base_ = cot.base();
    const int shift = ceil(log2(range_ / base_ / pow(2., width_)));
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::z0, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::zT, Process::zht> zT(setup);
    width_ = setup->tfpWidthZ0();
    range_ = 2. * setup->beamWindowZ();
    base_ = zT.base();
    const int shift = ceil(log2(range_ / base_ / pow(2., width_)));
    base_ *= pow(2., shift);
  }

}  // namespace trackerTFP
