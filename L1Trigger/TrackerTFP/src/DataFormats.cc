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

  // default constructor, trying to need space as proper constructed object
  DataFormats::DataFormats()
      : numDataFormats_(0),
        formats_(+Variable::end, std::vector<DataFormat*>(+Process::end, nullptr)),
        numUnusedBitsStubs_(+Process::end, TTBV::S_ - 1),
        numUnusedBitsTracks_(+Process::end, TTBV::S_ - 1),
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
    cout << "r     " << format(Variable::r, Process::kfin).width() << endl;
    cout << "phi   " << format(Variable::phi, Process::kfin).width() << endl;
    cout << "z     " << format(Variable::z, Process::kfin).width() << endl;
    cout << "dPhi  " << format(Variable::dPhi, Process::kfin).width() << endl;
    cout << "dZ    " << format(Variable::dZ, Process::kfin).width() << endl;
    cout << "inv2R " << format(Variable::inv2R, Process::kfin).width() << endl;
    cout << "phiT  " << format(Variable::phiT, Process::kfin).width() << endl;
    cout << "zT    " << format(Variable::zT, Process::kfin).width() << endl;
    for (const Process p : Processes)
      for (const Variable v : stubs_[+p])
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
    numChannel_[+Process::kfin] = setup_->kfNumWorker();
    numChannel_[+Process::kf] = setup_->kfNumWorker();
    transform(numChannel_.begin(), numChannel_.end(), back_inserter(numStreams_), [this](int channel) {
      return channel * setup_->numRegions();
    });
    numStreamsStubs_ = numStreams_;
    numStreamsTracks_ = vector<int>(+Process::end, 0);
    numStreamsTracks_[+Process::kfin] = numStreams_[+Process::kfin];
    numStreamsTracks_[+Process::kf] = numStreams_[+Process::kf];
    // Print digi data format of all variables of any specified algo step
    // (Look at DataFormat.h::tracks_ to see variable names).
    //for (const Variable v : tracks_[+Process::kf]) {
    //  const DataFormat& f = format(v, Process::kf);
    //  cout <<" KF "<< f.base() << " " << f.range() << " " << f.width() << endl;
    //}
    numStreamsStubs_[+Process::kfin] = setup_->numLayers() * numStreamsTracks_[+Process::kfin];
    numStreamsStubs_[+Process::kf] = setup_->numLayers() * numStreamsTracks_[+Process::kf];
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
        frame_(stub.ttStubRef(), Frame()),
        data_(data...),
        trackId_(0) {}

  // construct Stub from TTStubRef
  template <typename... Ts>
  Stub<Ts...>::Stub(const TTStubRef& ttStubRef, const DataFormats* dataFormats, Process p, Ts... data)
      : dataFormats_(dataFormats), p_(p), frame_(ttStubRef, Frame()), data_(data...), trackId_(0) {}

  // construct Stub from TTStubRef
  template <typename... Ts>
  Stub<Ts...>::Stub(const TTStubRef& ttStubRef, const DataFormats* dataFormats, Process p)
      : dataFormats_(dataFormats), p_(p), frame_(ttStubRef, Frame()), trackId_(0) {}

  // construct StubDTC from TTStubRef
  StubDTC::StubDTC(const TTStubRef& ttStubRef, const DataFormats* formats, int region, int layer, double r, double phi, double z, pair<double, double> phiT, pair<double, double> inv2R) : Stub(ttStubRef, formats, Process::dtc) {
    // stub position, phi w.r.t. processing region centre in rad
    const double offset = (region - .5) * format(Variable::phiT).range();
    get<0>(data_) = r;
    get<1>(data_) = phi - offset;
    get<2>(data_) = z;
    // convert unique layer id [1-6,11-15] into reduced layer id [0-6]
    // a fiducial track may not cross more then 7 detector layers, for stubs from a given track the reduced layer id is actually unique
    if (layer == 1)
      get<3>(data_) = 0;
    else if (layer == 2)
      get<3>(data_) = 1;
    else if (layer == 6 || layer == 11)
      get<3>(data_) = 2;
    else if (layer == 5 || layer == 12)
      get<3>(data_) = 3;
    else if (layer == 4 || layer == 13)
      get<3>(data_) = 4;
    else if (layer == 14)
      get<3>(data_) = 5;
    else if (layer == 3 || layer == 15)
      get<3>(data_) = 6;
    // assign stub to phiT sectors within a processing region
    const DataFormat& df = format(Variable::phiT);
    get<4>(data_) = df.integerF(phiT.first - offset);
    get<5>(data_) = df.integerF(phiT.second - offset);
    // inv2R min and max
    get<6>(data_) = format(Variable::inv2R).integer(inv2R.first);
    get<7>(data_) = format(Variable::inv2R).integer(inv2R.second);
  }

  // construct StubPP from Frame
  StubPP::StubPP(const FrameStub& frame, const DataFormats* formats) : Stub(frame, formats, Process::pp) {
    const double r = this->r() + setup()->chosenRofPhi();
    const double ratioRZ = setup()->chosenRofZ() / r;
    // extrapolated z at radius T assuming z0=0
    const double zT = this->z() * ratioRZ;
    // extrapolated z0 window at radius T
    const double dZT = setup()->beamWindowZ() * abs(1. - ratioRZ);
    const double zTMin = format(Variable::zT, Process::gp).integerF(zT - dZT);
    const double zTMax = format(Variable::zT, Process::gp).integerF(zT + dZT);
    for (int zT = zTMin; zT <= zTMax; zT++) {
      const int zTu = zT + setup()->numSectorsEta() / 2;
      for (int phiT = phiTMin(); phiT <= phiTMax(); phiT++) {
        const int phiTu = phiT + setup()->numSectorsPhi() / 2;;
        sectors_.set(zTu * setup()->numSectorsPhi() + phiTu);
      }
    }
  }

  // construct StubGP from Frame
  StubGP::StubGP(const FrameStub& frame, const DataFormats* formats, int channel) : Stub(frame, formats, Process::gp) {
    phiT_ = channel % setup()->numSectorsPhi() - setup()->numSectorsPhi() / 2;
    zT_ = channel / setup()->numSectorsPhi() - setup()->numSectorsEta() / 2;
    inv2RBins_ = TTBV(0, setup()->htNumBinsInv2R());
    for (int inv2R = inv2RMin(); inv2R <= inv2RMax(); inv2R++)
      inv2RBins_.set(inv2R + inv2RBins_.size() / 2);
  }

  // construct StubGP from StubPP
  StubGP::StubGP(const StubPP& stub, int phiT, int zT) : Stub(stub, stub.r(), stub.phi(), stub.z(), stub.layer(), TTBV(), stub.inv2RMin(), stub.inv2RMax()), phiT_(phiT), zT_(zT) {
    const double r = stub.r() + setup()->chosenRofPhi();
    const double cot = format(Variable::zT).floating(zT) / setup()->chosenRofZ();
    get<1>(data_) -= format(Variable::phiT).floating(phiT);
    get<2>(data_) -= r * cot;
    get<4>(data_) = setup()->module(r, stub.z());
  }

  // construct StubHT from Frame
  StubHT::StubHT(const FrameStub& frame, const DataFormats* formats, int channel) : Stub(frame, formats, Process::ht), valid_(false), array_(setup()->mhtNumBinsInv2R(), TTBV(0, setup()->mhtNumBinsPhiT())) {
    inv2R_ = format(Variable::inv2R).toSigned(channel);
    const double pitchRow = this->ps() ? setup()->pitchRowPS() : setup()->pitchRow2S();
    const double pitchCol = this->ps() ? setup()->pitchColPS() : setup()->pitchCol2S();
    const double pitchColR = this->barrel() ? (this->tilted() ? setup()->tiltUncertaintyR() : 0.0) : pitchCol;
    const double r = this->r() + setup()->chosenRofPhi();
    const double inv2Rf = format(Variable::inv2R).floating(inv2R_);
    const double sigma = format(Variable::phi).digi(pitchRow / r);
    const double dR = format(Variable::r).digi(setup()->scattering() + pitchColR);
    dPhi_ = format(Variable::phi).digi(sigma + dR * abs(inv2Rf) + format(Variable::phi).base());
    TTBV ttBV(bv());
    trackId_ = ttBV.extract(width(Variable::phiT) + width(Variable::zT));
  }

  // construct StubHT from StubGP and HT cell assignment
  StubHT::StubHT(const StubGP& stub, int inv2R, int phiT)
      : Stub(stub, false, stub.r(), stub.phi(), stub.z(), stub.layer(), stub.module(), phiT, stub.zT()), inv2R_(inv2R), phiT_(phiT) {
    get<2>(data_) -= format(Variable::inv2R).floating(inv2R) * r() + format(Variable::phiT).floating(phiT);
    get<6>(data_) += (stub.phiT() + .5) * setup()->htNumBinsPhiT();
    TTBV ttBV(bv());
    trackId_ = ttBV.extract(width(Variable::phiT) + width(Variable::zT));
  }

  // construct StubMHT from Frame
  StubMHT::StubMHT(const FrameStub& frame, const DataFormats* formats) : Stub(frame, formats, Process::mht), valid_(false), array_(setup()->zhtNumBinsCot(), TTBV(0, setup()->zhtNumBinsZT())) {
    const double cot = format(Variable::zT).floating(this->zT()) / setup()->chosenRofZ();
    const double pitchCol = this->ps() ? setup()->pitchColPS() : setup()->pitchCol2S();
    double sigma = abs(cot) * pitchCol;
    if (this->barrel())
      sigma = this->tilted() ? setup()->tiltApproxSlope() * abs(cot) + setup()->tiltApproxIntercept() : pitchCol;
    dZ_ = format(Variable::z).digi(sigma + format(Variable::z).base());
    TTBV ttBV(bv());
    trackId_ = ttBV.extract(width(Variable::phiT) + width(Variable::zT));
  }

  StubMHT::StubMHT(const StubHT& stub, bool last)
      : Stub(stub, last, stub.r(), stub.phi(), stub.z(), stub.layer(), stub.module(), stub.phiT(), stub.zT()), dPhi_(stub.dPhi()), valid_(stub.valid()) {
    trackId_ = stub.trackId();
  }

  // construct StubZHT from Frame
  StubZHT::StubZHT(const FrameStub& frame, const DataFormats* formats, int channel) : Stub(frame, formats, Process::zht) {
    const DataFormat& df = format(Variable::inv2R, Process::ht);
    inv2R_ = df.toSigned(channel);
    TTBV ttBV(bv());
    ttBV += df.ttBV(inv2R_);
    trackId_ = ttBV.extract(width(Variable::phiT) + width(Variable::zT) + df.width());
  }

  StubZHT::StubZHT(const StubMHT& stub, bool last)
      : Stub(stub, last, stub.r(), stub.phi(), stub.z(), stub.layer(), stub.module(), stub.phiT(), stub.zT()), dZ_(stub.dZ()), valid_(stub.valid()) {
    trackId_ = stub.trackId();
  }

  // construct StubKFin from Frame
  StubKFin::StubKFin(const FrameStub& frame, const DataFormats* formats, int layer)
      : Stub(frame, formats, Process::kfin), layer_(layer) {}

  // construct StubKFin from StubZHT
  StubKFin::StubKFin(const StubZHT& stub, int layer, int trackId, int stubId)
      : Stub(stub, stub.r(), stub.phi(), stub.z(), 0., 0.), layer_(layer), stubId_(stubId) {
    const double pitchRow = stub.ps() ? setup()->pitchRowPS() : setup()->pitchRow2S();
    const double pitchCol = stub.ps() ? setup()->pitchColPS() : setup()->pitchCol2S();
    const double pitchColR = stub.barrel() ? (stub.tilted() ? setup()->tiltUncertaintyR() : 0.0) : pitchCol;
    const double b = base(Variable::r) * pow(2, width(Variable::r) - setup()->kfinWidthAddrDPhi());
    const double r = (floor(this->r() / b + 1.e-12) + .5) * b + setup()->chosenRofPhi();
    const double inv2Rf = format(Variable::inv2R).floating(stub.inv2R());
    const double sigmaPhi = pitchRow / r;
    const double dR = setup()->scattering() + pitchColR;
    get<3>(data_) = format(Variable::phi).digi(sigmaPhi + dR * abs(inv2Rf) + format(Variable::phi).base());
    const double zT = format(Variable::zT).floating(stub.zT());
    const double cot = abs(zT) / setup()->chosenRofZ();
    double sigmaZ = cot * pitchCol;
    if (stub.barrel())
      sigmaZ = stub.tilted() ? setup()->tiltApproxSlope() * cot + setup()->tiltApproxIntercept() : pitchCol;
    get<4>(data_) = format(Variable::z).digi(sigmaZ + format(Variable::z).base());
    trackId_ = trackId;
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
      : Stub(ttStubRef, dataFormats, Process::kfin, r, phi, z, dPhi, dZ), layer_(layer) {}

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
      : dataFormats_(track.dataFormats()), p_(++track.p()), frame_(track.ttTrackRef(), Frame()), data_(data...) {}

  // construct Track from TTTrackRef
  template <typename... Ts>
  Track<Ts...>::Track(const TTTrackRef& ttTrackRef, const DataFormats* dataFormats, Process p, Ts... data)
      : dataFormats_(dataFormats), p_(p), frame_(ttTrackRef, Frame()), data_(data...) {}

  // construct TrackKFin from Frame
  TrackKFin::TrackKFin(const FrameTrack& frame, const DataFormats* dataFormats, const vector<StubKFin*>& stubs)
      : Track(frame, dataFormats, Process::kfin), stubs_(setup()->numLayers()), cot_(this->zT() / setup()->chosenRofZ()), hitPattern_(0, setup()->numLayers()) {
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

  // construct TrackKFin from TTTrackRef
  TrackKFin::TrackKFin(const TTTrackRef& ttTrackRef, const DataFormats* dataFormats, const TTBV& maybePattern, const TTBV& hitPattern, const vector<int>& lmap, int inv2R, int phiT, int zT, int trackId)
      : Track(ttTrackRef, dataFormats, Process::kfin, 0., 0., 0.), trackId_(trackId), maybePattern_(maybePattern), hitPattern_(hitPattern), lmap_(0, 0) {
    get<0>(data_) = format(Variable::inv2R, Process::ht).floating(inv2R);
    get<1>(data_) = format(Variable::phiT, Process::ht).floating(phiT);
    get<2>(data_) = format(Variable::zT, Process::gp).floating(zT);
    lmap_ = setup()->layerMap(hitPattern, lmap);
    size_ = *max_element(lmap.begin(), lmap.end());
    cot_ = this->zT() / setup()->chosenRofZ();
  }

  //
  void TrackKFin::kill(StubKFin* stub) {
    for (vector<StubKFin*>& layer : stubs_) {
      const auto s = find(layer.begin(), layer.end(), stub);
      if (s == layer.end())
        continue;
      *s = nullptr;
      break;
    }
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
  TrackKF::TrackKF(const FrameTrack& frame, const DataFormats* dataFormats) : Track(frame, dataFormats, Process::kf), cot_(this->cot()) {
    cot_ += format(Variable::zT, Process::gp).digi(this->zT()) / setup()->chosenRofZ();
  }

  // construct TrackKF from TrackKfin
  TrackKF::TrackKF(const TrackKFin& track, const vector<StubKF*> stubs, int trackId, double inv2R, double phiT, double cot, double zT)
      : Track(track, false, track.inv2R(), track.phiT(), 0., track.zT()),
        stubs_(stubs), trackId_(trackId), hitPattern_(0, setup()->numLayers()), numSkippedLayers_(0), numConsistentLayers_(0) {
    get<1>(data_) += inv2R;
    get<2>(data_) += phiT;
    get<3>(data_) += cot;
    get<4>(data_) += zT;
    const DataFormat& dfInv2R = format(Variable::inv2R, Process::ht);
    const DataFormat& dfphiT = format(Variable::phiT, Process::ht);
    get<0>(data_) = dfInv2R.integer(track.inv2R()) == dfInv2R.integer(this->inv2R()) &&
                    dfphiT.integer(track.phiT()) == dfphiT.integer(this->phiT());
    cot_ = cot + track.zT() / setup()->chosenRofZ();
    for (StubKF* stub : stubs_)
      hitPattern_.set(stub->layer());
    // Check stub consistent with helix
    auto consistent = [this](int& sum, const StubKF* stub) {
      static const DataFormat& phi = format(Variable::phi, Process::kf);
      static const DataFormat& z = format(Variable::z, Process::kf);
      const bool inRange0 = 2. * abs(stub->phi()) - stub->dPhi() < phi.base();
      const bool inRange1 = 2. * abs(stub->z()) - stub->dZ() < z.base();
      return sum += (inRange0 && inRange1 ? 1 : 0);
    };
    numConsistentLayers_ = accumulate(stubs.begin(), stubs.end(), 0, consistent);
    TTBV pattern = hitPattern_;
    pattern |= track.maybePattern();
    // Skipped layers before final stub on state
    numSkippedLayers_ = pattern.count(0, hitPattern_.pmEncode(), false);
  }

  // conversion to TTTrack with given stubs
  TTTrack<Ref_Phase2TrackerDigi_> TrackKF::ttTrack(int region, const vector<StubKF>& stubs) const {
    const double invR = -this->inv2R() * 2.;
    const double phi0 = deltaPhi(this->phiT() - this->inv2R() * setup()->chosenRofPhi() + region * format(Variable::phiT).range());
    const double cot = this->cotGlobal();
    const double z0 = this->zT() - cot * setup()->chosenRofZ();
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
    ttTrack.setPhiSector(region);
    ttTrack.setEtaSector(format(Variable::zT, Process::gp).integer(this->zT()));
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
    get<1>(data_) = track.phiT() + track.inv2R() * setup()->chosenRofPhi();
    const double cotSector = format(Variable::zT, Process::gp).digi(track.zT()) / setup()->chosenRofZ();
    get<2>(data_) = track.cot() + cotSector;
    get<3>(data_) = track.zT() - track.cot() * setup()->chosenRofZ();
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
  Format<Variable::r, Process::dtc>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phiT, Process::ht> phiT(setup);
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    width_ = setup->tmttWidthR();
    range_ = 2. * setup->maxRphi();
    base_ = phiT.base() / inv2R.base();
    const int shift = ceil(log2(range_ / base_)) - width_;
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::phi, Process::dtc>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phiT, Process::gp> phiT(setup);
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    const Format<Variable::r, Process::dtc> r(setup);
    width_ = setup->tmttWidthPhi();
    range_ = phiT.range() + inv2R.range() * setup->maxRphi();
    const int shift = ceil(log2(range_ / phiT.base())) - width_;
    base_ = phiT.base() * pow(2., shift);
  }
  template <>
  Format<Variable::z, Process::dtc>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::zT, Process::gp> zT(setup);
    width_ = setup->tmttWidthZ();
    range_ = 2. * setup->halfLength();
    const int shift = ceil(log2(range_ / zT.base())) - width_;
    base_ = zT.base() * pow(2., shift);
  }
  template <>
  Format<Variable::layer, Process::dtc>::Format(const Setup* setup) : DataFormat(false) {
    range_ = setup->numLayers();
    width_ = ceil(log2(range_));
  }

  template <>
  Format<Variable::phi, Process::gp>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi, Process::dtc> phi(setup);
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    const Format<Variable::phiT, Process::gp> phiT(setup);
    base_ = phi.base();
    range_ = phiT.base() + inv2R.range() * setup->maxRphi();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::z, Process::gp>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::dtc> z(setup);
    const Format<Variable::zT, Process::gp> zT(setup);
    base_ = z.base();
    range_ = zT.base() + setup->maxRz() * (zT.base() + 2. * setup->beamWindowZ()) / setup->chosenRofZ();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::module, Process::gp>::Format(const Setup* setup) : DataFormat(false) {
    width_ = setup->gpWidthModule();
  }
  template <>
  Format<Variable::phiT, Process::gp>::Format(const Setup* setup) : DataFormat(true) {
    range_ = 2. * M_PI / setup->numRegions();
    width_ = ceil(log2(setup->numSectorsPhi()));
    base_ = range_ / pow(2., width_);
  }
  template <>
  Format<Variable::zT, Process::gp>::Format(const Setup* setup) : DataFormat(true) {
    range_ = 2. * sinh(setup->maxEta()) * setup->chosenRofZ();
    base_ = range_ / setup->numSectorsEta();
    width_ = ceil(log2(setup->numSectorsEta()));
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
  Format<Variable::inv2R, Process::ht>::Format(const Setup* setup) : DataFormat(true) {
    range_ = 2. * setup->invPtToDphi() / setup->minPt();
    base_ = range_ / (double)setup->htNumBinsInv2R();
    width_ = ceil(log2(setup->htNumBinsInv2R()));
  }
  template <>
  Format<Variable::phiT, Process::ht>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phiT, Process::gp> phiT(setup);
    range_ = phiT.range();
    base_ = phiT.base() / (double)setup->htNumBinsPhiT();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::bit, Process::ht>::Format(const Setup* setup) : DataFormat(false) {
    width_ = 1;
  }

  template <>
  Format<Variable::phi, Process::kfin>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi, Process::ht> phi(setup);
    range_ = phi.range() * pow(2, setup->kfinShiftRangePhi());
    base_ = phi.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::z, Process::kfin>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::gp> z(setup);
    range_ = z.range() * pow(2, setup->kfinShiftRangeZ());
    base_ = z.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::dPhi, Process::kfin>::Format(const Setup* setup) : DataFormat(false) {
    const Format<Variable::phi, Process::dtc> phi(setup);
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    range_ = setup->pitchRowPS() / setup->innerRadius() + (setup->pitchCol2S() + setup->scattering()) * inv2R.limit() + phi.base();
    base_ = phi.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::dZ, Process::kfin>::Format(const Setup* setup) : DataFormat(false) {
    const Format<Variable::z, Process::dtc> z(setup);
    range_ = setup->pitchCol2S() * setup->maxCot() + z.base();
    base_ = z.base();
    width_ = ceil(log2(range_ / base_));
  }

  template <>
  Format<Variable::phi, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi, Process::kfin> phi(setup);
    const Format<Variable::inv2R, Process::kf> inv2R(setup);
    const Format<Variable::phiT, Process::kf> phiT(setup);
    const Format<Variable::dPhi, Process::kfin> dPhi(setup);
    range_ = (phiT.range() + setup->maxRphi() * inv2R.range()) * setup->kfRangeFactor() + dPhi.range();
    base_ = phi.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::z, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::kfin> z(setup);
    const Format<Variable::dZ, Process::kfin> dZ(setup);
    const Format<Variable::zT, Process::kf> zT(setup);
    const Format<Variable::cot, Process::kf> cot(setup);
    range_ = (zT.range() + cot.range() * setup->maxRz()) * setup->kfRangeFactor() + dZ.range();
    base_ = z.base();
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::inv2R, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::inv2R, Process::dr> dr(setup);
    const Format<Variable::inv2R, Process::ht> ht(setup);
    range_ = ht.range();
    base_ = dr.base() * pow(2., setup->kfBaseShift());
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::phiT, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi0, Process::dr> dr(setup);
    const Format<Variable::phiT, Process::ht> ht(setup);
    range_ = ht.range();
    base_ = dr.base() * pow(2., setup->kfBaseShift());
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::cot, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::cot, Process::dr> dr(setup);
    range_ = dr.range();
    base_ = dr.base() * pow(2., setup->kfBaseShift());
    width_ = ceil(log2(range_ / base_));
  }
  template <>
  Format<Variable::zT, Process::kf>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z0, Process::dr> z0(setup);
    const Format<Variable::zT, Process::gp> zT(setup);
    range_ = zT.range();
    base_ = z0.base() * pow(2., setup->kfBaseShift());
    width_ = ceil(log2(range_ / base_));
  }

  template <>
  Format<Variable::inv2R, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::inv2R, Process::ht> inv2R(setup);
    width_ = setup->tfpWidthInv2R();
    range_ = inv2R.range();
    base_ = inv2R.base();
    const int shift = ceil(log2(range_ / base_)) - width_;
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::phi0, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::phi, Process::dtc> phi(setup);
    width_ = setup->tfpWidthPhi0();
    range_ = phi.range();
    base_ = phi.base();
    const int shift = ceil(log2(range_ / base_)) - width_;
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::cot, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::zT, Process::gp> zT(setup);
    width_ = setup->tfpWidthCot();
    range_ = (zT.range() + 2. * setup->beamWindowZ()) / setup->chosenRofZ();
    base_ = (zT.base() + 2. * setup->beamWindowZ()) / setup->chosenRofZ();
    const int shift = ceil(log2(range_ / base_)) - width_;
    base_ *= pow(2., shift);
  }
  template <>
  Format<Variable::z0, Process::dr>::Format(const Setup* setup) : DataFormat(true) {
    const Format<Variable::z, Process::dtc> z(setup);
    width_ = setup->tfpWidthZ0();
    range_ = 2. * setup->beamWindowZ();
    base_ = z.base();
    const int shift = ceil(log2(range_ / base_)) - width_;
    base_ *= pow(2., shift);
  }

}  // namespace trackerTFP
