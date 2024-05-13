#include "L1Trigger/TrackerTFP/interface/KalmanFilterFormats.h"

#include <vector>
#include <deque>
#include <cmath>
#include <tuple>
#include <iterator>
#include <algorithm>
#include <limits>
#include <cstring>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  constexpr auto variableKFstrs_ = {
      "x0",           "x1",           "x2",        "x3",        "H00",    "H12",    "m0",       "m1",
      "v0",           "v1",           "r0",        "r1",        "S00",    "S01",    "S12",      "S13",
      "K00",          "K10",          "K21",       "K31",       "R00",    "R11",    "R00Rough", "R11Rough",
      "invR00Approx", "invR11Approx", "invR00Cor", "invR11Cor", "invR00", "invR11", "C00",      "C01",
      "C11",          "C22",          "C23",       "C33",       "r02",    "r12",    "chi20",    "chi21"};

  void KalmanFilterFormats::endJob() {
    const int wName =
        strlen(*max_element(variableKFstrs_.begin(), variableKFstrs_.end(), [](const auto& a, const auto& b) {
          return strlen(a) < strlen(b);
        }));
    for (VariableKF v = VariableKF::begin; v != VariableKF::end; v = VariableKF(+v + 1)) {
      double min = 1.e12;
      double abs = 1.e12;
      double max = -1.e12;
      vector<int> deltas;
      deltas.reserve(setup_->numLayers());
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        min = std::min(min, format(v, layer).min());
        abs = std::min(abs, format(v, layer).abs());
        max = std::max(max, format(v, layer).max());
        const double r = format(v, layer).twos()
                             ? std::max(std::abs(format(v, layer).min()), std::abs(format(v, layer).max())) * 2.
                             : format(v, layer).max();
        deltas.emplace_back(format(v, layer).width() - ceil(log2(r / format(v, layer).base())));
      }
      cout << setw(wName) << *next(variableKFstrs_.begin(), +v) << ": ";
      for (int delta : deltas)
        cout << setw(3) << (delta == -2147483648 ? "-" : to_string(delta)) << " ";
      cout << "| " << setw(14) << min << " " << setw(14) << max << " " << setw(14) << abs << " " << setw(14)
           << format(v, 0).base() << endl;
    }
  }

  KalmanFilterFormats::KalmanFilterFormats(const ParameterSet& iConfig, const DataFormats* dataFormats)
      : iConfig_(iConfig), dataFormats_(dataFormats), setup_(dataFormats_->setup()) {
    formats_.reserve(+VariableKF::end);
    fillFormats();
  }

  template <VariableKF it>
  void KalmanFilterFormats::fillFormats() {
    vector<DataFormatKF> v;
    v.reserve(8);
    fillFormats<it>(v);
    formats_.push_back(v);
    if constexpr (++it != VariableKF::end)
      fillFormats<++it>();
  }

  template <VariableKF it, int layer>
  void KalmanFilterFormats::fillFormats(vector<DataFormatKF>& v) {
    v.emplace_back(FormatKF<it>(dataFormats_, iConfig_, layer));
    if constexpr (layer < 7)
      fillFormats<it, layer + 1>(v);
  }

  DataFormatKF::DataFormatKF(const VariableKF& v, int layer, bool twos)
      : v_(v),
        layer_(layer),
        twos_(twos),
        width_(0),
        base_(1.),
        range_(0.),
        min_(numeric_limits<double>::max()),
        abs_(numeric_limits<double>::max()),
        max_(numeric_limits<double>::lowest()) {}

  // returns false if data format would oferflow for this double value
  bool DataFormatKF::inRange(double d) const {
    if (twos_)
      return d >= -range_ / 2. && d < range_ / 2.;
    return d >= 0 && d < range_;
  }

  void DataFormatKF::updateRangeActual(double d) {
    min_ = std::min(min_, d);
    abs_ = std::min(abs_, std::abs(d));
    max_ = std::max(max_, d);
    if (!inRange(d)) {
      string v = *next(variableKFstrs_.begin(), +v_);
      cms::Exception exception("out_of_range");
      exception.addContext("trackerTFP:DataFormatKF::updateRangeActual");
      exception << "Variable " << v << " = " << d << " in layer " << to_string(layer_) << " is out of range "
                << (twos_ ? -range_ / 2. : 0) << " to " << (twos_ ? range_ / 2. : range_) << "." << endl;
      if (twos_ || d >= 0.)
        exception.addAdditionalInfo("Consider raising BaseShift" + v + " for layer " + to_string(layer_) +
                                    " in KalmnaFilterFormats_cfi.py.");
      throw exception;
    }
  }

  template <>
  FormatKF<VariableKF::x0>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::x0, layer, true) {
    const DataFormat& input = dataFormats->format(Variable::inv2R, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftx0")[layer];
    base_ = pow(2, baseShift) * input.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::x1>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::x1, layer, true) {
    const DataFormat& input = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftx1")[layer];
    base_ = pow(2, baseShift) * input.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::x2>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::x2, layer, true) {
    const DataFormat& input = dataFormats->format(Variable::cot, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftx2")[layer];
    base_ = pow(2, baseShift) * input.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::x3>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::x3, layer, true) {
    const DataFormat& input = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftx3")[layer];
    base_ = pow(2, baseShift) * input.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::H00>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::H00, layer, true) {
    const DataFormat& ctb = dataFormats->format(Variable::r, Process::ctb);
    base_ = ctb.base();
    width_ = ctb.width();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::H12>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::H12, layer, true) {
    const Setup* setup = dataFormats->setup();
    const DataFormat& ctb = dataFormats->format(Variable::r, Process::ctb);
    base_ = ctb.base();
    range_ = 2. * setup->maxRz();
    width_ = ceil(log2(range_ / base_));
    calcRange();
  }

  template <>
  FormatKF<VariableKF::m0>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::m0, layer, true) {
    const DataFormat& ctb = dataFormats->format(Variable::phi, Process::ctb);
    base_ = ctb.base();
    width_ = ctb.width();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::m1>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::m1, layer, true) {
    const DataFormat& ctb = dataFormats->format(Variable::z, Process::ctb);
    base_ = ctb.base();
    width_ = ctb.width();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::v0>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::v0, layer, false) {
    const DataFormat& ctb = dataFormats->format(Variable::dPhi, Process::ctb);
    base_ = ctb.base() * ctb.base();
    width_ = ctb.width() + ctb.width();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::v1>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::v1, layer, false) {
    const DataFormat& ctb = dataFormats->format(Variable::dZ, Process::ctb);
    base_ = ctb.base() * ctb.base();
    width_ = ctb.width() + ctb.width();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::r0>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::r0, layer, true) {
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftr0")[layer];
    base_ = pow(2., baseShift) * x1.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::r1>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::r1, layer, true) {
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftr1")[layer];
    base_ = pow(2., baseShift) * x3.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::S00>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::S00, layer, true) {
    const DataFormat& x0 = dataFormats->format(Variable::inv2R, Process::kf);
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftS00")[layer];
    base_ = pow(2., baseShift) * x0.base() * x1.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::S01>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::S01, layer, true) {
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftS01")[layer];
    base_ = pow(2., baseShift) * x1.base() * x1.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::S12>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::S12, layer, true) {
    const DataFormat& x2 = dataFormats->format(Variable::cot, Process::kf);
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftS12")[layer];
    base_ = pow(2., baseShift) * x2.base() * x3.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::S13>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::S13, layer, true) {
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftS13")[layer];
    base_ = pow(2., baseShift) * x3.base() * x3.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::K00>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::K00, layer, true) {
    const DataFormat& x0 = dataFormats->format(Variable::inv2R, Process::kf);
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftK00")[layer];
    base_ = pow(2., baseShift) * x0.base() / x1.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::K10>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::K10, layer, true) {
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftK10")[layer];
    base_ = pow(2., baseShift);
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::K21>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::K21, layer, true) {
    const DataFormat& x2 = dataFormats->format(Variable::cot, Process::kf);
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftK21")[layer];
    base_ = pow(2., baseShift) * x2.base() / x3.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::K31>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::K31, layer, true) {
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftK31")[layer];
    base_ = pow(2., baseShift);
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::R00>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::R00, layer, false) {
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftR00")[layer];
    base_ = pow(2., baseShift) * x1.base() * x1.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::R11>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::R11, layer, false) {
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftR11")[layer];
    base_ = pow(2., baseShift) * x3.base() * x3.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::R00Rough>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::R00Rough, layer, false) {
    const FormatKF<VariableKF::R00> R00(dataFormats, iConfig, layer);
    width_ = dataFormats->setup()->widthAddrBRAM18();
    range_ = R00.range();
    const int baseShift = R00.width() - width_;
    base_ = pow(2., baseShift) * R00.base();
  }

  template <>
  FormatKF<VariableKF::R11Rough>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::R11Rough, layer, false) {
    const FormatKF<VariableKF::R11> R11(dataFormats, iConfig, layer);
    width_ = dataFormats->setup()->widthAddrBRAM18();
    range_ = R11.range();
    const int baseShift = R11.width() - width_;
    base_ = pow(2., baseShift) * R11.base();
  }

  template <>
  FormatKF<VariableKF::invR00Approx>::FormatKF(const DataFormats* dataFormats,
                                               const edm::ParameterSet& iConfig,
                                               int layer)
      : DataFormatKF(VariableKF::invR00Approx, layer, false) {
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftInvR00Approx")[layer];
    base_ = pow(2., baseShift) / x1.base() / x1.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::invR11Approx>::FormatKF(const DataFormats* dataFormats,
                                               const edm::ParameterSet& iConfig,
                                               int layer)
      : DataFormatKF(VariableKF::invR11Approx, layer, false) {
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftInvR11Approx")[layer];
    base_ = pow(2., baseShift) / x3.base() / x3.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::invR00Cor>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::invR00Cor, layer, false) {
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftInvR00Cor")[layer];
    base_ = pow(2., baseShift);
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::invR11Cor>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::invR11Cor, layer, false) {
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftInvR11Cor")[layer];
    base_ = pow(2., baseShift);
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::invR00>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::invR00, layer, false) {
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftInvR00")[layer];
    base_ = pow(2., baseShift) / x1.base() / x1.base();
    width_ = dataFormats->setup()->widthDSPau();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::invR11>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::invR11, layer, false) {
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftInvR11")[layer];
    base_ = pow(2., baseShift) / x3.base() / x3.base();
    width_ = dataFormats->setup()->widthDSPau();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::C00>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::C00, layer, false) {
    const DataFormat& x0 = dataFormats->format(Variable::inv2R, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftC00")[layer];
    base_ = pow(2., baseShift) * x0.base() * x0.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::C01>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::C01, layer, true) {
    const DataFormat& x0 = dataFormats->format(Variable::inv2R, Process::kf);
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftC01")[layer];
    base_ = pow(2., baseShift) * x0.base() * x1.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::C11>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::C11, layer, false) {
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftC11")[layer];
    base_ = pow(2., baseShift) * x1.base() * x1.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::C22>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::C22, layer, false) {
    const DataFormat& x2 = dataFormats->format(Variable::cot, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftC22")[layer];
    base_ = pow(2., baseShift) * x2.base() * x2.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::C23>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::C23, layer, true) {
    const DataFormat& x2 = dataFormats->format(Variable::cot, Process::kf);
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftC23")[layer];
    base_ = pow(2., baseShift) * x2.base() * x3.base();
    width_ = dataFormats->setup()->widthDSPbb();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::C33>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::C33, layer, false) {
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftC33")[layer];
    base_ = pow(2., baseShift) * x3.base() * x3.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::r02>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::r02, layer, false) {
    const DataFormat& x1 = dataFormats->format(Variable::phiT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftr02")[layer];
    base_ = pow(2., baseShift) * x1.base() * x1.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::r12>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::r12, layer, false) {
    const DataFormat& x3 = dataFormats->format(Variable::zT, Process::kf);
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftr12")[layer];
    base_ = pow(2., baseShift) * x3.base() * x3.base();
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::chi20>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::chi20, layer, false) {
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftchi20")[layer];
    base_ = pow(2., baseShift);
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

  template <>
  FormatKF<VariableKF::chi21>::FormatKF(const DataFormats* dataFormats, const edm::ParameterSet& iConfig, int layer)
      : DataFormatKF(VariableKF::chi21, layer, false) {
    const int baseShift = iConfig.getParameter<vector<int>>("BaseShiftchi21")[layer];
    base_ = pow(2., baseShift);
    width_ = dataFormats->setup()->widthDSPbu();
    calcRange();
  }

}  // namespace trackerTFP
