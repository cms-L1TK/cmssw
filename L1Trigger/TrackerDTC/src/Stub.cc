#include "L1Trigger/TrackerDTC/interface/Stub.h"

#include <cmath>
#include <iterator>
#include <algorithm>
#include <utility>

using namespace edm;
using namespace std;
using namespace tt;
using namespace trackerTFP;

namespace trackerDTC {

  Stub::Stub(const DataFormats* dataFormats, const SensorModule* sm, const TTStubRef& ttStubRef)
      : setup_(dataFormats->setup()),
        dataFormats_(dataFormats),
        layerEncoding_(nullptr),
        sm_(sm),
        ttStubRef_(ttStubRef),
        hybrid_(false),
        valid_(true),
        regions_(0, setup_->numOverlappingRegions()) {
    const DataFormat& dfR = dataFormats_->format(Variable::r, Process::dtc);
    const DataFormat& dfPhi = dataFormats_->format(Variable::phi, Process::dtc);
    const DataFormat& dfZ = dataFormats_->format(Variable::z, Process::dtc);
    const DataFormat& dfInv2R = dataFormats_->format(Variable::inv2R, Process::ht);
    // get stub local coordinates
    const MeasurementPoint& mp = ttStubRef->clusterRef(0)->findAverageLocalCoordinatesCentered();
    // convert to uniformed local coordinates
    // column number in pitch units
    col_ = (int)floor(pow(-1, sm_->signCol()) * (mp.y() - sm_->numColumns() / 2) / setup_->baseCol());
    // row number in half pitch units
    row_ = (int)floor(pow(-1, sm_->signRow()) * (mp.x() - sm_->numRows() / 2) / setup_->baseRow());
    // bend number in quarter pitch units
    bend_ = (int)floor(pow(-1, sm_->signBend()) * (ttStubRef->bendBE()) / setup_->baseBend());
    // reduced row number for look up
    rowLUT_ = (int)floor((double)row_ / pow(2., setup_->widthRow() - setup_->dtcWidthRowLUT()));
    // sub row number inside reduced row number
    rowSub_ = row_ - (rowLUT_ + .5) * pow(2, setup_->widthRow() - setup_->dtcWidthRowLUT());
    // convert local to global coordinates
    const double y = (col_ + .5) * setup_->baseCol() * sm_->pitchCol();
    // radius of a column of strips/pixel in cm
    d_ = sm_->r() + y * sm_->sinTilt();
    // stub z in cm
    z_ = dfZ.digi(sm_->z() + y * sm_->cosTilt());
    const double x = (rowLUT_ + .5) * setup_->baseRow() * setup_->dtcNumMergedRows() * sm_->pitchRow();
    // stub r wrt chosen RofPhi in cm
    r_ = dfR.digi(sqrt(d_ * d_ + x * x) - setup_->chosenRofPhi());
    const double x0 = rowLUT_ * setup_->baseRow() * setup_->dtcNumMergedRows() * sm_->pitchRow();
    const double x1 = (rowLUT_ + 1) * setup_->baseRow() * setup_->dtcNumMergedRows() * sm_->pitchRow();
    const double phi0 = sm_->phi() + atan2(x0, d_);
    const double phi1 = sm_->phi() + atan2(x1, d_);
    const double c = (phi0 + phi1) / 2.;
    const double m = (phi1 - phi0) / setup_->dtcNumMergedRows();
    // intercept of linearized stub phi in rad
    c_ = digi(c, dfPhi.base() / 2.);
    // slope of linearized stub phi in rad / strip
    m_ = digi(m, setup_->dtcBaseM());
    // stub phi w.r.t. detector region centre in rad
    phi_ = dfPhi.digi(c_ + rowSub_ * m_);
    // assaign stub to processing regions
    // radial (cylindrical) component of sensor separation
    const double dr = sm_->sep() / (sm_->cosTilt() - sm_->sinTilt() * z_ / d_);
    // converts bend into inv2R in 1/cm
    const double inv2ROverBend = sm_->pitchRow() / dr / d_;
    // inv2R in 1/cm
    const double inv2R = -bend_ * setup_->baseBend() * inv2ROverBend;
    // inv2R uncertainty in 1/cm
    const double dInv2R = setup_->bendCut() * inv2ROverBend;
    inv2R_.first = dfInv2R.digi(inv2R - dInv2R);
    inv2R_.second = dfInv2R.digi(inv2R + dInv2R);
    static const double maxInv2R = dfInv2R.limit();
    // cut on pt
    if (inv2R_.first > maxInv2R || inv2R_.second < -maxInv2R)
      valid_ = false;
    else {
      inv2R_.first = max(inv2R_.first, -maxInv2R);
      inv2R_.second = min(inv2R_.second, maxInv2R);
    }
    // range of stub extrapolated phi to radius chosenRofPhi in rad
    phiT_.first = phi_ - r_ * inv2R_.first;
    phiT_.second = phi_ - r_ * inv2R_.second;
    if (phiT_.first > phiT_.second)
      swap(phiT_.first, phiT_.second);
    if (phiT_.first < 0.)
      regions_.set(0);
    if (phiT_.second >= 0.)
      regions_.set(1);
  }

  Stub::Stub(const ParameterSet& iConfig,
             const Setup* setup,
             const DataFormats* dataFormats,
             const LayerEncoding* layerEncoding,
             const SensorModule* sm,
             const TTStubRef& ttStubRef)
      : setup_(setup),
        dataFormats_(dataFormats),
        layerEncoding_(layerEncoding),
        sm_(sm),
        ttStubRef_(ttStubRef),
        hybrid_(iConfig.getParameter<bool>("UseHybrid")) {
    const Stub stub(dataFormats, sm, ttStubRef);
    valid_ = stub.valid_;
    r_ = stub.r_;
    phi_ = stub.phi_;
    z_ = stub.z_;
    phiT_ = stub.phiT_;
    inv2R_ = stub.inv2R_;
    regions_ = stub.regions_;
    // apply "eta" cut 
    const DataFormat& dfZT = dataFormats->format(Variable::zT, Process::gp);
    const double r = r_ + setup->chosenRofPhi();
    if (hybrid_) {
      if (abs(z_ / r) > setup->hybridMaxCot())
        valid_ = false;
    } else {
      const double ratioRZ = setup->chosenRofZ() / r;
      // extrapolated z at radius T assuming z0=0
      const double zT = z_ * ratioRZ;
      // extrapolated z0 window at radius T
      const double dZT = setup->beamWindowZ() * abs(1. - ratioRZ);
      if (abs(zT) > dfZT.limit() + dZT)
        valid_ = false;
    }
    // apply data format specific manipulations
    if (!hybrid_)
      return;
    // stub r w.r.t. an offset in cm
    r_ -= sm->offsetR() - setup->chosenRofPhi();
    // stub z w.r.t. an offset in cm
    z_ -= sm->offsetZ();
    if (sm->type() == SensorModule::Disk2S) {
      // encoded r
      r_ = sm->encodedR() + (sm->side() ? -col_ : (col_ + sm->numColumns() / 2));
      r_ = (r_ + 0.5) * setup->hybridBaseR(sm->type());
    }
    // encode bend
    const vector<double>& encodingBend = setup->encodingBend(sm->windowSize(), sm->psModule());
    const auto pos = find(encodingBend.begin(), encodingBend.end(), abs(ttStubRef->bendBE()));
    const int uBend = distance(encodingBend.begin(), pos);
    bend_ = pow(-1, signbit(bend_)) * uBend;
  }

  // returns bit accurate representation of Stub
  Frame Stub::frame(int region) const { return hybrid_ ? formatHybrid(region) : formatTMTT(region); }

  // returns true if stub belongs to region
  bool Stub::inRegion(int region) const { return regions_[region]; }

  // truncates double precision to f/w integer equivalent
  double Stub::digi(double value, double precision) const { return (floor(value / precision + 1.e-12) + .5) * precision; }

  // returns 64 bit stub in hybrid data format
  Frame Stub::formatHybrid(int region) const {
    const SensorModule::Type type = sm_->type();
    // layer encoding
    const int decodedLayerId = layerEncoding_->decode(sm_);
    // stub phi w.r.t. processing region border in rad
    double phi = phi_ - (region - .5) * setup_->baseRegion() + setup_->hybridRangePhi() / 2.;
    // convert stub variables into bit vectors
    const bool twosR = type == SensorModule::BarrelPS || type == SensorModule::Barrel2S;
    const TTBV hwR(r_, setup_->hybridBaseR(type), setup_->hybridWidthR(type), twosR);
    const TTBV hwPhi(phi, setup_->hybridBasePhi(type), setup_->hybridWidthPhi(type));
    const TTBV hwZ(z_, setup_->hybridBaseZ(type), setup_->hybridWidthZ(type), true);
    const TTBV hwAlpha(row_, setup_->hybridBaseAlpha(type), setup_->hybridWidthAlpha(type), true);
    const TTBV hwBend(bend_, setup_->hybridWidthBend(type), true);
    const TTBV hwLayer(decodedLayerId, setup_->hybridWidthLayerId());
    const TTBV hwGap(0, setup_->hybridNumUnusedBits(type));
    const TTBV hwValid(1, 1);
    // assemble final bitset
    return Frame(hwGap.str() + hwR.str() + hwZ.str() + hwPhi.str() + hwAlpha.str() + hwBend.str() + hwLayer.str() +
                 hwValid.str());
  }

  Frame Stub::formatTMTT(int region) const {
    return StubDTC(ttStubRef_, dataFormats_, region, sm_->layerId(), r_, phi_, z_, phiT_, inv2R_).frame().second;
  }

}  // namespace trackerDTC