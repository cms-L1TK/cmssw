#include "L1Trigger/TrackFindingTracklet/interface/Setup.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/L1TrackTrigger/interface/TTBV.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <set>

namespace trklet {

  Setup::Setup(const Config& config, const trackerDTC::Setup* setup) : dtc_(setup), config_(config) {
    numFrames_ = dtc_->sysNumFrames() + dtc_->sysNumFramesInfra() - 1;
    tbMaxCot_ = std::sinh(dtc_->regMaxEta());
    const Settings settings;
    tbBaseInv2R_ = .5 * settings.krinvpars();
    tbBasePhi0_ = settings.kphi0pars();
    tbBaseCot_ = settings.ktpars();
    tbBaseZ0_ = settings.kz0pars();
    tbBaseR_ = settings.kr();
    tbBasePhi_ = settings.kphi1();
    tbBaseZ_ = settings.kz();
    tbBaseRs_.reserve(trackerDTC::SensorModule::Type::NumTypes);
    tbBasePhis_.reserve(trackerDTC::SensorModule::Type::NumTypes);
    tbBaseRZs_.reserve(trackerDTC::SensorModule::Type::NumTypes);
    for (trackerDTC::SensorModule::Type type : {trackerDTC::SensorModule::Type::BarrelPS,
                                                trackerDTC::SensorModule::Type::Barrel2S,
                                                trackerDTC::SensorModule::Type::DiskPS})
      tbBaseRs_.emplace_back(tbBaseR_ * std::pow(2, config_.tbBaseShiftsR[type]));
    tbBaseRs_.push_back(1.0);
    for (int baseShift : config_.tbBaseShiftsPhi)
      tbBasePhis_.emplace_back(tbBasePhi_ * std::pow(2, baseShift));
    for (trackerDTC::SensorModule::Type type :
         {trackerDTC::SensorModule::Type::BarrelPS, trackerDTC::SensorModule::Type::Barrel2S})
      tbBaseRZs_.emplace_back(tbBaseZ_ * std::pow(2, config_.tbBaseShiftsRZ[type]));
    for (trackerDTC::SensorModule::Type type :
         {trackerDTC::SensorModule::Type::DiskPS, trackerDTC::SensorModule::Type::Disk2S})
      tbBaseRZs_.emplace_back(tbBaseR_ * std::pow(2, config_.tbBaseShiftsRZ[type]));
    tbWidthSeedType_ = tt::ilog2(config_.tbNumSeedTypes);
    tmNumLayers_ = config_.tbNumLayers;
    int maxProjLayers(0);
    for (const std::vector<int>& projLayers : config_.tbSeedTypesProjectionLayers)
      maxProjLayers = std::max(maxProjLayers, static_cast<int>(projLayers.size()));
    const int baseShiftInvCot = tt::ilog2(config_.tbMaxR / config_.tbMinZ) - dtc_->widthDSPbu();
    drBaseInvCot_ = std::pow(2, baseShiftInvCot);
    drBaseCot_ = tbBaseCot_ * std::pow(2, config_.tbWidthCot - setup->widthAddrBRAM18() - 1);
    drBaseInv2R_ = tbBaseInv2R_ * std::pow(2, config_.tbWidthInv2R - setup->widthAddrBRAM18() - 1);
    const int tbWidthR = tt::ilog2(dtc_->sysOuterRadius() / tbBaseR_);
    drBaseR_ = tbBaseR_ * std::pow(2, tbWidthR - setup->widthAddrBRAM18());
    drNumLayers_ = tmNumLayers_;
    drWidthLayerId_ = tt::ilog2(drNumLayers_);
    kfNumLayers_ = config_.kfNumSeed + config_.kfNumProj;
    kfMinProj_ = config_.kfMinLayers - config_.kfNumSeed;
    tmMuxOrder_.reserve(config_.tmMuxOrder.size());
    for (const std::string& seedType : config_.tmMuxOrder) {
      const auto it = std::find(config_.tbSeedTypes.begin(), config_.tbSeedTypes.end(), seedType);
      tmMuxOrder_.push_back(std::distance(config_.tbSeedTypes.begin(), it));
    }
    tqScaleFactorZ0_ = std::pow(2, config_.tqScalePowZ0);
    tqScaleFactorCot_ = std::pow(2, config_.tqScalePowCot);
  }

  // returns radial position and phi, z residuals of a TB stub
  GlobalPoint Setup::stubPosTB(const tt::FrameStub& frame, double cot) const {
    // sennsor module for this seed stub
    const trackerDTC::SensorModule* sm = dtc_->sensorModule(frame.first);
    const trackerDTC::SensorModule::Type t = sm->type();
    // parse residuals from tt::Frame and take layerId from tt::TTStubRef
    TTBV bv(frame.second);
    const double z = bv.val(tbWidthRZ(t), 0, true) * tbBaseRZ(t) * (sm->barrel() ? 1. : -cot);
    bv >>= tbWidthRZ(t);
    const double phi = bv.val(tbWidthPhi(t), 0, true) * tbBasePhi(t);
    bv >>= tbWidthPhi(t);
    double r = bv.val(tbWidthR(t), 0, sm->barrel()) * tbBaseR(t) +
               (sm->barrel() ? tt::digiR(dtc_->stubLayerR(sm->layerIndex()), tbBaseR_) : 0.0);
    if (t == trackerDTC::SensorModule::Disk2S)
      r = tt::digiR(dtc_->stubDiskR(sm->layerIndex(), r), tbBaseR_);
    return GlobalPoint(GlobalPoint::Cylindrical(r, phi, z));
  }

  // returns radial position and phi, z residuals of a TB fake seed stub
  GlobalPoint Setup::stubPosTB(const TTStubRef& ttStubRef, double cot, double z0) const {
    // sennsor module for this seed stub
    const trackerDTC::SensorModule* sm = dtc_->sensorModule(ttStubRef);
    const int layerIndex = sm->layerIndex();
    double r;
    if (sm->barrel()) {
      r = tt::digiR(dtc_->stubLayerR(layerIndex), tbBaseR_);
      const double z = tt::digi(z0 + r * cot, tbBaseZ_);
      if (std::abs(z) > tt::digi(config_.tbMinZ, tbBaseZ_) && layerIndex == 0)
        r = tt::digiR(config_.tbInnerRadius, tbBaseR_);
    } else {
      const double z = tt::digi((sm->side() ? 1. : -1.) * dtc_->stubDiskZ(layerIndex), tbBaseZ0_);
      const double invCot = tt::digi(1. / std::abs(cot), drBaseInvCot_);
      r = tt::digiR((z - z0) * invCot, tbBaseR_);
    }
    return GlobalPoint(GlobalPoint::Cylindrical(r, 0, 0));
  }

}  // namespace trklet
