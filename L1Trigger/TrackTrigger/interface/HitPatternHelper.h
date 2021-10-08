// This is a helper function that can be used to decode hitpattern, which is a 7-bit integer produced by KF.
//
// There are three classes declared in HitPatternHelper (HPH) namesapce:
// 1)SensorModule: This is used to store important information about the sensor modules. For example r,z coordinates.
// 2)Setup: This is used to produce a collection of <SensorModule> needed by HPH.
// 3)HitPatternHelper: This is used to decode hitpattern with the help of the information from sensor modules and layermap.
//
// Two predictions on which layers particles will hit are made using different information:
// i)Loop over sensor modules and make predictions based on spatial coordinates of tracks. This prediction is considered more accurate.
// ii)Make predictions based on a hard-coded layermap. This prediction is considered less accurate and is used by Old KF to encode hitpattern.
//
//
//  Created by J.Li on 1/23/21.
//

#ifndef L1Trigger_TrackTrigger_interface_HitPatternHelper_h
#define L1Trigger_TrackTrigger_interface_HitPatternHelper_h

#include "FWCore/Framework/interface/data_default_record_trait.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "L1Trigger/TrackTrigger/interface/HitPatternHelperRcd.h"

#include <bitset>
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;

namespace HPH {

  class SensorModule {
  public:
    SensorModule() {}
    SensorModule(
        bool isbarrel, bool isPS, int numColumns, int layerid, double r, double z, double pitchCol, double tilt);
    ~SensorModule() {}

    bool isbarrel() const { return isbarrel_; }
    bool isPS() const { return isPS_; }
    bool isMaybe() const { return isMaybe_; }
    int numColumns() const { return numColumns_; }
    int layerid() const { return layerid_; }
    double r() const { return r_; }
    double z() const { return z_; }
    double pitchCol() const { return pitchCol_; }
    double tilt() const { return tilt_; }
    double sin() const { return sin_; }
    double cos() const { return cos_; }

    void setMaybe() { isMaybe_ = true; }

  private:
    bool isbarrel_;
    bool isPS_;
    bool isMaybe_;
    int numColumns_;
    int layerid_;
    double r_;
    double z_;
    double pitchCol_;
    double tilt_;
    double sin_;
    double cos_;
  };

  class Setup {
  public:
    Setup() {}
    Setup(const edm::ParameterSet& iConfig,
          const TrackerGeometry& trackerGeometry,
          const TrackerTopology& trackerTopology);
    ~Setup() {}

    static auto smallerR(SensorModule lhs, SensorModule rhs) { return lhs.r() < rhs.r(); }
    static auto smallerZ(SensorModule lhs, SensorModule rhs) { return lhs.z() < rhs.z(); }
    static auto equalRZ(SensorModule lhs, SensorModule rhs) {
      return abs(lhs.r() - rhs.r()) < delta_ && abs(lhs.z() - rhs.z()) < delta_;
    }
    std::vector<SensorModule> SensorModules() const { return SensorModules_; }

    bool HPHdebug() const { return iConfig_.getParameter<bool>("HPHdebug"); }
    bool useNewKF() const { return iConfig_.getParameter<bool>("useNewKF"); }
    double chosenRofZ() const { return iConfig_.getParameter<double>("chosenRofZ"); }
    double deltaTanL() const { return iConfig_.getParameter<double>("deltaTanL"); }

  private:
    edm::ParameterSet iConfig_;
    const TrackerGeometry* trackerGeometry_;
    const TrackerTopology* trackerTopology_;
    static constexpr double delta_ = 1.e-3;
    std::vector<SensorModule> SensorModules_;
  };

  class HitPatternHelper {
  public:
    HitPatternHelper() {}
    HitPatternHelper(const Setup* setup, int hitpattern, double cot, double z0);
    ~HitPatternHelper() {}

    int getetaSector() { return etaSector_; }        //Eta sectors defined in KF
    int getnumExpLayer() { return numExpLayer_; }    //The number of layers KF expects
    int getnumMissingPS() { return numMissingPS_; }  //The number of PS layers that are missing
    int getnumMissing2S() { return numMissing2S_; }  //The number of 2S layers that are missing
    int getnumPS() { return numPS_; }                //The number of PS layers are found in hitpattern
    int getnum2S() { return num2S_; }                //The number of 2S layers are found in hitpattern
    int getnumMissingInterior1() {
      return numMissingInterior1_;
    }  //The number of missing interior layers (using only hitpattern)
    int getnumMissingInterior2() {
      return numMissingInterior2_;
    }  //The number of missing interior layers (using hitpattern and sensor modules)
    std::vector<int> getbinary() { return binary_; }  //11-bit hitmask needed by TrackQuality.cc (0~5->L1~L6;6~10->D1~D5)
    static auto smallerID(SensorModule lhs, SensorModule rhs) { return lhs.layerid() < rhs.layerid(); }
    static auto equalID(SensorModule lhs, SensorModule rhs) { return lhs.layerid() == rhs.layerid(); }

    int ReducedId(
        int layerId);  //Converts layer id (1~6->L1~L6;11~15->D1~D5) to reduced layer id (0~5->L1~L6;6~10->D1~D5)
    int findLayer(int layerId);  //Search for a layer id from sensor modules

  private:
    int etaSector_;
    int hitpattern_;
    int numExpLayer_;
    int numMissingLayer_;
    int numMissingPS_;
    int numMissing2S_;
    int numPS_;
    int num2S_;
    int numMissingInterior1_;
    int numMissingInterior2_;
    double cot_;
    double z0_;
    const Setup* Setup_;
    std::vector<SensorModule> layers_;  //Sensor modules that particles are expected to hit
    std::vector<int> binary_;
    bool HPHdebug_;
    bool useNewKF_;
    float chosenRofZ_;
    float deltaTanL_;
    std::vector<float> etaRegions_ = {
        -2.4, -2.08, -1.68, -1.26, -0.90, -0.62, -0.41, -0.20, 0.0, 0.20, 0.41, 0.62, 0.90, 1.26, 1.68, 2.08, 2.4};

    //Layermap used in Old KF
    //Ultimate config is assumed (with maybe layer)
    //Index across is kalman layer
    //Index down is eta sector
    //Element is layer id where barrel layers=1,2,3,4,5,6 & endcap wheels=11,12,13,14,15; 0 is invalid.
    std::vector<int> hitmap_[8][7] = {
        {{1}, {2}, {3}, {4}, {5}, {6}, {0}},
        {{1}, {2}, {3}, {4}, {5}, {6}, {0}},
        {{1}, {2}, {3}, {4}, {5}, {6}, {0}},
        {{1}, {2}, {3}, {4}, {5}, {6}, {0}},
        {{1}, {2}, {3}, {4}, {5, 11}, {6, 12}, {13}},
        {{1}, {2}, {3, 4}, {11}, {12}, {13}, {14, 15}},
        {{1}, {2}, {11}, {12}, {13}, {14}, {15}},
        {{1}, {11}, {12}, {13}, {14}, {15}, {0}},
    };
  };

}  // namespace HPH

EVENTSETUP_DATA_DEFAULT_RECORD(HPH::Setup, HPH::SetupRcd);

#endif
