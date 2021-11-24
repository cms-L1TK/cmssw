// This is a helper function that can be used to decode hitpattern, which is a 7-bit integer produced by the Kalman filter.
//
// There are three classes declared in HitPatternHelper (hph) namesapce:
// 1)SensorModule: This is used to store important information about the sensor modules. (e.g. r,z coordinates of the tracker modules)
// 2)Setup: This is used to produce a collection of <SensorModule> needed by HitPatternHelper.
// 3)HitPatternHelper: This function returns more specific information (e.g. module type, layer id,...) about each stub on the TTTrack objects.
// This function needs three variables from TTTrack: hitPattern(),tanL() and z0().
// It makes predictions in two different ways depending on which version of the Kalman filter is deployed:
//
// Old KF (L1Trigger/TrackFindingTMTT/plugins/TMTrackProducer.cc)
// With this version of the KF, HitPatternHelper relys on a hard-coded layermap to deterimne layer id of each stub.
//
// New KF (L1Trigger/TrackerTFP/plugins/ProducerKF.cc)
// With this version of the KF, HitPatternHelper makes predictions based on the spatial coordinates of tracks and detector modules.
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
#include "L1Trigger/TrackTrigger/interface/Setup.h"

#include <bitset>
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;

namespace hph {

  //Class that stores configuration for HitPatternHelper
  class Setup {
  public:
    Setup() {}
    Setup(const edm::ParameterSet& iConfig, const tt::Setup& setupTT);
    ~Setup() {}

    bool hphDebug() const { return iConfig_.getParameter<bool>("hphDebug"); }
    bool useNewKF() const { return iConfig_.getParameter<bool>("useNewKF"); }
    double deltaTanL() const { return iConfig_.getParameter<double>("deltaTanL"); }
    double chosenRofZ() const { return setupTT_.chosenRofZ(); }
    std::vector<double> etaRegions() const { return setupTT_.boundarieEta(); }
    std::vector<tt::SensorModule> sensorModules() const { return setupTT_.sensorModules(); }

  private:
    edm::ParameterSet iConfig_;
    const tt::Setup setupTT_;
  };

  //Class that returns decoded information from hitpattern
  class HitPatternHelper {
  public:
    HitPatternHelper() {}
    HitPatternHelper(const Setup* setup, int hitpattern, double cot, double z0);
    ~HitPatternHelper() {}

    int etaSector() { return etaSector_; }        //Eta sectors defined in KF
    int numExpLayer() { return numExpLayer_; }    //The number of layers KF expects
    int numMissingPS() { return numMissingPS_; }  //The number of PS layers that are missing
    int numMissing2S() { return numMissing2S_; }  //The number of 2S layers that are missing
    int numPS() { return numPS_; }                //The number of PS layers are found in hitpattern
    int num2S() { return num2S_; }                //The number of 2S layers are found in hitpattern
    int numMissingInterior1() {
      return numMissingInterior1_;
    }  //The number of missing interior layers (using only hitpattern)
    int numMissingInterior2() {
      return numMissingInterior2_;
    }  //The number of missing interior layers (using hitpattern and sensor modules)
    std::vector<int> binary() { return binary_; }  //11-bit hitmask needed by TrackQuality.cc (0~5->L1~L6;6~10->D1~D5)
    static auto smallerID(tt::SensorModule lhs, tt::SensorModule rhs) { return lhs.layerId() < rhs.layerId(); }
    static auto equalID(tt::SensorModule lhs, tt::SensorModule rhs) { return lhs.layerId() == rhs.layerId(); }

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
    const Setup* setup_;
    std::vector<tt::SensorModule> layers_;  //Sensor modules that particles are expected to hit
    std::vector<int> binary_;
    bool hphDebug_;
    bool useNewKF_;
    float chosenRofZ_;
    float deltaTanL_;
    std::vector<double> etaRegions_;

    //Layermap used in Old KF
    //Ultimate config is assumed (with maybe layer)
    //Index across is kalman layer
    //Index down is eta sector
    //Element is layer id where barrel layers=1,2,3,4,5,6 & endcap wheels=11,12,13,14,15; 0 is invalid.
    std::vector<int> layermap_[8][7] = {
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

}  // namespace hph

EVENTSETUP_DATA_DEFAULT_RECORD(hph::Setup, hph::SetupRcd);

#endif
