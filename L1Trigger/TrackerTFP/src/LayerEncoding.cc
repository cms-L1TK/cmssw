#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"

#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>

using namespace std;
using namespace tt;

namespace trackerTFP {

  LayerEncoding::LayerEncoding(const DataFormats* dataFormats)
      : setup_(dataFormats->setup()),
        dataFormats_(dataFormats),
        zT_(&dataFormats->format(Variable::zT, Process::gp)),
        layerEncoding_(vector<vector<int>>(pow(2, zT_->width()))),
        maybePattern_(vector<TTBV>(pow(2, zT_->width()), TTBV(0, setup_->numLayers()))) {
    // number of boundaries of fiducial area in r-z plane for a given set of rough r-z track parameter
    static constexpr int boundaries = 2;
    // find unique sensor mouldes in r-z
    // allowed distance in r and z in cm between modules to consider them not unique
    static constexpr double delta = 1.e-3;
    vector<const SensorModule*> sensorModules;
    sensorModules.reserve(setup_->sensorModules().size());
    for (const SensorModule& sm : setup_->sensorModules())
      sensorModules.push_back(&sm);
    auto smallerR = [](const SensorModule* lhs, const SensorModule* rhs) { return lhs->r() < rhs->r(); };
    auto smallerZ = [](const SensorModule* lhs, const SensorModule* rhs) { return lhs->z() < rhs->z(); };
    auto equalRZ = [](const SensorModule* lhs, const SensorModule* rhs) {
      return abs(lhs->r() - rhs->r()) < delta && abs(lhs->z() - rhs->z()) < delta;
    };
    stable_sort(sensorModules.begin(), sensorModules.end(), smallerZ);
    stable_sort(sensorModules.begin(), sensorModules.end(), smallerR);
    sensorModules.erase(unique(sensorModules.begin(), sensorModules.end(), equalRZ), sensorModules.end());
    stable_sort(sensorModules.begin(), sensorModules.end(), smallerZ);
    sensorModules.erase(unique(sensorModules.begin(), sensorModules.end(), equalRZ), sensorModules.end());
    // find set of moudles for each set of rough r-z track parameter
    // loop over zT bins
    for (int binZT = 0; binZT < pow(2, zT_->width()); binZT++) {
      // z at radius chosenRofZ
      const double zT = zT_->floating(zT_->toSigned(binZT));
      // cotTheta of eta sector centre
      const double cot = zT / setup_->chosenRofZ();
      // cot uncertainty
      const double dCot = (zT_->base() / 2. + setup_->beamWindowZ()) / setup_->chosenRofZ();
      // z at radius chosenRofZ wrt zT of sectorZT of this bin boundaries
      const vector<double> zTs = {zT - zT_->base() / 2., zT + zT_->base() / 2.};
      // cotTheta wrt sectorCot of this bin boundaries
      const vector<double> cots = {cot - dCot, cot + dCot};
      // layer ids crossed by left and right rough r-z parameter shape boundaries
      vector<set<int>> layers(boundaries);
      // loop over all unique modules
      for (const SensorModule* sm : sensorModules) {
        // check if module is crossed by left and right rough r-z parameter shape boundaries
        for (int i = 0; i < boundaries; i++) {
          const double zTi = zTs[i];
          for (double coti : cots) {
            // distance between module and boundary in moudle tilt angle direction
            const double d =
                (zTi - sm->z() + (sm->r() - setup_->chosenRofZ()) * coti) / (sm->cosTilt() - sm->sinTilt() * coti);
            // compare distance with module size and add module layer id to layers if module is crossed
            if (abs(d) < sm->numColumns() * sm->pitchCol() / 2.)
              layers[i].insert(sm->layerId());
          }
        }
      }
      // mayber layers are given by layer ids crossed by only one booundary
      set<int> maybeLayer;
      set_symmetric_difference(layers[0].begin(),
                                layers[0].end(),
                                layers[1].begin(),
                                layers[1].end(),
                                inserter(maybeLayer, maybeLayer.end()));
      // layerEncoding is given by sorted layer ids crossed by any booundary
      set<int> layerEncoding;
      set_union(layers[0].begin(),
                layers[0].end(),
                layers[1].begin(),
                layers[1].end(),
                inserter(layerEncoding, layerEncoding.end()));
      // fill layerEncoding_
      vector<int>& le = layerEncoding_[binZT];
      le = vector<int>(layerEncoding.begin(), layerEncoding.end());
      // fill maybePattern_
      TTBV& mp = maybePattern_[binZT];
      for (int m : maybeLayer)
        mp.set(min((int)distance(le.begin(), find(le.begin(), le.end(), m)), setup_->numLayers() - 1));
      /*if (zT_->toSigned(binZT) == 4) {
        cout << mp << endl;
        for (int i : le)
          cout << i << " ";
        cout << endl;
      }*/
    }
    const bool print = false;
    if (!print)
      return;
    static constexpr int widthLayer = 3;
    stringstream ssLE;
    stringstream ssMP;
    for (int binZT = 0; binZT < pow(2, zT_->width()); binZT++) {
      const vector<int>& le = layerEncoding_[binZT];
      /*if (zT_->toSigned(binZT) == 11) {
        for (int i : le)
          cout << i << " ";
        cout << endl;
      }*/
      const TTBV& mp = maybePattern_[binZT];
      ssMP << mp << endl;
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        vector<int> layerIds;
        if (layer == 0)
          layerIds = {1, 16};
        else if (layer == 1)
          layerIds = {2, 16};
        else if (layer == 2)
          layerIds = {6, 11};
        else if (layer == 3)
          layerIds = {5, 12};
        else if (layer == 4)
          layerIds = {4, 13};
        else if (layer == 5)
          layerIds = {0, 14};
        else if (layer == 6)
          layerIds = {3, 15};
        for (int layerId : layerIds) {
          const auto it = find(le.begin(), le.end(), layerId);
          const bool valid = it != le.end();
          const int kfLayerId = min((int)distance(le.begin(), it), setup_->numLayers() - 1);
          /*if (zT_->toSigned(binZT) == 11) {
            cout << layer << " " << valid << " " << kfLayerId << endl;
          }*/
          if (valid)
            ssLE << "1" << TTBV(kfLayerId, widthLayer);
          else
            ssLE << "0" << string(widthLayer, '-');
        }
      }
      //if (zT_->toSigned(binZT) == 4)
        //throw cms::Exception("...");
      ssLE << endl;
    }
    fstream file;
    file.open("layerEncoding.mem", ios::out);
    file << ssLE.rdbuf();
    file.close();
    file.open("maybePatterns.mem", ios::out);
    file << ssMP.rdbuf();
    file.close();
  }

  // Set of layers for given bin in zT
  const vector<int>& LayerEncoding::layerEncoding(int zT) const {
    const int binZT = zT_->toUnsigned(zT);
    return layerEncoding_.at(binZT);
  }

  // pattern of maybe layers for given bin in zT
  const TTBV& LayerEncoding::maybePattern(int zT) const {
    const int binZT = zT_->toUnsigned(zT);
    return maybePattern_.at(binZT);
  }

}  // namespace trackerTFP
