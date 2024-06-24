#include "L1Trigger/TrackFindingTracklet/interface/ProjectionCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletConfigBuilder.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace std;
using namespace trklet;

ProjectionCalculator::ProjectionCalculator(string name, Settings const& settings, Globals* global)
    : ProcessBase(name, settings, global) {

      //Constants for coordinates and track parameter definitions (taken from TrackletCalculatorBase.cc)
      n_phi_ = 17;
      n_r_ = 12;
      n_z_ = 11;
      n_phi0_ = 16;
      n_rinv_ = 13;
      n_t_ = 9;
      n_phidisk_ = n_phi_-3;
      n_rdisk_ = n_r_-1;

      //Constants used for projectison to layers
      n_s_ = 12;
      n_s6_ = 14;

      //Constants used for projectison to disks
      n_tinv_ = 12;
      n_y_ = 14;
      n_x_ = 14;
      n_xx6_ = 14;

      phiHG_ = settings_.dphisectorHG();

      LUT_itinv_.resize(8192);

      for (int it = 0; it < 8192; it++) {
        if (it<100) {
          LUT_itinv_[it] =  0;
        } else {
          LUT_itinv_[it] = (1 << (n_t_ + n_tinv_)) / abs(it);
        }
      }

      for (unsigned int layerdisk = 0; layerdisk < N_LAYER + N_DISK; layerdisk++) {
        std::vector<TrackletProjectionsMemory*> tmp(settings_.nallstubs(layerdisk), nullptr);
        outputproj_.push_back(tmp);
      }

}

  // Project to layer (taken from TrackletCalculatorBase.cc)
void ProjectionCalculator::projLayer(int ir, int irinv, int iphi0, int it, int iz0, int &iz, int &iphi) {
  int irtilde = ir*phiHG_/sqrt(6.0)+0.5;
  int is =  (irtilde*irinv) >> n_s_;
  int is6 =  (1 << n_s6_) + ((is*is) >> (2 + 2*n_r_ + 2*n_rinv_ - 2*n_s_ - n_s6_));
  int iu = (ir*irinv) >> (n_r_ + n_rinv_ + 1 - n_phi_);
  iphi = (iphi0 << (n_phi_ - n_phi0_)) - ((iu*is6) >> n_s6_);
  int iv = (it*ir) >> (n_r_ + n_t_ - n_z_);
  iz = iz0 + ((iv*is6) >> n_s6_);
}
  // Project to disk (taken from TrackletCalculatorBase.cc)
void ProjectionCalculator::projDisk(int iz, int irinv, int iphi0, int it,int iz0, int &ir, int &iphi, int &iderphi, int &iderr) {

    int iz0_sign = (it>0)?iz0:-iz0;

    assert(abs(it)<LUT_itinv_.size());
    int itinv = LUT_itinv_[abs(it)];
    
    iderphi = (-irinv*itinv) >> 17;
    iderr = itinv >> 5;
    
    if (it<0) {
      iderphi = -iderphi;
      iderr = -iderr;
    }	   
    
    int iw = (((iz << (n_r_ - n_z_)) - (iz0_sign << (n_r_ - n_z_)))*itinv) >> n_tinv_;

    iphi = (iphi0 >> (n_phi0_ - n_phidisk_)) - ((iw*irinv) >> (1 + n_r_ + n_rinv_ - n_phidisk_));

    int ifact = (1 << n_y_)*phiHG_/sqrt(6.0);

    int iy = (ifact * irinv) >> n_y_;

    int ix = (iw*iy) >> n_x_;

    int ix6 = (1 << n_xx6_) - ( (ix*ix) >> (2 + 2*n_r_ + 2*n_rinv_ - 2*n_x_ - n_xx6_));

    ir = (iw*ix6) >> (n_r_ - n_rdisk_ + n_xx6_);
}

void ProjectionCalculator::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  if (output == "projout") {
    auto* tmp = dynamic_cast<TrackletProjectionsMemory*>(memory);
    unsigned int layerdisk = memory->getName()[memory->getName().size() - 5] - '1';   //layer or disk counting from 0
    unsigned int phiregion = memory->getName()[memory->getName().size() - 1] - 'A';  //phiregion counting from 0
    if (memory->getName()[memory->getName().size() - 6] == 'D') layerdisk += N_LAYER;
    projnames_.push_back(memory->getName());
    assert(layerdisk < N_LAYER + N_DISK);
    assert(phiregion < outputproj_[layerdisk].size());
    //check that phiregion not already initialized
    assert(outputproj_[layerdisk][phiregion] == nullptr);
    assert(tmp != nullptr);
    outputproj_[layerdisk][phiregion] = tmp;
    return;
  }

  if (output == "tparout") {
    auto* tmp = dynamic_cast<TrackletParametersMemory*>(memory);
    assert(tmp != nullptr);
    outputpars_.push_back(tmp);
    return;
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find output: " << output;
}

void ProjectionCalculator::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }

  if (input == "tparin" ) {
    auto* tmp = dynamic_cast<TrackletParametersMemory*>(memory);
    assert(tmp != nullptr);
    inputpars_.push_back(tmp);
    return;
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find input: " << input;
}

void ProjectionCalculator::execute() {

  for(unsigned int i = 0 ; i < inputpars_.size() ; i++) { // send copy of tpars to TB 
    int projPage = 0;
    std::string iname = inputpars_[i]->getName();

    for(unsigned int k = 0 ; k < outputpars_.size() ; k++) { // add copy of par to merged par output memory
      std::string oname = outputpars_[k]->getName();
      int parPage = iname[9]-oname[9];
      for(unsigned int j = 0; j < inputpars_[i]->nTracklets(); j++) {
	      outputpars_[k]->addTracklet(inputpars_[i]->getTracklet(j), parPage);
      }
    }
    
    std::string seed = iname.substr(5, 4); // extract seed from name 
    for(unsigned int k = 0; k < inputpars_[i]->nTracklets(); k++){
    auto tracklet = inputpars_[i]->getTracklet(k);
    //double phi0 = tracklet->phi0(); // non-digi track params, currently unneeded / unused 
    //double z0 = tracklet->z0();
    double t = tracklet->t(); // used later for t cut 
    //double rinv = tracklet->rinv();

    int irinv = tracklet->fpgarinv().value(); // digi track params
    int iphi0 = tracklet->fpgaphi0().value();
    int iz0 = tracklet->fpgaz0().value();
    int it = tracklet->fpgat().value();

    std::vector<int> izr_LD(N_LAYER + N_DISK, 0);
    std::vector<int> iphi_LD(N_LAYER + N_DISK, 0);
    std::vector<bool> valid_LD(N_LAYER + N_DISK, 0);
    std::vector<bool> addedLayer(N_LAYER, 0);
    std::vector<int> der_phi_LD(2, 0);
    std::vector<int> der_zr_LD(2, 0);
    
    Projection projs[N_LAYER + N_DISK];

    /////////////////////////////////
    // calculate layer projections //
    /////////////////////////////////
    int ir;
    bool valid_zmin, valid_zmax, valid_phimin, valid_phimax; // CODE ADAPTED FROM FIRMWARE-HLS PROJECTIONCALCULATOR - could be simplified to only consider layers/disks in wiring config. 
    const int zmin = -(1 << (settings_.nzbitsstub(0) - 1));
    const int zmax = (1 << (settings_.nzbitsstub(0) - 1));
    const int phimax = (1 << (settings_.nphibitsstub(3))) - 1;
    const int phimin = 0;
    // L1 
    ir = settings_.irmean(0);
    projLayer(ir, irinv, iphi0, it, iz0, izr_LD[0], iphi_LD[0]);
    valid_zmin = izr_LD[0] >= zmin;
    valid_zmax = izr_LD[0] < zmax;
    valid_phimax = iphi_LD[0] < phimax;
    valid_phimin = iphi_LD[0] > phimin;
    valid_LD[0] = valid_zmin & valid_zmax & valid_phimax & valid_phimin;
    iphi_LD[0] = iphi_LD[0] >> 3; //now done later 

    // L2 
    ir = settings_.irmean(1);
    projLayer(ir, irinv, iphi0, it, iz0, izr_LD[1], iphi_LD[1]);
    valid_zmin = izr_LD[1] >= zmin;
    valid_zmax = izr_LD[1] < zmax;
    valid_phimax = iphi_LD[1] < phimax;
    valid_phimin = iphi_LD[1] > phimin;
    valid_LD[1] = valid_zmin & valid_zmax & valid_phimax & valid_phimin;
    iphi_LD[1] = iphi_LD[1] >> 3;

    // L3
    ir = settings_.irmean(2);
    projLayer(ir, irinv, iphi0, it, iz0, izr_LD[2], iphi_LD[2]);
    valid_zmin = izr_LD[2] >= zmin;
    valid_zmax = izr_LD[2] < zmax;
    valid_phimax = iphi_LD[2] < phimax;
    valid_phimin = iphi_LD[2] > phimin;
    valid_LD[2] = valid_zmin & valid_zmax & valid_phimax & valid_phimin;
    iphi_LD[2] = iphi_LD[2] >> 3;

    // L4 
    ir = settings_.irmean(3);
    projLayer(ir, irinv, iphi0, it, iz0, izr_LD[3], iphi_LD[3]);
    valid_zmin = izr_LD[3] >= zmin;
    valid_zmax = izr_LD[3] < zmax;
    valid_phimax = iphi_LD[3] < phimax;
    valid_phimin = iphi_LD[3] > phimin;
    valid_LD[3] = valid_zmin & valid_zmax & valid_phimax & valid_phimin;
    izr_LD[3] = izr_LD[3] >> 4;

    // L5
    ir = settings_.irmean(4);
    projLayer(ir, irinv, iphi0, it, iz0, izr_LD[4], iphi_LD[4]);
    valid_zmin = izr_LD[4] >= zmin;
    valid_zmax = izr_LD[4] < zmax;
    valid_phimax = iphi_LD[4] < phimax;
    valid_phimin = iphi_LD[4] > phimin;
    valid_LD[4] = valid_zmin & valid_zmax & valid_phimax & valid_phimin;
    izr_LD[4] = izr_LD[4] >> 4;
    

    // L6
    ir = settings_.irmean(5);
    projLayer(ir, irinv, iphi0, it, iz0, izr_LD[5], iphi_LD[5]);
    valid_zmin = izr_LD[5] >= zmin;
    valid_zmax = izr_LD[5] < zmax;
    valid_phimax = iphi_LD[5] < phimax;
    valid_phimin = iphi_LD[5] > phimin;
    valid_LD[5] = valid_zmin & valid_zmax & valid_phimax & valid_phimin;
    izr_LD[5] = izr_LD[5] >> 4;

    // Derivatives
    der_phi_LD[0] = -(irinv >> (1+3)); 
    der_zr_LD[0] = it >> 3;

    ////////////////////////////////
    // calculate disk projections //
    ////////////////////////////////
    double irmindisk = settings_.rmindisk() / settings_.krprojshiftdisk(); 
    double irmaxdisk = settings_.rmaxdisk() / settings_.krprojshiftdisk();

    int tcut = 1.0/(settings_.ktpars());

    // D1 
    int izproj0 = settings_.izmean(0);
    projDisk(izproj0, irinv, iphi0, it, iz0, izr_LD[6], iphi_LD[6], der_phi_LD[1], der_zr_LD[1]);
    valid_LD[6] = izr_LD[6] >= irmindisk && izr_LD[6] < irmaxdisk && ((it > tcut) || (it < -tcut));

    // D2 
    int izproj1 = settings_.izmean(1);
    projDisk(izproj1, irinv, iphi0, it, iz0, izr_LD[7], iphi_LD[7], der_phi_LD[1], der_zr_LD[1]);
    valid_LD[7] = izr_LD[7] >= irmindisk && izr_LD[7] < irmaxdisk && ((it > tcut) || (it < -tcut));

    // D3
    int izproj2 = settings_.izmean(2);
    projDisk(izproj2, irinv, iphi0, it, iz0, izr_LD[8], iphi_LD[8], der_phi_LD[1], der_zr_LD[1]);
    valid_LD[8] = izr_LD[8] >= irmindisk && izr_LD[8] < irmaxdisk && ((it > tcut) || (it < -tcut));
    
    // D4
    int izproj3 = settings_.izmean(3);
    projDisk(izproj3, irinv, iphi0, it, iz0, izr_LD[9], iphi_LD[9], der_phi_LD[1], der_zr_LD[1]);
    valid_LD[9] = izr_LD[9] >= irmindisk && izr_LD[9] < irmaxdisk && ((it > tcut) || (it < -tcut));

    // D5
    int izproj4 = settings_.izmean(4);
    projDisk(izproj4, irinv, iphi0, it, iz0, izr_LD[10], iphi_LD[10], der_phi_LD[1], der_zr_LD[1]);
    valid_LD[10] = izr_LD[10] >= irmindisk && izr_LD[10] < irmaxdisk && ((t > tcut) || (t < -tcut));

    ///////////////////////////////////
    // Write projections to memories //
    ///////////////////////////////////
    std::vector<std::string> seedNames = {"L1L2", "L2L3", "L3L4", "L5L6", "D1D2", "D3D4", "L1D1", "L2D1"};

    for (int iSeed = 0; iSeed < 8; ++iSeed){
      bool overlap = (iSeed >= 6);
      if (seed == seedNames[iSeed]){ // FIXME find easier way to get iSeed (probably from seed name)
      unsigned int numTCs = nMergedTC[iSeed];
      for(unsigned int iTC = 0; iTC < numTCs; ++iTC){
        std::string tcStr = TrackletConfigBuilder::iMergedTCStr(iSeed, iTC); 
        size_t index = tcStr.find(iname[9]); // find index in merged TC string to find page
        if (index != std::string::npos) {
            projPage = static_cast<int>(index); // calculate projPage FIXME find better way of doing this
        } else {
            continue;
        }
      }
      

      for (unsigned int j = 0; j < settings_.projlayers()[iSeed].size(); ++j){
        unsigned int layer = settings_.projlayers()[iSeed][j]; // Loop through layers/disks projected to
        if (layer == 0) continue; // for seeds not projecting to any layers 
        if(valid_LD[layer - 1]){ // If projection to layer/disk is valid
            if ((izr_LD[layer - 1] == - (1 << (settings_.nzbitsstub(layer - 1) - 1))) || (izr_LD[layer - 1] == ((1 << (settings_.nzbitsstub(layer - 1) - 1)) - 1))) { // reject extreme z values 
              continue;
            }
            if (std::abs(izr_LD[layer - 1]) > 2048) {
              continue;
            }

          double phiprojlayer = iphi_LD[layer - 1] * settings_.kphi(layer - 1); // get un-digi projections 
          double zprojlayer = izr_LD[layer - 1] * settings_.kz(); // FIXME find better way to calculate these - but don't have stub coordinates to use exacttracklet function
          double phiderlayer = der_phi_LD[0] * settings_.kphider(); 
          double zderlayer = der_zr_LD[0] * settings_.kzder();
          
          projs[layer - 1].init(settings_,
                            layer - 1,
                            iphi_LD[layer - 1],
                            izr_LD[layer - 1],
                            der_phi_LD[0],
                            der_zr_LD[0],
                            phiprojlayer,
                            zprojlayer,
                            phiderlayer,
                            zderlayer,
                            phiprojlayer,
                            zprojlayer,
                            phiderlayer,
                            zderlayer,
                            overlap);
          addedLayer[layer - 1] = true;
        }
      }
      for (unsigned int j = 0; j < settings_.projdisks()[iSeed].size(); ++j){
        unsigned int disk = settings_.projdisks()[iSeed][j];
        if (disk == 0) continue; // for seeds not projecting to any disks  
        if(valid_LD[N_LAYER + disk - 1]){ // If projection to layer/disk is valid

        if (iphi_LD[N_LAYER + disk - 1] <= 0) // reject extreme phi values
          continue;
        if (iphi_LD[N_LAYER + disk - 1] >= (1 << settings_.nphibitsstub(0)) - 1)
          continue;

        if (iSeed <= 4){ // if barrel seed, need to check if haven't already projected to layer
          if (disk == 1 && addedLayer[5]) continue; 
          if (disk == 2 && addedLayer[4]) continue; 
          if (disk == 3 && addedLayer[3]) continue; 
          if (disk == 4 && addedLayer[2]) continue; 
        }
          
        double phiprojdisk = iphi_LD[N_LAYER + disk - 1]*settings_.kphi(N_LAYER); // un-digitize values using granularities for initializing projections
        double rprojdisk = izr_LD[N_LAYER + disk - 1] * settings_.kr();
        double phiderdisk = der_phi_LD[1]* settings_.kphiderdisk();
        double rderdisk = der_zr_LD[1] * settings_.krder();

        projs[N_LAYER + disk - 1].init(settings_,
                            N_LAYER + disk - 1,
                            iphi_LD[N_LAYER + disk - 1],
                            izr_LD[N_LAYER + disk - 1],
                            der_phi_LD[1],
                            der_zr_LD[1],
                            phiprojdisk,
                            rprojdisk,
                            phiderdisk,
                            rderdisk,
                            phiprojdisk,
                            rprojdisk,
                            phiderdisk,
                            rderdisk,
                            overlap);
        }
      }
      tracklet->addProjs(projs);
      for(unsigned int layerdisk = 0; layerdisk < N_LAYER + N_DISK; ++layerdisk){
      if (tracklet->validProj(layerdisk)) {
        if (layerdisk < N_LAYER){
          
          FPGAWord fpgaz = tracklet->proj(layerdisk).fpgarzproj();
          FPGAWord fpgaphi = tracklet->proj(layerdisk).fpgaphiproj();

          if (fpgaphi.atExtreme())
            edm::LogProblem("Tracklet") << "at extreme! " << fpgaphi.value();

          assert(!fpgaphi.atExtreme());

          if (fpgaz.atExtreme())
            continue;

          if (std::abs(fpgaz.value() * settings_.kz()) > settings_.zlength())
            continue;

          int iphivmRaw = fpgaphi.value() >> (fpgaphi.nbits() - 5);
          int iphi = iphivmRaw / (32 / settings_.nallstubs(layerdisk));

          if (outputproj_[layerdisk][iphi] == nullptr) continue; 
          outputproj_[layerdisk][iphi]->addProj(tracklet, projPage); // FIXME write to correct page - though doesn't affect emulation
        } else{

          FPGAWord fpgar = tracklet->proj(layerdisk).fpgarzproj();

          if (fpgar.value() * settings_.krprojshiftdisk() < settings_.rmindiskvm())
            continue;
          if (fpgar.value() * settings_.krprojshiftdisk() > settings_.rmaxdisk())
            continue;

          FPGAWord fpgaphi = tracklet->proj(layerdisk).fpgaphiproj();
          int iphivmRaw = fpgaphi.value() >> (fpgaphi.nbits() - 5);
          int iphi = iphivmRaw / (32 / settings_.nallstubs(layerdisk));//>> settings_.nbitsallstubs(layerdisk);
          
          if (outputproj_[layerdisk][iphi] == nullptr) continue; 

          outputproj_[layerdisk][iphi]->addProj(tracklet, projPage); // FIXME write to correct page 
        }
      }
    }
    }
  }
  }
  }

  return;

}
