#include "L1Trigger/TrackFindingTracklet/interface/ProjectionCalculator.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace std;
using namespace trklet;

ProjectionCalculator::ProjectionCalculator(string name, Settings const& settings, Globals* global)
    : ProcessBase(name, settings, global) {

}

void ProjectionCalculator::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  if (output == "projout") {
    auto* tmp = dynamic_cast<TrackletProjectionsMemory*>(memory);
    assert(tmp != nullptr);
    outputproj_.push_back(tmp);
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

  if (input == "projin" ) {
    auto* tmp = dynamic_cast<TrackletProjectionsMemory*>(memory);
    assert(tmp != nullptr);
    inputproj_.push_back(tmp);
    return;
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

  for(unsigned int i = 0 ; i < inputproj_.size() ; i++) {
    std::string iname = inputproj_[i]->getName();
    std::string ireg = iname.substr(iname.size()-6,iname.size());
    for(unsigned int k = 0 ; k < outputproj_.size() ; k++) {
      std::string oname = outputproj_[k]->getName();
      std::string oreg = oname.substr(oname.size()-6,oname.size());
      if (oreg == ireg) {
	int page = iname[10]-oname[10];
	for(unsigned int j = 0; j < inputproj_[i]->nTracklets(); j++) {
	  outputproj_[k]->addProj(inputproj_[i]->getTracklet(j), page);
	}
      }
    }
  }

  for(unsigned int i = 0 ; i < inputpars_.size() ; i++) {
    std::string iname = inputpars_[i]->getName();
    for(unsigned int k = 0 ; k < outputpars_.size() ; k++) {
      std::string oname = outputpars_[k]->getName();
      int page = iname[9]-oname[9];
      for(unsigned int j = 0; j < inputpars_[i]->nTracklets(); j++) {
	outputpars_[k]->addTracklet(inputpars_[i]->getTracklet(j), page);
      }
	//}
    }
  }

  return;

}
