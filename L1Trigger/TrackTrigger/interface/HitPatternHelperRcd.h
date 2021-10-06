//
//  Created by J.Li on 1/23/21.
//

#ifndef L1Trigger_TrackTrigger_interface_HitPatternHelperRcd_h
#define L1Trigger_TrackTrigger_interface_HitPatternHelperRcd_h

#include "FWCore/Framework/interface/DependentRecordImplementation.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "FWCore/Utilities/interface/mplVector.h"

namespace HPH {

  typedef edm::mpl::Vector<TrackerDigiGeometryRecord, TrackerTopologyRcd> Rcds;

  // record of HPH::SetupRcd
  class SetupRcd : public edm::eventsetup::DependentRecordImplementation<SetupRcd, Rcds> {};

}  // namespace HPH

#endif
