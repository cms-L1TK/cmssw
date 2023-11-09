#ifndef L1Trigger_TrackerDTC_SetupRcd_h
#define L1Trigger_TrackerDTC_SetupRcd_h

#include "FWCore/Framework/interface/DependentRecordImplementation.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "CondFormats/DataRecord/interface/TrackerDetToDTCELinkCablingMapRcd.h"
#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithmRecord.h"

#include "FWCore/Utilities/interface/mplVector.h"

namespace tt {

  typedef edm::mpl::
      Vector<TrackerDigiGeometryRecord, TrackerTopologyRcd, TrackerDetToDTCELinkCablingMapRcd, TTStubAlgorithmRecord>
          Rcds;

  // record of tt::Setup
  class SetupRcd : public edm::eventsetup::DependentRecordImplementation<SetupRcd, Rcds> {};

}  // namespace tt

#endif