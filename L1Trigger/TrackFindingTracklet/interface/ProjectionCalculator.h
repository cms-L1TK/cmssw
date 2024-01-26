#ifndef L1Trigger_TrackFindingTracklet_interface_ProjectionCalculator_h
#define L1Trigger_TrackFindingTracklet_interface_ProjectionCalculator_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletLUT.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"

namespace trklet {

  class Settings;
  class Globals;
  class MemoryBase;

  class ProjectionCalculator : public ProcessBase {
  public:
    ProjectionCalculator(std::string name, Settings const& settings, Globals* global);

    ~ProjectionCalculator() override = default;

    void addOutput(MemoryBase* memory, std::string output) override;
    void addInput(MemoryBase* memory, std::string input) override;

    void execute();

  private:

    std::vector<TrackletProjectionsMemory*> inputproj_;

    std::vector<TrackletProjectionsMemory*> outputproj_;

  };

};  // namespace trklet
#endif
