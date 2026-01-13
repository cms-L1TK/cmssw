#ifndef L1Trigger_TrackerTFP_Demonstrator_h
#define L1Trigger_TrackerTFP_Demonstrator_h

#include "FWCore/Framework/interface/data_default_record_trait.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"

#include <vector>
#include <string>

namespace trackerTFP {

  /*! \class  trackerTFP::Demonstrator
   *  \brief  ESProduct providing the algorithm to run input data through modelsim
   *          and to compares results with expected output data
   *  \date   2021, April
   */
  class Demonstrator {
  public:
    // configuration
    struct Config {
      std::string dirIPBB_;
      double runTime_;
      std::vector<int> linkMappingIn_;
      std::vector<int> linkMappingOut_;
    };
    Demonstrator() {}
    Demonstrator(const Config& iConfig, const tt::Setup* setup);
    ~Demonstrator() = default;
    // plays input through modelsim and compares result with output
    bool analyze(const std::vector<std::vector<tt::Frame>>& input,
                 const std::vector<std::vector<tt::Frame>>& output) const;
    // create testvector compatible with L1T from the TFP output - each tracker region gets its own set of links
    void l1tInput(const std::vector<std::vector<tt::Frame>>& output, const int nL1TEvent) const;
    // test on hardware
    int hw(std::string dirTGZ, std::string boardAdress, std::string boardType, std::string fpgaName, std::string dockerImage, int bufferOffset, int nEvents, bool saveAllFiles = false) const;
  private:
    // converts streams of bv into stringstream
    void convert(const std::vector<std::vector<tt::Frame>>& bits,
                 std::stringstream& ss,
                 const std::vector<int>& mapping,
                 const int nL1TEvent = -1) const;
    // plays stringstream through modelsim
    void sim(const std::stringstream& ss) const;
    // compares stringstream with modelsim output
    bool compare(std::stringstream& ss) const;
    // creates emp file header
    std::string header(const std::vector<int>& links) const;
    // creates 6 frame gap between packets
    std::string infraGap(int& nFrame, int numLinks) const;
    // creates frame number
    std::string frame(int& nFrame) const;
    // converts bv into hex
    std::string hex(const tt::Frame& bv, bool first = false) const;

    // path to ipbb proj area
    std::string dirIPBB_;
    // runtime in ms
    double runTime_;
    // which input links/channels to use
    std::vector<int> linkMappingIn_;
    // which ouput links/channels to use
    std::vector<int> linkMappingOut_;
    // path to input text file
    std::string dirIn_;
    // path to output text file
    std::string dirOut_;
    // path to expected output text file
    std::string dirPre_;
    // path to diff text file
    std::string dirDiff_;
    // number of frames per event (161-6)
    int numFrames_;
    // number of emp reset frames per event (6)
    int numFramesInfra_;
    // number of TFPs per time node (9)
    int numRegions_;
    // number of L1T testvectors
    mutable int nL1TFiles_ = 0;
  };

}  // namespace trackerTFP

EVENTSETUP_DATA_DEFAULT_RECORD(trackerTFP::Demonstrator, tt::SetupRcd);

#endif
