#include "L1Trigger/TrackerTFP/interface/Demonstrator.h"

#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <numeric>

namespace trackerTFP {

  Demonstrator::Demonstrator(const Config& iConfig, const tt::Setup* setup)
      : dirIPBB_(iConfig.dirIPBB_),
        runTime_(iConfig.runTime_),
        linkMappingIn_(iConfig.linkMappingIn_),
        linkMappingOut_(iConfig.linkMappingOut_),
        dirIn_(dirIPBB_ + "in.txt"),
        dirOut_(dirIPBB_ + "out.txt"),
        dirPre_(dirIPBB_ + "pre.txt"),
        dirDiff_(dirIPBB_ + "diff.txt"),
        numFrames_(setup->numFramesIOHigh()),
        numFramesInfra_(setup->numFramesInfra()),
        numRegions_(setup->numRegions()) {}

  // plays input through modelsim and compares result with output
  bool Demonstrator::analyze(const std::vector<std::vector<tt::Frame>>& input,
                             const std::vector<std::vector<tt::Frame>>& output) const {
    // default link mapping
    auto linkMapping =
        [this](const std::vector<int>& mapC, std::vector<int>& map, const std::vector<std::vector<tt::Frame>>& data) {
          if (mapC.empty()) {
            map.resize(data.size() / numRegions_);
            std::iota(map.begin(), map.end(), 0);
          } else
            map = mapC;
        };
    // converts input into stringstream
    std::stringstream ss;
    std::vector<int> map;
    linkMapping(linkMappingIn_, map, input);
    convert(input, ss, map);
    // play input through modelsim
    sim(ss);
    // converts output into stringstream
    map.clear();
    linkMapping(linkMappingOut_, map, output);
    convert(output, ss, map);
    // compares output with modelsim output
    auto result = compare(ss);
    return result;
  }

  // converts streams of bv into stringstream
  void Demonstrator::convert(const std::vector<std::vector<tt::Frame>>& bits,
                             std::stringstream& ss,
                             const std::vector<int>& mapping) const {
    // reset ss
    ss.str("");
    ss.clear();
    // number of tranceiver in a quad
    static constexpr int quad = 4;
    std::set<int> quads;
    for (int channel : mapping)
      quads.insert(channel / quad);
    std::vector<int> links;
    links.reserve(quads.size() * quad);
    for (int q : quads) {
      const int offset = q * quad;
      for (int c = 0; c < quad; c++)
        links.push_back(offset + c);
    }
    // start with header
    ss << header(links);
    int nFrame(0);
    // create one packet per region
    bool first = true;
    for (int region = 0; region < numRegions_; region++) {
      const int offset = region * mapping.size();
      // start with emp 6 frame gap
      ss << infraGap(nFrame, links.size());
      for (int frame = 0; frame < numFrames_; frame++) {
        // write one frame for all channel
        ss << this->frame(nFrame);
        for (int link : links) {
          const auto channel = find(mapping.begin(), mapping.end(), link);
          if (channel == mapping.end())
            ss << "  0000 " << std::string(TTBV::S_ / 4, '0');
          else {
            const std::vector<tt::Frame>& bvs = bits[offset + std::distance(mapping.begin(), channel)];
            ss << (frame < static_cast<int>(bvs.size()) ? hex(bvs[frame], first) : hex(tt::Frame(), first));
          }
        }
        ss << std::endl;
        first = false;
      }
    }
  }

  // plays stringstream through modelsim
  void Demonstrator::sim(const std::stringstream& ss) const {
    // write ss to disk
    std::fstream fs;
    fs.open(dirIn_.c_str(), std::fstream::out);
    fs << ss.rdbuf();
    fs.close();
    // run modelsim
    std::stringstream cmd;
    cmd << "cd " << dirIPBB_ << " && ./run_sim -quiet -c work.top -do 'run " << runTime_
        << "us' -do 'quit' &> /dev/null";
    std::system(cmd.str().c_str());
  }

  // compares stringstream with modelsim output
  bool Demonstrator::compare(std::stringstream& ss) const {
    // write ss to disk
    std::fstream fs;
    fs.open(dirPre_.c_str(), std::fstream::out);
    fs << ss.rdbuf();
    fs.close();
    // use linux diff on disk
    const std::string c = "diff " + dirPre_ + " " + dirOut_ + " &> " + dirDiff_;
    // read diff output
    std::system(c.c_str());
    ss.str("");
    ss.clear();
    fs.open(dirDiff_.c_str(), std::fstream::in);
    ss << fs.rdbuf();
    fs.close();
    // count lines, 4 are expected
    int n(0);
    std::string token;
    while (getline(ss, token))
      n++;
    return n == 4;
  }

  // creates emp file header
  std::string Demonstrator::header(const std::vector<int>& links) const {
    std::stringstream ss;
    // file header
    ss << "ID: CMSSW" << std::endl
       << "Metadata: (strobe,) start of orbit, start of packet, end of packet, valid" << std::endl
       << std::endl;
    // link header
    ss << "      Link  ";
    for (int link : links)
      ss << "            " << std::setfill('0') << std::setw(3) << link << "        ";
    ss << std::endl;
    return ss.str();
  }

  // creates 6 frame gap between packets
  std::string Demonstrator::infraGap(int& nFrame, int numLinks) const {
    std::stringstream ss;
    for (int gap = 0; gap < numFramesInfra_; gap++) {
      ss << frame(nFrame);
      for (int link = 0; link < numLinks; link++)
        ss << "  0000 " << std::string(TTBV::S_ / 4, '0');
      ss << std::endl;
    }
    return ss.str();
  }

  // creates frame number
  std::string Demonstrator::frame(int& nFrame) const {
    std::stringstream ss;
    ss << "Frame " << std::setfill('0') << std::setw(4) << nFrame++ << "  ";
    return ss.str();
  }

  // converts bv into hex
  std::string Demonstrator::hex(const tt::Frame& bv, bool first) const {
    std::stringstream ss;
    ss << (first ? "  1001 " : "  0001 ") << std::setfill('0') << std::setw(TTBV::S_ / 4) << std::hex << bv.to_ullong();
    return ss.str();
  }

  // converts streams of bv into stringstream
  void Demonstrator::convertL1T(const std::vector<std::vector<tt::Frame>>& bits,
                             std::stringstream& ss,
                             const std::vector<int>& mapping,
                             const int nL1TEvent,
                             const int region) const {
    // reset ss
    ss.str("");
    ss.clear();
    // number of tranceiver in a quad
    static constexpr int quad = 4;
    std::set<int> quads;
    for (int channel : mapping)
      quads.insert(channel / quad);
    std::vector<int> links;
    links.reserve(quads.size() * quad);
    for (int q : quads) {
      const int offset = q * quad;
      for (int c = 0; c < quad; c++)
        links.push_back(offset + c);
    }
    // start with header
    if (nL1TEvent == 0)
      ss << header(links);
    int nFrame = (numFrames_ + numFramesInfra_) * nL1TEvent;
    // create one packet per region
    bool first = nL1TEvent > 0 ? false : true;
    const int offset = region * mapping.size();
    // start with emp 6 frame gap
    ss << infraGap(nFrame, links.size());
    for (int frame = 0; frame < numFrames_; frame++) {
      // write one frame for all channel
      ss << this->frame(nFrame);
      for (int link : links) {
        const auto channel = find(mapping.begin(), mapping.end(), link);
        if (channel == mapping.end())
          ss << "  0000 " << std::string(TTBV::S_ / 4, '0');
        else {
          const std::vector<tt::Frame>& bvs = bits[offset + std::distance(mapping.begin(), channel)];
          ss << (frame < static_cast<int>(bvs.size()) ? hex(bvs[frame], first) : hex(tt::Frame(), first));
        }
      }
      ss << std::endl;
      first = false;
    }

  }

  // create testvector compatible with integration testing - tracker regions are not stacked in the same file
  void Demonstrator::makeInputFiles(const std::vector<std::vector<tt::Frame>>& input, const int nL1TEvent) const {
    // default link mapping
    auto linkMapping =
        [this](const std::vector<int>& mapC, std::vector<int>& map, const std::vector<std::vector<tt::Frame>>& data) {
          if (mapC.empty()) {
            map.resize(data.size()/numRegions_);
            std::iota(map.begin(), map.end(), 0);
          } else
            map = mapC;
        };
    // converts input into stringstream
    for (int region = 0; region < numRegions_; region++) {
      std::stringstream ss;
      std::vector<int> map;
      linkMapping(linkMappingIn_, map, input);
      convertL1T(input, ss, map, nL1TEvent, region);
      // write ss to disk
      std::fstream fs;
      if (nL1TEvent == 0) { // new file
        const std::string dirL1T = dirIPBB_ + "in_nonant" + std::to_string(region) + "_" + std::to_string(nInputFiles_) + ".txt";
        fs.open(dirL1T.c_str(), std::fstream::out);
        if (region == numRegions_ - 1)
          nInputFiles_++;
      } else { // append to existing file
        const std::string dirL1T = dirIPBB_ + "in_nonant" + std::to_string(region) + "_" + std::to_string(nInputFiles_ - 1) + ".txt";
        fs.open(dirL1T.c_str(), std::fstream::app);
      }
      fs << ss.rdbuf();
      fs.close();
    }
  }

  // create testvector compatible with integration testing - tracker regions are not stacked in the same file
  void Demonstrator::makeOutputFiles(const std::vector<std::vector<tt::Frame>>& output, const int nL1TEvent) const {
    // default link mapping
    auto linkMapping =
        [this](const std::vector<int>& mapC, std::vector<int>& map, const std::vector<std::vector<tt::Frame>>& data) {
          if (mapC.empty()) {
            map.resize(data.size()/numRegions_);
            std::iota(map.begin(), map.end(), 0);
          } else
            map = mapC;
        };
    // converts output into stringstream
    for (int region = 0; region < numRegions_; region++) {
      std::stringstream ss;
      std::vector<int> map;
      linkMapping(linkMappingOut_, map, output);
      convertL1T(output, ss, map, nL1TEvent, region);
      // write ss to disk
      std::fstream fs;
      if (nL1TEvent == 0) { // new file
        const std::string dirL1T = dirIPBB_ + "pre_nonant" + std::to_string(region) + "_" + std::to_string(nOutputFiles_) + ".txt";
        fs.open(dirL1T.c_str(), std::fstream::out);
        if (region == numRegions_ - 1)
          nOutputFiles_++;
      } else { // append to existing file
        const std::string dirL1T = dirIPBB_ + "pre_nonant" + std::to_string(region) + "_" + std::to_string(nOutputFiles_ - 1) + ".txt";
        fs.open(dirL1T.c_str(), std::fstream::app);
      }
      fs << ss.rdbuf();
      fs.close();
    }
  }

  // create testvector compatible with L1T from the TFP output - each tracker region gets its own set of links
  void Demonstrator::makeL1TFiles(const std::vector<std::vector<tt::Frame>>& output, const int nL1TEvent) const {
    // default link mapping
    auto linkMapping =
        [this](const std::vector<int>& mapC, std::vector<int>& map, const std::vector<std::vector<tt::Frame>>& data) {
          map.resize(data.size());
          std::iota(map.begin(), map.end(), 0);
        };
    // converts output into stringstream
    std::stringstream ss;
    std::vector<int> map;
    linkMapping(linkMappingOut_, map, output);
    convertL1T(output, ss, map, nL1TEvent);
    // write ss to disk
    std::fstream fs;
    if (nL1TEvent == 0) { // new file
      const std::string dirL1T = dirIPBB_ + "l1t_" + std::to_string(nL1TFiles_++) + ".txt";
      fs.open(dirL1T.c_str(), std::fstream::out);
    } else { // append to existing file
      const std::string dirL1T = dirIPBB_ + "l1t_" + std::to_string(nL1TFiles_ - 1) + ".txt";
      fs.open(dirL1T.c_str(), std::fstream::app);
    }
    fs << ss.rdbuf();
    fs.close();
  }

  // plays input files through hardware and compares with simulation output
  int Demonstrator::hw(std::string dirTGZ, std::string boardAdress, std::string boardType, std::string nameFPGA, std::string dockerImage, int bufferOffset, int nEvents, bool saveAllFiles) const {
    const std::string dirHWOut = dirIPBB_ + "/hw_out.txt";
    const std::string dirYAML = dirIPBB_ + "/default_" + boardType + ".yml";
    // make yaml file
    if (!nEvents) {
      std::stringstream ss;
      ss << R"(contexts:
  )" << boardType << R"(.processors:
    parameters:
      program:package:
        replace_with: firmware_package
      reset:source: internal
      powerUp:wait (s): "5"
      configureRxBuffers:startOffset: [0]
      configureTxBuffers:startOffset: [611]
      §configureRxBuffers:
        mode: [PlayOnce]
        ids: [replace_with: input_channels]
        dataURI: ['file']
        dataFile: [replace_with: input_data_file]
      §configureTxBuffers:
        mode: [Capture]
        ids: [replace_with: output_channels]
        dataURI: ['']
        dataFile: [type: file]
      saveRxBuffers:ids:
        replace_with: input_channels
      saveTxBuffers:ids:
        replace_with: output_channels)";
      std::ofstream ofs;
      ofs.open(dirYAML.c_str(), std::fstream::out);
      ofs << ss.rdbuf();
      ofs.close();
    }
    // extract link mapping from file
    auto getLinkMap = [](std::string& links, std::string fileName){
      std::ifstream ifs;
      ifs.open(fileName.c_str());
      std::string line;
      while (std::getline(ifs, line)) {
        if (line.find("Link") != std::string::npos) {
          std::stringstream ss(line);
          std::string token;
          while (ss >> token) {
            if (std::isdigit(token[0])) {
              if (!links.empty())
                  links += ",";
              links += std::to_string(std::stoi(token));
            }
          }
        }
      }
      ifs.close();
    };
    std::string inputLinks;
    std::string outputLinks;
    getLinkMap(inputLinks, dirIn_);
    getLinkMap(outputLinks, dirOut_);
    // run singularity
    std::stringstream cmd;
    cmd << "singularity run " << dockerImage << " " << boardAdress << " " << nameFPGA << " " << dirYAML << " " << dirTGZ << " " << dirIn_ << " " << dirHWOut << " --input-channels=" << inputLinks << " --output-channels=" << outputLinks << " --expected-output " << dirPre_ << " --max-frames 971";
    if (nEvents)
      cmd << " --continue";
    int result = std::system(cmd.str().c_str());
    // save input/output files
    if (saveAllFiles) {
      const std::string suffix = "_" + std::to_string(nEvents);
      // copy file function
      auto copyFile = [suffix](const std::string dir){
        std::string dirCopy = dir;
        dirCopy.insert(dirCopy.find_last_of('.'), suffix);
        // dirCopy.insert(dirCopy.find_last_of('/'), "/hwtest"); // TODO: CREATE THIS DIRECTORY
        std::ifstream ifs(dir.c_str());
        std::ofstream ofs(dirCopy.c_str());
        ofs << ifs.rdbuf();
        ifs.close();
        ofs.close();
      };
      copyFile(dirIn_);
      copyFile(dirOut_);
      copyFile(dirPre_);
      copyFile(dirDiff_);
      copyFile(dirHWOut);
    }
    return result;
  }

}  // namespace trackerTFP
