#ifndef L1Trigger_TrackerTFP_DataFormats_h
#define L1Trigger_TrackerTFP_DataFormats_h

/*----------------------------------------------------------------------
Classes to calculate and provide dataformats used by Track Trigger emulator
enabling automated conversions from frames to stubs/tracks and vice versa
In data members of classes Stub* & Track* below, the variables describing
stubs/tracks are stored both in digitial format as a 64b word in frame_,
and in undigitized format in an std::tuple. (This saves CPU)
----------------------------------------------------------------------*/

#include "FWCore/Framework/interface/data_default_record_trait.h"
#include "L1Trigger/TrackerTFP/interface/DataFormatsRcd.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "DataFormats/L1TrackTrigger/interface/TTBV.h"

#include <vector>
#include <cmath>
#include <initializer_list>
#include <tuple>
#include <iostream>
#include <string>

namespace trackerTFP {

  // track trigger processes
  enum class Process { begin, fe = begin, dtc, pp, gp, ht, mht, zht, kfin, kf, dr, end, x };
  // track trigger variables
  enum class Variable { begin, r = begin, phi, z, dPhi, dZ, inv2R, phiT, cot, zT, layer, module, hitPattern, match, phi0, z0, end, x };
  // track trigger process order
  constexpr std::initializer_list<Process> Processes = {Process::fe, Process::dtc, Process::pp, Process::gp, Process::ht, Process::mht, Process::zht, Process::kfin, Process::kf, Process::dr};
  // conversion: Process to int
  inline constexpr int operator+(Process p) { return static_cast<int>(p); }
  // conversion: Variable to int
  inline constexpr int operator+(Variable v) { return static_cast<int>(v); }
  // increment of Process
  inline constexpr Process operator++(Process p) { return Process(+p + 1); }
  // increment of Variable
  inline constexpr Variable operator++(Variable v) { return Variable(+v + 1); }

  //Base class representing format of a variable
  class DataFormat {
  public:
    DataFormat(bool twos, bool biased = true) : twos_(twos), width_(0), base_(1.), range_(0.) {}
    ~DataFormat() {}
    // converts int to bitvector
    TTBV ttBV(int i) const { return TTBV(i, width_, twos_); }
    // converts double to bitvector
    TTBV ttBV(double d) const { return TTBV(d, base_, width_, twos_); }
    // extracts int from bitvector, removing these bits from bitvector
    void extract(TTBV& in, int& out) const { out = in.extract(width_, twos_); }
    // extracts double from bitvector, removing these bits from bitvector
    void extract(TTBV& in, double& out) const { out = in.extract(base_, width_, twos_); }
    // extracts double from bitvector, removing these bits from bitvector
    void extract(TTBV& in, TTBV& out) const { out = in.slice(width_, twos_); }
    // extracts bool from bitvector, removing these bits from bitvector
    void extract(TTBV& in, bool& out) const { out = in.extract(); }
    // attaches integer to bitvector
    void attach(const int i, TTBV& ttBV) const { ttBV += TTBV(i, width_, twos_); }
    // attaches double to bitvector
    void attach(const double d, TTBV& ttBV) const { ttBV += TTBV(d, base_, width_, twos_); }
    // attaches bitvector to bitvector
    void attach(const TTBV bv, TTBV& ttBV) const { ttBV += bv; }
    // converts int to double
    double floating(int i) const { return (i + .5) * base_; }
    // converts double to int
    int integer(double d) const { return std::floor(d / base_ + 1.e-12); }
    // converts double to int, shifts value to fit into range
    int integerF(double d) const { return integer(std::max((twos_ ? -limit() : 0), std::min(d, limit()))); }
    // converts double to int and back to double
    double digi(double d) const { return floating(integer(d)); }
    // converts binary integer value to twos complement integer value
    int toSigned(int i) const { return i - std::pow(2, width_) / 2; }
    // converts twos complement integer value to binary integer value
    int toUnsigned(int i) const { return i + std::pow(2, width_) / 2; }
    // converts floating point value to binary integer value
    int toUnsigned(double d) const { return this->integer(d) + std::pow(2, width_) / 2; }
    // biggest representable floating point value
    double limit() const { return (range_ - base_) / (twos_ ? 2. : 1.); }
    // returns false if data format would oferflow for this double value
    bool inRange(double d, bool digi = false) const {
      const double range = digi ? base_ * pow(2, width_) : range_;
      return d >= -range / 2. && d < range / 2.;
    }
    // returns false if data format would oferflow for this int value
    bool inRange(int i) const { return inRange(floating(i)); }
    // true if twos'complement or false if binary representation is chosen
    bool twos() const { return twos_; }
    // number of used bits
    int width() const { return width_; }
    // precision
    double base() const { return base_; }
    // covered range
    double range() const { return range_; }
  protected:
    // true if twos'complement or false if binary representation is chosen
    bool twos_;
    // number of used bits
    int width_;
    // precision
    double base_;
    // covered range
    double range_;
  };

  // class representing format of a specific variable
  template <Variable v, Process p>
  class Format : public DataFormat {
  public:
    Format(const tt::Setup* setup);
    ~Format() {}
  };

  template <>
  Format<Variable::r, Process::dtc>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::phi, Process::dtc>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::z, Process::dtc>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::layer, Process::dtc>::Format(const tt::Setup* setup);
  
  template <>
  Format<Variable::phi, Process::gp>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::z, Process::gp>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::phiT, Process::gp>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::cot, Process::gp>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::zT, Process::gp>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::module, Process::gp>::Format(const tt::Setup* setup);

  template <>
  Format<Variable::phi, Process::ht>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::inv2R, Process::ht>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::phiT, Process::ht>::Format(const tt::Setup* setup);

  template <>
  Format<Variable::inv2R, Process::mht>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::phiT, Process::mht>::Format(const tt::Setup* setup);

  template <>
  Format<Variable::cot, Process::zht>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::zT, Process::zht>::Format(const tt::Setup* setup);

  template <>
  Format<Variable::phi, Process::kfin>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::z, Process::kfin>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::dPhi, Process::kfin>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::dZ, Process::kfin>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::hitPattern, Process::kfin>::Format(const tt::Setup* setup);

  template <>
  Format<Variable::phi, Process::kf>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::z, Process::kf>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::inv2R, Process::kf>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::phiT, Process::kf>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::cot, Process::kf>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::zT, Process::kf>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::match, Process::kf>::Format(const tt::Setup* setup);

  template <>
  Format<Variable::inv2R, Process::dr>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::phi0, Process::dr>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::cot, Process::dr>::Format(const tt::Setup* setup);
  template <>
  Format<Variable::z0, Process::dr>::Format(const tt::Setup* setup);

  /*! \class  trackerTFP::DataFormats
   *  \brief  Class to calculate and provide dataformats used by Track Trigger emulator
   *  \author Thomas Schuh
   *  \date   2020, June
   */
  class DataFormats {
  private:
    // variable flavour mapping, Each row below declares which processing steps use the variable named in the comment at the end of the row
    static constexpr std::array<std::array<Process, +Process::end>, +Variable::end> config_ = {{
    //  Process::fe  Process::dtc  Process::pp   Process::gp   Process::ht   Process::mht  Process::zht  Process::kfin  Process::kf    Process::dr
      {{Process::x,  Process::dtc, Process::dtc, Process::dtc, Process::dtc, Process::dtc, Process::dtc, Process::dtc,  Process::dtc,  Process::x }},  // Variable::r
      {{Process::x,  Process::dtc, Process::dtc, Process::gp,  Process::ht,  Process::ht,  Process::ht,  Process::kfin, Process::kf,   Process::x }},  // Variable::phi
      {{Process::x,  Process::dtc, Process::dtc, Process::gp,  Process::gp,  Process::gp,  Process::gp,  Process::kfin, Process::kf,   Process::x }},  // Variable::z
      {{Process::x,  Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::kfin, Process::kfin, Process::x }},  // Variable::dPhi
      {{Process::x,  Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::kfin, Process::kfin, Process::x }},  // Variable::dZ
      {{Process::x,  Process::ht,  Process::ht,  Process::ht,  Process::ht,  Process::mht, Process::mht, Process::ht,   Process::kf,   Process::dr}},  // Variable::inv2R
      {{Process::x,  Process::gp,  Process::gp,  Process::gp,  Process::ht,  Process::mht, Process::mht, Process::ht,   Process::kf,   Process::x }},  // Variable::phiT
      {{Process::x,  Process::x,   Process::x,   Process::gp,  Process::x,   Process::x,   Process::zht, Process::x,    Process::kf,   Process::dr}},  // Variable::cot
      {{Process::x,  Process::x,   Process::gp,  Process::gp,  Process::gp,  Process::gp,  Process::zht, Process::gp,   Process::kf,   Process::x }},  // Variable::zT
      {{Process::x,  Process::dtc, Process::dtc, Process::dtc, Process::dtc, Process::dtc, Process::dtc, Process::x,    Process::x,    Process::x }},  // Variable::layer
      {{Process::x,  Process::x,   Process::x,   Process::gp,  Process::gp,  Process::gp,  Process::gp,  Process::x,    Process::x,    Process::x }},  // Variable::module
      {{Process::x,  Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::kfin, Process::x,    Process::x }},  // Variable::hitPattern
      {{Process::x,  Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,    Process::kf,   Process::x }},  // Variable::match
      {{Process::x,  Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,    Process::x,    Process::dr}},  // Variable::phi0
      {{Process::x,  Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,   Process::x,    Process::x,    Process::dr}}   // Variable::z0
    }};
    // stub word assembly, shows which stub variables are used by each process
    static constexpr std::array<std::initializer_list<Variable>, +Process::end> stubs_ = {{
        {},  // Process::fe
        {Variable::r, Variable::phi, Variable::z, Variable::layer, Variable::phiT, Variable::phiT, Variable::inv2R, Variable::inv2R}, // Process::dtc
        {Variable::r, Variable::phi, Variable::z, Variable::layer, Variable::phiT, Variable::phiT, Variable::inv2R, Variable::inv2R}, // Process::pp
        {Variable::r, Variable::phi, Variable::z, Variable::layer, Variable::module, Variable::inv2R, Variable::inv2R},               // Process::gp
        {Variable::r, Variable::phi, Variable::z, Variable::layer, Variable::module, Variable::phiT, Variable::zT},                   // Process::ht
        {Variable::r, Variable::phi, Variable::z, Variable::layer, Variable::module, Variable::phiT, Variable::zT},                   // Process::mht
        {Variable::r, Variable::phi, Variable::z, Variable::layer, Variable::module, Variable::phiT, Variable::zT},                   // Process::zht
        {Variable::r, Variable::phi, Variable::z, Variable::dPhi, Variable::dZ},                                                      // Process::kfin
        {Variable::r, Variable::phi, Variable::z, Variable::dPhi, Variable::dZ},                                                      // Process::kf
        {}                                                                                                                            // Process::dr
    }};
    // track word assembly, shows which track variables are used by each process
    static constexpr std::array<std::initializer_list<Variable>, +Process::end> tracks_ = {{
        {},                                                                              // Process::fe
        {},                                                                              // Process::dtc
        {},                                                                              // Process::pp
        {},                                                                              // Process::gp
        {},                                                                              // Process::ht
        {},                                                                              // Process::mht
        {},                                                                              // Process::zht
        {Variable::hitPattern, Variable::inv2R, Variable::phiT, Variable::zT},           // Process::kfin
        {Variable::match, Variable::inv2R, Variable::phiT, Variable::cot, Variable::zT}, // Process::kf
        {Variable::inv2R, Variable::phi0, Variable::cot, Variable::z0}                   // Process::dr
    }};

  public:
    DataFormats();
    DataFormats(const tt::Setup* setup);
    ~DataFormats() {}
    // converts bits to ntuple of variables
    template <typename... Ts>
    void convertStub(Process p, const tt::Frame& bv, std::tuple<Ts...>& data) const;
    // converts ntuple of variables to bits
    template <typename... Ts>
    void convertStub(Process p, const std::tuple<Ts...>& data, tt::Frame& bv) const;
    // converts bits to ntuple of variables
    template <typename... Ts>
    void convertTrack(Process p, const tt::Frame& bv, std::tuple<Ts...>& data) const;
    // converts ntuple of variables to bits
    template <typename... Ts>
    void convertTrack(Process p, const std::tuple<Ts...>& data, tt::Frame& bv) const;
    // access to run-time constants
    const tt::Setup* setup() const { return setup_; }
    // number of bits being used for specific variable flavour
    int width(Variable v, Process p) const { return formats_[+v][+p]->width(); }
    // precision being used for specific variable flavour
    double base(Variable v, Process p) const { return formats_[+v][+p]->base(); }
    // number of unused frame bits for a given Stub flavour
    int numUnusedBitsStubs(Process p) const { return numUnusedBitsStubs_[+p]; }
    // number of unused frame bits for a given Track flavour
    int numUnusedBitsTracks(Process p) const { return numUnusedBitsTracks_[+p]; }
    // number of channels of a given process on a TFP
    int numChannel(Process p) const { return numChannel_[+p]; }
    // number of channels of a given process for whole system
    int numStreams(Process p) const { return numStreams_[+p]; }
    //
    int numStreamsStubs(Process p) const { return numStreamsStubs_[+p]; }
    //
    int numStreamsTracks(Process p) const { return numStreamsTracks_[+p]; }
    // access to spedific format
    const DataFormat& format(Variable v, Process p) const { return *formats_[+v][+p]; }
  private:
    // number of unique data formats
    int numDataFormats_;
    // method to count number of unique data formats
    template <Variable v = Variable::begin, Process p = Process::begin>
    void countFormats();
    // constructs data formats of all unique used variables and flavours
    template <Variable v = Variable::begin, Process p = Process::begin>
    void fillDataFormats();
    // helper (loop) data formats of all unique used variables and flavours
    template <Variable v, Process p, Process it = Process::begin>
    void fillFormats();
    // helper (loop) to convert bits to ntuple of variables
    template <int it = 0, typename... Ts>
    void extractStub(Process p, TTBV& ttBV, std::tuple<Ts...>& data) const;
    // helper (loop) to convert bits to ntuple of variables
    template <int it = 0, typename... Ts>
    void extractTrack(Process p, TTBV& ttBV, std::tuple<Ts...>& data) const;
    // helper (loop) to convert ntuple of variables to bits
    template <int it = 0, typename... Ts>
    void attachStub(Process p, const std::tuple<Ts...>& data, TTBV& ttBV) const;
    // helper (loop) to convert ntuple of variables to bits
    template <int it = 0, typename... Ts>
    void attachTrack(Process p, const std::tuple<Ts...>& data, TTBV& ttBV) const;
    // configuration during construction
    edm::ParameterSet iConfig_;
    // stored run-time constants
    const tt::Setup* setup_;
    // collection of unique formats
    std::vector<DataFormat> dataFormats_;
    // variable flavour mapping
    std::vector<std::vector<DataFormat*>> formats_;
    // number of unused frame bits for a all Stub flavours
    std::vector<int> numUnusedBitsStubs_;
    // number of unused frame bits for a all Track flavours
    std::vector<int> numUnusedBitsTracks_;
    // number of channels of all processes on a TFP
    std::vector<int> numChannel_;
    // number of channels of all processes for whole system
    std::vector<int> numStreams_;
    //
    std::vector<int> numStreamsStubs_;
    //
    std::vector<int> numStreamsTracks_;
  };

  // base class to represent stubs
  template <typename... Ts>
  class Stub {
  public:
    // construct Stub from Frame
    Stub(const tt::FrameStub& frame, const DataFormats* dataFormats, Process p);
    template <typename... Others>
    // construct Stub from other Stub
    Stub(const Stub<Others...>& stub, Ts... data);
    // construct Stub from TTStubRef
    Stub(const TTStubRef& ttStubRef, const DataFormats* dataFormats, Process p, Ts... data);
    Stub() {}
    ~Stub() {}
    // true if frame valid, false if gap in data stream
    explicit operator bool() const { return frame_.first.isNonnull(); }
    // access to DataFormats
    const DataFormats* dataFormats() const { return dataFormats_; }
    // stub flavour
    Process p() const { return p_; }
    // acess to frame
    const tt::FrameStub& frame() const { return frame_; }
    // access to TTStubRef
    const TTStubRef& ttStubRef() const { return frame_.first; }
    // access to bitvector
    const tt::Frame& bv() const { return frame_.second; }
    // id of collection this stub belongs to
    int trackId() const { return trackId_; }
  protected:
    // access to run-time constants
    const tt::Setup* setup() const { return dataFormats_->setup(); }
    // format of given variable
    const DataFormat& format(Variable v) const { return dataFormats_->format(v, p_); }
    // format of given variable and process
    const DataFormat& format(Variable v, Process p) const { return dataFormats_->format(v, p); }
    //number of bits uesd of given variable
    int width(Variable v) const { return dataFormats_->width(v, p_); }
    // precision of given variable
    double base(Variable v) const { return dataFormats_->base(v, p_); }
    // all dataformats
    const DataFormats* dataFormats_;
    // stub flavour
    Process p_;
    // underlying TTStubRef and bitvector
    tt::FrameStub frame_;
    // ntuple of variables this stub is assemled of
    std::tuple<Ts...> data_;
    // id of collection this stub belongs to
    int trackId_;
  };

  // class to represent stubs generated by process DTC
  class StubDTC : public Stub<double, double, double, int, int, int, int, int> {
  public:
    // construct StubDTC from TTStubRef
    StubDTC(const TTStubRef& ttStubRef, const DataFormats* dataFormats, int region, int layer, double r, double phi, double z, std::pair<double, double> phiT, std::pair<double, double> inv2R);
    ~StubDTC() {}
    // stub radius wrt chosenRofPhi
    double r() const { return std::get<0>(data_); }
    // stub phi wrt processing nonant centre
    double phi() const { return std::get<1>(data_); }
    // stub z
    double z() const { return std::get<2>(data_); }
    // reduced layer id
    int layer() const { return std::get<3>(data_); }
    // first phi sector this stub belongs to
    int phiTMin() const { return std::get<4>(data_); }
    // last phi sector this stub belongs to
    int phiTMax() const { return std::get<5>(data_); }
    // first inv2R bin this stub belongs to
    int inv2RMin() const { return std::get<6>(data_); }
    // last inv2R bin this stub belongs to
    int inv2RMax() const { return std::get<7>(data_); }
  };

  // class to represent stubs generated by process patch pannel
  class StubPP : public Stub<double, double, double, int, int, int, int, int> {
  public:
    // construct StubPP from Frame
    StubPP(const tt::FrameStub& frame, const DataFormats* dataFormats);
    ~StubPP() {}
    // true if stub belongs to given sector
    bool inSector(int sector) const { return sectors_[sector]; }
    // sectors this stub belongs to
    std::vector<int> sectors() const { return sectors_.ids(); }
    // stub radius wrt chosenRofPhi
    double r() const { return std::get<0>(data_); }
    // stub phi wrt processing nonant centre
    double phi() const { return std::get<1>(data_); }
    // stub z
    double z() const { return std::get<2>(data_); }
    // reduced layer id
    int layer() const { return std::get<3>(data_); }
    // first phi sector this stub belongs to
    int phiTMin() const { return std::get<4>(data_); }
    // last phi sector this stub belongs to
    int phiTMax() const { return std::get<5>(data_); }
    // first inv2R bin this stub belongs to
    int inv2RMin() const { return std::get<6>(data_); }
    // last inv2R bin this stub belongs to
    int inv2RMax() const { return std::get<7>(data_); }
  private:
    // sectors this stub belongs to
    TTBV sectors_;
  };

  // class to represent stubs generated by process geometric processor
  class StubGP : public Stub<double, double, double, int, TTBV, int, int> {
  public:
    // construct StubGP from Frame
    StubGP(const tt::FrameStub& frame, const DataFormats* dataFormats, int channel);
    // construct StubGP from StubPP
    StubGP(const StubPP& stub, int phiT, int zT);
    ~StubGP() {}
    // true if stub belongs to given inv2R bin
    bool inInv2RBin(int inv2RBin) const { return inv2RBins_[inv2RBin]; }
    // inv2R bins  this stub belongs to
    std::vector<int> inv2RBins() const { return inv2RBins_.ids(); }
    // stub phi sector
    int phiT() const { return phiT_; }
    // stub eta sector
    int zT() const { return zT_; }
    // stub radius wrt chosenRofPhi
    double r() const { return std::get<0>(data_); }
    // stub phi wrt phi sector centre
    double phi() const { return std::get<1>(data_); }
    // stub z residual wrt eta sector
    double z() const { return std::get<2>(data_); }
    // reduced layer id
    int layer() const { return std::get<3>(data_); }
    // module type
    const TTBV& module() const { return std::get<4>(data_); }
    // first inv2R bin this stub belongs to
    int inv2RMin() const { return std::get<5>(data_); }
    // last inv2R bin this stub belongs to
    int inv2RMax() const { return std::get<6>(data_); }
  private:
    // inv2R bins this stub belongs to
    TTBV inv2RBins_;
    // stub phi sector
    int phiT_;
    // stub eta sector
    int zT_;
  };

  // class to represent stubs generated by process hough transform
  class StubHT : public Stub<double, double, double, int, TTBV, int, int> {
  public:
    // construct StubHT from Frame
    StubHT(const tt::FrameStub& frame, const DataFormats* dataFormats, int channel);
    // construct StubHT from StubGP and HT cell assignment
    StubHT(const StubGP& stub, int inv2R, int phiT);
    ~StubHT() {}
    // stub radius wrt chosenRofPhi
    double r() const { return std::get<0>(data_); };
    // stub phi residual wrt track parameter
    double phi() const { return std::get<1>(data_); };
    // stub z residual wrt eta sector
    double z() const { return std::get<2>(data_); };
    // reduced layer id
    int layer() const { return std::get<3>(data_); };
    // module type
    const TTBV& module() const { return std::get<4>(data_); }
    // stub phi at radius chosenRofPhi wrt processing nonant centre
    int phiT() const { return std::get<5>(data_); };
    // eta sector
    int zT() const { return std::get<6>(data_); };
    // stub inv2R wrt HT inv2R bin center
    int inv2R() const { return inv2R_; };
    // stub phiT wrt sector center
    int phiTlocal() const { return phiT_; }
    //
    bool ps() const { return setup()->ps(this->module()); }
    //
    bool barrel() const { return setup()->barrel(this->module()); }
    //
    bool tilted() const { return setup()->tilted(this->module()); }
    // stub phi uncertainty in rad
    double dPhi() const { return dPhi_; }
    //
    void mark(int inv2R, int phiT) { array_[inv2R].set(phiT); }
    //
    void hit(int inv2R, int phiT) { if (!valid_ && array_[inv2R][phiT]) valid_ = true; }
    //
    bool valid() const { return valid_; }
  private:
    // stub inv2R wrt HT inv2R bin center
    int inv2R_;
    // stub phiT wrt sector center
    int phiT_;
    // stub phi uncertainty in rad
    double dPhi_;
    //
    bool valid_;
    //
    std::vector<TTBV> array_;
  };

  // class to represent stubs generated by process mini hough transform
  class StubMHT : public Stub<double, double, double, int, TTBV, int, int> {
  public:
    // construct StubMHT from Frame
    StubMHT(const tt::FrameStub& frame, const DataFormats* dataFormats);
    //
    StubMHT(const StubHT& stub, bool last);
    ~StubMHT() {}
    // stub radius wrt choenRofPhi
    double r() const { return std::get<0>(data_); }
    // stub phi residual wrt finer track parameter
    double phi() const { return std::get<1>(data_); }
    // stub z rsidual wrt eta sector
    double z() const { return std::get<2>(data_); }
    // reduced layer id
    int layer() const { return std::get<3>(data_); }
    // module type
    const TTBV& module() const { return std::get<4>(data_); }
    // stub phi at radius chosenRofPhi wrt processing nonant centre
    int phiT() const { return std::get<5>(data_); }
    // eta sector
    int zT() const { return std::get<6>(data_); }
    //
    bool ps() const { return setup()->ps(this->module()); }
    //
    bool barrel() const { return setup()->barrel(this->module()); }
    //
    bool tilted() const { return setup()->tilted(this->module()); }
    // stub phi uncertainty in rad
    double dPhi() const { return dPhi_; }
    // stub z uncertainty in cm
    double dZ() const { return dZ_; }
    //
    bool valid() const { return valid_; }
    //
    bool last() const { return last_; }
    //
    void mark(int cot, int zT) { array_[cot].set(zT); }
    //
    void hit(int cot, int zT) { if (!valid_ && array_[cot][zT]) valid_ = true; }
  private:
    // stub phi uncertainty in rad
    double dPhi_;
    // stub z uncertainty in cm
    double dZ_;
    //
    bool valid_;
    //
    bool last_;
    //
    std::vector<TTBV> array_;
  };

  // class to represent stubs generated by process z hough transform
  class StubZHT : public Stub<double, double, double, int, TTBV, int, int> {
  public:
    // construct StubZHT from Frame
    StubZHT(const tt::FrameStub& frame, const DataFormats* dataFormats, int channel);
    //
    StubZHT(const StubMHT& stub, bool last);
    ~StubZHT() {}
    // stub radius wrt chonseRofPhi
    double r() const { return std::get<0>(data_); }
    // stub phiresiudal wrt finer track parameter
    double phi() const { return std::get<1>(data_); }
    // stub z residual to track parameter
    double z() const { return std::get<2>(data_); }
    // reduced layer id
    int layer() const { return std::get<3>(data_); }
    // module type
    const TTBV& module() const { return std::get<4>(data_); }
    // stub phi at radius chosenRofPhi wrt processing nonant centre
    int phiT() const { return std::get<5>(data_); }
    // stub z at radius chosenRofZ
    int zT() const { return std::get<6>(data_); }
    //
    bool ps() const { return setup()->ps(this->module()); }
    //
    bool barrel() const { return setup()->barrel(this->module()); }
    //
    bool tilted() const { return setup()->tilted(this->module()); }
    // stub z uncertainty in cm
    double dZ() const { return dZ_; }
    // stub inv2R
    int inv2R() const { return inv2R_; }
    //
    bool valid() const { return valid_; }
    //
    bool last() const { return last_; }
  private:
    // stub z uncertainty in cm
    double dZ_;
    // stub inv2R
    int inv2R_;
    //
    bool valid_;
    //
    bool last_;
  };

  // class to represent stubs generated by process kfin
  class StubKFin : public Stub<double, double, double, double, double> {
  public:
    // construct StubKFin from Frame
    StubKFin(const tt::FrameStub& frame, const DataFormats* dataFormats, int layer);
    // construct StubKFin from StubZHT
    StubKFin(const StubZHT& stub, int layer);
    // construct StubKFin from TTStubRef
    StubKFin(const TTStubRef& ttStubRef,
             const DataFormats* dataFormats,
             double r,
             double phi,
             double z,
             double dPhi,
             double dZ,
             int layer);
    ~StubKFin() {}
    // kf layer id
    int layer() const { return layer_; }
    // stub radius wrt chosenRofPhi
    double r() const { return std::get<0>(data_); }
    // stub phi residual wrt finer track parameter
    double phi() const { return std::get<1>(data_); }
    // stub z residual wrt track parameter
    double z() const { return std::get<2>(data_); }
    // stub phi uncertainty
    double dPhi() const { return std::get<3>(data_); }
    // stub z uncertainty
    double dZ() const { return std::get<4>(data_); }
  private:
    // kf layer id
    int layer_;
  };

  // class to represent stubs generated by process kalman filter
  class StubKF : public Stub<double, double, double, double, double> {
  public:
    // construct StubKF from Frame
    StubKF(const tt::FrameStub& frame, const DataFormats* dataFormats, int layer);
    // construct StubKF from StubKFin
    StubKF(const StubKFin& stub, double inv2R, double phiT, double cot, double zT);
    ~StubKF() {}
    // kf layer id
    int layer() const { return layer_; }
    // stub radius wrt choenRofPhi
    double r() const { return std::get<0>(data_); }
    // stub phi residual wrt fitted parameter
    double phi() const { return std::get<1>(data_); }
    // stub z residual wrt fitted parameter
    double z() const { return std::get<2>(data_); }
    // stub phi uncertainty
    double dPhi() const { return std::get<3>(data_); }
    // stub z uncertainty
    double dZ() const { return std::get<4>(data_); }
  private:
    // kf layer id
    int layer_;
  };

  // base class to represent tracks
  template <typename... Ts>
  class Track {
  public:
    // construct Track from Frame
    Track(const tt::FrameTrack& frame, const DataFormats* dataFormats, Process p);
    // construct Track from other Track
    template <typename... Others>
    Track(const Track<Others...>& track, Ts... data);
    // construct Track from TTTrackRef
    Track(const TTTrackRef& ttTrackRef, const DataFormats* dataFormats, Process p, Ts... data);
    ~Track() {}
    // true if frame valid, false if gap in data stream
    explicit operator bool() const { return frame_.first.isNonnull(); }
    // access to DataFormats
    const DataFormats* dataFormats() const { return dataFormats_; }
    // track flavour
    Process p() const { return p_; }
    // acces to frame
    const tt::FrameTrack& frame() const { return frame_; }
    // access to TTTrackRef
    const TTTrackRef& ttTrackRef() const { return frame_.first; }
    // access to bitvector
    const tt::Frame& bv() const { return frame_.second; }
  protected:
    // access to run-time constants
    const tt::Setup* setup() const { return dataFormats_->setup(); }
    // format of given variable
    const DataFormat& format(Variable v) const { return dataFormats_->format(v, p_); }
    // format of given variable and process
    const DataFormat& format(Variable v, Process p) const { return dataFormats_->format(v, p); }
    //number of bits uesd of given variable
    int width(Variable v) const { return dataFormats_->width(v, p_); }
    // precision of given variable
    double base(Variable v) const { return dataFormats_->base(v, p_); }
    // all data formats
    const DataFormats* dataFormats_;
    // track flavour
    Process p_;
    // underlying TTTrackRef and bitvector
    tt::FrameTrack frame_;
    // ntuple of variables this track is assemled of
    std::tuple<Ts...> data_;
  };

  class TrackKFin : public Track<TTBV, double, double, double> {
  public:
    // construct TrackKFin from Frame
    TrackKFin(const tt::FrameTrack& frame, const DataFormats* dataFormats, const std::vector<StubKFin*>& stubs);
    // construct TrackKFin from TTTrackRef
    TrackKFin(const TTTrackRef& ttTrackRef, const DataFormats* dataFormats, const TTBV& maybePattern, int inv2R, int phiT, int zT);
    ~TrackKFin() {}
    // pattern of layers which are only maybe crossed by found candidate
    const TTBV& maybePattern() const { return std::get<0>(data_); }
    // track inv2R
    double inv2R() const { return std::get<1>(data_); }
    // track phi at radius chosenRofPhi wrt pprocessing centre
    double phiT() const { return std::get<2>(data_); }
    // track z at radius chosenRofZ
    double zT() const { return std::get<3>(data_); }
    //
    const TTBV& hitPattern() const { return hitPattern_; }
    // true if given layer has a hit
    bool hitPattern(int index) const { return hitPattern_[index]; }
    // true if given layer has a hit or is a maybe layer
    bool maybePattern(int index) const { return hitPattern_[index] || maybePattern()[index]; }
    // stubs on a given layer
    const std::vector<StubKFin*>& layerStubs(int layer) const { return stubs_[layer]; }
    // firts stub on a given layer
    StubKFin* layerStub(int layer) const { return stubs_[layer].front(); }
    // selection of ttStubRefs for given hit ids on given layers
    std::vector<TTStubRef> ttStubRefs(const TTBV& hitPattern, const std::vector<int>& layerMap) const;
    // stubs organized in layer
    const std::vector<std::vector<StubKFin*>>& stubs() const { return stubs_; }
    // cotTheta
    double cot() const { return cot_; }
  private:
    // stubs organized in layer
    std::vector<std::vector<StubKFin*>> stubs_;
    //
    TTBV hitPattern_;
    // cotTheta
    double cot_;
  };

  // class to represent tracks generated by process kalman filter
  class TrackKF : public Track<bool, double, double, double, double> {
  public:
    // construct TrackKF from Frame
    TrackKF(const tt::FrameTrack& frame, const DataFormats* dataFormats);
    // construct TrackKF from TrackKFKFin
    TrackKF(const TrackKFin& track, const std::vector<StubKF*> stubs, int trackId, double inv2R, double phiT, double cot, double zT);
    ~TrackKF() {}
    // true if kf prameter consistent with mht parameter
    bool match() const { return std::get<0>(data_); }
    // track qOver pt
    double inv2R() const { return std::get<1>(data_); }
    // track phi at radius chosenRofPhi wrt processing nonant centre
    double phiT() const { return std::get<2>(data_); }
    // track cotTheta wrt eta sector cotTheta
    double cot() const { return std::get<3>(data_); }
    // track z at radius chosenRofZ
    double zT() const { return std::get<4>(data_); }
    // conversion to TTTrack with given stubs
    TTTrack<Ref_Phase2TrackerDigi_> ttTrack(int region, const std::vector<StubKF>& stubs) const;
    // stubs forming this track
    const std::vector<StubKF*>& stubs() const { return stubs_; }
    // Track id
    int trackId() const { return trackId_; }
    // hitPattern
    const TTBV& hitPattern() const { return hitPattern_; }
    // number of skipped layers
    int numSkippedLayers() const { return numSkippedLayers_; }
    // number of consitent layers
    int numConsistentLayers() const { return numConsistentLayers_; }
    // global cot
    double cotGlobal() const { return cot_; }
  private:
    // stubs forming this track
    std::vector<StubKF*> stubs_;
    // Track id
    int trackId_;
    // hitPattern
    TTBV hitPattern_;
    // number of skipped layers
    int numSkippedLayers_;
    // number of consitent layers
    int numConsistentLayers_;
    // global cot
    double cot_;
  };

  //Class to represent KFout 96-bit track for use in distribution server
  class TrackKFOut {
  public:
    TrackKFOut() : TrackKFOut(0, 0, 0, 0, tt::FrameTrack(), 0, 0, false) {}
    // construct TrackKF from Partial Tracks
    TrackKFOut(TTBV PartialTrack1,
               TTBV PartialTrack2,
               TTBV PartialTrack3,
               int sortKey,
               const tt::FrameTrack& track,
               int trackID,
               int linkID,
               bool valid)
        : PartialTrack1_(PartialTrack1),
          PartialTrack2_(PartialTrack2),
          PartialTrack3_(PartialTrack3),
          sortKey_(sortKey),
          track_(track),
          trackID_(trackID),
          linkID_(linkID),
          valid_(valid){};
    ~TrackKFOut() {}
    int sortKey() const { return sortKey_; }
    bool dataValid() const { return valid_; }
    int trackID() const { return trackID_; }
    int linkID() const { return linkID_; }
    TTBV PartialTrack1() const { return PartialTrack1_; }
    TTBV PartialTrack2() const { return PartialTrack2_; }
    TTBV PartialTrack3() const { return PartialTrack3_; }
    tt::FrameTrack track() const { return track_; }
  private:
    TTBV PartialTrack1_;
    TTBV PartialTrack2_;
    TTBV PartialTrack3_;
    int sortKey_;
    tt::FrameTrack track_;
    int trackID_;
    int linkID_;
    bool valid_;
  };

  typedef std::vector<TrackKFOut> TrackKFOutSACollection;
  typedef std::shared_ptr<TrackKFOut> TrackKFOutSAPtr;
  typedef std::vector<TrackKFOutSAPtr> TrackKFOutSAPtrCollection;
  typedef std::vector<std::vector<std::shared_ptr<TrackKFOut>>> TrackKFOutSAPtrCollections;
  typedef std::vector<std::vector<std::vector<std::shared_ptr<TrackKFOut>>>> TrackKFOutSAPtrCollectionss;

  // class to represent tracks generated by process duplicate removal
  class TrackDR : public Track<double, double, double, double> {
  public:
    // construct TrackDR from Frame
    TrackDR(const tt::FrameTrack& frame, const DataFormats* dataFormats);
    // construct TrackDR from TrackKF
    TrackDR(const TrackKF& track);
    ~TrackDR() {}
    // track inv2R
    double inv2R() const { return std::get<0>(data_); }
    // track phi at radius 0 wrt processing nonant centre
    double phi0() const { return std::get<1>(data_); }
    // track cotThea
    double cot() const { return std::get<2>(data_); }
    // track z at radius 0
    double z0() const { return std::get<3>(data_); }
    // conversion to TTTrack
    TTTrack<Ref_Phase2TrackerDigi_> ttTrack() const;
  private:
  };

}  // namespace trackerTFP

EVENTSETUP_DATA_DEFAULT_RECORD(trackerTFP::DataFormats, trackerTFP::DataFormatsRcd);

#endif