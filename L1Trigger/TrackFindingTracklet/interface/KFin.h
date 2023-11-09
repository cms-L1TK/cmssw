#ifndef L1Trigger_TrackFindingTracklet_KFin_h
#define L1Trigger_TrackFindingTracklet_KFin_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "L1Trigger/TrackFindingTracklet/interface/ChannelAssignment.h"

#include <vector>

namespace trklet {

  /*! \class  trklet::KFin
   *  \brief  Class to emulate the data transformation happening betwwen DR and KF
   *  \author Thomas Schuh
   *  \date   2023, Feb
   */
  class KFin {
  public:
    KFin(const edm::ParameterSet& iConfig,
         const tt::Setup* setup_,
         const trackerTFP::DataFormats* dataFormats,
         const trackerTFP::LayerEncoding* layerEncoding,
         const ChannelAssignment* channelAssignment,
         int region);
    ~KFin() {}
    // read in and organize input tracks and stubs
    void consume(const tt::StreamsTrack& streamsTrack, const tt::StreamsStub& streamsStub);
    // fill output products
    void produce(tt::StreamsStub& accpetedStubs,
                 tt::StreamsTrack& acceptedTracks,
                 tt::StreamsStub& lostStubs,
                 tt::StreamsTrack& lostTracks);

  private:
    // truncates double precision of val into base precision, +1.e-12 restores robustness of addition of 2 digitised values
    double digi(double val, double base) const { return (floor(val / base + 1.e-12) + .5) * base; }
    struct Stub {
      Stub(const TTStubRef& ttStubRef, double r, double phi, double z, int layerId, bool psTilt, int channel)
          : ttStubRef_(ttStubRef), r_(r), phi_(phi), z_(z), layerId_(layerId), psTilt_(psTilt), channel_(channel) {}
      TTStubRef ttStubRef_;
      double r_;
      double phi_;
      double z_;
      int layerId_;
      bool psTilt_;
      int channel_;
      // phi uncertainty * sqrt(12) + additional terms in rad
      double dPhi_;
      // z uncertainty * sqrt(12) + additional terms in cm
      double dZ_;
    };
    struct Track {
      static constexpr int max_ = 7;
      Track() { stubs_.reserve(max_); }
      Track(const tt::FrameTrack& frame,
            const std::vector<Stub*>& stubs,
            double cot,
            double zT,
            const std::vector<Stub*>& stubs)
          : ttTrackRef_(ttTrackRef), valid_(valid), inv2R_(inv2R), phiT_(phiT), cot_(cot), zT_(zT), stubs_(stubs) {}
      TTTrackRef ttTrackRef_;
      bool valid_;
      TTBV maybe_;
      double inv2R_;
      double phiT_;
      double cot_;
      double zT_;
      double inv2R_;
      int sectorEta_;
    };
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;
    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // provides dataformats
    const trackerTFP::DataFormats* dataFormats_;
    // helper class to encode layer
    const trackerTFP::LayerEncoding* layerEncoding_;
    // helper class to assign tracks to channel
    const ChannelAssignment* channelAssignment_;
    // processing region (0 - 8) aka processing phi nonant
    const int region_;
    // storage of input tracks
    std::vector<Track> tracks_;
    // storage of input stubs
    std::vector<Stub> stubs_;
    // h/w liked organized pointer to input tracks
    std::vector<std::vector<Track*>> input_;
    // unified tracklet digitisation granularity
    double baseUinv2R_;
    double baseUphiT_;
    double baseUcot_;
    double baseUzT_;
    double baseUr_;
    double baseUphi_;
    double baseUz_;
    // KF input format digitisation granularity (identical to TMTT)
    double baseLinv2R_;
    double baseLphiT_;
    double baseLzT_;
    double baseLr_;
    double baseLphi_;
    double baseLz_;
    // Finer granularity (by powers of 2) than the TMTT one. Used to transform from Tracklet to TMTT base.
    double baseHinv2R_;
    double baseHphiT_;
    double baseHzT_;
    double baseHr_;
    double baseHphi_;
    double baseHz_;
    // digitisation granularity used for inverted cot(theta)
    double baseInvCot_;
  };

}  // namespace trklet

#endif
