#ifndef L1Trigger_TrackFindingTracklet_TrackQuality_h
#define L1Trigger_TrackFindingTracklet_TrackQuality_h

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackFindingTracklet/interface/DataFormats.h"

#include <vector>

namespace trklet {

  /*! \class  trklet::TrackQuality
   *  \brief  Bit accurate emulation of the track quality BDT including calculation of chi2s.
   *  \author Thomas Schuh
   *  \date   Aug 2025
   */
  class TrackQuality {
  public:
    // various internal dataformats
    struct InternalFormats {
      DataFormat m02_;
      DataFormat m12_;
      DataFormat invV0_;
      DataFormat invV1_;
      DataFormat chi20_;
      DataFormat chi21_;
      DataFormat chi20BDT_;
      DataFormat chi21BDT_;
    };
    TrackQuality(const tt::Setup* setup, const DataFormats* df, const InternalFormats& internal, int region)
        : setup_(setup), dataFormats_(df), internalFormats_(&internal), region_(region) {}
    ~TrackQuality() = default;
    // read in and organize input tracks and stubs
    void consume(const tt::StreamsTrack&, const tt::StreamsStub&);
    // fills output products
    void produce(tt::StreamTrack&) const;

  private:
    // representation of an input Frame
    struct Frame {
      Frame() : track_(nullptr), stubs_() {}
      Frame(TrackKF* track, int n) : track_(track), stubs_(n, nullptr) {}
      TrackKF* track_;
      std::vector<StubKF*> stubs_;
    };
    // helper class to store configurations
    const tt::Setup* setup_;
    // provides dataformats
    const DataFormats* dataFormats_;
    // collection of internal formats
    const InternalFormats* internalFormats_;
    // region id
    int region_;
    // KF tracks
    std::vector<TrackKF> tracks_;
    // KF stubs
    std::vector<StubKF> stubs_;
    // input data
    std::vector<Frame> input_;
  };

}  // namespace trklet

#endif
