#ifndef L1Trigger_TrackFindingTracklet_KFin_h
#define L1Trigger_TrackFindingTracklet_KFin_h

#include "L1Trigger/TrackTrigger/interface/Setup.h"
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
    struct Track {
      static constexpr int max_ = 7;
      Track(const tt::FrameTrack& frame, const std::vector<tt::FrameStub>& stubs) : frame_(frame), stubs_(stubs) {}
      tt::FrameTrack frame_;
      std::vector<tt::FrameStub> stubs_ = std::vector<tt::FrameStub>(max_);
    };
    // remove and return first element of deque, returns nullptr if empty
    template <class T>
    T* pop_front(std::deque<T*>& ts) const;
    // true if truncation is enbaled
    bool enableTruncation_;
    // provides run-time constants
    const tt::Setup* setup_;
    // helper class to assign tracks to channel
    const ChannelAssignment* channelAssignment_;
    // processing region (0 - 8) aka processing phi nonant
    const int region_;
    // storage of input tracks
    std::vector<Track> tracks_;
    // h/w liked organized pointer to input tracks
    std::vector<std::vector<Track*>> input_;
  };

}  // namespace trklet

#endif
