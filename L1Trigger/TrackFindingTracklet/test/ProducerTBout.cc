#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackBuilderChannel.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <numeric>

using namespace std;
using namespace edm;
using namespace trackerTFP;
using namespace tt;

namespace trackFindingTracklet {

  /*! \class  trackFindingTracklet::ProducerTBout
   *  \brief  Transforms TTTracks from Tracklet pattern reco. into f/w comparable format
   *  \author Thomas Schuh
   *  \date   2021, Oct
   */
  class ProducerTBout : public stream::EDProducer<> {
  public:
    explicit ProducerTBout(const ParameterSet&);
    ~ProducerTBout() override {}

  private:
    virtual void beginRun(const Run&, const EventSetup&) override;
    virtual void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}

    // dtc stub
    struct Stub {
      Stub(const FrameStub& frame, int region, int channel, const Setup* setup) {
        region_ = region;
        ttStubRef_ = frame.first;
        frame_ = frame.second;
        gp_ = setup->stubPos(true, frame, region, channel);
        barrel_ = setup->barrel(ttStubRef_);
        const bool ps = setup->psModule(ttStubRef_);
        const TTBV ttBV(frame.second);
        r_ = TTBV(ttBV, 39, 39 - (barrel_ || ps ? 7 : 12), true);
        if (!barrel_ && !ps)
          r_.resize(12);
      }
      int region_;
      TTStubRef ttStubRef_;
      Frame frame_;
      GlobalPoint gp_;
      bool barrel_;
      TTBV r_;
    };
    // return h/w bits for given ttTrackRef
    Frame conv(const TTTrackRef& ttTrackRef) const;
    // return h/w bits for given dtc stub and origin ttTrackRef
    Frame conv(const TTTrackRef& ttTrackRef, const Stub& stub) const;

    // ED input token of TTTracks
    EDGetTokenT<TTTracks> edGetTokenTTTracks_;
    // ED input token of DTC Stubs
    EDGetTokenT<TTDTC> edGetTokenTTDTC_;
    // ED output token for stubs
    EDPutTokenT<StreamsStub> edPutTokenAcceptedStubs_;
    EDPutTokenT<StreamsStub> edPutTokenLostStubs_;
    // ED output token for tracks
    EDPutTokenT<StreamsTrack> edPutTokenAcceptedTracks_;
    EDPutTokenT<StreamsTrack> edPutTokenLostTracks_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // DataFormats token
    ESGetToken<DataFormats, DataFormatsRcd> esGetTokenDataFormats_;
    // TrackBuilderChannel token
    ESGetToken<TrackBuilderChannel, TrackBuilderChannelRcd> esGetTokenTrackBuilderChannel_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_;
    // helper class to extract structured data from TTDTC::Frames
    const DataFormats* dataFormats_;
    // helper class to assign tracks to channel
    TrackBuilderChannel* trackBuilderChannel_;
    //
    bool enableTruncation_;
    //
    trklet::Settings settings_;
  };

  ProducerTBout::ProducerTBout(const ParameterSet& iConfig) :
    iConfig_(iConfig)
  {
    const InputTag& inputTag = iConfig.getParameter<InputTag>("InputTag");
    const InputTag& inputTagDTC = iConfig.getParameter<InputTag>("InputTagDTC");
    const string& branchAcceptedStubs = iConfig.getParameter<string>("BranchAcceptedStubs");
    const string& branchAcceptedTracks = iConfig.getParameter<string>("BranchAcceptedTracks");
    const string& branchLostStubs = iConfig.getParameter<string>("BranchLostStubs");
    const string& branchLostTracks = iConfig.getParameter<string>("BranchLostTracks");
    // book in- and output ED products
    edGetTokenTTTracks_ = consumes<TTTracks>(inputTag);
    edGetTokenTTDTC_ = consumes<TTDTC>(inputTagDTC);
    edPutTokenAcceptedStubs_ = produces<StreamsStub>(branchAcceptedStubs);
    edPutTokenAcceptedTracks_ = produces<StreamsTrack>(branchAcceptedTracks);
    edPutTokenLostStubs_ = produces<StreamsStub>(branchLostStubs);
    edPutTokenLostTracks_ = produces<StreamsTrack>(branchLostTracks);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    esGetTokenTrackBuilderChannel_ = esConsumes<TrackBuilderChannel, TrackBuilderChannelRcd, Transition::BeginRun>();
    // initial ES products
    setup_ = nullptr;
    dataFormats_ = nullptr;
    trackBuilderChannel_ = nullptr;
    //
    enableTruncation_ = iConfig.getParameter<bool>("EnableTruncation");
  }

  void ProducerTBout::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    if (!setup_->configurationSupported())
      return;
    // check process history if desired
    if (iConfig_.getParameter<bool>("CheckHistory"))
      setup_->checkHistory(iRun.processHistory());
    // helper class to extract structured data from TTDTC::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    // helper class to assign tracks to channel
    trackBuilderChannel_ = const_cast<TrackBuilderChannel*>(&iSetup.getData(esGetTokenTrackBuilderChannel_));
  }

  void ProducerTBout::produce(Event& iEvent, const EventSetup& iSetup) {
    const int numStreamsTracks = setup_->numRegions() * trackBuilderChannel_->numChannels();
    const int numStreamsStubs = numStreamsTracks * trackBuilderChannel_->maxNumProjectionLayers();
    // empty KFin products
    StreamsStub streamAcceptedStubs(numStreamsStubs);
    StreamsTrack streamAcceptedTracks(numStreamsTracks);
    StreamsStub streamLostStubs(numStreamsStubs);
    StreamsTrack streamLostTracks(numStreamsTracks);
    // read in hybrid track finding product and produce KFin product
    if (setup_->configurationSupported()) {
      // create DTC stub frames
      Handle<TTDTC> handleTTDTC;
      iEvent.getByToken<TTDTC>(edGetTokenTTDTC_, handleTTDTC);
      const TTDTC& ttDTC = *handleTTDTC;
      vector<Stub> stubsDTC;
      stubsDTC.reserve(ttDTC.nStubs());
      for (int region : ttDTC.tfpRegions())
        for (int channel : ttDTC.tfpChannels())
          for (const FrameStub& frame : ttDTC.stream(region, channel))
            if (frame.first.isNonnull())
              stubsDTC.emplace_back(frame, region, channel, setup_);
      // create TTrackRefs
      Handle<TTTracks> handleTTTracks;
      iEvent.getByToken<TTTracks>(edGetTokenTTTracks_, handleTTTracks);
      vector<TTTrackRef> ttTrackRefs;
      ttTrackRefs.reserve(handleTTTracks->size());
      for (int i = 0; i < (int)handleTTTracks->size(); i++)
        ttTrackRefs.emplace_back(TTTrackRef(handleTTTracks, i));
      // count tracks per channel and size output products
      vector<int> nTTTracksStreams(numStreamsTracks, 0);
      int channelId;
      for (const TTTrackRef& ttTrackRef : ttTrackRefs)
        if (trackBuilderChannel_->channelId(ttTrackRef, channelId))
          nTTTracksStreams[channelId]++;
      for (int channelTrack = 0; channelTrack < numStreamsTracks; channelTrack++) {
        const int num = nTTTracksStreams[channelTrack];
        const int lost = enableTruncation_ && num > setup_->numFrames() ? num - setup_->numFrames() : 0;
        const int accepted = lost == 0 ? num : setup_->numFrames();
        streamAcceptedTracks[channelTrack].reserve(accepted);
        streamLostTracks[channelTrack].reserve(lost);
        for (int projection = 0; projection < trackBuilderChannel_->maxNumProjectionLayers(); projection++) {
          const int channelStub = channelTrack * trackBuilderChannel_->maxNumProjectionLayers() + projection;
          streamAcceptedStubs[channelStub].reserve(accepted);
          streamLostStubs[channelStub].reserve(lost);
        }
      }
      // fill output products
      for (const TTTrackRef& ttTrackRef : ttTrackRefs) {
        const bool valid = trackBuilderChannel_->channelId(ttTrackRef, channelId);
        const bool truncate = enableTruncation_ && (int)streamAcceptedTracks[channelId].size() > setup_->numFrames();
        StreamTrack& tracks = truncate ? streamLostTracks[channelId] : streamAcceptedTracks[channelId];
        if (!valid && !truncate) { // fill gap
          tracks.emplace_back(FrameTrack());
          for (int projection = 0; projection < trackBuilderChannel_->maxNumProjectionLayers(); projection++) {
            const int channelStub = channelId * trackBuilderChannel_->maxNumProjectionLayers() + projection;
            streamAcceptedStubs[channelStub].emplace_back(FrameStub());
          }
          continue;
        }
        // conv track word
        tracks.emplace_back(ttTrackRef, conv(ttTrackRef));
        // conv stub words
        StreamsStub& streams = truncate ? streamLostStubs : streamAcceptedStubs;
        TTBV pattern(0, trackBuilderChannel_->maxNumProjectionLayers());
        int layerId;
        for (const TTStubRef& ttStubRef : ttTrackRef->getStubRefs()) {
          if (!trackBuilderChannel_->layerId(ttTrackRef, ttStubRef, layerId))
            continue;
          pattern.set(layerId);
          StreamStub& stubs = streams[channelId * trackBuilderChannel_->maxNumProjectionLayers() + layerId];
          // find dtc stub
          auto found = [ttTrackRef, ttStubRef](const Stub& stub) {
            return stub.ttStubRef_ == ttStubRef && stub.region_ == (int)ttTrackRef->phiSector();
          };
          const Stub& stub = *find_if(stubsDTC.begin(), stubsDTC.end(), found);
          // conv stub word
          stubs.emplace_back(ttStubRef, conv(ttTrackRef, stub));
        }
        // add gaps to layer w/o stub
        for (int layerId : pattern.ids(false))
          streams[channelId * trackBuilderChannel_->maxNumProjectionLayers() + layerId].emplace_back(FrameStub());
      }
    }
    // store products
    iEvent.emplace(edPutTokenAcceptedStubs_, move(streamAcceptedStubs));
    iEvent.emplace(edPutTokenAcceptedTracks_, move(streamAcceptedTracks));
    iEvent.emplace(edPutTokenLostStubs_, move(streamLostStubs));
    iEvent.emplace(edPutTokenLostTracks_, move(streamLostTracks));
  }

  // return h/w bits for given ttTrackRef
  Frame ProducerTBout::conv(const TTTrackRef& ttTrackRef) const {
    static constexpr int widthSeedType = 3;
    static constexpr int widthInvR = 14;
    static constexpr int widthPhi0 = 18;
    static constexpr int widthZ0 = 10;
    static constexpr int widthTanL = 14;
    static const double baseInvR = settings_.kphi1() / settings_.kr() * pow(2, settings_.rinv_shift());
    static const double basePhi0 = settings_.kphi1() * pow(2, settings_.phi0_shift());
    static const double baseZ0 = settings_.kz() * pow(2, settings_.z0_shift());
    static const double baseTanL = settings_.kz() / settings_.kr() * pow(2, settings_.t_shift());
    // sub words
    // phi0 w.r.t. processing region border in rad
    double phi0 = deltaPhi(ttTrackRef->phi() - ttTrackRef->phiSector() * setup_->baseRegion() + setup_->hybridRangePhi() / 2.);
    if (phi0 < 0.)
      phi0 += 2. * M_PI;
    const TTBV hwValid(1, 1);
    const TTBV hwSeedType((int)ttTrackRef->trackSeedType(), widthSeedType);
    const TTBV hwInvR(ttTrackRef->rInv(), baseInvR, widthInvR, true);
    const TTBV hwPhi0(phi0, basePhi0, widthPhi0, false);
    const TTBV hwZ0(ttTrackRef->z0(), baseZ0, widthZ0, true);
    const TTBV hwTanL(ttTrackRef->tanL(), baseTanL, widthTanL, true);
    const TTBV hw(hwValid.str() + hwSeedType.str() + hwInvR.str() + hwZ0.str() + hwTanL.str());
    return hw.bs();
  }

  // return h/w bits for given dtc stub and origin ttTrackRef
  Frame ProducerTBout::conv(const TTTrackRef& ttTrackRef, const Stub& stub) const {
    static constexpr int widthPhi = 12;
    static constexpr int widthZ = 9;
    static constexpr int widthR = 7;
    static const double basePhi = settings_.kphi1();
    static const double baseR = settings_.kr();
    static const double baseZ = settings_.kz();
    static const double baseInvR = settings_.kphi1() / settings_.kr() * pow(2, settings_.rinv_shift());
    static const double basePhi0 = settings_.kphi1() * pow(2, settings_.phi0_shift());
    static const double baseZ0 = settings_.kz() * pow(2, settings_.z0_shift());
    static const double baseTanL = settings_.kz() / settings_.kr() * pow(2, settings_.t_shift());
    const int widthRZ = stub.barrel_ ? widthZ : widthR;
    const double baseRZ = stub.barrel_ ? baseZ : baseR;
    // calc residuals
    const double rInv = (ttTrackRef->rInv() / baseInvR + .5) * baseInvR;
    const double phi0 = (ttTrackRef->phi() / basePhi0 + .5) * basePhi0;
    const double z0 = (ttTrackRef->z0() / baseZ0 + .5) * baseZ0;
    const double tanL = (ttTrackRef->tanL() / baseTanL + .5) * baseTanL;
    const double phi = deltaPhi(phi0 - rInv * stub.gp_.perp() / 2. - stub.gp_.phi());
    const double r = (stub.gp_.z() - z0) / tanL - stub.gp_.perp();
    const double z = z0 + tanL * stub.gp_.perp() - stub.gp_.z();
    const double rz = stub.barrel_ ? r : z;
    // sub words
    const TTBV hwValid(1, 1);
    const TTBV hwR(stub.r_);
    const TTBV hwPhi(phi, basePhi, widthPhi, true);
    const TTBV hwRZ(rz, baseRZ, widthRZ, true);
    const TTBV hw(hwValid.str() + hwR.str() + hwPhi.str() + hwRZ.str());
    return hw.bs();
  }

} // namespace trackFindingTracklet

DEFINE_FWK_MODULE(trackFindingTracklet::ProducerTBout); 