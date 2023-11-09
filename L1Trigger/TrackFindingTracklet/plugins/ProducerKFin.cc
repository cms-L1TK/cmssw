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
#include "L1Trigger/TrackFindingTracklet/interface/ChannelAssignment.h"
#include "L1Trigger/TrackFindingTracklet/interface/KFin.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTypes.h"

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

namespace trklet {

  /*! \class  trklet::ProducerKFin
   *  \brief  Transforms format of DR into that expected by KF input.
   *  \author Thomas Schuh
   *  \date   2023, Feb
   */
  class ProducerKFin : public stream::EDProducer<> {
  public:
    explicit ProducerKFin(const ParameterSet&);
    ~ProducerKFin() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override;
    void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}
    // ED input token of Tracks
    EDGetTokenT<StreamsTrack> edGetTokenTracks_;
    // ED input token of Stubs
    EDGetTokenT<StreamsStub> edGetTokenStubs_;
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
    // LayerEncoding token
    ESGetToken<LayerEncoding, LayerEncodingRcd> esGetTokenLayerEncoding_;
    // ChannelAssignment token
    ESGetToken<ChannelAssignment, ChannelAssignmentRcd> esGetTokenChannelAssignment_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_ = nullptr;
    // helper class to extract structured data from tt::Frames
    const DataFormats* dataFormats_ = nullptr;
    // helper class to encode layer
    const LayerEncoding* layerEncoding_ = nullptr;
    // helper class to assign tracks to channel
    const ChannelAssignment* channelAssignment_ = nullptr;
  };

  ProducerKFin::ProducerKFin(const ParameterSet& iConfig) : iConfig_(iConfig) {
    const string& label = iConfig.getParameter<string>("LabelTBout");
    const string& branchAcceptedStubs = iConfig.getParameter<string>("BranchStubsAccepted");
    const string& branchAcceptedTracks = iConfig.getParameter<string>("BranchTracksAccepted");
    const string& branchLostStubs = iConfig.getParameter<string>("BranchStubsTruncated");
    const string& branchLostTracks = iConfig.getParameter<string>("BranchTracksTruncated");
    // book in- and output ED products
    edGetTokenTracks_ = consumes<StreamsTrack>(InputTag(label, branchAcceptedTracks));
    edGetTokenStubs_ = consumes<StreamsStub>(InputTag(label, branchAcceptedStubs));
    edPutTokenAcceptedStubs_ = produces<StreamsStub>(branchAcceptedStubs);
    edPutTokenAcceptedTracks_ = produces<StreamsTrack>(branchAcceptedTracks);
    edPutTokenLostStubs_ = produces<StreamsStub>(branchLostStubs);
    edPutTokenLostTracks_ = produces<StreamsTrack>(branchLostTracks);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    esGetTokenLayerEncoding_ = esConsumes<LayerEncoding, LayerEncodingRcd, Transition::BeginRun>();
    esGetTokenChannelAssignment_ = esConsumes<ChannelAssignment, ChannelAssignmentRcd, Transition::BeginRun>();
  }

  void ProducerKFin::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    // helper class to extract structured data from tt::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    // helper class to encode layer
    layerEncoding_ = &iSetup.getData(esGetTokenLayerEncoding_);
    // helper class to assign tracks to channel
    channelAssignment_ = &iSetup.getData(esGetTokenChannelAssignment_);
  }

  void ProducerKFin::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty KFin products
    StreamsStub acceptedStubs(dataFormats_->numStreamsStubs(Process::ctb));
    StreamsTrack acceptedTracks(dataFormats_->numStreamsTracks(Process::ctb));
    StreamsStub lostStubs(dataFormats_->numStreamsStubs(Process::ctb));
    StreamsTrack lostTracks(dataFormats_->numStreamsTracks(Process::ctb));
    // read in TBout Product and produce KFin product
    Handle<StreamsStub> handleStubs;
    iEvent.getByToken<StreamsStub>(edGetTokenStubs_, handleStubs);
    const StreamsStub& stubs = *handleStubs;
    Handle<StreamsTrack> handleTracks;
    iEvent.getByToken<StreamsTrack>(edGetTokenTracks_, handleTracks);
    const StreamsTrack& tracks = *handleTracks;
    for (int region = 0; region < setup_->numRegions(); region++) {
      // object to reformat tracks from tracklet fromat to TMTT format in a processing region
      KFin kfin(iConfig_, setup_, dataFormats_, layerEncoding_, channelAssignment_, &settings_, region);
      // read in and organize input tracks and stubs
      kfin.consume(tracks, stubs);
      // fill output products
      kfin.produce(acceptedStubs, acceptedTracks, lostStubs, lostTracks);
    }
    // store products
    iEvent.emplace(edPutTokenAcceptedStubs_, std::move(acceptedStubs));
    iEvent.emplace(edPutTokenAcceptedTracks_, std::move(acceptedTracks));
    iEvent.emplace(edPutTokenLostStubs_, std::move(lostStubs));
    iEvent.emplace(edPutTokenLostTracks_, std::move(lostTracks));
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerKFin);
