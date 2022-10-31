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

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "L1Trigger/TrackerTFP/interface/KFin.h"

#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  /*! \class  trackerTFP::ProducerKFin
   *  \brief  transforms TTTracks into KF input
   *  \author Thomas Schuh
   *  \date   2020, July
   */
  class ProducerKFin : public stream::EDProducer<> {
  public:
    explicit ProducerKFin(const ParameterSet&);
    ~ProducerKFin() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override;
    void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}

    // ED input token of TTTracks
    EDGetTokenT<vector<TTTrack<Ref_Phase2TrackerDigi_>>> edGetTokenTTTracks_;
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
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_ = nullptr;
    // helper class to extract structured data from tt::Frames
    const DataFormats* dataFormats_ = nullptr;
    // helper class to encode layer
    const LayerEncoding* layerEncoding_ = nullptr;
    //
    bool enableTruncation_;
  };

  ProducerKFin::ProducerKFin(const ParameterSet& iConfig) : iConfig_(iConfig) {
    const string& labelTTTracks = iConfig.getParameter<string>("LabelTTFound");
    const string& labelStubs = iConfig.getParameter<string>("LabelZHT");
    const string& branchAcceptedStubs = iConfig.getParameter<string>("BranchAcceptedStubs");
    const string& branchAcceptedTracks = iConfig.getParameter<string>("BranchAcceptedTracks");
    const string& branchLostStubs = iConfig.getParameter<string>("BranchLostStubs");
    const string& branchLostTracks = iConfig.getParameter<string>("BranchLostTracks");
    // book in- and output ED products
    edGetTokenTTTracks_ =
        consumes<vector<TTTrack<Ref_Phase2TrackerDigi_>>>(InputTag(labelTTTracks, branchAcceptedTracks));
    edGetTokenStubs_ = consumes<StreamsStub>(InputTag(labelStubs, branchAcceptedStubs));
    edPutTokenAcceptedStubs_ = produces<StreamsStub>(branchAcceptedStubs);
    edPutTokenAcceptedTracks_ = produces<StreamsTrack>(branchAcceptedTracks);
    edPutTokenLostStubs_ = produces<StreamsStub>(branchLostStubs);
    edPutTokenLostTracks_ = produces<StreamsTrack>(branchLostTracks);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    esGetTokenLayerEncoding_ = esConsumes<LayerEncoding, LayerEncodingRcd, Transition::BeginRun>();
    //
    enableTruncation_ = iConfig.getParameter<bool>("EnableTruncation");
  }

  void ProducerKFin::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    if (!setup_->configurationSupported())
      return;
    // check process history if desired
    if (iConfig_.getParameter<bool>("CheckHistory"))
      setup_->checkHistory(iRun.processHistory());
    // helper class to extract structured data from tt::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    // helper class to encode layer
    layerEncoding_ = &iSetup.getData(esGetTokenLayerEncoding_);
  }

  void ProducerKFin::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty KFin products
    StreamsStub acceptedStubs(dataFormats_->numStreamsStubs(Process::kfin));
    StreamsTrack acceptedTracks(dataFormats_->numStreamsTracks(Process::kfin));
    StreamsStub lostStubs(dataFormats_->numStreamsStubs(Process::kfin));
    StreamsTrack lostTracks(dataFormats_->numStreamsTracks(Process::kfin));
    // read in SFout Product and produce KFin product
    if (setup_->configurationSupported()) {
      Handle<StreamsStub> handleStubs;
      iEvent.getByToken<StreamsStub>(edGetTokenStubs_, handleStubs);
      const StreamsStub& streams = *handleStubs.product();
      Handle<vector<TTTrack<Ref_Phase2TrackerDigi_>>> handleTTTracks;
      iEvent.getByToken<vector<TTTrack<Ref_Phase2TrackerDigi_>>>(edGetTokenTTTracks_, handleTTTracks);
      vector<TTTrackRef> ttTrackRefs;
      const int nTTTracks = handleTTTracks->size();
      ttTrackRefs.reserve(nTTTracks);
      for (int i = nTTTracks - 1; i > -1; i--)
        ttTrackRefs.emplace_back(handleTTTracks, i);
      for (int region = 0; region < setup_->numRegions(); region++) {
        // object to connect track candidates with track fit
        KFin kfin(iConfig_, setup_, dataFormats_, layerEncoding_, region);
        // read in and organize input product
        kfin.consume(streams);
        // read in TTTrack Refs
        kfin.consume(ttTrackRefs);
        // fill output products
        kfin.produce(acceptedTracks, acceptedStubs, lostTracks, lostStubs);
      }
    }
    // store products
    iEvent.emplace(edPutTokenAcceptedStubs_, move(acceptedStubs));
    iEvent.emplace(edPutTokenAcceptedTracks_, move(acceptedTracks));
    iEvent.emplace(edPutTokenLostStubs_, move(lostStubs));
    iEvent.emplace(edPutTokenLostTracks_, move(lostTracks));
  }

}  // namespace trackerTFP

DEFINE_FWK_MODULE(trackerTFP::ProducerKFin);
