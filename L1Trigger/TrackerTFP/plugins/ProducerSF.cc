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
#include "L1Trigger/TrackerTFP/interface/SeedFilter.h"

#include <string>
#include <numeric>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  /*! \class  trackerTFP::ProducerSF
   *  \brief  L1TrackTrigger Seed Filter emulator
   *  \author Thomas Schuh
   *  \date   2020, June
   */
  class ProducerSF : public stream::EDProducer<> {
  public:
    explicit ProducerSF(const ParameterSet&);
    ~ProducerSF() override {}

  private:
    virtual void beginRun(const Run&, const EventSetup&) override;
    virtual void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}

    // ED input token of mht stubs
    EDGetTokenT<StreamsStub> edGetToken_;
    // ED output token for accepted stubs
    EDPutTokenT<StreamsStub> edPutTokenAccepted_;
    // ED output token for lost stubs
    EDPutTokenT<StreamsStub> edPutTokenLost_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // DataFormats token
    ESGetToken<DataFormats, DataFormatsRcd> esGetTokenDataFormats_;
    // LayerEncoding token
    ESGetToken<LayerEncoding, LayerEncodingRcd> esGetTokenLayerEncoding_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_;
    // helper class to extract structured data from tt::Frames
    const DataFormats* dataFormats_;
    //
    const LayerEncoding* layerEncoding_;
  };

  ProducerSF::ProducerSF(const ParameterSet& iConfig) :
    iConfig_(iConfig)
  {
    const string& label = iConfig.getParameter<string>("LabelMHT");
    const string& branchAccepted = iConfig.getParameter<string>("BranchAcceptedStubs");
    const string& branchLost = iConfig.getParameter<string>("BranchLostStubs");
    // book in- and output ED products
    edGetToken_ = consumes<StreamsStub>(InputTag(label, branchAccepted));
    edPutTokenAccepted_ = produces<StreamsStub>(branchAccepted);
    edPutTokenLost_ = produces<StreamsStub>(branchLost);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    esGetTokenLayerEncoding_ = esConsumes<LayerEncoding, LayerEncodingRcd, Transition::BeginRun>();
    // initial ES products
    setup_ = nullptr;
    dataFormats_ = nullptr;
    layerEncoding_ = nullptr;
  }

  void ProducerSF::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    if (!setup_->configurationSupported())
      return;
    // check process history if desired
    if (iConfig_.getParameter<bool>("CheckHistory"))
      setup_->checkHistory(iRun.processHistory());
    // helper class to extract structured data from tt::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    //
    layerEncoding_ = &iSetup.getData(esGetTokenLayerEncoding_);
  }

  void ProducerSF::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty SF products
    StreamsStub accepted(dataFormats_->numStreams(Process::mht));
    StreamsStub lost(dataFormats_->numStreams(Process::mht));
    // read in MHT Product and produce SF product
    if (setup_->configurationSupported()) {
      Handle<StreamsStub> handle;
      iEvent.getByToken<StreamsStub>(edGetToken_, handle);
      const StreamsStub& streams = *handle.product();
      for (int region = 0; region < setup_->numRegions(); region++) {
        // object to find in a region finer rough candidates in r-z
        SeedFilter sf(iConfig_, setup_, dataFormats_, layerEncoding_, region);
        // read in and organize input product
        sf.consume(streams);
        // fill output products
        sf.produce(accepted, lost);
      }
    }
    // store products
    iEvent.emplace(edPutTokenAccepted_, move(accepted));
    iEvent.emplace(edPutTokenLost_, move(lost));
  }

} // namespace trackerTFP

DEFINE_FWK_MODULE(trackerTFP::ProducerSF);