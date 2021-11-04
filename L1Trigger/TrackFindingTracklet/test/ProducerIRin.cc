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

#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <numeric>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackFindingTracklet {

  /*! \class  trackFindingTracklet::ProducerIRin
   *  \brief  Transforms TTTDCinto f/w comparable format for summer chain configuratiotn
   *  \author Thomas Schuh
   *  \date   2021, Oct
   */
  class ProducerIRin : public stream::EDProducer<> {
  public:
    explicit ProducerIRin(const ParameterSet&);
    ~ProducerIRin() override {}

  private:
    virtual void beginRun(const Run&, const EventSetup&) override;
    virtual void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}
    // ED input token of DTC Stubs
    EDGetTokenT<TTDTC> edGetTokenTTDTC_;
    // ED output token for stubs
    EDPutTokenT<StreamsStub> edPutTokenStubs_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_;
    // reduce l1 tracking to summer chain configuration
    bool summerChain_;
    // map of used tfp channels in summer chain config
    vector<int> channelEncoding_;
  };

  ProducerIRin::ProducerIRin(const ParameterSet& iConfig) : iConfig_(iConfig) {
    const InputTag& inputTag = iConfig.getParameter<InputTag>("InputTagDTC");
    const string& branchStubs = iConfig.getParameter<string>("BranchAcceptedStubs");
    // book in- and output ED products
    edGetTokenTTDTC_ = consumes<TTDTC>(inputTag);
    edPutTokenStubs_ = produces<StreamsStub>(branchStubs);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    // initial ES products
    setup_ = nullptr;
  }

  void ProducerIRin::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    if (!setup_->configurationSupported())
      return;
    // check process history if desired
    if (iConfig_.getParameter<bool>("CheckHistory"))
      setup_->checkHistory(iRun.processHistory());
    // reduce l1 tracking to summer chain configuration
    summerChain_ = iConfig_.getParameter<bool>("SummerChain");
    // map of used tfp channels in summer chain config
    channelEncoding_ = iConfig_.getParameter<vector<int>>("SummerChainChannels");
  }

  void ProducerIRin::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty IRin product
    StreamsStub streamStubs;
    // read in hybrid track finding product and produce KFin product
    if (setup_->configurationSupported()) {
      Handle<TTDTC> handleTTDTC;
      iEvent.getByToken<TTDTC>(edGetTokenTTDTC_, handleTTDTC);
      const int numChannel =
          summerChain_ ? channelEncoding_.size() : handleTTDTC->tfpRegions().size() * handleTTDTC->tfpChannels().size();
      streamStubs.reserve(numChannel);
      for (int tfpRegion : handleTTDTC->tfpRegions())
        for (int tfpChannel : summerChain_ ? channelEncoding_ : handleTTDTC->tfpChannels())
          streamStubs.emplace_back(handleTTDTC->stream(tfpRegion, tfpChannel));
    }
    // store products
    iEvent.emplace(edPutTokenStubs_, move(streamStubs));
  }

}  // namespace trackFindingTracklet

DEFINE_FWK_MODULE(trackFindingTracklet::ProducerIRin);