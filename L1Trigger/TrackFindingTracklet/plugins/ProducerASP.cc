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

#include "L1Trigger/TrackFindingTracklet/interface/Setup.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubProcessor.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <numeric>

namespace trklet {

  /*! \class  trklet::ProducerASP
   *  \brief  Processes AllStubs
   *  \author Thomas Schuh
   *  \date   2026, Nov
   */
  class ProducerASP : public edm::stream::EDProducer<> {
  public:
    explicit ProducerASP(const edm::ParameterSet&);
    ~ProducerASP() override = default;

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;

    // ED input token of Stubs
    edm::EDGetTokenT<tt::StreamsStub> edGetToken_;
    // ED output token for stubs
    edm::EDPutTokenT<tt::StreamsStub> edPutToken_;
    // Setup token
    edm::ESGetToken<Setup, trackerDTC::SetupRcd> esGetTokenSetup_;
    // DataFormats token
    edm::ESGetToken<DataFormats, trackerDTC::SetupRcd> esGetTokenDataFormats_;
  };

  ProducerASP::ProducerASP(const edm::ParameterSet& iConfig) {
    const std::string& label = iConfig.getParameter<std::string>("InputLabelASP");
    const std::string& branch = iConfig.getParameter<std::string>("BranchStubs");
    // book in- and output ED products
    edGetToken_ = consumes(edm::InputTag(label, branch));
    edPutToken_ = produces(branch);
    // book ES products
    esGetTokenSetup_ = esConsumes();
    esGetTokenDataFormats_ = esConsumes();
  }

  void ProducerASP::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // helper class to store configurations
    const Setup* setup = &iSetup.getData(esGetTokenSetup_);
    const DataFormats* df = &iSetup.getData(esGetTokenDataFormats_);
    tt::StreamsStub output(setup->sysNumRegion() * setup->aspMux().size());
    const tt::StreamsStub& input = iEvent.get(edGetToken_);
    for (int region = 0; region < setup->sysNumRegion(); region++) {
      // object to reformat tracks from tracklet fromat to TMTT format in a processing region
      AllStubProcessor asp(setup, df, region);
      // read in and organize input stubs
      asp.consume(input);
      // fill output products
      asp.produce(output);
    }
    // store products
    iEvent.emplace(edPutToken_, std::move(output));
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerASP);
