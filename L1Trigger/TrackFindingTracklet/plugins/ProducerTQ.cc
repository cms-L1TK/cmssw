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

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackFindingTracklet/interface/ChannelAssignment.h"
#include "L1Trigger/TrackFindingTracklet/interface/DataFormats.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackQuality.h"

#include <string>

namespace trklet {

  /*! \class  trklet::ProducerTQ
   *  \brief  Bit accurate emulation of the track quality BDT
   *  \author Thomas Schuh
   *  \date   20245, Aug
   */
  class ProducerTQ : public edm::stream::EDProducer<> {
  public:
    explicit ProducerTQ(const edm::ParameterSet&);
    ~ProducerTQ() override = default;
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    // ED input token of kf stubs
    edm::EDGetTokenT<tt::StreamsStub> edGetTokenStubs_;
    // ED input token of kf tracks
    edm::EDGetTokenT<tt::StreamsTrack> edGetTokenTracks_;
    // ED output token for additional track variables created by TQ
    edm::EDPutTokenT<tt::StreamsTrack> edPutToken_;
    // Setup token
    edm::ESGetToken<tt::Setup, tt::SetupRcd> esGetTokenSetup_;
    // ChannelAssignment token
    edm::ESGetToken<ChannelAssignment, ChannelAssignmentRcd> esGetTokenChannelAssignment_;
    // DataFormats token
    edm::ESGetToken<DataFormats, ChannelAssignmentRcd> esGetTokenDataFormats_;
    // helper class to extract structured data from tt::Frames
    const DataFormats* dataFormats_ = nullptr;
    // Internal data formats
    TrackQuality::InternalFormats internalFormats_;
  };

  ProducerTQ::ProducerTQ(const edm::ParameterSet& iConfig) {
    const std::string& label = iConfig.getParameter<std::string>("InputLabelTQ");
    const std::string& branchStubs = iConfig.getParameter<std::string>("BranchStubs");
    const std::string& branchTracks = iConfig.getParameter<std::string>("BranchTracks");
    // book in- and output ED products
    edGetTokenStubs_ = consumes(edm::InputTag(label, branchStubs));
    edGetTokenTracks_ = consumes(edm::InputTag(label, branchTracks));
    edPutToken_ = produces(branchTracks);
    // book ES products
    esGetTokenSetup_ = esConsumes();
    esGetTokenChannelAssignment_ = esConsumes<edm::Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<edm::Transition::BeginRun>();
  }

  void ProducerTQ::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    // helper class to extract structured data from tt::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    const ChannelAssignment* channelAssignment = &iSetup.getData(esGetTokenChannelAssignment_);
    int width;
    int shift;
    double base;
    double range;
    // m02
    const DataFormat& phi = dataFormats_->format(Variable::phi, Process::kf);
    width = 2 * phi.width();
    base = std::pow(phi.base(), 2);
    range = std::pow(phi.range(), 2) / 4.;
    internalFormats_.m02_ = DataFormat(false, width, base, range);
    // m12
    const DataFormat& z = dataFormats_->format(Variable::z, Process::kf);
    width = 2 * z.width();
    base = std::pow(z.base(), 2);
    range = std::pow(z.range(), 2) / 4.;
    internalFormats_.m12_ = DataFormat(false, width, base, range);
    // invV0
    const DataFormat& dPhi = dataFormats_->format(Variable::dPhi, Process::kf);
    width = channelAssignment->tqWidthInvV0();
    base = std::pow(dPhi.base(), -2);
    range = base * std::pow(2, width) / (std::pow(2, width) - 1);
    shift = std::ceil(std::log2(range / base)) - width;
    base *= std::pow(2, shift);
    internalFormats_.invV0_ = DataFormat(false, width, base, range);
    // invV1
    const DataFormat& dZ = dataFormats_->format(Variable::dZ, Process::kf);
    width = channelAssignment->tqWidthInvV1();
    base = std::pow(dZ.base(), -2);
    range = base * std::pow(2, width) / (std::pow(2, width) - 1);
    shift = std::ceil(std::log2(range / base)) - width;
    base *= std::pow(2, shift);
    internalFormats_.invV1_ = DataFormat(false, width, base, range);
    // chi20BDT
    shift = channelAssignment->tqBaseShiftBDTchi20();
    width = channelAssignment->tqWidthBDTchi20();
    base = std::pow(2., shift);
    range = base * std::pow(2, width);
    internalFormats_.chi20BDT_ = DataFormat(false, width, base, range);
    // chi21BDT
    shift = channelAssignment->tqBaseShiftBDTchi21();
    width = channelAssignment->tqWidthBDTchi21();
    base = std::pow(2., shift);
    range = base * std::pow(2, width);
    internalFormats_.chi21BDT_ = DataFormat(false, width, base, range);
  }

  void ProducerTQ::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // helper class to store configurations
    const tt::Setup* setup = &iSetup.getData(esGetTokenSetup_);
    // read in KF Product
    const tt::StreamsTrack& streamsTracks = iEvent.get(edGetTokenTracks_);
    const tt::StreamsStub& streamsStubs = iEvent.get(edGetTokenStubs_);
    // empty TQ product
    tt::StreamsTrack output(setup->numRegions());
    //produce TQ product
    for (int region = 0; region < setup->numRegions(); region++) {
      // object emulating tq algorithm
      TrackQuality tq(setup, dataFormats_, internalFormats_, region);
      // read in and organize input tracks and stubs
      tq.consume(streamsTracks, streamsStubs);
      // fills output products
      tq.produce(output[region]);
    }
    // store TQ product
    iEvent.emplace(edPutToken_, std::move(output));
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerTQ);
