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
#include "L1Trigger/TrackFindingTracklet/interface/Setup.h"
#include "L1Trigger/TrackFindingTracklet/interface/DataFormats.h"
#include "L1Trigger/TrackFindingTracklet/interface/KalmanFilter.h"

#include <string>
#include <vector>

namespace trklet {

  /*! \class  trklet::ProducerKF
   *  \brief  L1TrackTrigger Kamlan Filter emulator
   *  \author Thomas Schuh
   *  \date   2020, July
   */
  class ProducerKF : public edm::stream::EDProducer<> {
  public:
    explicit ProducerKF(const edm::ParameterSet&);
    ~ProducerKF() override = default;

  private:
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override {
      std::stringstream ss;
      if (printDebug_)
        kalmanFilterFormats_.endJob(ss);
      edm::LogPrint(moduleDescription().moduleName()) << ss.str();
    }
    // merge woker output
    void merge(
        const tt::StreamsStub&, const tt::StreamsTrack&, tt::StreamsStub&, tt::StreamsTrack&, int, const Setup*) const;
    // ED input token of sf stubs and tracks
    edm::EDGetTokenT<tt::StreamsStub> edGetTokenStubs_;
    edm::EDGetTokenT<tt::StreamsTrack> edGetTokenTracks_;
    // ED output token for accepted stubs and tracks
    edm::EDPutTokenT<tt::TTTracks> edPutTokenTTTracks_;
    edm::EDPutTokenT<tt::StreamsStub> edPutTokenStubs_;
    edm::EDPutTokenT<tt::StreamsTrack> edPutTokenTracks_;
    // Setup token
    edm::ESGetToken<Setup, trackerDTC::SetupRcd> esGetTokenSetup_;
    // DataFormats token
    edm::ESGetToken<DataFormats, trackerDTC::SetupRcd> esGetTokenDataFormats_;
    // helper class to extract structured data from tt::Frames
    const DataFormats* dataFormats_;
    // provides dataformats of Kalman filter internals
    KalmanFilterFormats kalmanFilterFormats_;
    //
    ConfigKF iConfig_;
    // print end job internal unused MSB
    bool printDebug_;
  };

  ProducerKF::ProducerKF(const edm::ParameterSet& iConfig) {
    iConfig_.enableIntegerEmulation_ = iConfig.getParameter<bool>("EnableIntegerEmulation");
    iConfig_.widthR00_ = iConfig.getParameter<int>("WidthR00");
    iConfig_.widthR11_ = iConfig.getParameter<int>("WidthR11");
    iConfig_.widthC00_ = iConfig.getParameter<int>("WidthC00");
    iConfig_.widthC01_ = iConfig.getParameter<int>("WidthC01");
    iConfig_.widthC11_ = iConfig.getParameter<int>("WidthC11");
    iConfig_.widthC22_ = iConfig.getParameter<int>("WidthC22");
    iConfig_.widthC23_ = iConfig.getParameter<int>("WidthC23");
    iConfig_.widthC33_ = iConfig.getParameter<int>("WidthC33");
    iConfig_.baseShiftx0_ = iConfig.getParameter<int>("BaseShiftx0");
    iConfig_.baseShiftx1_ = iConfig.getParameter<int>("BaseShiftx1");
    iConfig_.baseShiftx2_ = iConfig.getParameter<int>("BaseShiftx2");
    iConfig_.baseShiftx3_ = iConfig.getParameter<int>("BaseShiftx3");
    iConfig_.baseShiftr0_ = iConfig.getParameter<int>("BaseShiftr0");
    iConfig_.baseShiftr1_ = iConfig.getParameter<int>("BaseShiftr1");
    iConfig_.baseShiftS00_ = iConfig.getParameter<int>("BaseShiftS00");
    iConfig_.baseShiftS01_ = iConfig.getParameter<int>("BaseShiftS01");
    iConfig_.baseShiftS12_ = iConfig.getParameter<int>("BaseShiftS12");
    iConfig_.baseShiftS13_ = iConfig.getParameter<int>("BaseShiftS13");
    iConfig_.baseShiftR00_ = iConfig.getParameter<int>("BaseShiftR00");
    iConfig_.baseShiftR11_ = iConfig.getParameter<int>("BaseShiftR11");
    iConfig_.baseShiftInvR00Approx_ = iConfig.getParameter<int>("BaseShiftInvR00Approx");
    iConfig_.baseShiftInvR11Approx_ = iConfig.getParameter<int>("BaseShiftInvR11Approx");
    iConfig_.baseShiftInvR00Cor_ = iConfig.getParameter<int>("BaseShiftInvR00Cor");
    iConfig_.baseShiftInvR11Cor_ = iConfig.getParameter<int>("BaseShiftInvR11Cor");
    iConfig_.baseShiftInvR00_ = iConfig.getParameter<int>("BaseShiftInvR00");
    iConfig_.baseShiftInvR11_ = iConfig.getParameter<int>("BaseShiftInvR11");
    iConfig_.baseShiftS00Shifted_ = iConfig.getParameter<int>("BaseShiftS00Shifted");
    iConfig_.baseShiftS01Shifted_ = iConfig.getParameter<int>("BaseShiftS01Shifted");
    iConfig_.baseShiftS12Shifted_ = iConfig.getParameter<int>("BaseShiftS12Shifted");
    iConfig_.baseShiftS13Shifted_ = iConfig.getParameter<int>("BaseShiftS13Shifted");
    iConfig_.baseShiftK00_ = iConfig.getParameter<int>("BaseShiftK00");
    iConfig_.baseShiftK10_ = iConfig.getParameter<int>("BaseShiftK10");
    iConfig_.baseShiftK21_ = iConfig.getParameter<int>("BaseShiftK21");
    iConfig_.baseShiftK31_ = iConfig.getParameter<int>("BaseShiftK31");
    iConfig_.baseShiftC00_ = iConfig.getParameter<int>("BaseShiftC00");
    iConfig_.baseShiftC01_ = iConfig.getParameter<int>("BaseShiftC01");
    iConfig_.baseShiftC11_ = iConfig.getParameter<int>("BaseShiftC11");
    iConfig_.baseShiftC22_ = iConfig.getParameter<int>("BaseShiftC22");
    iConfig_.baseShiftC23_ = iConfig.getParameter<int>("BaseShiftC23");
    iConfig_.baseShiftC33_ = iConfig.getParameter<int>("BaseShiftC33");
    iConfig_.baseShiftInvDH_ = iConfig.getParameter<int>("BaseShiftInvDH");
    iConfig_.baseShiftInvDH2_ = iConfig.getParameter<int>("BaseShiftInvDH2");
    iConfig_.baseShiftHv0_ = iConfig.getParameter<int>("BaseShiftHv0");
    iConfig_.baseShiftHv1_ = iConfig.getParameter<int>("BaseShiftHv1");
    iConfig_.baseShiftH2v0_ = iConfig.getParameter<int>("BaseShiftH2v0");
    iConfig_.baseShiftH2v1_ = iConfig.getParameter<int>("BaseShiftH2v1");
    printDebug_ = iConfig.getParameter<bool>("PrintKFDebug");
    const std::string& label = iConfig.getParameter<std::string>("InputLabelKF");
    const std::string& branchStubs = iConfig.getParameter<std::string>("BranchStubs");
    const std::string& branchTracks = iConfig.getParameter<std::string>("BranchTracks");
    // book in- and output ED products
    edGetTokenStubs_ = consumes(edm::InputTag(label, branchStubs));
    edGetTokenTracks_ = consumes(edm::InputTag(label, branchTracks));
    edPutTokenStubs_ = produces(branchStubs);
    edPutTokenTracks_ = produces(branchTracks);
    edPutTokenTTTracks_ = produces(branchTracks);
    // book ES products
    esGetTokenSetup_ = esConsumes();
    esGetTokenDataFormats_ = esConsumes<edm::Transition::BeginRun>();
  }

  void ProducerKF::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    // helper class to extract structured data from tt::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    kalmanFilterFormats_.consume(dataFormats_, iConfig_);
  }

  void ProducerKF::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // helper class to store configurations
    const Setup* setup = &iSetup.getData(esGetTokenSetup_);
    // empty KF products
    tt::StreamsStub streamsStub(setup->sysNumRegion() * setup->kfNumLayers());
    tt::StreamsTrack streamsTrack(setup->sysNumRegion());
    // read in DR Product and produce KF product
    const tt::StreamsStub& stubs = iEvent.get(edGetTokenStubs_);
    const tt::StreamsTrack& tracks = iEvent.get(edGetTokenTracks_);
    // prep TTTracks
    tt::TTTracks ttTracks;
    std::vector<TTTrackRef> ttTrackRefs;
    for (int region = 0; region < setup->sysNumRegion(); region++) {
      // object to fit tracks in a processing region
      KalmanFilter kf(setup, dataFormats_, &kalmanFilterFormats_, region, ttTracks);
      // read in and organize input tracks and stubs
      kf.consume(tracks, stubs);
      // fill output products
      kf.produce(streamsStub, streamsTrack);
    }
    // store products
    iEvent.emplace(edPutTokenStubs_, std::move(streamsStub));
    iEvent.emplace(edPutTokenTracks_, std::move(streamsTrack));
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerKF);
