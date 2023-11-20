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
#include "L1Trigger/TrackTrigger/interface/StubPtConsistency.h"

#include <string>
#include <numeric>

using namespace std;
using namespace edm;
using namespace trackerTFP;
using namespace tt;

namespace trklet {

  /*! \class  trklet::ProducerTT
   *  \brief  Converts KF output into TTTracks
   *  \author Thomas Schuh
   *  \date   2021, Aug
   */
  class ProducerTT : public stream::EDProducer<> {
  public:
    explicit ProducerTT(const ParameterSet&);
    ~ProducerTT() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override;
    void produce(Event&, const EventSetup&) override;
    void endJob() {}

    // ED input token of kf stubs
    EDGetTokenT<StreamsStub> edGetTokenStubs_;
    // ED input token of kf tracks
    EDGetTokenT<StreamsTrack> edGetTokenTracks_;
    // ED output token for TTTracks
    EDPutTokenT<TTTracks> edPutToken_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // DataFormats token
    ESGetToken<DataFormats, DataFormatsRcd> esGetTokenDataFormats_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_ = nullptr;
    // helper class to extract structured data from tt::Frames
    const DataFormats* dataFormats_ = nullptr;
  };

  ProducerTT::ProducerTT(const ParameterSet& iConfig) : iConfig_(iConfig) {
    const string& label = iConfig.getParameter<string>("LabelKF");
    const string& branchStubs = iConfig.getParameter<string>("BranchStubsAccepted");
    const string& branchTracks = iConfig.getParameter<string>("BranchTracksAccepted");
    // book in- and output ED products
    edGetTokenStubs_ = consumes<StreamsStub>(InputTag(label, branchStubs));
    edGetTokenTracks_ = consumes<StreamsTrack>(InputTag(label, branchTracks));
    edPutToken_ = produces<TTTracks>(branchTracks);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
  }

  void ProducerTT::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    // helper class to extract structured data from tt::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
  }

  void ProducerTT::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty KFout product
    TTTracks ttTracks;
    // read in KF Product and produce KFout product
    Handle<StreamsStub> handleStubs;
    iEvent.getByToken<StreamsStub>(edGetTokenStubs_, handleStubs);
    const StreamsStub& streamsStubs = *handleStubs.product();
    Handle<StreamsTrack> handleTracks;
    iEvent.getByToken<StreamsTrack>(edGetTokenTracks_, handleTracks);
    const StreamsTrack& streamsTracks = *handleTracks.product();
    // count number of kf tracks
    int nTracks(0);
    for (const StreamTrack& stream : streamsTracks)
      nTracks += accumulate(stream.begin(), stream.end(), 0, [](int sum, const FrameTrack& frame) {
        return sum + (frame.first.isNonnull() ? 1 : 0);
      });
    ttTracks.reserve(nTracks);
    auto put = [&ttTracks, this](const TrackKF& track, const vector<StubKF*>& stubs, int region) {
      const double zT = dataFormats_->format(Variable::zT, Process::gp).digi(track.zT());
      const double cot = zT / setup_->chosenRofZ() + track.cot();
      const double z0 = track.zT() - setup_->chosenRofZ() * cot;
      const double inv2R = track.inv2R();
      const double phi0 = deltaPhi(track.phiT() - setup_->chosenRofPhi() * inv2R +
                                   region * dataFormats_->format(Variable::phiT, Process::kf).range());
      const double invR = -2. * track.inv2R();
      TTBV hitVector(0, setup_->numLayers());
      for (int layer = 0; layer < (int)stubs.size(); layer++)
        if (stubs[layer])
          hitVector.set(layer);
      double chi2phi(0.);
      double chi2z(0.);
      vector<TTStubRef> ttStubRefs;
      ttStubRefs.reserve(hitVector.count());
      for (int layer : hitVector.ids()) {
        StubKF* stub = stubs[layer];
        chi2phi += pow(stub->phi() / stub->dPhi(), 2);
        chi2z += pow(stub->z() / stub->dZ(), 2);
        ttStubRefs.push_back(stub->frame().first);
      }
      static constexpr int nPar = 4;
      static constexpr double d0 = 0.;
      static constexpr double trkMVA1 = 0.;
      static constexpr double trkMVA2 = 0.;
      static constexpr double trkMVA3 = 0.;
      const int hitPattern = hitVector.val();
      const double bField = setup_->bField();
      TTTrack<Ref_Phase2TrackerDigi_> ttTrack(
          invR, phi0, cot, z0, d0, chi2phi, chi2z, trkMVA1, trkMVA2, trkMVA3, hitPattern, nPar, bField);
      ttTrack.setStubRefs(ttStubRefs);
      ttTrack.setPhiSector(track.frame().first->phiSector());
      ttTrack.setEtaSector(track.frame().first->etaSector());
      ttTrack.setTrackSeedType(track.frame().first->trackSeedType());
      ttTrack.setStubPtConsistency(StubPtConsistency::getConsistency(
          ttTrack, setup_->trackerGeometry(), setup_->trackerTopology(), bField, nPar));
      return ttTrack;
    };
    // convert kf track frames per channel and stub frames per channel and layer to TTTracks
    for (int channel = 0; channel < (int)streamsTracks.size(); channel++) {
      const int region = channel / dataFormats_->numChannel(Process::kf);
      const int offset = channel * setup_->numLayers();
      int iTrk(0);
      for (const FrameTrack& frameTrack : streamsTracks[channel]) {
        if (frameTrack.first.isNull())
          continue;
        // convert stub frames to kf stubs
        vector<StubKF> stubs;
        vector<StubKF*> stubsPtr;
        stubs.reserve(setup_->numLayers());
        for (int layer = 0; layer < setup_->numLayers(); layer++) {
          const FrameStub& frameStub = streamsStubs[offset + layer][iTrk];
          StubKF* stub = nullptr;
          if (frameStub.first.isNonnull()) {
            stubs.emplace_back(frameStub, dataFormats_);
            stub = &stubs.back();
          }
          stubsPtr.push_back(stub);
        }
        // convert track frame to kf track
        TrackKF track(frameTrack, dataFormats_);
        // convert kf track and kf stubs to TTTrack
        put(track, stubsPtr, region);
        iTrk++;
      }
    }
    // store products
    iEvent.emplace(edPutToken_, std::move(ttTracks));
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerTT);
