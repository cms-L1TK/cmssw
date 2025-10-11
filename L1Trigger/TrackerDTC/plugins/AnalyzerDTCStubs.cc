#define EDM_ML_DEBUG

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

// ROOT file service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT histos
#include <TH1F.h>
#include <TH2D.h>

// TTStubDetSetVec + related typedefs
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
// DetId
#include "DataFormats/DetId/interface/DetId.h"
// tt::Setup & SensorModule
#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackTrigger/interface/SensorModule.h"

#include <vector>

namespace trackerDTC 
{

  class AnalyzerDTCStubs : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> 
  {
    public:
    
        explicit AnalyzerDTCStubs(const edm::ParameterSet& iConfig);
        ~AnalyzerDTCStubs() override = default;

        void beginJob() override;
        void beginRun(const edm::Run&, const edm::EventSetup&) override;
        void analyze(const edm::Event&, const edm::EventSetup&) override;
        void endRun(const edm::Run&, const edm::EventSetup&) override;
        void endJob() override;

    private:

        const edm::InputTag inputTag_;
        edm::EDGetTokenT<TTStubDetSetVec> edGetTokenTTStubs_;

        // Setup ES tokens
        edm::ESGetToken<tt::Setup, tt::SetupRcd> esGetTokenSetup_;             // Event (used in analyze)
        edm::ESGetToken<tt::Setup, tt::SetupRcd> esGetTokenSetupBeginRun_;     // BeginRun (used to book TH2D / throughput bins)

        // Histogram: total number of stubs across ALL DTCs per event
        TH1F* hisAllDTCStubs_{nullptr};

        // TH2D: x = DTC id, y = number of stubs (per event)
        TH2D* h2DTCvsStubs_{nullptr};

        // NEW: Throughput in Gbps per event (sum over all DTCs)
        TH1F* hisThroughputGbps_{nullptr};

        static constexpr int kMaxOccY_ = 500;
        static constexpr int kNBinsY_  = 100;

        static constexpr double kEventRateHz_ = 750e3;  // 750 kHz
        static constexpr int    kBitsPerStub_ = 64;     // 64 bits per stub
  };

  AnalyzerDTCStubs::AnalyzerDTCStubs(const edm::ParameterSet& iConfig)
      : inputTag_(iConfig.getParameter<edm::InputTag>("InputTag")) 
  {
    usesResource("TFileService");
    edGetTokenTTStubs_       = consumes<TTStubDetSetVec>(inputTag_);
    esGetTokenSetup_         = esConsumes();
    esGetTokenSetupBeginRun_ = esConsumes<edm::Transition::BeginRun>();

    LogDebug("AnalyzerDTCStubs") << "Constructed with InputTag: " << inputTag_;
  }

  void AnalyzerDTCStubs::beginJob() 
  {
    LogDebug("AnalyzerDTCStubs") << "beginJob()";

    edm::Service<TFileService> fs;

    const double lo = -0.5;
    const double hi = kMaxOccY_ - 0.5;
    hisAllDTCStubs_ = fs->make<TH1F>("HisAllDTCStubOccupancy",
                                     "All DTCs Stub Occupancy;# stubs (sum over all DTCs) per event;Events",
                                     kNBinsY_, lo, hi);
  }

  void AnalyzerDTCStubs::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) 
  {
    LogDebug("AnalyzerDTCStubs") << "beginRun(): run=" << iRun.run();

    const tt::Setup& setup = iSetup.getData(esGetTokenSetupBeginRun_);
    const int numDTCs = setup.numDTCs();

    edm::Service<TFileService> fs;

    const double xlo = -0.5;
    const double xhi = numDTCs - 0.5;
    const double ylo = -0.5;
    const double yhi = kMaxOccY_ - 0.5;

    h2DTCvsStubs_ = fs->make<TH2D>("DTCvsStubOccupancy",
                                   "Stub Occupancy per DTC;DTC id;# stubs per event",
                                   numDTCs, xlo, xhi,
                                   kNBinsY_, ylo, yhi);

    const double maxGbps         = 25.0;
    const int    nBinsGbps       = 100;
    const double loGbps          = 0.0;
    const double hiGbps          = maxGbps;

    hisThroughputGbps_ = fs->make<TH1F>("HisThroughputGbps",
                                        "Throughput;Gbps;Events",
                                        nBinsGbps, loGbps, hiGbps);

    (void)iRun;
  }

  void AnalyzerDTCStubs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
  {
    LogDebug("AnalyzerDTCStubs") << "analyze(): event=" << iEvent.id().event()
                                 << " run=" << iEvent.id().run()
                                 << " lumi=" << iEvent.id().luminosityBlock();

    const tt::Setup& setup      = iSetup.getData(esGetTokenSetup_);
    const int number_of_dtcs    = setup.numDTCs();
    const int number_of_modules = setup.numModulesPerDTC();

    edm::Handle<TTStubDetSetVec> handle;
    iEvent.getByToken(edGetTokenTTStubs_, handle);

    if (!handle.isValid()) 
    {
      edm::LogError("AnalyzerDTCStubs")
          << "Failed to get TTStubDetSetVec with InputTag " << inputTag_;
      return;
    }

    std::vector<std::vector<std::vector<TTStubRef>>> stubsDTCs(
        number_of_dtcs, std::vector<std::vector<TTStubRef>>(number_of_modules));

    for (auto module = handle->begin(); module != handle->end(); ++module) 
    {
        const DetId detId = module->detId() + setup.offsetDetIdDSV();
        tt::SensorModule* sm = setup.sensorModule(detId);

        std::vector<TTStubRef>& stubsModule = stubsDTCs[sm->dtcId()][sm->modId()];
        stubsModule.reserve(module->size());

        for (TTStubDetSet::const_iterator ttStub = module->begin(); ttStub != module->end(); ++ttStub)
            stubsModule.emplace_back(makeRefTo(handle, ttStub));
    }

    if (h2DTCvsStubs_) 
    {
      for (int dtcId = 0; dtcId < number_of_dtcs; ++dtcId) 
      {
        int nStubsThisDTC = 0;
        for (const auto& modVec : stubsDTCs[dtcId])
        {
          nStubsThisDTC += static_cast<int>(modVec.size());
        }
        h2DTCvsStubs_->Fill(dtcId, nStubsThisDTC);
      }
    }

    for (int dtcId = 0; dtcId < number_of_dtcs; ++dtcId)
    {
        int nStubsAllDTCs = 0;
        for (const auto& modVec : stubsDTCs[dtcId])
        {
            nStubsAllDTCs += static_cast<int>(modVec.size());
        }
        if (hisAllDTCStubs_) hisAllDTCStubs_->Fill(nStubsAllDTCs);
        if (hisThroughputGbps_) 
        {
            const double gbps = (static_cast<double>(nStubsAllDTCs) * kBitsPerStub_ * kEventRateHz_) / 1e9;
            hisThroughputGbps_->Fill(gbps);
        }
    }

  }

  void AnalyzerDTCStubs::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) 
  {
    LogDebug("AnalyzerDTCStubs") << "endRun(): run=" << iRun.run();
    (void)iRun;
    (void)iSetup;
  }

  void AnalyzerDTCStubs::endJob() 
  {
    LogDebug("AnalyzerDTCStubs") << "endJob()";
  }

}  // namespace trackerDTC

DEFINE_FWK_MODULE(trackerDTC::AnalyzerDTCStubs);
