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

        // Setup ES token (Event transition, used in analyze)
        edm::ESGetToken<tt::Setup, tt::SetupRcd> esGetTokenSetup_;

        // Histogram: number of stubs on DTC 0 per event
        TH1F* hisDTC0Stubs_{nullptr};
  };

  AnalyzerDTCStubs::AnalyzerDTCStubs(const edm::ParameterSet& iConfig)
      : inputTag_(iConfig.getParameter<edm::InputTag>("InputTag")) 
  {
    usesResource("TFileService");  // we’ll book a histogram
    edGetTokenTTStubs_ = consumes<TTStubDetSetVec>(inputTag_);
    esGetTokenSetup_   = esConsumes();  // default = Event (used in analyze)

    LogDebug("AnalyzerDTCStubs") << "Constructed with InputTag: " << inputTag_;
  }

  void AnalyzerDTCStubs::beginJob() 
  {
    LogDebug("AnalyzerDTCStubs") << "beginJob()";

    // Book histogram once
    edm::Service<TFileService> fs;
    // Choose reasonable range; adjust later if you see saturation.
    // Here: range 0..4096 with 256 bins (like example style: maxOcc/16 bins).
    const int    maxOcc = 2000;
    const int    nBins  = 100;
    const double lo     = -0.5;
    const double hi     = maxOcc - 0.5;

    hisDTC0Stubs_ = fs->make<TH1F>("HisDTC0StubOccupancy",
                                   "DTC 0 Stub Occupancy;# stubs in DTC 0 per event;Events",
                                   nBins, lo, hi);
  }

  void AnalyzerDTCStubs::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) 
  {
    LogDebug("AnalyzerDTCStubs") << "beginRun(): run=" << iRun.run();
    (void)iRun;
    (void)iSetup;
  }

  void AnalyzerDTCStubs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
  {
    LogDebug("AnalyzerDTCStubs") << "analyze(): event=" << iEvent.id().event()
                                 << " run=" << iEvent.id().run()
                                 << " lumi=" << iEvent.id().luminosityBlock();

    // Fetch setup (Event transition)
    const tt::Setup& setup      = iSetup.getData(esGetTokenSetup_);
    const int number_of_dtcs    = setup.numDTCs();
    const int number_of_modules = setup.numModulesPerDTC();

    // Fetch TTStubDetSetVec
    edm::Handle<TTStubDetSetVec> handle;
    iEvent.getByToken(edGetTokenTTStubs_, handle);

    if (!handle.isValid()) 
    {
      edm::LogError("AnalyzerDTCStubs")
          << "Failed to get TTStubDetSetVec with InputTag " << inputTag_;
      return;
    }

    // Reorganize stubs into [DTC][module] containers
    std::vector<std::vector<std::vector<TTStubRef>>> stubsDTCs(
        number_of_dtcs, std::vector<std::vector<TTStubRef>>(number_of_modules));

    for (auto module = handle->begin(); module != handle->end(); ++module) 
    {
        // DetSetVec->detId + 1 = tk layout det id (per ProducerDTC pattern)
        const DetId detId = module->detId() + setup.offsetDetIdDSV();
        tt::SensorModule* sm = setup.sensorModule(detId);

        std::vector<TTStubRef>& stubsModule = stubsDTCs[sm->dtcId()][sm->modId()];
        stubsModule.reserve(module->size());

        for (TTStubDetSet::const_iterator ttStub = module->begin(); ttStub != module->end(); ++ttStub)
            stubsModule.emplace_back(makeRefTo(handle, ttStub));
    }

    // Count stubs in DTC 0 across all its modules
    int nStubsDTC0 = 0;
    if (!stubsDTCs.empty()) {
      const auto& dtc0 = stubsDTCs.at(0);
      for (const auto& modVec : dtc0)
        nStubsDTC0 += static_cast<int>(modVec.size());
    }

    if (hisDTC0Stubs_) hisDTC0Stubs_->Fill(nStubsDTC0);

    LogDebug("AnalyzerDTCStubs") << "DTC 0 has " << nStubsDTC0
                                 << " stubs in this event; DetSets in input = " << handle->size();
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
