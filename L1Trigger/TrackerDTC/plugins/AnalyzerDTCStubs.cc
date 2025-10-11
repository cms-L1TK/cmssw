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

// TTStubDetSetVec + related typedefs
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
// tt::Setup
#include "L1Trigger/TrackTrigger/interface/Setup.h"

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

        // Setup ES token
        edm::ESGetToken<tt::Setup, tt::SetupRcd> esGetTokenSetup_;
  };

  AnalyzerDTCStubs::AnalyzerDTCStubs(const edm::ParameterSet& iConfig)
      : inputTag_(iConfig.getParameter<edm::InputTag>("InputTag")) 
  {
    edGetTokenTTStubs_ = consumes<TTStubDetSetVec>(inputTag_);
    esGetTokenSetup_   = esConsumes();

    LogDebug("AnalyzerDTCStubs") << "Constructed with InputTag: " << inputTag_;
  }

  void AnalyzerDTCStubs::beginJob() 
  {
    LogDebug("AnalyzerDTCStubs") << "beginJob()";
  }

  void AnalyzerDTCStubs::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) 
  {
    LogDebug("AnalyzerDTCStubs") << "beginRun(): run=" << iRun.run();

    // Fetch Setup from EventSetup
    // LogDebug("AnalyzerDTCStubs") << "Setup retrieved. numDTCs=" << setup.numDTCs();

    (void)iRun;
  }

  void AnalyzerDTCStubs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
  {
    LogDebug("AnalyzerDTCStubs") << "analyze(): event=" << iEvent.id().event()
                                 << " run=" << iEvent.id().run()
                                 << " lumi=" << iEvent.id().luminosityBlock();
    (void)iSetup;

    edm::Handle<TTStubDetSetVec> handleTTStubs;
    iEvent.getByToken(edGetTokenTTStubs_, handleTTStubs);

    if (!handleTTStubs.isValid()) 
    {
      edm::LogError("AnalyzerDTCStubs")
          << "Failed to get TTStubDetSetVec with InputTag " << inputTag_;
      return;
    }

    const tt::Setup& setup      = iSetup.getData(esGetTokenSetup_);
    const int number_of_dtcs    = setup.numDTCs();
    const int number_of_modules = setup.numModulesPerDTC();

    std::vector<std::vector<std::vector<TTStubRef>>> stubsDTCs(
    number_of_dtcs, std::vector<std::vector<TTStubRef>>(number_of_modules));

    edm::Handle<TTStubDetSetVec> handle;
    iEvent.getByToken(edGetTokenTTStubs_, handle);

    for (auto module = handle->begin(); module != handle->end(); module++) 
    {
        const DetId detId = module->detId() + setup.offsetDetIdDSV();
        tt::SensorModule* sm = setup.sensorModule(detId);
        std::vector<TTStubRef>& stubsModule = stubsDTCs[sm->dtcId()][sm->modId()];
        stubsModule.reserve(module->size());
        for (TTStubDetSet::const_iterator ttStub = module->begin(); ttStub != module->end(); ttStub++)
            stubsModule.emplace_back(makeRefTo(handle, ttStub));
    }

    LogDebug("AnalyzerDTCStubs") << "Retrieved TTStubDetSetVec with " << handleTTStubs->size() << " DetSets.";
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