#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/StubAssociation.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"

#include <TProfile.h>

#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

namespace tt {

  class AnalyzerMC : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
  public:
    // Constructor
    explicit AnalyzerMC(const edm::ParameterSet& iConfig);
    // Destructor
    ~AnalyzerMC() override {}
    // Mandatory methods
    void beginJob() override {}
    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    void beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override;
    void endRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override {}
    void endJob() override;

  private:
    // ED input token of TTStubs
    edm::EDGetTokenT<TTStubDetSetVec> edGetTokenTTStubs_;
    // ED input token of TTClusterAssociation
    edm::EDGetTokenT<TTClusterAssMap> getTokenTTClusterAssMap_;
    // ED input token of StubAssociation with reconstructable TPs
    edm::EDGetTokenT<StubAssociation> getTokenFake_;
    // ED input token of StubAssociation with selected TPs
    edm::EDGetTokenT<StubAssociation> getTokenEff_;
    // Setup token
    edm::ESGetToken<tt::Setup, tt::SetupRcd> esGetTokenSetup_;
    // Histograms
    TProfile* prof_;
    // printout
    std::stringstream log_;
  };

  AnalyzerMC::AnalyzerMC(edm::ParameterSet const& iConfig) {
    const std::string& label = iConfig.getParameter<std::string>("StubAssociation");
    const std::string& branchFake = iConfig.getParameter<std::string>("BranchFake");
    const std::string& branchEff = iConfig.getParameter<std::string>("BranchEff");
    const edm::InputTag& inputTagTTStubs = iConfig.getParameter<edm::InputTag>("InputTagTTStubs");
    const edm::InputTag& inputTagTTClusterAssMap = iConfig.getParameter<edm::InputTag>("InputTagTTClusterAssMap");
    edGetTokenTTStubs_ = consumes(inputTagTTStubs);
    getTokenTTClusterAssMap_ = consumes(inputTagTTClusterAssMap);
    getTokenFake_ = consumes(edm::InputTag(label, branchFake));
    getTokenEff_ = consumes(edm::InputTag(label, branchEff));
    // book ES product
    esGetTokenSetup_ = esConsumes();
    // log config
    log_.setf(std::ios::fixed, std::ios::floatfield);
    log_.precision(4);
  }

  void AnalyzerMC::beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) {
    // book histograms
    edm::Service<TFileService> fs;
    TFileDirectory dir;
    dir = fs->mkdir("MC");
    prof_ = dir.make<TProfile>("Counts", ";", 6, 0.5, 6.5);
    prof_->GetXaxis()->SetBinLabel(1, "Stubs");
    prof_->GetXaxis()->SetBinLabel(2, "matched Stubs");
    prof_->GetXaxis()->SetBinLabel(3, "reco Stubs");
    prof_->GetXaxis()->SetBinLabel(4, "any TPs");
    prof_->GetXaxis()->SetBinLabel(5, "reco TPs");
    prof_->GetXaxis()->SetBinLabel(6, "eff TPs");
  }

  void AnalyzerMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // helper class to store configurations
    const Setup* setup = &iSetup.getData(esGetTokenSetup_);
    // count stubs & matched stubs
    const TTStubDetSetVec& ttStubDetSetVec = iEvent.get(edGetTokenTTStubs_);
    const TTClusterAssMap& ttClusterAssMap = iEvent.get(getTokenTTClusterAssMap_);
    int nStubs(0);
    int nMatched(0);
    for (const auto& module : ttStubDetSetVec) {
      nStubs += module.size();
      for (const auto& ttStub : module) {
        bool matched(false);
        for (unsigned int iClus = 0; iClus < 2; iClus++)
          for (const TPPtr& tpPtr : ttClusterAssMap.findTrackingParticlePtrs(ttStub.clusterRef(iClus)))
            if (tpPtr.isNonnull())
              matched = true;
        if (matched)
          nMatched++;
      }
    }
    // get number of TPs
    const StubAssociation& forFake = iEvent.get(getTokenFake_);
    const StubAssociation& forEff = iEvent.get(getTokenEff_);
    const double numRegions = setup->numRegions();
    // store
    prof_->Fill(1, nStubs / numRegions);
    prof_->Fill(2, nMatched / numRegions);
    prof_->Fill(3, forFake.numStubs() / numRegions);
    prof_->Fill(4, ttClusterAssMap.getTrackingParticleToTTClustersMap().size());
    prof_->Fill(5, forFake.numTPs() / numRegions);
    prof_->Fill(6, forEff.numTPs() / numRegions);
  }

  // prints out Monte Carlo summary
  void AnalyzerMC::endJob() {
    const double numStubs = prof_->GetBinContent(1);
    const double numStubsMatched = prof_->GetBinContent(2);
    const double numStubsReco = prof_->GetBinContent(3);
    const double numTPsAny = prof_->GetBinContent(4);
    const double numTPsReco = prof_->GetBinContent(5);
    const double numTPsEff = prof_->GetBinContent(6);
    const double errStubs = prof_->GetBinError(1);
    const double errStubsMatched = prof_->GetBinError(2);
    const double errStubsReco = prof_->GetBinError(3);
    const double errTPsAny = prof_->GetBinError(4);
    const double errTPsReco = prof_->GetBinError(5);
    const double errTPsEff = prof_->GetBinError(6);
    const std::vector<double> nums = {numStubs, numStubsMatched, numStubsReco, numTPsAny, numTPsReco, numTPsEff};
    const std::vector<double> errs = {errStubs, errStubsMatched, errStubsReco, errTPsAny, errTPsReco, errTPsEff};
    const int wNums = std::ceil(std::log10(*std::max_element(nums.begin(), nums.end()))) + 5;
    const int wErrs = std::ceil(std::log10(*std::max_element(errs.begin(), errs.end()))) + 5;
    log_ << "=============================================================" << std::endl;
    log_ << "                         Monte Carlo  SUMMARY                         " << std::endl;
    log_ << "number of stubs         per TFP = " << std::setw(wNums) << numStubs << " +- " << std::setw(wErrs)
         << errStubs << std::endl;
    log_ << "number of matched stubs per TFP = " << std::setw(wNums) << numStubsMatched << " +- " << std::setw(wErrs)
         << errStubsMatched << std::endl;
    log_ << "number of    reco stubs per TFP = " << std::setw(wNums) << numStubsReco << " +- " << std::setw(wErrs)
         << errStubsReco << std::endl;
    log_ << "number of any TPs       per TFP = " << std::setw(wNums) << numTPsAny << " +- " << std::setw(wErrs)
         << errTPsAny << std::endl;
    log_ << "number of reco TPs      per TFP = " << std::setw(wNums) << numTPsReco << " +- " << std::setw(wErrs)
         << errTPsReco << std::endl;
    log_ << "number of TPs for eff   per TFP = " << std::setw(wNums) << numTPsEff << " +- " << std::setw(wErrs)
         << errTPsEff << std::endl;
    log_ << "=============================================================";
    edm::LogPrint(moduleDescription().moduleName()) << log_.str();
  }

}  // namespace tt

DEFINE_FWK_MODULE(tt::AnalyzerMC);
