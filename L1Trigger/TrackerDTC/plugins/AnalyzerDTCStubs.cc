#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

namespace trackerDTC {

  /*! \class  trackerDTC::AnalyzerDTCStubs
   *  \brief  Empty analyzer skeleton
   *  \author Andrew Mastronikolis
   *  \date   2025, October
   */

  class AnalyzerDTCStubs : public edm::stream::EDAnalyzer<> {
  public:

    explicit AnalyzerDTCStubs(const edm::ParameterSet&);
    ~AnalyzerDTCStubs() override = default;

  private:

    void beginRun   (const edm::Run&, const edm::EventSetup&)   override {}
    void analyze    (const edm::Event&, const edm::EventSetup&) override {}
    void endRun     (const edm::Run&, const edm::EventSetup&)   override {}

  };

  AnalyzerDTCStubs::AnalyzerDTCStubs(const edm::ParameterSet& iConfig) 
  {
    
  }

}  // namespace trackerDTC

DEFINE_FWK_MODULE(trackerDTC::AnalyzerDTCStubs);