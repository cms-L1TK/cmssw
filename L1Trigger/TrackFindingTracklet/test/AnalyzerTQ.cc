#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/StubAssociation.h"
#include "L1Trigger/TrackTrigger/interface/Associator.h"
#include "L1Trigger/TrackFindingTracklet/interface/DataFormats.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"

#include <TProfile.h>
#include <TTree.h>
#include <TH1F.h>

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <cmath>
#include <numeric>
#include <sstream>

namespace trklet {

  /*! \class  trklet::AnalyzerTQ
   *  \brief  Class to analyze emulated Track Quality 
   *  \author Thomas Schuh
   *  \date   2025, Aug
   */
  class AnalyzerTQ : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
  public:
    AnalyzerTQ(const edm::ParameterSet& iConfig);
    void beginJob() override {}
    void beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override;
    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    void endRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) override {}
    void endJob() override;

  private:
    // ED input token of stubs
    edm::EDGetTokenT<tt::StreamsStub> edGetTokenStubs_;
    // ED input token of tracks
    edm::EDGetTokenT<tt::StreamsTrack> edGetTokenTracks_;
    // ED output token for stub association for fake rate
    edm::EDGetTokenT<tt::StubAssociation> edGetTokenFake_;
    // ED output token for stub association for tracking efficiency
    edm::EDGetTokenT<tt::StubAssociation> edGetTokenEff_;
    // Associator token
    edm::ESGetToken<tt::Associator, trackerDTC::SetupRcd> esGetTokenAssociator_;
    // Setup token
    edm::ESGetToken<Setup, trackerDTC::SetupRcd> esGetTokenSetup_;
    // DataFormats token
    edm::ESGetToken<DataFormats, trackerDTC::SetupRcd> esGetTokenDataFormats_;
    // TTTrackAssociationMap Token (Failed Attempt to Work with Those below)
    edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> ttTrackMCTruthToken_;
    // TTTracks Token
    edm::EDGetTokenT<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> ttTrackToken_;
    //
    int numMVA_ = 8;
    //
    bool looseMatching_;
    //
    int nEvents_ = 0;
    // Histograms
    std::vector<TProfile*> prof_;
    TH1D* BDTScoreHist_Real;
    TH1D* BDTScoreHist_Fake;
    // printout
    std::stringstream log_;
    // private
    edm::InputTag MCTruthTrackInputTag;
    edm::InputTag L1TrackInputTag;
    // tree
    bool training_run;
    bool evaluation_run;
    TTree* tree_;
    std::vector<int> f_nstubs;
    std::vector<double> f_zT;
    std::vector<double> f_cot;
    std::vector<double> f_chi20;
    std::vector<double> f_chi21;
    std::vector<int> f_ngaps;

    // private helper methods
    std::vector<TTStubRef> getStubRefs(int region, int frame, int numRegions, const tt::StreamsStub& streamsStub);
    void clearBranches();
  };

  AnalyzerTQ::AnalyzerTQ(const edm::ParameterSet& iConfig)
      : looseMatching_(iConfig.getParameter<bool>("LooseMatching")) {
    L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
    MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
    training_run = iConfig.getParameter<bool>("TrainingMode");
    evaluation_run = iConfig.getParameter<bool>("EvaluationMode");
    usesResource("TFileService");
    // book in- and output ED products
    const std::string& labelKF = iConfig.getParameter<std::string>("OutputLabelKF");
    const std::string& labelTQ = iConfig.getParameter<std::string>("OutputLabelTQ");
    const std::string& branchStubs = iConfig.getParameter<std::string>("BranchStubs");
    const std::string& branchTracks = iConfig.getParameter<std::string>("BranchTracks");
    edGetTokenStubs_ = consumes(edm::InputTag(labelKF, branchStubs));
    edGetTokenTracks_ = consumes(edm::InputTag(labelTQ, branchTracks));
    const std::string& labelMC = iConfig.getParameter<std::string>("LabelMC");
    const std::string& branchFake = iConfig.getParameter<std::string>("BranchFake");
    const std::string& branchEff = iConfig.getParameter<std::string>("BranchEff");
    edGetTokenFake_ = consumes(edm::InputTag(labelMC, branchFake));
    edGetTokenEff_ = consumes(edm::InputTag(labelMC, branchEff));
    ttTrackToken_ = consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>(L1TrackInputTag);
    ttTrackMCTruthToken_ = consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>>(MCTruthTrackInputTag);
    // book ES products
    esGetTokenAssociator_ = esConsumes();
    esGetTokenSetup_ = esConsumes();
    esGetTokenDataFormats_ = esConsumes();
    // log config
    log_.setf(std::ios::fixed, std::ios::floatfield);
    log_.precision(4);
  }

  void AnalyzerTQ::beginRun(const edm::Run& iEvent, const edm::EventSetup& iSetup) {
    // book histograms
    edm::Service<TFileService> fs;
    TFileDirectory dir;
    dir = fs->mkdir("TQ");
    BDTScoreHist_Real = dir.make<TH1D>("BDTScoreHist_Real", ";BDT Score;Entries", 7, 0, numMVA_ - 1);
    BDTScoreHist_Fake = dir.make<TH1D>("BDTScoreHist_Fake", ";BDT Score;Entries", 7, 0, numMVA_ - 1);
    prof_ = std::vector<TProfile*>(numMVA_);
    for (int mva = 0; mva < numMVA_; mva++) {
      prof_[mva] = dir.make<TProfile>(("Counts for MVA" + std::to_string(mva)).c_str(), ";", 4, 0.5, 4.5);
      prof_[mva]->GetXaxis()->SetBinLabel(1, "All TPs");
      prof_[mva]->GetXaxis()->SetBinLabel(2, "All Tracks");
      prof_[mva]->GetXaxis()->SetBinLabel(3, "Matched to any Tracks");
      prof_[mva]->GetXaxis()->SetBinLabel(4, "Found Perfect TPs");
    }

    if (training_run) {
      tree_ = fs->make<TTree>("TrackQualityAttributes", "Track Quality Analysis");
      // Branches to Train Against
      tree_->Branch("trk_feature_1", &f_nstubs);
      tree_->Branch("trk_feature_2", &f_zT);
      tree_->Branch("trk_feature_3", &f_cot);
      tree_->Branch("trk_feature_4", &f_chi20);
      tree_->Branch("trk_feature_5", &f_chi21);
      tree_->Branch("trk_feature_6", &f_ngaps);
    }
  }

  void AnalyzerTQ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // TSCHUH PART
    // read in tracks and stubs products
    const tt::StreamsStub& streamsStub = iEvent.get(edGetTokenStubs_);
    const tt::StreamsTrack& streamsTrack = iEvent.get(edGetTokenTracks_);
    const Setup* setup = &iSetup.getData(esGetTokenSetup_);
    const DataFormats* df = &iSetup.getData(esGetTokenDataFormats_);
    // read in MCTruth
    tt::Associator forFake = iSetup.getData(esGetTokenAssociator_);
    tt::Associator forEff = iSetup.getData(esGetTokenAssociator_);
    forFake.consume(iEvent.get(edGetTokenFake_));
    forEff.consume(iEvent.get(edGetTokenEff_));
    for (TProfile* prof : prof_)
      prof->Fill(1, forEff.numTPs());
    // analyze and associate tracks with TrackingParticles per mva categorie
    for (int mva = 0; mva < numMVA_; mva++) {
      std::set<TPPtr> tpPtrsPerfect;
      int nTracks(0);
      int nMatched(0);
      for (int region = 0; region < setup->sysNumRegion(); region++) {
        const tt::StreamTrack& streamTrack = streamsTrack[region * 2 + 1];
        const int numFrames = streamTrack.size();
        for (int frame = 0; frame < numFrames; frame++) {
          if (streamTrack[frame].first.isNull())
            continue;
          const TrackTQ trackTQ(streamTrack[frame], df);
          if (trackTQ.mva() < mva)
            continue;
          nTracks++;
          const int offset = region * setup->kfNumLayers();
          std::vector<TTStubRef> ttStubRefs;
          ttStubRefs.reserve(setup->kfNumLayers());
          for (int layer = 0; layer < setup->kfNumLayers(); layer++) {
            const TTStubRef& ttStubRef = streamsStub[offset + layer][frame].first;
            if (ttStubRef.isNonnull())
              ttStubRefs.push_back(ttStubRef);
          }
          const std::vector<TPPtr>& any =
              looseMatching_ ? forFake.associate(ttStubRefs) : forFake.associateFinal(ttStubRefs);
          if (any.empty())
            continue;
          nMatched++;
          const std::vector<TPPtr> perfect = forEff.associateFinal(ttStubRefs);
          tpPtrsPerfect.insert(perfect.begin(), perfect.end());
        }
      }
      prof_[mva]->Fill(2, nTracks);
      prof_[mva]->Fill(3, nMatched);
      prof_[mva]->Fill(4, tpPtrsPerfect.size());
    }
    // END OF TSCHUH PART

    // AMASTRON PART
    if (training_run || evaluation_run) {
      // Initialize MCTruthTTTrackHandle
      edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTTrackHandle;
      iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

      // Initialize TTTrackHandle
      edm::Handle<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> TTTrackHandle;
      iEvent.getByToken(ttTrackToken_, TTTrackHandle);

      if (training_run)
        clearBranches();

      // Initialize Dummy Variable for L1 TTTrack Counting and TQ BDT Features
      int this_l1track = 0;
      double feature_1 = 0;
      double feature_2 = 0;
      double feature_3 = 0;
      double feature_4 = 0;
      double feature_5 = 0;
      double feature_6 = 0;

      std::vector<TTTrack<Ref_Phase2TrackerDigi_>>::const_iterator iterL1Track;
      for (iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++) {
        // Fetch L1 TTTrack from Collection.
        edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> l1track_ptr(TTTrackHandle, this_l1track);

        // Get StubRefs for this TTTrack.
        std::vector<TTStubRef> l1track_stubrefs = l1track_ptr->getStubRefs();

        for (int region = 0; region < setup->sysNumRegion(); region++) {
          // Fetch Stream with Offset #0, which corresponds to TQ Input (= KF Output) Tracks.
          const tt::StreamTrack& streamTrack0 = streamsTrack[region * 2 + 0];

          // Fetch Stream with Offset #1, which corresponds to TQ Converted Tracks.
          const tt::StreamTrack& streamTrack1 = streamsTrack[region * 2 + 1];

          // Fetch Num Frames = Num Track Objects
          const int numFrames1 = streamTrack1.size();
          const int numFrames0 = streamTrack0.size();
          int numFrames = 0;

          // Throw Exception when Said Frames are mismatched.
          if (numFrames1 != numFrames0) {
            throw cms::Exception("FrameMismatch")
                << "numFrames1 (" << numFrames1 << ") != numFrames0 (" << numFrames0 << "). This should not happen.";
          } else {
            numFrames = numFrames1;
          }

          // Iterate Over Common
          for (int frame = 0; frame < numFrames; frame++) {
            if (streamTrack1[frame].first.isNull() || streamTrack0[frame].first.isNull())
              continue;

            // Fetch TQ Frame for the Chi Squared Values.
            const TrackTQ trackTQ(streamTrack1[frame], df);

            // Fetch KF Frame for the z0, cot and hitpattern.
            const TrackKF trackKF(streamTrack0[frame], df);

            // Fetch DataFormats for said variables.
            const DataFormat& dfChi20 = df->format(Variable::chi20, Process::tq);
            const DataFormat& dfChi21 = df->format(Variable::chi21, Process::tq);
            const DataFormat& dfZT = df->format(Variable::z0, Process::tq);
            const DataFormat& dfCot = df->format(Variable::cot, Process::tq);

            std::vector<TTStubRef> stubRefs = getStubRefs(region, frame, setup->kfNumLayers(), streamsStub);
            bool matched = (stubRefs == l1track_stubrefs);

            // If there is no L1 TTTrack Stub Refs that Match TQ/KF Track StubRefs, Simply Continue.
            if (!matched) {
              continue;
            } else {
              // Compute BDT Input Features
              const double z0_scale = std::pow(2, 3);
              const double tanL_scale = std::pow(2, 7);
              feature_1 = stubRefs.size();
              feature_2 = dfZT.integer(trackKF.z0()) / z0_scale;
              feature_3 = dfCot.integer(trackKF.cot()) / tanL_scale;
              feature_4 = dfChi20.integer(trackTQ.chi20());
              feature_5 = dfChi21.integer(trackTQ.chi21());
              const TTBV hitPattern = trackTQ.hitPattern();
              feature_6 = hitPattern.count(hitPattern.plEncode(), hitPattern.pmEncode(), false);

              // Push to Event Vector.
              f_nstubs.push_back(feature_1);
              f_zT.push_back(feature_2);
              f_cot.push_back(feature_3);
              f_chi20.push_back(feature_4);
              f_chi21.push_back(feature_5);
              f_ngaps.push_back(feature_6);

              // Break out of the loop when a match is found.
              // Oddly enough, it certain cases, two or more matches are found!
              break;
            }
          }
        }
        this_l1track++;
      }

      if (training_run)
        tree_->Fill();

      nEvents_++;
    }
  }

  void AnalyzerTQ::endJob() {
    if (nEvents_ == 0)
      return;
    // printout summary
    log_ << "                         TQ  SUMMARY                         " << std::endl;
    for (int mva = 0; mva < numMVA_; mva++) {
      const double allTracks = prof_[mva]->GetBinContent(2);
      const double allMatched = prof_[mva]->GetBinContent(3);
      const double numPerfect = prof_[mva]->GetBinContent(4);
      const double allTPs = prof_[mva]->GetBinContent(1);
      const double fracFake = (allTracks - allMatched) / allTracks;
      const double effPerfect = numPerfect / allTPs;
      log_ << "mva " << mva << " (effi: " << effPerfect << ", fake rate: " << fracFake << ") numMatched "
           << std::setw(3) << (int)std::round(allMatched) << " numFake " << std::setw(3)
           << (int)std::round(allTracks - allMatched) << std::endl;
    }
    log_ << "=============================================================";
    edm::LogPrint(moduleDescription().moduleName()) << log_.str();
  }
  std::vector<TTStubRef> AnalyzerTQ::getStubRefs(int region,
                                                 int frame,
                                                 int numKFLayers,
                                                 const tt::StreamsStub& streamsStub) {
    const int offset = region * numKFLayers;
    std::vector<TTStubRef> ttStubRefs;
    ttStubRefs.reserve(numKFLayers);
    for (int layer = 0; layer < numKFLayers; layer++) {
      const TTStubRef& ttStubRef = streamsStub[offset + layer][frame].first;
      if (ttStubRef.isNonnull())
        ttStubRefs.push_back(ttStubRef);
    }
    return ttStubRefs;
  }
  void AnalyzerTQ::clearBranches() {
    f_zT.clear();
    f_cot.clear();
    f_chi20.clear();
    f_chi21.clear();
    f_nstubs.clear();
    f_ngaps.clear();
  }
}  // namespace trklet

DEFINE_FWK_MODULE(trklet::AnalyzerTQ);
