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
    // ED input token of KF metadata
    edm::EDGetTokenT<tt::StreamsTrack> edGetTokenMetadata_;
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
    // TQ Features [1 tp 6]
    std::vector<int>    trk_tq_feature_1;
    std::vector<double> trk_tq_feature_2;
    std::vector<double> trk_tq_feature_3;
    std::vector<double> trk_tq_feature_4;
    std::vector<double> trk_tq_feature_5;
    std::vector<int>    trk_tq_feature_6;
    // L1 Track Attributes
    std::vector<int>    trk_hitpattern;
    std::vector<double> trk_pt;
    std::vector<double> trk_z0;
    std::vector<double> trk_eta;
    std::vector<double> trk_tanL;
    std::vector<double> trk_phi;
    std::vector<double> trk_chi2bend;
    std::vector<int>    trk_matched_to_tp;
    std::vector<int>    trk_matchedtp_pdgid;
    std::vector<int>    trk_num_iterations;
    std::vector<double> trk_matchtp_pt;
    std::vector<double> trk_matchtp_eta;
    std::vector<double> trk_matchtp_phi;
    std::vector<double> trk_matchtp_z0;
    std::vector<double> trk_matchtp_tanL;
    // Per-layer chi2 values
    std::vector<double> trk_chi20_layer_0;
    std::vector<double> trk_chi20_layer_1;
    std::vector<double> trk_chi20_layer_2;
    std::vector<double> trk_chi20_layer_3;
    std::vector<double> trk_chi20_layer_4;
    std::vector<double> trk_chi20_layer_5;
    std::vector<double> trk_chi20_layer_6;
    // Per-layer chi2 values
    std::vector<double> trk_chi21_layer_0;
    std::vector<double> trk_chi21_layer_1;
    std::vector<double> trk_chi21_layer_2;
    std::vector<double> trk_chi21_layer_3;
    std::vector<double> trk_chi21_layer_4;
    std::vector<double> trk_chi21_layer_5;
    std::vector<double> trk_chi21_layer_6;
    // Per-layer bendFE values
    std::vector<double> trk_bendFE_layer_0;
    std::vector<double> trk_bendFE_layer_1;
    std::vector<double> trk_bendFE_layer_2;
    std::vector<double> trk_bendFE_layer_3;
    std::vector<double> trk_bendFE_layer_4;
    std::vector<double> trk_bendFE_layer_5;
    std::vector<double> trk_bendFE_layer_6;

    // private helper methods
    std::vector<TTStubRef> getStubRefs(int region, int frame, int numRegions, const tt::StreamsStub& streamsStub);
    void clearTrackBranches();
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
    const std::string& branchMetadata = iConfig.getParameter<std::string>("BranchMetadata");
    edGetTokenStubs_ = consumes(edm::InputTag(labelKF, branchStubs));
    edGetTokenTracks_ = consumes(edm::InputTag(labelTQ, branchTracks));
    edGetTokenMetadata_ = consumes(edm::InputTag(labelKF, branchMetadata));
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
      // TQ Attributes
      tree_->Branch("trk_feature_1", &trk_tq_feature_1);
      tree_->Branch("trk_feature_2", &trk_tq_feature_2);
      tree_->Branch("trk_feature_3", &trk_tq_feature_3);
      tree_->Branch("trk_feature_4", &trk_tq_feature_4);
      tree_->Branch("trk_feature_5", &trk_tq_feature_5);
      tree_->Branch("trk_feature_6", &trk_tq_feature_6);
      // Other L1 Track Attributes
      tree_->Branch("trk_hitpattern", &trk_hitpattern);
      tree_->Branch("trk_pt", &trk_pt);
      tree_->Branch("trk_z0", &trk_z0);
      tree_->Branch("trk_eta", &trk_eta);
      tree_->Branch("trk_tanL", &trk_tanL);
      tree_->Branch("trk_phi", &trk_phi);
      tree_->Branch("trk_chi2bend", &trk_chi2bend);
      tree_->Branch("trk_num_iterations", &trk_num_iterations);
      tree_->Branch("trk_matchtp_eventtype", &trk_matched_to_tp);
      tree_->Branch("trk_matchtp_pdgid", &trk_matchedtp_pdgid);
      tree_->Branch("trk_matchtp_pt", &trk_matchtp_pt);
      tree_->Branch("trk_matchtp_eta", &trk_matchtp_eta);
      tree_->Branch("trk_matchtp_phi", &trk_matchtp_phi);
      tree_->Branch("trk_matchtp_z0", &trk_matchtp_z0);
      tree_->Branch("trk_matchtp_tanL", &trk_matchtp_tanL);
      // Per-layer chi2 branches
      tree_->Branch("trk_chi20_layer_0", &trk_chi20_layer_0);
      tree_->Branch("trk_chi20_layer_1", &trk_chi20_layer_1);
      tree_->Branch("trk_chi20_layer_2", &trk_chi20_layer_2);
      tree_->Branch("trk_chi20_layer_3", &trk_chi20_layer_3);
      tree_->Branch("trk_chi20_layer_4", &trk_chi20_layer_4);
      tree_->Branch("trk_chi20_layer_5", &trk_chi20_layer_5);
      tree_->Branch("trk_chi20_layer_6", &trk_chi20_layer_6);
      // Per-layer chi2 branches
      tree_->Branch("trk_chi21_layer_0", &trk_chi21_layer_0);
      tree_->Branch("trk_chi21_layer_1", &trk_chi21_layer_1);
      tree_->Branch("trk_chi21_layer_2", &trk_chi21_layer_2);
      tree_->Branch("trk_chi21_layer_3", &trk_chi21_layer_3);
      tree_->Branch("trk_chi21_layer_4", &trk_chi21_layer_4);
      tree_->Branch("trk_chi21_layer_5", &trk_chi21_layer_5);
      tree_->Branch("trk_chi21_layer_6", &trk_chi21_layer_6);
      // Per-layer bendBE() branches
      tree_->Branch("trk_bendbe_0", &trk_bendFE_layer_0);
      tree_->Branch("trk_bendbe_1", &trk_bendFE_layer_1);
      tree_->Branch("trk_bendbe_2", &trk_bendFE_layer_2);
      tree_->Branch("trk_bendbe_3", &trk_bendFE_layer_3);
      tree_->Branch("trk_bendbe_4", &trk_bendFE_layer_4);
      tree_->Branch("trk_bendbe_5", &trk_bendFE_layer_5);
      tree_->Branch("trk_bendbe_6", &trk_bendFE_layer_6);
    }
  }

  void AnalyzerTQ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // TSCHUH PART
    // read in tracks and stubs products
    const tt::StreamsStub& streamsStub = iEvent.get(edGetTokenStubs_);
    const tt::StreamsTrack& streamsTrack = iEvent.get(edGetTokenTracks_);
    const tt::StreamsTrack& streamsMetadata = iEvent.get(edGetTokenMetadata_);
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

      edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTTrackHandle;
      iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);
      edm::Handle<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> TTTrackHandle;
      iEvent.getByToken(ttTrackToken_, TTTrackHandle);

      // Build a map from stub set to L1Track index for fast lookup
      std::map<std::set<TTStubRef>, int> stubSetToTrackIndex;

      if (TTTrackHandle.isValid()) {
          // std::cout << "TTTrackHandle is valid, size = " << TTTrackHandle->size() << std::endl;
          for (size_t idx = 0; idx < TTTrackHandle->size(); ++idx) {
              const auto& ttTrack = (*TTTrackHandle)[idx];
              std::set<TTStubRef> trackStubSet;
              for (const auto& stub : ttTrack.getStubRefs()) {
                  if (stub.isNonnull())
                      trackStubSet.insert(stub);
              }
              stubSetToTrackIndex[trackStubSet] = idx;
              // std::cout << "L1Track " << idx << " has " << trackStubSet.size() << " stubs" << std::endl;
          }
      } else {
          std::cout << "[AnalyzerTQ] TTTrackHandle is NOT valid" << std::endl;
      }

      if (training_run) {
        clearTrackBranches();
      }

      for (int region = 0; region < setup->sysNumRegion(); region++) {

        const tt::StreamTrack& streamTrack = streamsTrack[region * 2 + 1];
        const tt::StreamTrack& streamMetadata = streamsMetadata[region];
        const int numFrames = streamTrack.size();

        for (int frame = 0; frame < numFrames; frame++) {

          if (streamTrack[frame].first.isNull())
            continue;

          const TrackTQ trackTQ (streamTrack[frame], df);
          const DataFormat& dfChi20 = df->format(Variable::chi20, Process::tq);
          const DataFormat& dfChi21 = df->format(Variable::chi21, Process::tq);
          const DataFormat& dfZT = df->format(Variable::z0, Process::tq);
          const DataFormat& dfCot = df->format(Variable::cot, Process::tq);
          const int numIterationsKF = streamMetadata[frame].second.to_ullong();

          auto stubRefs = getStubRefs(region, frame, setup->kfNumLayers(), streamsStub);
          std::set<TTStubRef> tqStubSet;
          for (const auto& stub : stubRefs) {
            if (stub.isNonnull())
              tqStubSet.insert(stub);
          }
          
          bool matched = false;
          int matchedIndex = -1;
          auto it = stubSetToTrackIndex.find(tqStubSet);
          if (it != stubSetToTrackIndex.end()) {
            matched = true;
            matchedIndex = it->second;
          }
          
          // Print results
          // std::cout << "TQ Track " << frame << ": " 
          //           << tqStubSet.size() << " stubs, "
          //           << "matched = " << (matched ? "YES" : "NO");
          // if (matched) {
          //     std::cout << " (L1Track index " << matchedIndex << ")";
          // }
          // std::cout << std::endl;

          if (!matched) continue;

          const auto& matchedTrack = (*TTTrackHandle)[matchedIndex];

          edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> l1track_ptr(TTTrackHandle, matchedIndex);

          float tmp_trk_pt        = l1track_ptr->pt();
          float tmp_trk_eta       = l1track_ptr->eta();
          float tmp_trk_phi       = l1track_ptr->phi();
          float tmp_trk_z0        = l1track_ptr->z0();
          float tmp_trk_tanL      = l1track_ptr->tanL();
          float tmp_trk_bendchi2  = l1track_ptr->chi2BendRed();
          int   tmp_trk_genuine   = 0;
          

          int   tmp_matchtp_pdgid = -999;
          float tmp_matchtp_pt    = -999;
          float tmp_matchtp_eta   = -999;
          float tmp_matchtp_phi   = -999;
          float tmp_matchtp_z0    = -999;
          float tmp_matchtp_tanL  = -999;

          if (MCTruthTTTrackHandle->isGenuine(l1track_ptr))
            tmp_trk_genuine = 1;
          int tmp_matchtp_eventtype = -999;
          if (tmp_trk_genuine) {
            edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);
            if (my_tp.isNull()) assert(false);
            int tmp_eventid = my_tp->eventId().event();
            if (tmp_eventid > 0)
              tmp_matchtp_eventtype = 2;
            else
              tmp_matchtp_eventtype = 1;
            tmp_matchtp_pdgid = my_tp->pdgId();
            tmp_matchtp_pt    = my_tp->pt();
            tmp_matchtp_eta   = my_tp->eta();
            tmp_matchtp_phi   = my_tp->phi();
            tmp_matchtp_z0    = my_tp->z0();
            tmp_matchtp_tanL  = my_tp->tanl();
          }

          if (evaluation_run) {
            const float mva = trackTQ.mva();
            if (tmp_matchtp_eventtype != -999) {
              BDTScoreHist_Real->Fill(mva);
            } else {
              BDTScoreHist_Fake->Fill(mva);
            }
          }

          const int offset = region * setup->kfNumLayers();

          double Chi20F = 0.0;
          std::vector<double> chi20_layers(setup->kfNumLayers(), 0.0);

          double Chi21F = 0.0;
          std::vector<double> chi21_layers(setup->kfNumLayers(), 0.0);

          std::vector<double> bend_layers(setup->kfNumLayers(), 0.0);

          for (int layer = 0; layer < setup->kfNumLayers(); layer++) {
            const TTStubRef& ttStubRef = streamsStub[offset + layer][frame].first;
            if (ttStubRef.isNonnull()) {
              // Create StubKF from the frame stub
              const tt::FrameStub& frameStub = streamsStub[offset + layer][frame];
              StubKF stub(frameStub, df);
              if (!stub)
                continue;
              
              const double m02 = std::pow(stub.phi(), 2);
              const double m12 = std::pow(stub.z(), 2);
              const double invV0 = 3. / std::pow(stub.dPhi(), 2);
              const double invV1 = 3. / std::pow(stub.dZ(), 2);
              Chi20F += dfChi20.limit(dfChi20.digi(m02 * invV0));
              Chi21F += dfChi21.limit(dfChi21.digi(m12 * invV1));
              
              // Use stub.layerId() to index into the vector
              int layerId = stub.layerId() - 1; // [0 to 6]
              if (layerId >= 0 && layerId < setup->kfNumLayers()) {
                chi20_layers[layerId] = dfChi20.limit(dfChi20.digi(m02 * invV0));
                chi21_layers[layerId] = dfChi21.limit(dfChi21.digi(m12 * invV1));
                bend_layers[layerId] = (stub.frame().first)->bendBE();
              }
            }
          }
          if (training_run) {

            // fetch attributes (as they are on TrackQuality.cc)
            const double z0_scale = std::pow(2, 3);
            const double tanL_scale = std::pow(2, 7);
            const double feature_1 = stubRefs.size();
            const double feature_2 = dfZT.integer(trackTQ.z0()) / z0_scale;
            const double feature_3 = dfCot.integer(trackTQ.cot()) / tanL_scale;
            const double feature_4 = dfChi20.integer(trackTQ.chi20());
            const double feature_5 = dfChi21.integer(trackTQ.chi21());
            const TTBV hitPattern = trackTQ.hitPattern();
            const double feature_6 = hitPattern.count(hitPattern.plEncode(), hitPattern.pmEncode(), false);

            // Push TQ Attributes
            trk_tq_feature_1.push_back(feature_1);
            trk_tq_feature_2.push_back(feature_2);
            trk_tq_feature_3.push_back(feature_3);
            trk_tq_feature_4.push_back(feature_4);
            trk_tq_feature_5.push_back(feature_5);
            trk_tq_feature_6.push_back(feature_6);
            // Push Other L1 Track Attributes
            trk_hitpattern.push_back(matchedTrack.hitPattern());
            trk_pt.push_back(tmp_trk_pt);
            trk_z0.push_back(tmp_trk_z0);
            trk_eta.push_back(tmp_trk_eta);
            trk_tanL.push_back(tmp_trk_tanL);
            trk_phi.push_back(tmp_trk_phi);
            trk_chi2bend.push_back(tmp_trk_bendchi2);            
            trk_num_iterations.push_back(numIterationsKF);
            trk_matched_to_tp.push_back(tmp_matchtp_eventtype);
            trk_matchedtp_pdgid.push_back(tmp_matchtp_pdgid);
            trk_matchtp_pt.push_back(tmp_matchtp_pt);
            trk_matchtp_eta.push_back(tmp_matchtp_eta);
            trk_matchtp_phi.push_back(tmp_matchtp_phi);
            trk_matchtp_z0.push_back(tmp_matchtp_z0);
            trk_matchtp_tanL.push_back(tmp_matchtp_tanL);
            // Fill per-layer chi2 values (using layerId as index)
            trk_chi20_layer_0.push_back(chi20_layers[0]);
            trk_chi20_layer_1.push_back(chi20_layers[1]);
            trk_chi20_layer_2.push_back(chi20_layers[2]);
            trk_chi20_layer_3.push_back(chi20_layers[3]);
            trk_chi20_layer_4.push_back(chi20_layers[4]);
            trk_chi20_layer_5.push_back(chi20_layers[5]);
            trk_chi20_layer_6.push_back(chi20_layers[6]);
            // Fill per-layer chi2 values (using layerId as index)
            trk_chi21_layer_0.push_back(chi21_layers[0]);
            trk_chi21_layer_1.push_back(chi21_layers[1]);
            trk_chi21_layer_2.push_back(chi21_layers[2]);
            trk_chi21_layer_3.push_back(chi21_layers[3]);
            trk_chi21_layer_4.push_back(chi21_layers[4]);
            trk_chi21_layer_5.push_back(chi21_layers[5]);
            trk_chi21_layer_6.push_back(chi21_layers[6]);
            // Fill per-layer chi2 values (using layerId as index)
            trk_bendFE_layer_0.push_back(bend_layers[0]);
            trk_bendFE_layer_1.push_back(bend_layers[1]);
            trk_bendFE_layer_2.push_back(bend_layers[2]);
            trk_bendFE_layer_3.push_back(bend_layers[3]);
            trk_bendFE_layer_4.push_back(bend_layers[4]);
            trk_bendFE_layer_5.push_back(bend_layers[5]);
            trk_bendFE_layer_6.push_back(bend_layers[6]);
          }
        }
      }
    }
    // END OF AMASTRON PART
    if (training_run) tree_->Fill();
    nEvents_++;
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

  std::vector<TTStubRef> AnalyzerTQ::getStubRefs(int region, int frame, int numKFLayers, const tt::StreamsStub& streamsStub) {
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

  void AnalyzerTQ::clearTrackBranches() {
    // TQ Attributes
    trk_tq_feature_1.clear();
    trk_tq_feature_2.clear();
    trk_tq_feature_3.clear();
    trk_tq_feature_4.clear();
    trk_tq_feature_5.clear();
    trk_tq_feature_6.clear();
    // L1 Track Attributes
    trk_hitpattern.clear();
    trk_pt.clear();
    trk_z0.clear();
    trk_eta.clear();
    trk_tanL.clear();
    trk_phi.clear();
    trk_chi2bend.clear();
    trk_matched_to_tp.clear();
    trk_matchedtp_pdgid.clear();
    trk_num_iterations.clear();
    trk_matchtp_pt.clear();
    trk_matchtp_eta.clear();
    trk_matchtp_phi.clear();
    trk_matchtp_z0.clear();
    trk_matchtp_tanL.clear();
    // Per-layer chi2 values (iteration 0)
    trk_chi20_layer_0.clear();
    trk_chi20_layer_1.clear();
    trk_chi20_layer_2.clear();
    trk_chi20_layer_3.clear();
    trk_chi20_layer_4.clear();
    trk_chi20_layer_5.clear();
    trk_chi20_layer_6.clear();
    // Per-layer chi2 values (iteration 1)
    trk_chi21_layer_0.clear();
    trk_chi21_layer_1.clear();
    trk_chi21_layer_2.clear();
    trk_chi21_layer_3.clear();
    trk_chi21_layer_4.clear();
    trk_chi21_layer_5.clear();
    trk_chi21_layer_6.clear();
    // Per-layer bendFE values
    trk_bendFE_layer_0.clear();
    trk_bendFE_layer_1.clear();
    trk_bendFE_layer_2.clear();
    trk_bendFE_layer_3.clear();
    trk_bendFE_layer_4.clear();
    trk_bendFE_layer_5.clear();
    trk_bendFE_layer_6.clear();
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::AnalyzerTQ);
