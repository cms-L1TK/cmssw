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

#include <string>
#include <numeric>

using namespace std;
using namespace edm;
using namespace trackerTFP;
using namespace tt;

namespace trackFindingTracklet {

  /*! \class  trackFindingTracklet::ProducerKFout
   *  \brief  Converts KF output into TFP output
   *  \author Thomas Schuh
   *  \date   2021, Aug
   */
  class ProducerKFout : public stream::EDProducer<> {
  public:
    explicit ProducerKFout(const ParameterSet&);
    ~ProducerKFout() override {}
    template<typename T>
    int digitise(const vector<T> Bins, T Value, T factor = 1 );

  private:
    void beginRun(const Run&, const EventSetup&) override;
    void produce(Event&, const EventSetup&) override;
    void endJob() {}

    // ED input token of kf stubs
    EDGetTokenT<StreamsStub> edGetTokenStubs_;
    // ED input token of kf tracks
    EDGetTokenT<StreamsTrack> edGetTokenTracks_;
    // ED input token of kf input to kf output TTTrack map
    EDGetTokenT<TTTrackRefMap> edGetTokenTTTrackRefMap_;
    // ED output token for accepted kfout tracks
    EDPutTokenT<StreamsTrack> edPutTokenAccepted_;
    // ED output token for truncated kfout tracks
    EDPutTokenT<StreamsTrack> edPutTokenLost_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // DataFormats token
    ESGetToken<DataFormats, DataFormatsRcd> esGetTokenDataFormats_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_;
    // helper class to extract structured data from TTDTC::Frames
    const DataFormats* dataFormats_;
    // Cot bins used to convert from cot to TanL
    vector<double> cotBins_;
    // Factors used to convert between phi/Z at a radius T, modified to include scale conversions needed in calculation
    double modChosenRofZ_ = 0;
    double modChosenRofPhi_ = 0;
    // Corrections to Phi depending on which phi sector
    double BaseSectorCorr_ = 0;
    double UnsignedBaseSector_ = 0;
    // Bins for dPhi/dZ use to access weight LUT below
    vector<double> dPhiBins_;
    vector<double> dZBins_;
    // LUT for weighting functions for chi2 calculation
    vector<double> v0Bins_;
    vector<double> v1Bins_;
    // Bins for final Chi2 Packing
    vector<double> chi2rphiBins_;
    vector<double> chi2rzBins_;

    int maxTracksPerEvent_;
  };

  ProducerKFout::ProducerKFout(const ParameterSet& iConfig) :
    iConfig_(iConfig)
  {
    const string& labelKF = iConfig.getParameter<string>("LabelKF");
    const string& labelAS = iConfig.getParameter<string>("LabelAS");
    const string& branchStubs = iConfig.getParameter<string>("BranchAcceptedStubs");
    const string& branchTracks = iConfig.getParameter<string>("BranchAcceptedTracks");
    const string& branchLost = iConfig.getParameter<string>("BranchLostTracks");
    // book in- and output ED products
    edGetTokenStubs_ = consumes<StreamsStub>(InputTag(labelKF, branchStubs));
    edGetTokenTracks_ = consumes<StreamsTrack>(InputTag(labelKF, branchTracks));
    edGetTokenTTTrackRefMap_ = consumes<TTTrackRefMap>(InputTag(labelAS, branchTracks));
    edPutTokenAccepted_ = produces<StreamsTrack>(branchTracks);
    edPutTokenLost_ = produces<StreamsTrack>(branchLost);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    // initial ES products
    setup_ = nullptr;
    dataFormats_ = nullptr;
  }

  void ProducerKFout::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    if (!setup_->configurationSupported())
      return;
    // check process history if desired
    if (iConfig_.getParameter<bool>("CheckHistory"))
      setup_->checkHistory(iRun.processHistory());
    // helper class to extract structured data from TTDTC::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);

    for (int i = 0; i < setup_->numSectorsEta(); i++) cotBins_.push_back(setup_->sectorCot(i));

    modChosenRofZ_ = setup_->chosenRofZ();
    modChosenRofPhi_ = setup_->hybridChosenRofPhi();

    UnsignedBaseSector_ = (M_PI / (double)(setup_->numRegions() * setup_->numSectorsPhi()) );
    
    // Convert Integer bins to doubles for internal calculation
    for(size_t i = 0; i < setup_->kfoutdPhiBins().size(); i++)
      dPhiBins_.push_back((double)setup_->kfoutdPhiBins()[i] * dataFormats_->base(Variable::dPhi, Process::kfin));
    for(size_t i = 0; i < setup_->kfoutdZBins().size(); i++) 
      dZBins_.push_back((double)setup_->kfoutdZBins()[i] * dataFormats_->base(Variable::dZ, Process::kfin));
    for(size_t i = 0; i < setup_->kfoutv0Bins().size(); i++) 
      v0Bins_.push_back(setup_->kfoutv0Bins()[i]*dataFormats_->base(Variable::dPhi, Process::kfin));
    for(size_t i = 0; i < setup_->kfoutv1Bins().size(); i++) 
      v1Bins_.push_back(setup_->kfoutv1Bins()[i]*dataFormats_->base(Variable::dZ, Process::kfin));
    for(size_t i = 0; i < setup_->kfoutchi2rphiBins().size(); i++)
      chi2rphiBins_.push_back((double)setup_->kfoutchi2rphiBins()[i] *setup_->kfoutchi2rphiConv() * pow(dataFormats_->base(Variable::dPhi, Process::kfin),3));
    for(size_t i = 0; i < setup_->kfoutchi2rzBins().size(); i++)
      chi2rzBins_.push_back((double)setup_->kfoutchi2rzBins()[i] *setup_->kfoutchi2rzConv() * pow(dataFormats_->base(Variable::dZ, Process::kfin),3));
    maxTracksPerEvent_ = setup_->numFramesIO() * setup_->clockRatio();
  }


  template<typename T>
  int ProducerKFout::digitise(const vector<T> Bins, T Value, T factor ) {
    for (int i = 0; i < (int)Bins.size(); i++){
      if (Value > Bins[i]*factor && Value <= Bins[i+1]*factor) {return i;}
    }
    return -1;
  }

  void ProducerKFout::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty KFout product
    StreamsTrack accepted(setup_->numRegions() * setup_->tfpNumChannel());
    StreamsTrack lost(setup_->numRegions() * setup_->tfpNumChannel());
    // read in KF Product and produce KFout product
    if (setup_->configurationSupported()) {
      Handle<StreamsStub> handleStubs;
      iEvent.getByToken<StreamsStub>(edGetTokenStubs_, handleStubs);
      const StreamsStub& streamsStubs = *handleStubs.product();
      Handle<StreamsTrack> handleTracks;
      iEvent.getByToken<StreamsTrack>(edGetTokenTracks_, handleTracks);
      const StreamsTrack& streamsTracks = *handleTracks.product();
      Handle<TTTrackRefMap> handleTTTrackRefMap;
      iEvent.getByToken<TTTrackRefMap>(edGetTokenTTTrackRefMap_, handleTTTrackRefMap);
      const TTTrackRefMap& ttTrackRefMap = *handleTTTrackRefMap.product();
      // 18 Output Links (First Vector) each has a vector of tracks per event (second vector) each track is 3 32 bit TTBV partial tracks 
      vector<vector<TTBV>> SortedPartialTracks(setup_->numRegions() * setup_->tfpNumChannel(),vector<TTBV>(0));   
      StreamsTrack OutputStreamsTracks(setup_->numRegions() * setup_->tfpNumChannel());

      for (int iLink = 0; iLink < (int)streamsTracks.size(); iLink++ ){

        for (int iTrack = 0; iTrack < (int)streamsTracks[iLink].size(); iTrack++ ){
          auto track = streamsTracks[iLink].at(iTrack);
          TrackKF InTrack(track,dataFormats_);
          
          double temp_z0 = InTrack.zT() - ((InTrack.cot() * setup_->chosenRofZ()));

          // Correction to Phi calcuation depending if +ve/-ve phi sector
          if (InTrack.sectorPhi() == 0) BaseSectorCorr_ = -UnsignedBaseSector_;
          else BaseSectorCorr_ = UnsignedBaseSector_;

          double temp_phi0 = InTrack.phiT() - ((InTrack.inv2R()) * setup_->hybridChosenRofPhi()) + BaseSectorCorr_;
          
          double temp_tanL = cotBins_[InTrack.sectorEta()] + InTrack.cot();
        
          TTBV HitPattern(0,setup_->numLayers());

          double tempchi2rphi = 0;
          double tempchi2rz   = 0;

          for (int iStub = 0; iStub < setup_->numLayers() - 1; iStub++ ){
            auto stub = streamsStubs[setup_->numLayers()*iLink+iStub].at(iTrack);
            StubKF InStub(stub,dataFormats_,iStub);

            if (stub.first.isNonnull()){
              HitPattern.set(iStub);

              double phiSquared = InStub.phi() * InStub.phi();
              double zSquared   = InStub.z() * InStub.z();

              double tempv0 = (v0Bins_[digitise(dPhiBins_, InStub.dPhi())]);
              double tempv1 = (v1Bins_[digitise(dZBins_  , InStub.dZ())]);

              double tempRphi = phiSquared * tempv0;
              double tempRz   = zSquared * tempv1;

              tempchi2rphi += tempRphi;
              tempchi2rz   += tempRz;
            }
          }
          // TODO extract TTTrack bit widths from TTTrack word pending update to the TTTrack_word class
          TTBV TrackValid(1,1,false);
          TTBV extraMVA(0,6,false);
          TTBV TQMVA(0,3,false);
          TTBV BendChi2(0,3,false); 
          TTBV Chi2rphi(digitise(chi2rphiBins_,tempchi2rphi),4,false);
          TTBV Chi2rz(digitise(chi2rzBins_,tempchi2rz),4,false);
          TTBV D0(0,13,false);
          TTBV z0(temp_z0 ,dataFormats_->base(Variable::zT,Process::kf) ,12,true);
          TTBV TanL(temp_tanL,dataFormats_->base(Variable::cot,Process::kf),16,true);
          TTBV phi0(temp_phi0,dataFormats_->base(Variable::phiT,Process::kf),12,true);
          TTBV InvR(-1*InTrack.inv2R(),dataFormats_->base(Variable::inv2R,Process::kf) ,16,true );

                              // 6      +   3   +   7        +  3       + 13
          TTBV PartialTrack1((extraMVA + TQMVA + HitPattern + BendChi2 + D0),32,false);
                             // 4       + 12    + 16
          TTBV PartialTrack2((Chi2rz   + z0    + TanL),32,false);
                             // 4       + 12    +  15  +    1
          TTBV PartialTrack3((Chi2rphi + phi0  + InvR + TrackValid),32,false);

          // Sort Tracks based on eta
          if (iLink % 2 == 0){
            if (InTrack.sectorEta() < (int)(setup_->numSectorsEta()/2)){
                SortedPartialTracks[iLink].push_back(PartialTrack1);
                SortedPartialTracks[iLink].push_back(PartialTrack2);
                SortedPartialTracks[iLink].push_back(PartialTrack3);
                OutputStreamsTracks[iLink].emplace_back(track);
              }
              else{
                SortedPartialTracks[iLink+1].push_back(PartialTrack1);
                SortedPartialTracks[iLink+1].push_back(PartialTrack2);
                SortedPartialTracks[iLink+1].push_back(PartialTrack3);
                OutputStreamsTracks[iLink+1].emplace_back(track);
              }
            }
          else{
            if (InTrack.sectorEta() < (int)(setup_->numSectorsEta()/2)){
              SortedPartialTracks[iLink-1].push_back(PartialTrack1);
              SortedPartialTracks[iLink-1].push_back(PartialTrack2);
              SortedPartialTracks[iLink-1].push_back(PartialTrack3);
              OutputStreamsTracks[iLink-1].emplace_back(track);
            }
            else{
              SortedPartialTracks[iLink].push_back(PartialTrack1);
              SortedPartialTracks[iLink].push_back(PartialTrack2);
              SortedPartialTracks[iLink].push_back(PartialTrack3);
              OutputStreamsTracks[iLink].emplace_back(track);
            }
 
          }
        } // Tracks
      } // Links

      // Fill products and match up tracks
      TTBV NullBitTrack(0,32,false);
      for (int iLink = 0; iLink < (int)OutputStreamsTracks.size(); iLink++ ){
        // Iterate through partial tracks
        int numLinkTracks = (int)OutputStreamsTracks[iLink].size();
        if (numLinkTracks > 0){
          if ((numLinkTracks % 2 != 0)) { //If there is an odd number of tracks 
            SortedPartialTracks[iLink].push_back(NullBitTrack);  //Pad out final set of bits
            OutputStreamsTracks[iLink].emplace_back(OutputStreamsTracks[iLink][numLinkTracks]); //Pad out with final repeated track
            numLinkTracks++;
            } //If there is an odd number of tracks 
          for (int iTrack = 0; iTrack < (int)(SortedPartialTracks[iLink].size()); iTrack++ ){  
            if (iTrack % 2 == 1){
              TTTrackRef TrackRef;
              for (auto &it : ttTrackRefMap) {   //Iterate throguh ttTrackRefMap to find TTTrackRef Key by a TTTrack Value
                if(it.second == OutputStreamsTracks[iLink][(int)(iTrack-1)/3].first) { 
                  TrackRef = it.first;
                } 
              }
              if ((int)iTrack/3 <= maxTracksPerEvent_){
                accepted[iLink].emplace_back(std::make_pair(TrackRef,(SortedPartialTracks[iLink][iTrack].slice(32) + SortedPartialTracks[iLink][iTrack-1].slice(32)).bs()));
              }
              else{
                lost[iLink].emplace_back(std::make_pair(TrackRef,(SortedPartialTracks[iLink][iTrack].slice(32) + SortedPartialTracks[iLink][iTrack-1].slice(32)).bs()));
              }
            }
          }
        }
      }

    } // Config Supported
    // store products
    iEvent.emplace(edPutTokenAccepted_, move(accepted));
    iEvent.emplace(edPutTokenLost_, move(lost));
  }
} // namespace trackFindingTracklet

DEFINE_FWK_MODULE(trackFindingTracklet::ProducerKFout);
