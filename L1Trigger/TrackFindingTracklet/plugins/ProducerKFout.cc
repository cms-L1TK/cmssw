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
    vector<int> cotBins;
    // Factors used to convert between phi/Z at a radius T, modified to include scale conversions needed in calculation
    double modChosenRofZ = 0;
    double modChosenRofPhi = 0;
    // Corrections to Phi depending on which phi sector
    int BaseSectorCorr = 0;
    int UnsignedBaseSector = 0;
    // Redefine these here to avoid impossible to read code later
    // Widths used for each part of data stored in the KF track and KF stubs
    int wTm = 0;
    int wTSp = 0;
    int wTSe = 0;
    int wTphi = 0;
    int wTinvR = 0;
    int wTcot = 0;
    int wTz = 0;

    int wSr = 0;
    int wSphi = 0;
    int wSz  = 0;
    int wSdPhi = 0;
    int wSdZ = 0;

    int chi2rphiConv = 0;
    int chi2rzConv = 0;
    // Bins for dPhi/dZ use to access weight LUT below
    vector<int> dPhiBins;
    vector<int> dZBins;
    // LUT for weighting functions for chi2 calculation
    vector<int> v0Bins;
    vector<int> v1Bins;
    // Bins for final Chi2 Packing
    vector<double> chi2rphiBins;
    vector<double> chi2rzBins;

    int chi2ScaleFactor;

    int maxTracksPerEvent;
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

    for (int i = 0; i < setup_->numSectorsEta(); i++){
      cotBins.push_back(setup_->sectorCot(i)*(1/dataFormats_->base(Variable::cot,Process::kf)));
    };
    modChosenRofZ = setup_->chosenRofZ() * (dataFormats_->base(Variable::cot,Process::kf) / dataFormats_->base(Variable::zT,Process::kf));
    modChosenRofPhi = setup_->hybridChosenRofPhi() * dataFormats_->base(Variable::inv2R,Process::kf) / dataFormats_->base(Variable::phiT,Process::kf);

    UnsignedBaseSector = (M_PI / (double)(setup_->numRegions() * setup_->numSectorsPhi()) ) / (dataFormats_->base(Variable::phiT,Process::kf)); 

    wTm    = dataFormats_->width(Variable::match,Process::kf);
    wTSp   = dataFormats_->width(Variable::sectorPhi,Process::kf);
    wTSe   = dataFormats_->width(Variable::sectorEta,Process::kf);
    wTphi  = dataFormats_->width(Variable::phiT,Process::kf);
    wTinvR = dataFormats_->width(Variable::inv2R,Process::kf);
    wTcot  = dataFormats_->width(Variable::cot,Process::kf);
    wTz    = dataFormats_->width(Variable::zT,Process::kf);

    wSr    = dataFormats_->width(Variable::r,Process::kf);
    wSphi  = dataFormats_->width(Variable::phi,Process::kf);
    wSz    = dataFormats_->width(Variable::z,Process::kf);
    wSdPhi = dataFormats_->width(Variable::dPhi,Process::kf);
    wSdZ   = dataFormats_->width(Variable::dZ,Process::kf);

    dPhiBins = setup_->kfoutdPhiBins();
    dZBins   = setup_->kfoutdZBins();
    v0Bins   = setup_->kfoutv0Bins();
    v1Bins   = setup_->kfoutv1Bins();

    chi2rphiBins = setup_->kfoutchi2rphiBins();
    chi2rzBins   = setup_->kfoutchi2rzBins();

    chi2rphiConv = setup_->kfoutchi2rphiConv();
    chi2rzConv   = setup_->kfoutchi2rzConv();

    chi2ScaleFactor = setup_->kfoutchi2ScaleFactor();

    maxTracksPerEvent = setup_->kfoutmaxTracksPerEvent();

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

          std::vector<TTBV> Svalids;
          std::vector<TTBV> Srs;
          std::vector<TTBV> Sphis;
          std::vector<TTBV> Szs;
          std::vector<TTBV> SdPhis;
          std::vector<TTBV> SdZs;

          // Create bit vector (gives us all the slicing etc. )
          TTBV InputLinkTrack(track.second);
          // Slice Constructor (const TTBV& ttBV, int begin, int end = 0, bool twos =  false)
          // Unpack KF track 
          TTBV Tvalid(     InputLinkTrack , wTz + wTcot + wTinvR + wTphi + wTSe + wTSp + wTm + 1, wTz + wTcot + wTinvR + wTphi + wTSe + wTSp + wTm , false );
          TTBV Tmatch(     InputLinkTrack , wTz + wTcot + wTinvR + wTphi + wTSe + wTSp + wTm    , wTz + wTcot + wTinvR + wTphi + wTSe + wTSp , false );
          TTBV TsectorPhi( InputLinkTrack , wTz + wTcot + wTinvR + wTphi + wTSe + wTSp          , wTz + wTcot + wTinvR + wTphi + wTSe , false );
          TTBV TsectorEta( InputLinkTrack , wTz + wTcot + wTinvR + wTphi + wTSe                 , wTz + wTcot + wTinvR + wTphi , false );
          TTBV TphiT(      InputLinkTrack , wTz + wTcot + wTinvR + wTphi                        , wTz + wTcot + wTinvR , true );
          TTBV TinvR(      InputLinkTrack , wTz + wTcot + wTinvR                                , wTz + wTcot , true );
          TTBV Tcot(       InputLinkTrack , wTz + wTcot                                         , wTz , true );
          TTBV TzT(        InputLinkTrack , wTz                                                 , 0 , true );

          int temp_z0 = TzT.val() - ((Tcot.val() * modChosenRofZ));

          // Correction to Phi calcuation depending if +ve/-ve phi sector
          if (TsectorPhi.val() == 0) {
            BaseSectorCorr = -UnsignedBaseSector;
          }
          else{
            BaseSectorCorr = UnsignedBaseSector;
          }

          int temp_phi0 = TphiT.val() - ((TinvR.val()) * modChosenRofPhi) + BaseSectorCorr;
          
          int temp_tanL = cotBins[TsectorEta.val()] + Tcot.val();
        
          TTBV HitPattern(0,setup_->numLayers());

          int tempchi2rphi = 0;
          int tempchi2rz   = 0;

          for (int iStub = 0; iStub < setup_->numLayers() - 1; iStub++ ){
            auto stub = streamsStubs[setup_->numLayers()*iLink+iStub].at(iTrack);
            TTBV InputLinkStub(stub.second);

            // Unpack and store stub bits  
            TTBV Svalid( InputLinkStub , wSdZ + wSdPhi + wSz + wSphi + wSr + 1 , wSdZ + wSdPhi + wSz + wSphi + wSr , false );
            TTBV Sr(     InputLinkStub , wSdZ + wSdPhi + wSz + wSphi + wSr     , wSdZ + wSdPhi + wSz+ wSphi , true );
            TTBV Sphi(   InputLinkStub , wSdZ + wSdPhi + wSz + wSphi           , wSdZ + wSdPhi + wSz , true );
            TTBV Sz(     InputLinkStub , wSdZ + wSdPhi + wSz                   , wSdZ + wSdPhi , true );
            TTBV SdPhi(  InputLinkStub , wSdZ + wSdPhi                         , wSdZ , true );
            TTBV SdZ(    InputLinkStub , wSdZ                                  , 0 , true );

            if (Svalid.val() == 1){
              HitPattern.set(iStub);
              Svalids.push_back(Svalid);
              Srs.push_back(Sr);
              Sphis.push_back(Sphi);
              Szs.push_back(Sz);
              SdPhis.push_back(SdPhi);
              SdZs.push_back(SdZ);

              int phiSquared = Sphi.val() * Sphi.val();
              int zSquared   = Sz.val() * Sz.val();

              int tempv0 = (v0Bins[digitise(dPhiBins, SdPhi.val())]);
              int tempv1 = (v1Bins[digitise(dZBins  , SdZ.val())]);

              int tempRphi = phiSquared * tempv0;
              int tempRz   = zSquared * tempv1;

              tempchi2rphi += tempRphi / chi2ScaleFactor;
              tempchi2rz   += tempRz / chi2ScaleFactor;
            }
          }
          // TODO extract TTTrack bit widths from TTTrack word pending update to the TTTrack_word class
          TTBV TrackValid(Tvalid,1,false);
          TTBV extraMVA(0,6,false);
          TTBV TQMVA(0,3,false);
          TTBV BendChi2(0,3,false); 
          TTBV Chi2rphi(digitise(chi2rphiBins,(double)tempchi2rphi,(double)chi2rphiConv),4,false);
          TTBV Chi2rz(digitise(chi2rzBins,(double)tempchi2rz,(double)chi2rzConv),4,false);
          TTBV D0(0,13,false);
          TTBV z0(temp_z0 ,12,true);
          TTBV TanL(temp_tanL,16,true);
          TTBV phi0(temp_phi0 ,12,true);
          TTBV InvR(-1*TinvR.val() ,16,true );

                              // 6      +   3   +   7        +  3       + 13
          TTBV PartialTrack1((extraMVA + TQMVA + HitPattern + BendChi2 + D0),32,false);
                             // 4       + 12    + 16
          TTBV PartialTrack2((Chi2rz   + z0    + TanL),32,false);
                             // 4       + 12    +  15  +    1
          TTBV PartialTrack3((Chi2rphi + phi0  + InvR + TrackValid),32,false);

          // Sort Tracks based on eta
          if (iLink % 2 == 0){
            if (TsectorEta.val() < 8){
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
            if (TsectorEta.val() < 8){
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
        std::cout << SortedPartialTracks[iLink].size() << " " << OutputStreamsTracks[iLink].size() << std::endl;

        if (numLinkTracks > 0){
          if ((numLinkTracks % 2 != 0)) { 
            SortedPartialTracks[iLink].push_back(NullBitTrack);
            OutputStreamsTracks[iLink].emplace_back(OutputStreamsTracks[iLink][numLinkTracks]);
            numLinkTracks++;
            } //If there is an odd number of tracks 
          std::cout << SortedPartialTracks[iLink].size() << " " << OutputStreamsTracks[iLink].size() << std::endl;
          for (int iTrack = 0; iTrack < (int)(SortedPartialTracks[iLink].size()); iTrack++ ){  
            std::cout << iLink << " " << iTrack << " PT - " << SortedPartialTracks[iLink][iTrack].resize(32).bs().to_string() << std::endl;
            if (iTrack % 2 == 1){
              if ((int)iTrack/3 <= maxTracksPerEvent){
                accepted[iLink].emplace_back(std::make_pair(ttTrackRefMap.at(OutputStreamsTracks[iLink][(int)(iTrack-1)/3].first),(SortedPartialTracks[iLink][iTrack].slice(32) + SortedPartialTracks[iLink][iTrack-1].slice(32)).bs()));
              }
              else{
                lost[iLink].emplace_back(std::make_pair(ttTrackRefMap.at(OutputStreamsTracks[iLink][(int)(iTrack-1)/3].first),(SortedPartialTracks[iLink][iTrack].slice(32) + SortedPartialTracks[iLink][iTrack-1].slice(32)).bs()));
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
