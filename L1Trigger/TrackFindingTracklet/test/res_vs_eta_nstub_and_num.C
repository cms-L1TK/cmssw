// -----------------------------------------------------------------------------------------------------------
// This code creates a resolution vs eta plot separated by number of stubs in a track (4-7 stubs).
// To run, it just needs the ntuple made by "L1TrackNtupleMaker.C".
// There is a legend and options to put labels on the plot.
//
// There is a section titled "Things to Change" just below. 
// Follow the intructions in that section to create and save the plot.
//
// Info: 
// tp        - all tracking particles
// matchtrk  - *L1 track* properties, for tracking particles matched to an L1 track
// trk       - all L1 tracks
// tp_nmatch - number of matched tracks for a tracking particle
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan Dec 2023
// -----------------------------------------------------------------------------------------------------------

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

#include "plotstyle.h"

void res_vs_eta_nstub_and_num() {

        SetPlotStyle();

        /////////////////////////////////////////// Things to Change //////////////////////////////////////////////////

        // Deciding whether you want to check the resolution of all tracks or only the matched tracks
        bool matchtrk = true; // if matchtrk is true, you check matched tracks. If false, you check all tracks
        // WARNING: SETTING MATCHTRK TO false ONLY WORKS FOR SINGLE MUON AND SINGLE ELECTRON. DO NOT USE FOR TTBAR PU200 or ANYTHING WITH MORE THAN 2 PARTICLES PER EVENT

        // File and directory to get Ntuple from 
        TString file = "SingleMuon_PU0_D88_oldkf_DRon_HOon_all"; // without ".root" at the end
        TString fileDir = "/eos/user/d/dabraham/tschuh_comparisons/"; // must have "/" at end of directory

        // Folder to save to
        TString saveDir = "/eos/user/d/dabraham/tschuh_comparisons/hybrid_vs_newkf_plots/"; // must have "/" at end of directory
        TString saveFile = "nstub_res_vs_eta_trk_oldkf_DRon_HOon.root"; // save file name in save directory

        // Labels
        TString DataSetLabel = "Single Muon PU=0"; // Enter Data set name
        TString label2 = "Old KF"; // Enter label with extra information (in my case it was the Kalman Filter versions)

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // read ntuples ---------------------------------------------------------------------------------------------------
        TChain* tree = new TChain("L1TrackNtuple/eventTree");
        tree->Add(fileDir + file + ".root");

        if (tree->GetEntries() == 0) {
        cout << "File doesn't exist or is empty, returning..."
        << endl;
        return;
        }

        // define leaves and branches ----------------------------------------------------------------------------------------------------------------
        vector<float>* tp_z0;
        vector<float>* matchtrk_z0;
        vector<float>* trk_z0;
        vector<float>* tp_eta;
        vector<float>* trk_eta;
        vector<int>* tp_nmatch;
        vector<int>* matchtrk_nstub;
        vector<int>* trk_nstub;
        vector<float>* tp_phi;
        vector<float>* matchtrk_phi;
        vector<float>* trk_phi;
        vector<int>* matchtrk_seed;
        vector<int>* trk_seed;
        
        TBranch* b_tp_z0;
        TBranch* b_matchtrk_z0;
        TBranch* b_trk_z0;
        TBranch* b_tp_eta;
        TBranch* b_trk_eta;
        TBranch* b_tp_nmatch;
        TBranch* b_matchtrk_nstub;
        TBranch* b_trk_nstub;
        TBranch* b_tp_phi;
        TBranch* b_matchtrk_phi;
        TBranch* b_trk_phi;
        TBranch* b_matchtrk_seed;
        TBranch* b_trk_seed;

        tp_z0 = 0;
        matchtrk_z0 = 0;
        trk_z0 = 0;
        tp_phi = 0;
        matchtrk_phi = 0;
        trk_phi = 0;
        tp_eta = 0;
        trk_eta = 0;
        tp_nmatch = 0;
        matchtrk_nstub = 0;
        trk_nstub = 0;
        matchtrk_seed = 0;
        trk_seed = 0;

        tree->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
        tree->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
        tree->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
        tree->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
        tree->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
        tree->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
        tree->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
        tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
        tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
        tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
        tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
        tree->SetBranchAddress("matchtrk_seed", &matchtrk_seed, &b_matchtrk_seed);
        tree->SetBranchAddress("trk_seed", &trk_seed, &b_trk_seed);


        vector<vector<float>> resEta;
        vector<float> tempResVec;
        const int nETARANGE = 25;
        TString etarange[nETARANGE] = {"0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9",
                                 "1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8",
                                 "1.9", "2.0", "2.1", "2.2", "2.3", "2.4", "2.5"};
        
        int nevt = tree->GetEntries();

        vector<int> tempNstubVec;
        vector<vector<int>> nstubVec;

        float eta_max = 2.5;
        TCanvas c;
        TCanvas c2;
        TH1F* h_resEtaNstub[8];
        for (int i = 0; i<8; i++) {
                h_resEtaNstub[i] =
                        new TH1F("nmatchVsEta_" + i, ";Tracking particle |#eta|; z0 Resolution", nETARANGE, 0, eta_max);
        }

        vector<float> tempBinContent;
        int tempBinContentCounter = 0;
        float quantile68;
        float maxRes = 0;

        cout << "sup before events" << endl; 
        // event loop -------------------------------------------------------------------------------------------
        for (int im = 0; im < nETARANGE; im++) {
                for (int iNstub = 4; iNstub <= 7; iNstub++) {
                        //cout << "------------------------------------------------" << endl;
                        //cout << "eta Range: " << im << "   Nstub: " << iNstub << endl;
                        tempResVec.clear();
                        for (int i = 0; i < nevt; i++) { 
                                tree->GetEntry(i, 0);

                                // if you want to compare the matchtrk vs tp resolution (as opposed to trk)
                                if (matchtrk) {
                                        // tracking particle loop -----------------------------------------------------------------------
                                        for (int it = 0; it < (int)tp_eta->size(); it++) {
                                                //cout << "tp_eta: " << tp_eta->at(it) << "   matchtrk_nstub: " << matchtrk_nstub->at(it);

                                                if ((std::abs(tp_eta->at(it)) > (float)im * 0.1) && (std::abs(tp_eta->at(it)) < (float)(im + 1) * 0.1)) {

                                                        if (matchtrk_nstub->at(it) == iNstub){
                                                                //cout << "   PASSED";
                                                                tempResVec.push_back(std::abs(tp_z0->at(it)-matchtrk_z0->at(it)));
                                                        }
                                                }
                                                //cout << endl;
                                        }
                                }
                                // If you want to compare the trk vs tp resolution (as opposed to matchtrk): ONLY WORKS FOR SINGLE MUON OR SINGLE ELECTRON, DO NOT USE FOR TTBAR
                                else {
                                        // trk loop -----------------------------------------------------------------------
                                        for (int itrk = 0; itrk < (int)trk_eta->size(); itrk++) {
                                                //cout << "tp_eta: " << tp_eta->at(it) << "   matchtrk_nstub: " << matchtrk_nstub->at(it);

                                                if ((std::abs(trk_eta->at(itrk)) > (float)im * 0.1) && (std::abs(trk_eta->at(itrk)) < (float)(im + 1) * 0.1)) {

                                                        if (trk_nstub->at(itrk) == iNstub){
                                                                //cout << "   PASSED";
                                                                if (std::abs(trk_phi->at(itrk)-tp_phi->at(0)) < 1) {
                                                                        tempResVec.push_back(std::abs(tp_z0->at(0)-trk_z0->at(itrk)));
                                                                }
                                                                else {
                                                                        tempResVec.push_back(std::abs(tp_z0->at(1)-trk_z0->at(itrk)));
                                                                }
                                                        }
                                                }
                                                //cout << endl;
                                        }
                                }
                        }
                        //tempBinContent.push_back(std::accumulate(tempResVec.begin(), tempResVec.end(),0) / float(tempResVec.size()));
                        std::sort(tempResVec.begin(), tempResVec.end());//Sorting the vector
                        //cout << "before if" << endl;
                        if (tempResVec.empty()) {
                                quantile68 = 0;
                                //cout << "after quantile68 if" << endl;
                        }
                        else {
                                quantile68 = tempResVec[ceil(tempResVec.size() * 0.68)];
                                //cout << "after quantile68 else" << endl;
                        }
                        
                        if (quantile68 > maxRes) {
                                maxRes = quantile68;
                        }
                        //tempBinContent.push_back(quantile68);
                        //cout << "bedore h_..." << endl;
                        h_resEtaNstub[iNstub]->SetBinContent(im+1, quantile68);
                        //cout << "after h_..." << endl;
                        ////////////// DEBUG CODE ///////////////
                        // cout << "tempBinContent: " << quantile68 << endl;
                        // cout << "tempResVec Size: " << tempResVec.size() << endl;
                        // cout << "tempResVec Content: ";
                        for(float resVecContent : tempResVec) {
                                //cout << resVecContent << "   ";
                        }
                        // cout << endl;
                        // cout << "quantile68: " << quantile68 << endl;
                        /////////////////////////////////////////
                        //tempBinContentCounter++;
                }
        }
        cout << "yo after" << endl;

        // Nstub Hist for each eta value
        TH1F* h_numEtaNstub[8];
        for (int i = 4; i<8; i++) {
                h_numEtaNstub[i] =
                        new TH1F("numVsEtaNstub_" + i, ";Tracking particle |#eta|; # of tracks", nETARANGE, 0, eta_max);
        }


        for (int i = 0; i < nevt; i++) { 
        tree->GetEntry(i, 0);
                if (matchtrk) {
                        for (int it = 0; it < (int)tp_eta->size(); it++) {
                                for (int iNstub = 4; iNstub <= 7; iNstub++) {
                                        if (matchtrk_nstub->at(it) == iNstub) {
                                                h_numEtaNstub[iNstub]->Fill(std::abs(tp_eta->at(it)));
                                        }
                                }
                        }
                }
                else {
                        for (int itrk = 0; itrk < (int)trk_eta->size(); itrk++) {
                                for (int iNstub = 4; iNstub <= 7; iNstub++) {
                                        if (trk_nstub->at(itrk) == iNstub) {
                                                if (std::abs(trk_phi->at(itrk)-tp_phi->at(0)) < 1) {
                                                        h_numEtaNstub[iNstub]->Fill(std::abs(tp_eta->at(0)));
                                                }
                                                else {
                                                        h_numEtaNstub[iNstub]->Fill(std::abs(tp_eta->at(1)));
                                                }
                                        }
                                }
                        }
                }
        }
                                



        //vector <int> colors = {kBlack, kRed-3, kOrange-3, kBlue, kGreen+3, kCyan-2, kRed-5, kMagenta-3};
        vector <int> colors = {0, 0, 0, 0, kBlack, kGreen+3, kOrange-3, kBlue, kCyan-2, kRed-5, kMagenta-3};
        
        c.cd();

        h_resEtaNstub[4]->Draw("p"); 
        h_resEtaNstub[4]->SetMarkerColor(colors[4]);
        h_resEtaNstub[4]->SetMaximum(maxRes*1.1);
        for (int i = 5; i<8; i++) {
                h_resEtaNstub[i]->Draw("p same");
                h_resEtaNstub[i]->SetMarkerColor(colors[i]);
        }

        c2.cd();
        int maxNum = h_numEtaNstub[4]->GetMaximum();
        h_numEtaNstub[4]->Draw("p"); 
        h_numEtaNstub[4]->SetMarkerColor(colors[4]);
        for (int i = 5; i<8; i++) {
                h_numEtaNstub[i]->Draw("p same");
                h_numEtaNstub[i]->SetMarkerColor(colors[i]);
                if (h_numEtaNstub[i]->GetMaximum() > maxNum) {
                        maxNum = h_numEtaNstub[i]->GetMaximum();
                }
        }
        h_numEtaNstub[4]->SetMaximum(maxNum);


        c.cd();
        // make legend
        TLegend* leg = new TLegend(0.2, 0.63, 0.33, 0.93);
        TString nstubNum;
        for (int i = 4; i<8; i++) {
                nstubNum = "Nstub " + to_string(i);
                leg->AddEntry(h_resEtaNstub[i], nstubNum, "p");
        }

        char ctxt[500];
        char ctxt2[500];
        char ctxt3[500];
        sprintf(ctxt, DataSetLabel); // Add label saying 
        mySmallText(0.35, 0.88, 1, ctxt); // which data set it is
        sprintf(ctxt2, label2);
        mySmallText(0.35, 0.83, 1, ctxt2);
        if (matchtrk) {
                sprintf(ctxt3, "matchtrk");
        }
        else {
                sprintf(ctxt3, "trk");
        }
        mySmallText(0.35, 0.78, 1, ctxt3);

        leg->Draw();
        
        c.SaveAs(saveDir + saveFile);

        c2.cd();

        TLegend* leg2 = new TLegend(0.2, 0.63, 0.33, 0.93);
        for (int i = 4; i<8; i++) {
                nstubNum = "Nstub " + to_string(i);
                leg2->AddEntry(h_numEtaNstub[i], nstubNum, "p");
        }

        sprintf(ctxt, DataSetLabel); // Add label saying 
        mySmallText(0.35, 0.88, 1, ctxt); // which data set it is
        sprintf(ctxt2, label2);
        mySmallText(0.35, 0.83, 1, ctxt2);
        if (matchtrk) {
                sprintf(ctxt3, "matchtrk");
        }
        else {
                sprintf(ctxt3, "trk");
        }
        mySmallText(0.35, 0.78, 1, ctxt3);

        leg2->Draw();

        c2.SaveAs(saveDir + "num_" + saveFile);

}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.050;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}