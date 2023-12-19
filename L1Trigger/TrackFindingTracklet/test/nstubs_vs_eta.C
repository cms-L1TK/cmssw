////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This code plots average number or RMS of stubs in tracks against eta 
// (can toggle RMS or mean by commenting out in code lines 139-149)
// To run, it just needs an Ntuple made by L1TrackNtupleMaker.C
// There are options at the beginning to enter the filename.
//
// By David Abrahamyan Dec 2023
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "plotstyle.h"

void nstubs_vs_eta() {

        SetPlotStyle();

        ///////////////////////////////////////// Things to Change //////////////////////////////////////////////////////

        // File and directory to get Ntuple from
        TString file = "[ENTER FILENAME]"; // without ".root" at the end
        TString fileDir = "[ENTER DIRECTORY]"; // must have "/" at end of directory

        // Folder to save to
        TString saveDir = "[ENTER DIRECTORY]"; // must have "/" at end of directory
        TString saveFile = "[ENTER FILENAME]"; // save file name in save directory

        // Do you want to do stub mean or stub RMS?
        bool meanValue = true; // set to true if mean. false if RMS.

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        vector<int>* tp_nmatch;
        vector<int>* matchtrk_nstub;
        vector<int>* trk_nstub;
        vector<float>* tp_phi;
        vector<float>* matchtrk_phi;
        vector<float>* trk_phi;
        
        TBranch* b_tp_z0;
        TBranch* b_matchtrk_z0;
        TBranch* b_trk_z0;
        TBranch* b_tp_eta;
        TBranch* b_tp_nmatch;
        TBranch* b_matchtrk_nstub;
        TBranch* b_trk_nstub;
        TBranch* b_tp_phi;
        TBranch* b_matchtrk_phi;
        TBranch* b_trk_phi;

        tp_z0 = 0;
        matchtrk_z0 = 0;
        trk_z0 = 0;
        tp_phi = 0;
        matchtrk_phi = 0;
        trk_phi = 0;
        tp_eta = 0;
        tp_nmatch = 0;
        matchtrk_nstub = 0;
        trk_nstub = 0;

        tree->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
        tree->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
        tree->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
        tree->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
        tree->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
        tree->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
        tree->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
        tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
        tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
        tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);


        vector<vector<float>> nstubEta;
        vector<float> tempNstubVec;
        const int nETARANGE = 25;
        TString etarange[nETARANGE] = {"0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9",
                                 "1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8",
                                 "1.9", "2.0", "2.1", "2.2", "2.3", "2.4", "2.5"};
        
        int nevt = tree->GetEntries();

        // for 25 bins in hist
        
        cout << "sup before events" << endl; 
        // event loop -------------------------------------------------------------------------------------------
        for (int im = 0; im < nETARANGE; im++) {
                tempNstubVec.clear();
                for (int i = 0; i < nevt; i++) {
                        tree->GetEntry(i, 0);

                        // track loop -----------------------------------------------------------------------
                        for (int it = 0; it < (int)tp_eta->size(); it++) {
                                for (int itrk = 0; itrk < (int)trk_nstub->size(); itrk++) {
                                        if (std::abs(trk_phi->at(itrk) - tp_phi->at(it)) > 1)
                                                continue;

                                        if ((std::abs(tp_eta->at(it)) > (float)im * 0.1) && (std::abs(tp_eta->at(it)) < (float)(im + 1) * 0.1)) {
                                                tempNstubVec.push_back(trk_nstub->at(itrk));
                                                //nmatch[im].push_back(tp_nmatch->at(it)); // fill nstubEta[i] with all nmatch values for a range of eta values (i.e. 1.3-1.4)
                                        }  
                                }
                        }
                }
                nstubEta.push_back(tempNstubVec);
        }
        cout << "yo after" << endl;

        // Debug Code ///////////////////////////////////////////////////////
        // cout << "Size of nstubEta: " << nstubEta.size() << endl;
        // for (int i = 0; i < nETARANGE; i++) {
        //         cout << "Size of nstubEta[" << i << "]: " << nstubEta[i].size() << endl;
        // }


        /////////////////////////////////////////////////////////////////////
        
        

        float binContentList[nETARANGE];
        float eta_max = 2.5;

        TCanvas c;

        if (meanValue) {
                TH1F* h_nstubVsEta =
                new TH1F("nstubVsEta", ";Tracking particle |#eta|; trk n_{stub} Mean", nETARANGE, 0, eta_max);
        }
        else {
                TH1F* h_nstubVsEta =
                new TH1F("nstubVsEta", ";Tracking particle |#eta|; trk n_{stub} RMS", nETARANGE, 0, eta_max);
        }

        

        // take averages of each nstubEta[i] list
        float squaredSum;
        for (int i = 0; i < nETARANGE; i++) {

                if (meanValue) {
                        // Mean
                        binContentList[i] = std::accumulate(nstubEta[i].begin(), nstubEta[i].end(),0) / float(nstubEta[i].size());
                        h_nstubVsEta->SetBinContent(i+1, binContentList[i]);
                }
                else {
                        //RMS
                        squaredSum = 0;
                        for (int ii = 0; ii < nstubEta[i].size(); ii++) {
                                squaredSum = squaredSum + pow((nstubEta[i][ii]), 2);
                        }
                        binContentList[i] = sqrt(squaredSum / float(nstubEta[i].size())); 
                        h_nstubVsEta->SetBinContent(i+1, binContentList[i]);
                }

                
        }

        gPad->SetGridx();
        gPad->SetGridy();

        
        h_nstubVsEta->SetMinimum(0.9);
        h_nstubVsEta->SetMaximum(6.1);
        h_nstubVsEta->Draw("p");

        c.SaveAs(saveDir + saveFile);
}



// double x[100], y[100];
//    int n = 20;
//    for (int i=0;i<n;i++) {
//      x[i] = i*0.1;
//      y[i] = 10*sin(x[i]+0.2);
//    }
//    auto g = new TGraph(n,x,y);
//    g->SetTitle("Graph title;X title;Y title");
//    g->Draw("AC*");