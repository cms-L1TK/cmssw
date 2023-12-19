////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This code plots stubs on the Z and R axes and then draws a line showing the tracking particles path and
// the reconstructed matchtrk path both based on track parameters.
// To use you need to have a specific event and tracking particle that you want to plot with. 
//
// ******************* : Currently, it chooses a trk based on the lowest eta residual and assigns that trk to matchtrk. 
// ****** NOTE ******* : So the matchtrk used is not the matchtrk made by the L1 Track algo with the current code.
// ******************* : That can be reverted by commenting out line 163
//
// ******************* : This code is written so that it can only be used for Single Particle events (i.e. Single Muon or Single Electron)
// ***** WARNING ***** : See lines 169-200 to see why it is limited (essentially, it picks stubs assuming only two particles fly through detector and they fly in opposite directions)
// ******************* : This can be fixed with some creativity probably
//
// It needs an Ntuple made by L1TrackNtupleMaker.C
// There are file and save settings at the start of the code in a section labelled "Things to Change"
//
//
// Code definitions: 
// tp        - all tracking particles
// matchtrk  - L1 Tracks that are matched to a tracking particle
// trk       - all L1 tracks
// tp_nmatch - number of matched tracks for a tracking particle
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

void ZRplot() {
        
        /////////////////////////////////// Things to Change /////////////////////////////////////////////////
        
        // File and directory to get Ntuple from 
        TString file = "[FILE NAME]"; // without ".root" at the end
        TString fileDir = "[ENTER DIRECTORY]"; // must have "/" at end of directory

        // save settings
        TString saveDir = "[SAVE DIRECTORY]"; // must have "/" at end of directory
        TString saveName = "ZRplot.pdf"; // save file name in save directory

        // choose event and tp that you want to look at
        int eventNum = 2064; // choose event number that contains the tp in question
        int tpNum = 0; // in the ntuple, tp's are listed in a vector, so this is their index in that array

        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // read ntuples ---------------------------------------------------------------------------------------------------
        TChain* tree = new TChain("L1TrackNtuple/eventTree");
        tree->Add(fileDir + file + ".root");

        if (tree->GetEntries() == 0) {
        cout << "File doesn't exist or is empty, returning..."
        << endl;
        return; 
        }

        // -----------------------------------------------------------------------------------------------
        // define leaves and branches (reading from Ntuple)

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
        vector<float>* trk_eta;
        vector<float>* matchtrk_eta;
        vector<float>* allstub_z;
        vector<float>* allstub_x;
        vector<float>* allstub_y;
        vector<float>* allstub_isPSmodule;
        vector<float>* allstub_isTiltedBarrel;
        vector<float>* allstub_isBarrel;
        
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
        TBranch* b_trk_eta;
        TBranch* b_matchtrk_eta;
        TBranch* b_allstub_z;
        TBranch* b_allstub_x;
        TBranch* b_allstub_y;
        TBranch* b_allstub_isPSmodule;
        TBranch* b_allstub_isTiltedBarrel;
        TBranch* b_allstub_isBarrel;

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
        trk_eta = 0;
        matchtrk_eta = 0;
        allstub_z = 0;
        allstub_x = 0;
        allstub_y = 0;
        allstub_isPSmodule = 0;
        allstub_isTiltedBarrel = 0;
        allstub_isBarrel = 0;

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
        tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
        tree->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
        tree->SetBranchAddress("allstub_z", &allstub_z, &b_allstub_z);
        tree->SetBranchAddress("allstub_x", &allstub_x, &b_allstub_x);
        tree->SetBranchAddress("allstub_y", &allstub_y, &b_allstub_y);
        tree->SetBranchAddress("allstub_isPSmodule", &allstub_isPSmodule, &b_allstub_isPSmodule);
        tree->SetBranchAddress("allstub_isTiltedBarrel", &allstub_isTiltedBarrel, &b_allstub_isTiltedBarrel);
        tree->SetBranchAddress("allstub_isBarrel", &allstub_isBarrel, &b_allstub_isBarrel);

        TCanvas c;
        char ctxt[500];
        char ctxt2[500];
        char ctxt3[500];
        char ctxt4[500];
        

        tree->GetEntry(eventNum,0);

        // -----------------------------------------------------------------------------------------------
        // My own Duplicate Removal

        // Choosing best trk based on eta and assigning it to matchtrk_eta
        float etaResidual = 100; // set high because we are finding minimum
        int minResidualInd = 100; // put high so that error shows up if not reassigned
        float tempResidual;
        float residual;
        for (int itrk = 0; itrk < (int)trk_eta->size(); itrk++) {
                tempResidual = std::abs(tp_eta->at(tpNum) - trk_eta->at(itrk));

                if (tempResidual < etaResidual) {
                etaResidual = tempResidual;
                minResidualInd = itrk;
                }
        }

        // effectively, if minResidualInd was assigned to anything at all
        if (minResidualInd < 99) {
                matchtrk_eta->at(tpNum) = trk_eta->at(minResidualInd);
        };

        // --------------------------------------------------------------------------------------------
        // getting stub and line data

        // is tp_eta positive or negative (we can then know which stubs to use)
        bool zpos = false;
        bool zneg = false;
        bool zunclear = false;

        if (tp_eta->at(tpNum) > 0.1)
                zpos = true;
        else if (tp_eta->at(tpNum) < -0.1)
                zneg = true;
        else
                zunclear = true; // if its less than 0.1, some of the stubs might 
                                 // be back and forth between positive and negative z for one tp

        // get stub indicies
        vector<int> stubInd;
        if (zpos) {
                for (int istub = 0; istub < (int)allstub_z->size(); istub++) {
                        if (allstub_z->at(istub) > 0) {
                                stubInd.push_back(istub);
                        }
                }
        }
        else if (zneg) {
                for (int istub = 0; istub < (int)allstub_z->size(); istub++) {
                        if (allstub_z->at(istub) < 0) {
                                stubInd.push_back(istub);
                        }
                }
        }
        else if (zunclear) {
                throw std::invalid_argument("too close to eta=0 so it's unclear which stubs to use");
        }

        // get stub coordinates
        int nstub = stubInd.size();
        float z_stub[nstub], r_stub[nstub];
        for (int i = 0; i < nstub; i++) {
                z_stub[i] = allstub_z->at(stubInd[i]);
                r_stub[i] = sqrt(pow(allstub_x->at(stubInd[i]),2) + pow(allstub_y->at(stubInd[i]),2));
                // Debug stuff
                // cout << "z: " << z_stub[i] << endl;
                // cout << "r: " << r_stub[i] << endl;
        }

        // Get tp and matchmatchtrk lines
        float tp_x[2], tp_y[2], matchtrk_x[2], matchtrk_y[2];
        float matchtrk_theta = 2 * atan(exp(-1 * matchtrk_eta->at(tpNum))); // theta is in radians
        float tp_theta = 2 * atan(exp(-1 * tp_eta->at(tpNum)));
        float matchtrk_slope = cos(matchtrk_theta)/sin(matchtrk_theta);
        float tp_slope = cos(tp_theta)/sin(tp_theta);
        // first point
        tp_x[0] = 0;
        matchtrk_x[0] = 0;
        tp_y[0] = tp_z0->at(tpNum);
        matchtrk_y[0] = matchtrk_z0->at(tpNum);
        // second point
        tp_x[1] = 500;
        matchtrk_x[1] = 500;
        tp_y[1] = tp_y[0] + tp_x[1] * tp_slope;
        matchtrk_y[1] = matchtrk_y[0] + matchtrk_x[1] * matchtrk_slope;

        cout << "tp_slope: " << tp_slope << endl;
        cout << "matchtrk_slope: " << matchtrk_slope << endl;
        cout << "tp_z0: " << tp_z0->at(tpNum) << endl;
        cout << "matchtrk_z0: " << matchtrk_z0->at(tpNum) << endl;

        // ----------------------------------------------------------------------------------------
        // dZ error for stubs
        // some constants
        static constexpr double m = 0.884;
        static constexpr double b = 0.507;
        static constexpr double lengthPS = 0.01;
        static constexpr double length2S = 5.025;

        // calculating dZ in cm
        float dZ_stub[nstub];
        float dR_stub[nstub];
        double length;
        double cot;
        for (int i = 0; i < nstub; i++) {
                length = allstub_isPSmodule->at(stubInd[i]) ? lengthPS : length2S;
                cot = std::abs(sinh(tp_eta->at(tpNum))); //cot(theta) is same as sinh(eta)

                dZ_stub[i] = (std::abs((allstub_isBarrel->at(stubInd[i]) ? (allstub_isTiltedBarrel->at(stubInd[i]) ? m * cot + b : 1.0) : cot) * length)) / 2; // divide by two because the formula gives full Z strip length
                dR_stub[i] = 0;
                // debug code
                cout << "dZ " << i << ": " << dZ_stub[i]*2<< "   isBarrel: " << allstub_isBarrel->at(stubInd[i]) << "   is Tilted: " << allstub_isTiltedBarrel->at(stubInd[i]) << "   length: " << length << "   is PS: " << allstub_isPSmodule->at(stubInd[i])<< "    cot: " << cot << "   z: "<< z_stub[i] << "    r: " << r_stub[i] << endl;
        }

        // ---------------------------------------------------------------------------------------------
        // Draw things 

        // Draw Stubs
        TGraphErrors *stubGraph = new TGraphErrors(nstub, r_stub, z_stub, dR_stub, dZ_stub);
        stubGraph->SetMarkerStyle(kDot);
        stubGraph->SetTitle(";R;Z");
        stubGraph->Draw("AP");
        //stubGraph->SetMinimum(0); FAIL: that sets minimum of Y-axis. I want to set minimum of x-axis

        // Draw tp and matchtrk
        TGraph *tpGraph = new TGraph(2, tp_x, tp_y);
        TGraph *matchtrkGraph = new TGraph(2, matchtrk_x, matchtrk_y);
        tpGraph->SetLineColor(kGreen);
        matchtrkGraph->SetLineColor(kBlue);
        tpGraph->SetLineWidth(1);
        matchtrkGraph->SetLineWidth(1);
        tpGraph->Draw("same");
        matchtrkGraph->Draw("same");

        // ----------------------------------------------------------------------------------------
        // legend and labels

        // make legend
        TLegend* leg = new TLegend(0.72, 0.4, 0.87, 0.55);
        leg->AddEntry(stubGraph, "stubs"); 
        leg->AddEntry(tpGraph, "tp");
        leg->AddEntry(matchtrkGraph, "matchtrk"); 
        leg->Draw();

        // make labels 
        TString eventLabel = "Event: " + to_string(eventNum);
        sprintf(ctxt, eventLabel);
        mySmallText(0.15, 0.57, 1, ctxt);
        TString tpLabel = "tp: " + to_string(tpNum);
        sprintf(ctxt2, tpLabel);
        mySmallText(0.15, 0.52, 1, ctxt2);
        TString etaLabel = "#eta: " + to_string(tp_eta->at(tpNum));
        sprintf(ctxt4, etaLabel);
        mySmallText(0.15, 0.47, 1, ctxt4);
        TString residualLabel = "#eta Residual: " + to_string(etaResidual); 
        sprintf(ctxt3, residualLabel);
        mySmallText(0.15, 0.42, 1, ctxt3);

        c.SaveAs(saveDir + saveName);

        cout << "isBarrel Size: " << allstub_isBarrel->size() << endl;
        cout << "isTiltedBarrel Size: " << allstub_isTiltedBarrel->size() << endl;
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.04;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}