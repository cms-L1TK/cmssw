//--------------------------------------------------------------------------------------------------------
// This code overlays 3 resolution vs eta plots over each other.
// There is a legend, and option to put labels on the plot.
// To use, you must have already created the output root files from L1TrackNtuplePlot.C
// 
// There is a section titled "Things to Change" just below. 
// Follow the intructions in that section to create and save the plot.
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan December 2023
//--------------------------------------------------------------------------------------------------------

#include "plotstyle.h"

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

double GetMaxHists(vector<TH1*> hists);

//void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

void overlayHists_resVsEta_3plots (){
    
    SetPlotStyle(); // Make sure "plotstyle.h" is in your folder where you run this from

    ///////////////////////////// Things to Change ////////////////////////////////////////////////////////////

    // Load in the directory where the output ROOT files are stored, and name the output file names
    TString dir = "[ENTER DIRECTORY]"; // must have "/" at end of directory
    TString file1name = "[ENTER OUTPUT FILE]"; // enter file names without ".root" at the end
    TString file2name = "[ENTER OUTPUT FILE]";
    TString file3name = "[ENTER OUTPUT FILE]";

    // save directory and identifying file name (data set and "resVsEta" will be appended onto the end)
    TString saveDir = "[ENTER SAVE DIRECTORY]"; // must have "/" at end of directory
    TString saveFileStart = "[ENTER IDENTIFYING FILE NAME]"; // Final file name shown at the bottom of the code
    
    // Property you want to plot
    TString prop = "resVsEta_";

    // Data set to compare and the paramters you wil check
    vector <TString> dataSets = {"SingleMuon_PU0_D88"}; // Enter data set in a format that can be appended to a filename
    vector <TString> params = {"eta", "ptRel", "phi", "z0"}; // keep those parameters which you wish to keep plots for (must be spelled like that as it is part of the histogram names from the output ROOT file)
    vector <TString> labels = {"Single Muon PU=0"}; // Enter label for the data set which will be on the plot
    
    vector <float> maxVals = {0.035, 0.22, 0.0012, 2.4}; // enter max y-axis values for plots in order of "params" list above

    // if you would like add extra labels, set the boolean to true and 
    bool label1On = true;
    string label1 = "Label 1";
    bool label2On = true;
    string label2 = "Label 2";

    // Add legend labels
    string legLabel1 = "Data Set 1, Res=68%";
    string legLabel2 = "Data Set 2, Res=68%";
    string legLabel3 = "Data Set 3, Res=68%";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas c;
    char ctxt[500];
    char ctxt2[500];
    char ctxt3[500];
    double max;

    int maxCounter = 0;


    // For data sets like SingleMuon, SingleElectron
    for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {

        // For parameters like eta, phi, ...
        for (int iParam = 0; iParam < params.size(); iParam++) {
            // Load hybrid and newkf root files
            TFile *file1= new TFile(dir + file1name + ".root");
            TFile *file2= new TFile(dir + file2name + ".root");
            TFile *file3= new TFile(dir + file3name + ".root");

            // copy hybrid and newkf hists for param
            TH1F *hist1= (TH1F*)file1->Get(prop + params[iParam] + "_68");
            TH1F *hist2= (TH1F*)file2->Get(prop + params[iParam] + "_68");
            TH1F *hist3= (TH1F*)file3->Get(prop + params[iParam] + "_68");

            // Set colors, markers, etc.
            hist1->SetMarkerStyle(kOpenSquare);
            hist2->SetMarkerStyle(kOpenCircle);
            hist3->SetMarkerStyle(kStar);

            hist1->SetMarkerColor(kBlue);
            hist2->SetMarkerColor(kOrange-3);
            hist3->SetMarkerColor(kGreen+2);

            // Set max y-axis values
            hist1->GetYaxis()->SetRangeUser(0, maxVals[maxCounter]);
            maxCounter++;

            // make legend
            TLegend* leg = new TLegend(0.2, 0.65, 0.6, 0.93);
            leg->AddEntry(hist1, legLabel1, "p"); 
            leg->AddEntry(hist2, legLabel2, "p");
            leg->AddEntry(hist3, legLabel3, "p");  
            
            // Draw and save histos
            hist1->Draw("p");
            hist2->Draw("p same");
            hist3->Draw("p same");
            sprintf(ctxt, labels[iDataSet]); // Add label saying 
            mySmallText(0.2, 0.57, 1, ctxt); // which data set it is
            if (label1On) {
                sprintf(ctxt2, label1);
                mySmallText(0.2, 0.52, 1, ctxt2);
            }
            if (label2On) {
                sprintf(ctxt3, label2);
                mySmallText(0.2, 0.47, 1, ctxt3);
            }

            leg->Draw();
            gStyle->SetOptStat(1);
            c.SaveAs(saveDir + saveFileStart + "_" + dataSets[iDataSet] + "_" + prop + params[iParam] + ".pdf");
            gStyle->SetOptStat(0);

            // delete pointers you want to remake
            delete hist1, hist3, hist2, leg;
            delete file1, file2;
        }
    }
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.050;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}

double GetMaxHists(vector<TH1*> hists) {
  double maxValue = 0;
  double currMax = 0;

  for(int i=0; i<hists.size(); i++) {
    currMax = hists[i]->GetMaximum();

    if(currMax > maxValue) {
      maxValue = currMax;
    }
  }

  return maxValue;
}