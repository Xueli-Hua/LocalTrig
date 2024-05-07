/*
Input: Folder of evt info ntuples
Output: A plot of evt plane phi
*/

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TChain.h"

#include "TMath.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <regex>
#include <map>


using namespace std;


void GetFiles(char const* input, vector<string>& files) {
    TSystemDirectory dir(input, input);
    TList *list = dir.GetListOfFiles();

    if (list) {
        TSystemFile *file;
        string fname;
        TIter next(list);
        while ((file = (TSystemFile*) next())) {
            fname = file->GetName();

            if (file->IsDirectory() && (fname.find(".") == string::npos)) {
                string newDir = string(input) + fname + "/";
                GetFiles(newDir.c_str(), files);
            }
            else if ((fname.find(".root") != string::npos)) {
                files.push_back(string(input) + fname);
                cout << files.back() << endl;
            }
        }
    }

    return;
}

//

void FillChain(TChain& chain, vector<string>& files) {
    for (auto file : files) {
        chain.Add(file.c_str());
    }
}
//
int EvtPlaneAng(char const* input) {
    // read in all files in the input folder 
    vector<string> files;
    GetFiles(input, files);

    // read in evtinfo tree 
    TChain evtinfoChain("eventinfoana/EventInfoNtuple");
    FillChain(evtinfoChain, files);
    TTreeReader evtinfoReader(&evtinfoChain);

    TTreeReaderValue<double> trkQx(evtinfoReader, "trkQx");
    TTreeReaderValue<double> trkQy(evtinfoReader, "trkQy");
    TTreeReaderValue<double> twQx(evtinfoReader, "twQx");
    TTreeReaderValue<double> twQy(evtinfoReader, "twQy");

    // create histograms
    int nbins = 68;
    float min = -1.7;
    float max = 1.7;
    TH1F trkphiHist("trkPsi2", "", nbins, min, max);
    TH1F twphiHist("twPsi2", "", nbins, min, max);
    TH1F trkphiRHist("trkPsi2R", "", nbins, min, max);
    TH1F twphiRHist("twPsi2R", "", nbins, min, max);

    double trkphi,twphi,trkphiRecenter,twphiRecenter;
    double totrkQx=0,totrkQy=0,totwQx=0,totwQy=0;
    int ntrkQ=0,ntwQ=0;

    Long64_t totalEvents = evtinfoReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        evtinfoReader.Next();
        if (i % 20000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }
        if (!isnan(*trkQx)) {
            totrkQx += *trkQx;
            totrkQy += *trkQy;
            ntrkQ++;
        }
        if (!isnan(*twQx)) {
            totwQx += *twQx;
            totwQy += *twQy;
            ntwQ++;
        }
    }
    cout << "ntrkQ: " << ntrkQ << " ,ntwQ: " << ntwQ << endl;
    cout << "totrkQx: " << totrkQx << " ,totrkQy: " << totrkQy << endl;
    cout << "totwQx: " << totwQx << " ,totwQy: " << totwQy << endl;
    double AvetrkQx = totrkQx/ntrkQ;
    double AvetrkQy = totrkQy/ntrkQ;
    double AvetwQx = totwQx/ntwQ;
    double AvetwQy = totwQy/ntwQ;
    cout << "AvetrkQx: " << AvetrkQx << " ,AvetrkQy: " << AvetrkQy << endl;
    cout << "AvetwQx: " << AvetwQx << " ,AvetwQy: " << AvetwQy << endl;

    TTree * evtinfotree = (&evtinfoChain)->GetTree();
    for (Long64_t i = 0; i < totalEvents; i++) {
        evtinfotree->GetEntry(i);
        if (i % 20000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }
        trkphi = TMath::ATan2(*trkQy,*trkQx)/2;
        twphi = TMath::ATan2(*twQy,*twQx)/2;
        trkphiHist.Fill(trkphi);
        twphiHist.Fill(twphi);
    
        trkphiRecenter = TMath::ATan2(*trkQy-AvetrkQy,*trkQx-AvetrkQx)/2;
        twphiRecenter = TMath::ATan2(*twQy-AvetwQy,*twQx-AvetwQx)/2;
        trkphiRHist.Fill(trkphiRecenter);
        twphiRHist.Fill(twphiRecenter);
    }
    
    // plot phi 
    TCanvas phiCanvas("phiCanvas", "", 0, 0, 800, 600);
    phiCanvas.cd();

    trkphiHist.GetXaxis()->SetTitle("#Psi_{2} from track");
    trkphiHist.GetXaxis()->CenterTitle(true);
    trkphiHist.GetYaxis()->SetTitle("Entries");
    trkphiHist.GetYaxis()->CenterTitle(true);
    trkphiHist.SetMinimum(0);

    trkphiHist.SetLineColor(kBlack);
    trkphiHist.Draw();
    trkphiRHist.SetLineColor(kGreen+2);
    trkphiRHist.SetLineWidth(2);
    trkphiRHist.Draw("same");

    TLegend phiLegend(0.33, 0.22 ,0.43, 0.32);
    phiLegend.SetBorderSize(0);
    phiLegend.SetFillStyle(0);
    phiLegend.SetTextSize(0.04);
    //phiLegend.SetHeader("");
    phiLegend.AddEntry(&trkphiHist, "Raw", "l");
    phiLegend.AddEntry(&trkphiRHist, "Recenter", "l");
    phiLegend.Draw("same");

    phiCanvas.SaveAs("results/epangtrk.png");

    phiCanvas.cd();
    twphiHist.GetXaxis()->SetTitle("#Psi_{2} from HF tower");
    twphiHist.GetXaxis()->CenterTitle(true);
    twphiHist.GetYaxis()->SetTitle("Entries");
    twphiHist.GetYaxis()->CenterTitle(true);
    twphiHist.SetMinimum(0);

    twphiHist.SetLineColor(kBlack);
    twphiHist.Draw();
    twphiRHist.SetLineColor(kGreen+2);
    twphiRHist.SetLineWidth(2);
    twphiRHist.Draw("same");

    phiLegend.Clear();
    phiLegend.AddEntry(&twphiHist, "Raw", "l");
    phiLegend.AddEntry(&twphiRHist, "Recenter", "l");
    phiLegend.Draw("same");

    phiCanvas.SaveAs("results/epangtw.png");

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        return EvtPlaneAng(argv[1]);
    else {
        cout << "ERROR" << endl;
        return -1;
    }
}
