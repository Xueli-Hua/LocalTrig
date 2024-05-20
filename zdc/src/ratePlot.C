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
#include "TH2F.h"
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

void ratePlot()
{

    vector<string> seeds={"L1_ZeroBias_copy","L1_ZDC1n_Bkp1_OR","L1_SingleJet8_ZDC1n_AsymXOR","L1_SingleJet12_ZDC1n_AsymXOR","L1_SingleJet16_ZDC1n_AsymXOR","L1_SingleJet20_ZDC1n_AsymXOR","L1_SingleJet24_ZDC1n_AsymXOR","L1_SingleJet28_ZDC1n_AsymXOR","L1_SingleJet8_ZDC1n_OR","L1_SingleJet12_ZDC1n_OR","L1_SingleJet16_ZDC1n_OR","L1_SingleJet20_ZDC1n_OR","L1_SingleJet24_ZDC1n_OR","L1_SingleJet28_ZDC1n_OR","L1_ZDC1n_AsymXOR","L1_ZDC1n_OR"};
    
    TFile *f = TFile::Open("rootfiles/runNb.root");
    
    TH1F *hlumi_itrig[16];
    TH1F *hlumi = (TH1F *)f->Get("hlumi");
	hlumi->Rebin(1);

    // plot the rates vs lumi
    TCanvas lumiCanvas("lumiCanvas", "", 0, 0, 800, 600);
	lumiCanvas.SetLeftMargin(0.15);
    TLegend lumiLegend(0.55, 0.77 ,0.77, 0.88);
    TGraphAsymmErrors ZDCRate[16];
    for (int it=0;it<16;it++) {
		hlumi_itrig[it] = (TH1F *)f->Get(("hlumi_itrig"+to_string(it)).c_str());
		hlumi_itrig[it]->Rebin(1);
		ZDCRate[it].Divide(hlumi_itrig[it], hlumi);
    	lumiCanvas.cd();

		ZDCRate[it].Scale(11245.6*880);
		ZDCRate[it].GetXaxis()->SetTitle("Lumi");
    	ZDCRate[it].GetXaxis()->CenterTitle(true);
    	ZDCRate[it].GetYaxis()->SetTitle("pre-DT rate (Hz)");
    	ZDCRate[it].GetYaxis()->CenterTitle(true);
    	//ZDCRate[it].GetXaxis()->SetLimits(0,11);
    	//ZDCRate[it].SetMaximum(150);

    	ZDCRate[it].SetMarkerColor(46);
    	ZDCRate[it].SetLineColor(46);
    	ZDCRate[it].SetMarkerSize(0.5);
    	ZDCRate[it].SetMarkerStyle(20);
    	ZDCRate[it].Draw("AP");

    	lumiLegend.Clear();
    	lumiLegend.SetBorderSize(0);
   		lumiLegend.SetFillStyle(0);
    	lumiLegend.SetTextSize(0.03);
    	lumiLegend.SetHeader(seeds[it].c_str());
    	lumiLegend.Draw();
		lumiCanvas.SaveAs(("zdcTrigPlots/Rate_v0_"+seeds[it]+".png").c_str());
    }
}
