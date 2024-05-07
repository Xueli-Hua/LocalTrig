/*
Input: Folder of L1Ntuples
Output: A plot of the jet turn-ons with and with out L1 dR matching vs calo jet pT
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

#include "../include/GlobalAlgBlk.h"

using namespace std;
GlobalAlgBlk *l1uGT_;

unsigned int ParseAlias(std::string alias)
{
  std::smatch base_match;
  std::regex integer("L1uGT\\.m_algoDecisionInitial\\[([0-9]+)\\]");
  unsigned int nbit = 0;

  if (std::regex_match(alias, base_match, integer))
  {
    nbit = std::stoi(base_match[1].str(), nullptr);
  }

  return nbit;
}

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
int Efficiency(char const* input) {
    // read in all files in the input folder 
    vector<string> files;
    GetFiles(input, files);

    // read in reco mu and trk information 
    /*TChain recoMuChain("muonAnalyzer/MuonTree");
    TChain trkChain("PbPbTracks/trackTree");
    FillChain(recoMuChain, files);
    FillChain(trkChain, files);
    TTreeReader recoMuReader(&recoMuChain);
    TTreeReader trkReader(&trkChain);

    TTreeReaderValue<int> nTrk(trkReader, "nTrk");
    TTreeReaderArray<bool> isFake(trkReader, "isFakeVtx");
    TTreeReaderArray<bool> trkHP(trkReader, "highPurity");
    TTreeReaderArray<float> xVtx(trkReader, "xVtx");
    TTreeReaderArray<float> yVtx(trkReader, "yVtx");
    TTreeReaderArray<float> zVtx(trkReader, "zVtx");
    TTreeReaderValue<int> recomuN(recoMuReader, "nReco");
    TTreeReaderArray<float> recomuP(recoMuReader, "recoP");
    TTreeReaderArray<float> recomuPt(recoMuReader, "recoPt");
    TTreeReaderArray<float> recomuEta(recoMuReader, "recoEta");
    TTreeReaderArray<bool> recomuIsTrk(recoMuReader, "recoIsTracker");
    TTreeReaderArray<bool> recomuIDSoft(recoMuReader, "recoIDSoft"); // mu.passed(reco::Muon::SoftCutBasedId))
    TTreeReaderArray<bool> recomuIDHySoft(recoMuReader, "recoIDHybridSoft");
    TTreeReaderValue<int> innermuN(recoMuReader, "nInner");
    TTreeReaderArray<int> innermuTrkL(recoMuReader, "innerTrkLayers");
    TTreeReaderArray<int> innermuPixL(recoMuReader, "innerPixelLayers");
    TTreeReaderArray<float> innerDxy(recoMuReader, "innerDxy");
    TTreeReaderArray<float> innerDz(recoMuReader, "innerDz");
    TTreeReaderArray<bool> innerIsHPTrk(recoMuReader, "innerIsHighPurityTrack");
    
    // read in emulated mu information 
    TChain l1Chain("l1object/L1UpgradeFlatTree");
    FillChain(l1Chain, files);
    TTreeReader l1Reader(&l1Chain);
    TTreeReaderValue<vector<float>> l1muEt(l1Reader, "muonEt");
    TTreeReaderValue<vector<float>> l1muEta(l1Reader, "muonEta");
    TTreeReaderValue<vector<unsigned short>> l1muQual(l1Reader, "muonQual");*/

    // read in L1uGT information 
    TChain l1uGTChainForBit("l1uGTTree/L1uGTTree");
    FillChain(l1uGTChainForBit, files);

    TChain l1uGTChain("l1uGTTree/L1uGTTree");
    FillChain(l1uGTChain, files);
    TTreeReader l1uGTReader(&l1uGTChain);
    TTreeReaderArray<bool> m_algoDecisionInitial(l1uGTReader, "m_algoDecisionInitial");
    
    (&l1uGTChainForBit)->GetEntry(1);
    TTree * ugtree = (&l1uGTChainForBit)->GetTree();
    TList * aliases = ugtree->GetListOfAliases();
    TIter iter(aliases);
    std::vector<std::string> names;
    std::for_each(iter.Begin(), TIter::End(), [&](TObject* alias){ names.push_back(alias->GetName()); } );
    std::map<std::string, std::string> SeedAlias;
    for (auto const & name: names) {
      SeedAlias[name] = l1uGTChainForBit.GetAlias(name.c_str());
    }
    
    std::map<std::string, std::string> XMLConv; 
    std::map<std::string, unsigned int> SeedBit;
    for (auto const & name: SeedAlias) {
      if (XMLConv.find(name.first) != XMLConv.end())
        SeedBit[XMLConv[name.first]] = ParseAlias(name.second);
      else
        SeedBit[name.first] = ParseAlias(name.second);
    }

    ofstream trignames;
    trignames.open("results/trigs.txt");
    for (auto const & name: names) trignames << name.c_str() << endl;
    trignames.close();
    
    //string seed = "L1_SingleMuOpen_NotMinimumBiasHF2_AND_BptxAND"; 
    string seedzdc = "L1_ZDC1n_Bkp1_OR"; 
    string seedtrue = "L1_AlwaysTrue";
    string seedsgmo = "L1_SingleMuOpen";
    if (SeedBit.find(seedzdc.c_str()) == SeedBit.end()) return false;
    bool l1uGTdecision1;
    bool l1uGTdecision2;
    bool l1uGTdecision3;

    // read in l1UpgradeTree 
    TChain l1UpgChain("l1UpgradeTree/L1UpgradeTree");
    FillChain(l1UpgChain, files);
    TTreeReader l1UpgReader(&l1UpgChain);
    TTreeReaderArray<float> sumZDCEt(l1UpgReader, "sumZDCEt");

    // create histograms for efficiency plots 
    int nbins = 25;
    float min = 0;
    float max = 10;
    TH1F zVtxHist("zVtxHist", "", 200, 0, 200);
    TH1F l1muHist("l1muHist", "", nbins, min, max);
    TH1F recomuHist("recomuHist", "", nbins, min, max);

    Double_t zdcnum=0;
    Double_t truenum=0;
    Double_t sgmonum=0;
    Long64_t totalEvents = l1uGTReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        l1uGTReader.Next();l1UpgReader.Next();
        //l1Reader.Next(); recoMuReader.Next(); trkReader.Next(); 
        if (i % 20000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }

        if (SeedBit[seedzdc.c_str()]>=m_algoDecisionInitial.GetSize()) continue;  
        l1uGTdecision1 = m_algoDecisionInitial.At(SeedBit[seedzdc.c_str()]);
        l1uGTdecision2 = m_algoDecisionInitial.At(SeedBit[seedtrue.c_str()]);
        l1uGTdecision3 = m_algoDecisionInitial.At(SeedBit[seedsgmo.c_str()]);
        cout << "iEvt: " << i << ", sumZDCEt.size: " << sumZDCEt.GetSize() << endl;
        if (l1uGTdecision1) {
          zdcnum++;
        }
        if (l1uGTdecision2) truenum++;
        if (l1uGTdecision3) sgmonum++;
      
        /*int NtrkHP = 0;
        bool primaryVertex=false;

        // iterate through trks and do selection 
        for (int i = 0; i < *nTrk; ++i) {
            //rho = TMath::Sqrt(xVtx[i]*xVtx[i]+yVtx[i]*yVtx[i]);
            zVtxHist.Fill(TMath::Abs(zVtx[i]));
            if (!isFake[i] && TMath::Abs(zVtx[i])<25 && !(TMath::Sqrt(xVtx[i]*xVtx[i]+yVtx[i]*yVtx[i])>2)) primaryVertex = true;
            if (trkHP[i]) NtrkHP++;
        }
        if (!primaryVertex) continue;
        if (NtrkHP!=2) continue;

        // iterate through l1object muons and find max et 
        for (int i = 0; i < *recomuN; ++i) {
            if (recomuP[i]>2.5 && TMath::Abs(recomuEta[i]) < 2.4 && recomuIsTrk[i] && recomuIDHySoft[i]) { 
                recomuHist.Fill(recomuPt[i]); 
                if (l1uGTdecision) l1muHist.Fill(recomuPt[i]);
            }
        }
        */
    }
    cout << "L1_ZDC1n_Bkp1_OR rate: " << zdcnum << "/" << totalEvents << " = " << zdcnum/totalEvents << endl;
    cout << "L1_AlwaysTrue rate: " << truenum << "/" << totalEvents << " = " << truenum/totalEvents << endl;
    cout << "L1_SingleMuOpen rate: " << sgmonum << "/" << totalEvents << " = " << sgmonum/totalEvents << endl;
    TGraphAsymmErrors RecoMuEff(&l1muHist, &recomuHist);
    
    // plot the turn ons vs reco mu pt 
    TCanvas recoCanvas("recoCanvas", "", 0, 0, 500, 500);
    recoCanvas.cd();

    RecoMuEff.GetXaxis()->SetTitle("Reco #mu pT (GeV)");
    RecoMuEff.GetXaxis()->CenterTitle(true);
    RecoMuEff.GetYaxis()->SetTitle("Efficiency");
    RecoMuEff.GetYaxis()->CenterTitle(true);
    RecoMuEff.GetXaxis()->SetLimits(0,11);
    RecoMuEff.SetMinimum(0);

    RecoMuEff.SetMarkerColor(46);
    RecoMuEff.SetLineColor(46);
    RecoMuEff.SetMarkerSize(0.5);
    RecoMuEff.SetMarkerStyle(20);
    RecoMuEff.Draw("AP");

    TLegend recoLegend(0.13, 0.12 ,0.88, 0.2);
    recoLegend.SetBorderSize(0);
    recoLegend.SetFillStyle(0);
    recoLegend.SetTextSize(0.03);
    recoLegend.SetHeader(seedzdc.c_str());
    //recoLegend.AddEntry(&RecoMuEff, "#DeltaR Matched", "lep");
    recoLegend.Draw();

    TLine recoline(0,1,11,1);
    recoline.SetLineColor(46);
    recoline.Draw();

    recoCanvas.SaveAs("results/non.pdf");

    // save histograms to file so I can look at them 
    TFile* fout = new TFile("results/non.root", "recreate");

    zVtxHist.Write();
    l1muHist.Write();
    recomuHist.Write();

    fout->Close();
   
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        return Efficiency(argv[1]);
    else {
        cout << "ERROR" << endl;
        return -1;
    }
}
