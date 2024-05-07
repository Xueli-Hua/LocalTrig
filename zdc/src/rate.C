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

using namespace std;

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
int rate(char const* input) {
    // read in all files in the input folder 
    vector<string> files;
    GetFiles(input, files);

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
    TTreeReaderValue<int> nSumsZDC(l1UpgReader, "nSumsZDC");

    // create histograms for efficiency plots 
    int nbins = 160;
    float min = 0;
    float max = 1600;
    TH1F sumZDCEtHist[8]("sumZDCEt", "", nbins, min, max);

    Double_t zdcnum=0;
    Double_t truenum=0;
    Double_t sgmonum=0;
    Long64_t totalEvents = l1uGTReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        l1uGTReader.Next();l1UpgReader.Next();
        if (i % 20000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }

        if (SeedBit[seedzdc.c_str()]>=m_algoDecisionInitial.GetSize()) continue;  
        l1uGTdecision1 = m_algoDecisionInitial.At(SeedBit[seedzdc.c_str()]);
        if (l1uGTdecision1) zdcnum++;
        for (int izdc = 0; i < *nSumsZDC; ++izdc) {
            sumZDCEtHist[izdc].Fill(sumZDCEt[izdc]);
        }

        l1uGTdecision2 = m_algoDecisionInitial.At(SeedBit[seedtrue.c_str()]);
        l1uGTdecision3 = m_algoDecisionInitial.At(SeedBit[seedsgmo.c_str()]);
        if (l1uGTdecision2) truenum++;
        if (l1uGTdecision3) sgmonum++;
    }
    cout << "L1_ZDC1n_Bkp1_OR rate: " << zdcnum << "/" << totalEvents << " = " << zdcnum/totalEvents << endl;
    cout << "L1_AlwaysTrue rate: " << truenum << "/" << totalEvents << " = " << truenum/totalEvents << endl;
    cout << "L1_SingleMuOpen rate: " << sgmonum << "/" << totalEvents << " = " << sgmonum/totalEvents << endl;

    // save histograms to file so I can look at them 
    TFile* fout = new TFile("results/sumZDCEt.root", "recreate");
    for (int i = 0; i < 8; ++i) {sumZDCEtHist[i].Write();}
    fout->Close();
   
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        return rate(argv[1]);
    else {
        cout << "ERROR" << endl;
        return -1;
    }
}
