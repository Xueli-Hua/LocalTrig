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
    
    string seedmb30To100 = "L1_Centrality_30_100_MinimumBiasHF1_AND_BptxAND"; 
    string seedtrue = "L1_AlwaysTrue";
    string seedsgmo = "L1_SingleMuOpen";
    string seed80to100 = "L1_Centrality_80_100_BptxAND";
    string seed30to50 = "L1_Centrality_30_50_BptxAND";
    string seedmb = "L1_MinimumBiasHF1_AND_BptxAND";
    if (SeedBit.find(seedmb.c_str()) == SeedBit.end()) return false;
    bool l1uGTdecisionmb30To100;
    bool l1uGTdecisiontrue;
    bool l1uGTdecisionsgmo;
    bool l1uGTdecision80to100;
    bool l1uGTdecision30to50;
    bool l1uGTdecisionmb;

    // read in l1EventTree
    TChain l1EvtChain("l1EventTree/L1EventTree");
    FillChain(l1EvtChain, files);
    TTreeReader l1EvtReader(&l1EvtChain);
    TTreeReaderValue<UInt_t> runNb(l1EvtReader, "run");

    Int_t mb30To100num=0;
    Int_t truenum=0;
    Int_t sgmonum=0;
    Int_t num80to100=0;
    Int_t num30to50=0;
    Int_t mbnum=0;

    Int_t mb30To100ANDsgmonum=0;
    Int_t cen80to100ANDsgmonum=0;
    Int_t cen80to100ANDmbnum=0;
    Int_t cen30to50ANDsgmonum=0;
    Int_t cen30to50ANDmbnum=0;


    Int_t NEvts=0;
    Long64_t totalEvents = l1uGTReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        l1uGTReader.Next();l1EvtReader.Next();
        if (i % 20000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }

        if (*runNb!=375703) continue;
	NEvts++;

        if (SeedBit[seedmb.c_str()]>=m_algoDecisionInitial.GetSize()) continue;  
        l1uGTdecisionmb30To100 = m_algoDecisionInitial.At(SeedBit[seedmb30To100.c_str()]);
        l1uGTdecisiontrue = m_algoDecisionInitial.At(SeedBit[seedtrue.c_str()]);
        l1uGTdecisionsgmo = m_algoDecisionInitial.At(SeedBit[seedsgmo.c_str()]);
        l1uGTdecision30to50 = m_algoDecisionInitial.At(SeedBit[seed30to50.c_str()]);
        l1uGTdecision80to100 = m_algoDecisionInitial.At(SeedBit[seed80to100.c_str()]);
        l1uGTdecisionmb = m_algoDecisionInitial.At(SeedBit[seedmb.c_str()]);
        if (l1uGTdecisionmb30To100) {
            mb30To100num++;
            if (l1uGTdecisionsgmo) mb30To100ANDsgmonum++;
        }
        if (l1uGTdecisiontrue) truenum++;
        if (l1uGTdecisionsgmo) {
	   sgmonum++;
	   if (l1uGTdecision30to50) cen30to50ANDsgmonum++;
	   if (l1uGTdecision80to100) cen80to100ANDsgmonum++;
	}
	if (l1uGTdecisionmb) {
	   mbnum++;
	   if (l1uGTdecision30to50) cen30to50ANDmbnum++;
	   if (l1uGTdecision80to100) cen80to100ANDmbnum++;
	}
        if (l1uGTdecision30to50) num30to50++;
        if (l1uGTdecision80to100) num80to100++;
    }
    cout << "L1_Centrality_30_100_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << mb30To100num << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mb30To100num*11245.6*880/NEvts << endl;
    cout << "L1_AlwaysTrue rate: " << setw(20) << truenum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << truenum*11245.6*880/NEvts << endl;
    cout << "L1_SingleMuOpen rate: " << setw(20) << sgmonum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << sgmonum*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_30_100_MinimumBiasHF1_AND_BptxAND + L1_SingleMuOpen rate: " << setw(20) << mb30To100ANDsgmonum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mb30To100ANDsgmonum*11245.6*880/NEvts << endl;
   
    cout << "L1_Centrality_80_100_BptxAND rate: " << setw(20) << num80to100 << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << num80to100*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_30_50_BptxAND rate: " << setw(20) << num30to50 << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << num30to50*11245.6*880/NEvts << endl;
    cout << "L1_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << mbnum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mbnum*11245.6*880/NEvts << endl;
    
    cout << "L1_Centrality_80_100_BptxAND + L1_SingleMuOpen rate: " << setw(20) << cen80to100ANDsgmonum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << cen80to100ANDsgmonum*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_80_100_BptxAND + L1_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << cen80to100ANDmbnum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << cen80to100ANDmbnum*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_30_50_BptxAND + L1_SingleMuOpen rate: " << setw(20) << cen30to50ANDsgmonum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << cen30to50ANDsgmonum*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_30_50_BptxAND + L1_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << cen30to50ANDmbnum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << cen30to50ANDmbnum*11245.6*880/NEvts << endl;
   
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
