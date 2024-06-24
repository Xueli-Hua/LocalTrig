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

    string seedmumb30To100 = "L1_SIngleMuOpen_Centrality_30_100_MinimumBiasHF1_AND_BptxAND";
    string seedmb30To100 = "L1_Centrality_30_100_MinimumBiasHF1_AND_BptxAND"; 
    string seedmu30To100 = "L1_SIngleMuOpen_Centrality_30_100_BptxAND";
    string seedtrue = "L1_AlwaysTrue";
    string seed30to100 = "L1_Centrality_30_100_BptxAND";
    string seedmu = "L1_SingleMuOpen";
    string seedmb = "L1_MinimumBiasHF1_AND_BptxAND";

    string seedmumb80To100 = "L1_SingleMuOpen_Centrality_80_100_MinimumBiasHF1_AND_BptxAND";
    string seedmb80To100 = "L1_Centrality_80_100_MinimumBiasHF1_AND_BptxAND"; 
    string seedmu80To100 = "L1_SIngleMuOpen_Centrality_80_100_BptxAND";
    string seed80to100 = "L1_Centrality_80_100_BptxAND";
    
    if (SeedBit.find(seedmb.c_str()) == SeedBit.end()) return false;
    bool l1uGTmb30To100;
    bool l1uGTmumb30To100;
    bool l1uGTmu30To100;
    bool l1uGTtrue;
    bool l1uGTmu;
    bool l1uGT30to100;
    bool l1uGTmb;

    //bool l1uGTmb80To100;
    //bool l1uGTmumb80To100;
    //bool l1uGTmu80To100;
    bool l1uGT80to100;

    // read in l1EventTree
    TChain l1EvtChain("l1EventTree/L1EventTree");
    FillChain(l1EvtChain, files);
    TTreeReader l1EvtReader(&l1EvtChain);
    TTreeReaderValue<UInt_t> runNb(l1EvtReader, "run");

    Int_t mb30To100num=0;
    Int_t mu30To100num=0;
    Int_t mumb30To100num=0;
    Int_t truenum=0;
    Int_t munum=0;
    Int_t mbnum=0;
    Int_t num30T0100=0;

    Int_t mb80To100num=0;
    Int_t mu80To100num=0;
    Int_t mumb80To100num=0;
    Int_t num80T0100=0;

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
        l1uGTmb30To100 = m_algoDecisionInitial.At(SeedBit[seedmb30To100.c_str()]);
        l1uGTmu30To100 = m_algoDecisionInitial.At(SeedBit[seedmu30To100.c_str()]);
        l1uGTmumb30To100 = m_algoDecisionInitial.At(SeedBit[seedmumb30To100.c_str()]);
        l1uGTtrue = m_algoDecisionInitial.At(SeedBit[seedtrue.c_str()]);
        l1uGTmu = m_algoDecisionInitial.At(SeedBit[seedmu.c_str()]);
        l1uGT30to100 = m_algoDecisionInitial.At(SeedBit[seed30to100.c_str()]);
        l1uGTmb = m_algoDecisionInitial.At(SeedBit[seedmb.c_str()]);

	//l1uGTmb80To100 = m_algoDecisionInitial.At(SeedBit[seedmb80To100.c_str()]);
        //l1uGTmu80To100 = m_algoDecisionInitial.At(SeedBit[seedmu80To100.c_str()]);
        //l1uGTmumb80To100 = m_algoDecisionInitial.At(SeedBit[seedmumb80To100.c_str()]);
	l1uGT80to100 = m_algoDecisionInitial.At(SeedBit[seed80to100.c_str()]);

        if (l1uGTmumb30To100) mumb30To100num++;
        if (l1uGTmb30To100) mb30To100num++;
        if (l1uGTmu30To100) mu30To100num++;
        if (l1uGTtrue) truenum++;
        if (l1uGTmu) munum++;
        if (l1uGT30to100) num30T0100++;
        if (l1uGTmb) mbnum++;

	if (l1uGT80to100 && l1uGTmb && l1uGTmu) mumb80To100num++;
        if (l1uGT80to100 && l1uGTmb) mb80To100num++;
        if (l1uGT80to100 && l1uGTmu) mu80To100num++;
	if (l1uGT80to100) num80T0100++;
    }
    cout << "L1_SingleMuOpen_Centrality_30_100_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << mumb30To100num << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mumb30To100num*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_30_100_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << mb30To100num << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mb30To100num*11245.6*880/NEvts << endl;
    cout << "L1_SIngleMuOpen_Centrality_30_100_BptxAND rate: " << setw(20) << mu30To100num << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mu30To100num*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_30_100_BptxAND rate: " << setw(20) << num30T0100 << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << num30T0100*11245.6*880/NEvts << endl;
    cout << "L1_AlwaysTrue rate: " << setw(20) << truenum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << truenum*11245.6*880/NEvts << endl;
    cout << "L1_SingleMuOpen rate: " << setw(20) << munum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << munum*11245.6*880/NEvts << endl;
    cout << "L1_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << mbnum << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mbnum*11245.6*880/NEvts << endl;

    cout << "L1_SingleMuOpen_Centrality_80_100_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << mumb80To100num << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mumb80To100num*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_80_100_MinimumBiasHF1_AND_BptxAND rate: " << setw(20) << mb80To100num << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mb80To100num*11245.6*880/NEvts << endl;
    cout << "L1_SIngleMuOpen_Centrality_80_100_BptxAND rate: " << setw(20) << mu80To100num << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << mu80To100num*11245.6*880/NEvts << endl;
    cout << "L1_Centrality_80_100_BptxAND rate: " << setw(20) << num80T0100 << "*11245.6*880" << "/" << NEvts << " = "  << setw(20) << num80T0100*11245.6*880/NEvts << endl;
    
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
