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
int checkTrigFire(char const* input) {
    // read in all files in the input folder 
    vector<string> files;
    GetFiles(input, files);

    // read in L1uGTEmu information 
    TChain l1uGTChainForBit("l1uGTEmuTree/L1uGTTree");
    FillChain(l1uGTChainForBit, files);

    TChain l1uGTEmuChain("l1uGTEmuTree/L1uGTTree");
    FillChain(l1uGTEmuChain, files);
    TTreeReader l1uGTEmuReader(&l1uGTEmuChain);
    TTreeReaderArray<bool> m_algoDecisionInitial_Emu(l1uGTEmuReader, "m_algoDecisionInitial");

    // read in L1uGT information 
    TChain l1uGTChain("l1uGTTree/L1uGTTree");
    FillChain(l1uGTChain, files);
    TTreeReader l1uGTReader(&l1uGTChain);
    TTreeReaderArray<bool> m_algoDecisionInitial_unpacker(l1uGTReader, "m_algoDecisionInitial");
    
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

    /*ofstream trignames;
    trignames.open("results/trigs.txt");
    for (auto const & name: names) trignames << name.c_str() << endl;
    trignames.close();*/
    
    string seedzdc = "L1_ZDC1n_Bkp1_OR"; 
    if (SeedBit.find(seedzdc.c_str()) == SeedBit.end()) return false;
    bool l1uGTdecision1;

    // read in emulated information
    TChain emuChain("l1UpgradeEmuTree/L1UpgradeTree");
    FillChain(emuChain, files);
    TTreeReader emuReader(&emuChain);
    TTreeReaderValue<vector<float> > emuSum(emuReader, "sumZDCEt");
	 
    // read in l1UpgradeTree 
    TChain unpackerChain("l1UpgradeTree/L1UpgradeTree");
    FillChain(unpackerChain, files);
    TTreeReader unpackerReader(&unpackerChain);
    TTreeReaderValue<vector<float> > unpackerSum(unpackerReader, "sumZDCEt");
	 
    // read in l1EventTree
    TChain l1EvtChain("l1EventTree/L1EventTree");
    FillChain(l1EvtChain, files);
    TTreeReader l1EvtReader(&l1EvtChain);
    //TTreeReaderValue<UInt_t> runNb(l1EvtReader, "run");
	 
    Long64_t totalEvents = l1uGTReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        l1uGTReader.Next();l1uGTEmuReader.Next();unpackerReader.Next();emuReader.Next();l1EvtReader.Next();
        if (i % 200000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }
	
	
        if (SeedBit[seedzdc.c_str()]>=m_algoDecisionInitial_Emu.GetSize()) continue;  
        l1uGTdecision1 = m_algoDecisionInitial_Emu.At(SeedBit[seedzdc.c_str()]);
        if (l1uGTdecision1) zdcnum++; 

    }
    }

    // save histograms to file so I can look at them 
    TFile* fout = new TFile("results/runNb.root", "recreate");
    runNbHist.Write(); 
    hTrigvsSumMinus_unpacker.Write();
    hTrigvsSumPlus_unpacker.Write();
    hTrigvsSumMinus_Emu.Write();
    hTrigvsSumPlus_Emu.Write();
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
