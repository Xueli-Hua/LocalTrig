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
int rate(char const* input) {
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
    string seedzb = "L1_ZeroBias_copy";
    string seedsgmo = "L1_SingleMuOpen";
    if (SeedBit.find(seedzdc.c_str()) == SeedBit.end()) return false;
    bool l1uGTdecision1;
    bool l1uGTdecision2;
    bool l1uGTdecision3;

    vector<string> seeds={"L1_ZeroBias_copy","L1_ZDC1n_Bkp1_OR","L1_SingleJet8_ZDC1n_AsymXOR","L1_SingleJet12_ZDC1n_AsymXOR","L1_SingleJet16_ZDC1n_AsymXOR","L1_SingleJet20_ZDC1n_AsymXOR","L1_SingleJet24_ZDC1n_AsymXOR","L1_SingleJet28_ZDC1n_AsymXOR","L1_SingleJet8_ZDC1n_OR","L1_SingleJet12_ZDC1n_OR","L1_SingleJet16_ZDC1n_OR","L1_SingleJet20_ZDC1n_OR","L1_SingleJet24_ZDC1n_OR","L1_SingleJet28_ZDC1n_OR","L1_ZDC1n_AsymXOR","L1_ZDC1n_OR"};
    bool l1uGTEmudecisions[16];
    bool l1uGTdecisions[16];
    Double_t num[16];
    vector<UInt_t> runRange = {375245,375252,375256,375259,375300,375317,375413,375441,375448,375455,375463,375658,375665,375697,375703};
    
    // read in emulated information
    TChain emuChain("l1UpgradeEmuTree/L1UpgradeTree");
    FillChain(emuChain, files);
    TTreeReader emuReader(&emuChain);
    TTreeReaderValue<vector<float> > emuSum(emuReader, "sumZDCEt");
    TTreeReaderValue<vector<short>> emuType(emuReader, "sumZDCType"); 
    TTreeReaderValue<vector<float>>	emuBx(emuReader, "sumZDCBx");
	 
    // read in l1UpgradeTree 
    TChain unpackerChain("l1UpgradeTree/L1UpgradeTree");
    FillChain(unpackerChain, files);
    TTreeReader unpackerReader(&unpackerChain);
    TTreeReaderValue<vector<float> > unpackerSum(unpackerReader, "sumZDCEt");
    TTreeReaderValue<vector<short>>	unpackerType(unpackerReader, "sumZDCType");
    TTreeReaderValue<vector<float>>	unpackerBx(unpackerReader, "sumZDCBx");
	 
    // read in l1EventTree
    TChain l1EvtChain("l1EventTree/L1EventTree");
    FillChain(l1EvtChain, files);
    TTreeReader l1EvtReader(&l1EvtChain);
    TTreeReaderValue<UInt_t> runNb(l1EvtReader, "run");
    TTreeReaderValue<UInt_t> lumi(l1EvtReader, "lumi");
	 
    // create histograms for efficiency plots 
    //int nbins = 200;
    //float min = 0;
    //float max = 2000;
    TH1F runNbHist("runNb","runNb",900,374900,375800);
    TH2F hTrigvsSumMinus_unpacker("hTrigvsSumMinus_unpacker","hTrigvsSumMinus_unpacker",16,0,16,5000, 0, 500);
    TH2F hTrigvsSumPlus_unpacker("hTrigvsSumPlus_unpacker","hTrigvsSumPlus_unpacker",16,0,16,5000, 0, 500);
    TH2F hTrigvsSumMinus_Emu("hTrigvsSumMinus_Emu","hTrigvsSumMinus_Emu",16,0,16,5000, 0, 500);
    TH2F hTrigvsSumPlus_Emu("hTrigvsSumPlus_Emu","hTrigvsSumPlus_Emu",16,0,16,5000, 0, 500);
    TH1D hlumi;
    TH1D hlumi[16];
    TH1D hlumi_itrig[16];
	 
    Double_t zdcnum=0;
    Double_t zbnum=0;
    Double_t sgmonum=0;
    Long64_t totalEvents = l1uGTReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < 1000000; i++) {
        l1uGTReader.Next();l1uGTEmuReader.Next();unpackerReader.Next();emuReader.Next();l1EvtReader.Next();
        if (i % 200000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }
	
	runNbHist.Fill(*runNb);
	hlumi.Fill(*lumi);
	
        if (SeedBit[seedzdc.c_str()]>=m_algoDecisionInitial_Emu.GetSize()) continue;  
        l1uGTdecision1 = m_algoDecisionInitial_Emu.At(SeedBit[seedzdc.c_str()]);
        l1uGTdecision2 = m_algoDecisionInitial_Emu.At(SeedBit[seedzb.c_str()]);
        l1uGTdecision3 = m_algoDecisionInitial_Emu.At(SeedBit[seedsgmo.c_str()]);
        if (l1uGTdecision1) zdcnum++; 
        if (l1uGTdecision2) zbnum++;
        if (l1uGTdecision3) sgmonum++;
	
	for (int it=0;it<16;it++) {
	    hlumi[it].Fill(*lumi);
	    l1uGTEmudecisions[it]=m_algoDecisionInitial_Emu.At(SeedBit[seeds[it].c_str()]);
	    l1uGTdecisions[it]=m_algoDecisionInitial_unpacker.At(SeedBit[seeds[it].c_str()]);
	    if (l1uGTEmudecisions[it]) {
		num[it]++;
		hlumi_itrig[it].Fill(*lumi);
	    }
	}
	
	for (size_t j=0; j<(*unpackerSum).size(); j++){
	    int unpackedBx   = (*unpackerBx)[j];
      	    int unpackedType = (*unpackerType)[j];
            int unpackedSum  = (*unpackerSum)[j]*2;
	    if(unpackedBx == -2) continue; 
	    for (int it=0;it<16;it++) {
	        if (l1uGTdecisions[it]) {
			if (unpackedType == 28) hTrigvsSumMinus_unpacker.Fill(it,unpackedSum);
			else if (unpackedType == 27) hTrigvsSumPlus_unpacker.Fill(it,unpackedSum);
		}
	    }
	}
	for (size_t j = 0; j<(*emuSum).size(); j++){
      	    //int emBx       = (*emuBx)[j];
            int emType     = (*emuType)[j];
            int emSum      = ((*emuSum)[j])*2;
	    for (int it=0;it<16;it++) {
	        if (l1uGTEmudecisions[it]) {
			if (emType == 28) hTrigvsSumMinus_Emu.Fill(it,emSum);
			else if (emType == 27) hTrigvsSumPlus_Emu.Fill(it,emSum);
		}
	    }
        }

    }
    cout << "L1_ZDC1n_Bkp1_OR rate: " << zdcnum << "/" << totalEvents << " = " << zdcnum/totalEvents << endl;
    cout << "L1_ZeroBias_copy rate: " << zbnum << "/" << totalEvents << " = " << zbnum/totalEvents << endl;
    cout << "L1_SingleMuOpen rate: " << sgmonum << "/" << totalEvents << " = " << sgmonum/totalEvents << endl;

    const std::map<uint, int> BrNb_ = {    {0, 1088},    {1, 880},    {2, 960},    {3, 394},    {4, 204}    };

    for (int i=0;i<16;i++) {
	cout << "Nb of Branches: " << BrNb_.at(1) << ", " << seeds[i].c_str() << " rate: " << num[i] << "/" << totalEvents << "*11245.6*" << BrNb_.at(1) << " = " << num[i]/totalEvents*11245.6*BrNb_.at(1) << endl;
    }

    // save histograms to file so I can look at them 
    TFile* fout = new TFile("results/runNb.root", "recreate");
    runNbHist.Write(); 
    hTrigvsSumMinus_unpacker.Write();
    hTrigvsSumPlus_unpacker.Write();
    hTrigvsSumMinus_Emu.Write();
    hTrigvsSumPlus_Emu.Write();
    hlumi.SetName(("hlumi");
    hlumi.Write();
    for (int it=0;it<16;it++) {
	hlumi[it].SetName(("hlumi_"+to_string(it)).c_str());
	hlumi_itrig[it].SetName(("hlumi_itrig"+to_string(it)).c_str());
	hlumi[it].Write();
	hlumi_itrig[it].Write();
    }
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
