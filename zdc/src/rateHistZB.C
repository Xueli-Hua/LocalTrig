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
#include <iomanip>

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
int rate(char const* input, char const* output) {
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
    
    /*string seedzdc = "L1_SingleJet8_ZDC1n_AsymXOR"; 
    string seedzb = "L1_ZeroBias_copy";
    string seedsgmo = "L1_SingleMuOpen";
    if (SeedBit.find(seedzdc.c_str()) == SeedBit.end()) return false;
    bool l1uGTdecision1;
    bool l1uGTdecision2;
    bool l1uGTdecision3;
    Double_t zdcnum=0;
    Double_t zbnum=0;
    Double_t sgmonum=0;

    vector<string> seeds={"L1_ZeroBias_copy","L1_SingleMu5","L1_SingleMu3","L1_SingleJet8_BptxAND","L1_SingleJet16","L1_SingleJet28","L1_SingleEG5","L1_SingleEG15","L1_SingleMuOpen_SingleEG15","L1_ZDC1n_Bkp1_OR","L1_SingleJet8_ZDC1n_AsymXOR","L1_SingleJet12_ZDC1n_AsymXOR","L1_SingleJet16_ZDC1n_AsymXOR","L1_SingleJet20_ZDC1n_AsymXOR","L1_SingleJet24_ZDC1n_AsymXOR","L1_SingleJet28_ZDC1n_AsymXOR","L1_SingleJet8_ZDC1n_OR","L1_SingleJet12_ZDC1n_OR","L1_SingleJet16_ZDC1n_OR","L1_SingleJet20_ZDC1n_OR","L1_SingleJet24_ZDC1n_OR","L1_SingleJet28_ZDC1n_OR","L1_ZDC1n_AsymXOR","L1_ZDC1n_OR"};
    bool l1uGTEmudecisions[24];
    Double_t num[24];*/

    //vector<UInt_t> runRange880 = {375245,375252,375256,375259,375300,375317,375413,375441,375448,375455,375463,375658,375665,375697,375703};
    //vector<UInt_t> runRange1088 = {374925,374970,375007,375013,375055,375058,375064,375110,375145,375164,375195,375202};

    for (auto const & name: names) {
	if (SeedBit.find(name.c_str()) == SeedBit.end()) return false;
    }
    Int_t Nseeds = names.size();
    bool l1uGTEmu[Nseeds];
    Long64_t npass[Nseeds];
	 
    // read in l1EventTree
    TChain l1EvtChain("l1EventTree/L1EventTree");
    FillChain(l1EvtChain, files);
    TTreeReader l1EvtReader(&l1EvtChain);
    TTreeReaderValue<UInt_t> runNb(l1EvtReader, "run");
    TTreeReaderValue<UInt_t> lumi(l1EvtReader, "lumi");
	 
    // create histograms
    TH1F runNbHist("runNb","runNb",900,374900,375800);
    TH1F hlumi("lumiNb","lumiNb",1000,0,1000);
    /*TH1F hlumi_itrig[24];
    for (int i=0;i<24;i++) {
    	hlumi_itrig[i].SetName("lumiNb");
	hlumi_itrig[i].SetTitle("lumiNb");
	hlumi_itrig[i].SetBins(1000,0,1000);
    }*/
	 
    Int_t NEvts=0;
    Long64_t totalEvents = l1uGTReader.GetEntries(true);

    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        l1uGTReader.Next();l1uGTEmuReader.Next();l1EvtReader.Next();
        if (i % 200000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }

	runNbHist.Fill(*runNb);
	hlumi.Fill(*lumi);
	
	//if (std::find(runRange1088.begin(), runRange1088.end(), *runNb) == runRange1088.end()) continue;
	if (*runNb!=375703) continue;
	NEvts++;

        /*if (SeedBit[seedzdc.c_str()]>=m_algoDecisionInitial_Emu.GetSize()) continue;  
        l1uGTdecision1 = m_algoDecisionInitial_Emu.At(SeedBit[seedzdc.c_str()]);
        l1uGTdecision2 = m_algoDecisionInitial_Emu.At(SeedBit[seedzb.c_str()]);
        l1uGTdecision3 = m_algoDecisionInitial_Emu.At(SeedBit[seedsgmo.c_str()]);
        if (l1uGTdecision1) zdcnum++; 
        if (l1uGTdecision2) zbnum++;
        if (l1uGTdecision3) sgmonum++;
	
	for (int it=0;it<24;it++) {
	    l1uGTEmudecisions[it]=m_algoDecisionInitial_Emu.At(SeedBit[seeds[it].c_str()]);
	    //l1uGTdecisions[it]=m_algoDecisionInitial_unpacker.At(SeedBit[seeds[it].c_str()]);
	    if (l1uGTEmudecisions[it]) {
		num[it]++;
		hlumi_itrig[it].Fill(*lumi);
	    }
	}*/

	for (unsigned int is=0;is<names.size();is++) {
	    l1uGTEmu[is]=m_algoDecisionInitial_Emu.At(SeedBit[names[is].c_str()]);
	    if (l1uGTEmu[is]) npass[is]++;
	}

    }
    const std::map<uint, int> BrNb_ = {    {0, 1088},    {1, 880},    {2, 960},    {3, 394},    {4, 204},	{5,875} };

    /*ofstream Parttrig;
    Parttrig.open("results/Parttrig_"+string(output)+".txt");
    for (int i=0;i<24;i++) {
	Parttrig << "Nb of Branches: " << BrNb_.at(1) << ", " << seeds[i].c_str() << " rate: " << num[i] << "/" << NEvts << "*11245.6*" << BrNb_.at(1) << " = " << num[i]/NEvts*11245.6*BrNb_.at(1) << endl;
    }*/

    ofstream trigrates;
    trigrates.open("results/trigRates_"+string(output)+".txt");
    for (unsigned int j=0;j<names.size();j++){
	trigrates << names[j].c_str() << setw(20) << npass[j] << "/" << NEvts << "*11245.6*" << BrNb_.at(1) << " = "  << setw(20) << npass[j]*11245.6*BrNb_.at(1)/NEvts << endl;
    }
    trigrates.close();

    // save histograms to file so I can look at them 
    TFile* fout = new TFile(("results/runNb_"+string(output)+".root").c_str(), "recreate");
    runNbHist.Write(); 
    hlumi.SetName("hlumi");
    hlumi.Write();
    /*for (int it=0;it<24;it++) {
	hlumi_itrig[it].SetName(("hlumi_itrig"+to_string(it)).c_str());
	hlumi_itrig[it].Write();
    }*/
    fout->Close();
   
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return rate(argv[1],argv[2]);
    else {
        cout << "ERROR" << endl;
        return -1;
    }
}
