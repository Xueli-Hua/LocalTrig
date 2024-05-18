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
    TChain l1uGTChainForBit("l1uGTEmuTree/L1uGTTree");
    FillChain(l1uGTChainForBit, files);

    TChain l1uGTChain("l1uGTEmuTree/L1uGTTree");
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
    Double_t num[16][5];
    vector<Double_t> ZBrate = {11245.6*1088,11245.6*880,11245.6*960,11245.6*394,11245.6*204};
    vector<TTreeReaderValue<UInt_t>> runRange880 = {375245,375252,375256,375259,375300,375317,375413,375441,375448,375455,375463,375658,375665,375697,375703};
    vector<TTreeReaderValue<UInt_t>> runRange960 = {375483,375491,375507,375513,375531,375545,375549};
    vector<TTreeReaderValue<UInt_t>> runRange394 = {375391};
    vector<TTreeReaderValue<UInt_t>> runRange1088 = {374925,374970,375007,375013,375055,375058,375064,375110,375145,375164,375195,375202};
    vector<TTreeReaderValue<UInt_t>> runRange204 = {374950,374961};


    // read in l1UpgradeTree 
    TChain l1UpgChain("l1UpgradeTree/L1UpgradeTree");
    FillChain(l1UpgChain, files);
    TTreeReader l1UpgReader(&l1UpgChain);
    TTreeReaderArray<float> sumZDCEt(l1UpgReader, "sumZDCEt");
    TTreeReaderValue<unsigned short> nSumsZDC(l1UpgReader, "nSumsZDC");

    // read in l1EventTree
    TChain l1EvtChain("l1EventTree/L1EventTree");
    FillChain(l1EvtChain, files);
    TTreeReader l1EvtReader(&l1EvtChain);
    TTreeReaderValue<UInt_t> runNb(l1EvtReader, "run");

    // create histograms for efficiency plots 
    //int nbins = 200;
    float min = 0;
    //float max = 2000;
    TH1F sumZDCEtHist("sumZDCEt", "", 140, min, 1400);
    TH1D *sumPlusEmuHist = new TH1D("sumPlusEmu", "sumPlusEmu", 5000, 0, 10000);
    TH1D *sumMinusEmuHist = new TH1D("sumMinusEmu", "sumMinusEmu", 5000, 0, 10000);

    Double_t zdcnum=0;
    Double_t zbnum=0;
    Double_t sgmonum=0;
    Long64_t totalEvents = l1uGTReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        l1uGTReader.Next();l1UpgReader.Next();l1EvtReader.Next();
        if (i % 20000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }

        if (SeedBit[seedzdc.c_str()]>=m_algoDecisionInitial.GetSize()) continue;  
        l1uGTdecision1 = m_algoDecisionInitial.At(SeedBit[seedzdc.c_str()]);
        if (l1uGTdecision1) {
            zdcnum++;
            for (int izdc = 0; izdc < *nSumsZDC; ++izdc) {
            sumZDCEtHist.Fill(sumZDCEt[izdc]);
            }
        }

        int sumPInt = sumZDCEt[3]*2;
        int sumNInt = sumZDCEt[4]*2;
        sumMinusEmuHist->Fill(sumNInt); 
        sumPlusEmuHist->Fill(sumPInt); 
        
        l1uGTdecision2 = m_algoDecisionInitial.At(SeedBit[seedzb.c_str()]);
        l1uGTdecision3 = m_algoDecisionInitial.At(SeedBit[seedsgmo.c_str()]);
        if (l1uGTdecision2) zbnum++;
        if (l1uGTdecision3) sgmonum++;

	for (int i=0;i<16;i++) {
	    l1uGTEmudecisions[i]=m_algoDecisionInitial.At(SeedBit[seeds[i].c_str()]);
	    if (l1uGTEmudecisions[i]) {
		if (std::find(runRange1088.begin(), runRange1088.end(), runNb) != runRange1088.end()) num[i][0]++;
		if (std::find(runRange880.begin(), runRange880.end(), runNb) != runRange880.end()) num[i][1]++;
		if (std::find(runRange960.begin(), runRange960.end(), runNb) != runRange960.end()) num[i][2]++;
		if (std::find(runRange394.begin(), runRange394.end(), runNb) != runRange394.end()) num[i][3]++;
		if (std::find(runRange204.begin(), runRange204.end(), runNb) != runRange204.end()) num[i][4]++;
	    }
	}

    }
    cout << "L1_ZDC1n_Bkp1_OR rate: " << zdcnum << "/" << totalEvents << " = " << zdcnum/totalEvents << endl;
    cout << "L1_ZeroBias_copy rate: " << zbnum << "/" << totalEvents << " = " << zbnum/totalEvents << endl;
    cout << "L1_SingleMuOpen rate: " << sgmonum << "/" << totalEvents << " = " << sgmonum/totalEvents << endl;

    const std::map<uint, int> BrNb_ = {    {0, 1088},    {1, 880},    {2, 960},    {3, 394},    {4, 204}    };

    for (int i=0;i<16;i++) {
	for (int j=0;j<5;j++) {
	cout << "Nb of Branches: " << BrNb_.at(j) << ", " << seeds[i].c_str() << " rate: " << num[i][j] << "/" << totalEvents << "*11245.6*" << BrNb_.at(j) << " = " << num[i][j]/totalEvents*ZBrate[j] << endl;
    	}
    }

    // save histograms to file so I can look at them 
    TFile* fout = new TFile("results/sumZDCEt2.root", "recreate");
    sumPlusEmuHist->Write(); 
    sumMinusEmuHist->Write(); 
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
