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
#include "TH1D.h"
#include "TH2D.h"
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
    
    if (SeedBit.find("L1_ZDCM22") == SeedBit.end()) return false;
    bool l1uGTZDCP22_emu;
    bool l1uGTZDCM22_emu;
    bool l1uGTZDCP22_unpacker;
    bool l1uGTZDCM22_unpacker;
    bool l1uGTJet28_emu;
    bool l1uGTJet28_unpacker;

    // read in emulated information
    TChain emuChain("l1UpgradeEmuTree/L1UpgradeTree");
    FillChain(emuChain, files);
    TTreeReader emuReader(&emuChain);
    TTreeReaderValue<vector<float>> emuSum(emuReader, "sumZDCEt");
    TTreeReaderValue<vector<float>> emuJetPt(emuReader, "jetEt");
	 
    // read in l1UpgradeTree 
    TChain unpackerChain("l1UpgradeTree/L1UpgradeTree");
    FillChain(unpackerChain, files);
    TTreeReader unpackerReader(&unpackerChain);
    TTreeReaderValue<vector<float>> unpackerSum(unpackerReader, "sumZDCEt");
    TTreeReaderValue<vector<float>> unpackerJetPt(unpackerReader, "jetEt");
    

    /* create histograms */
    TH1D hZDCP22_emu("hZDCP22_emu", "ZDC Plus 22", 1024, -0.5, 1023.5);
    TH1D hZDCM22_emu("hZDCM22_emu", "ZDC Minus 22", 1024, -0.5, 1023.5);
    TH1D hZDCP22_unpacker("hZDCP22_unpacker", "ZDC Plus 22", 1024, -0.5, 1023.5);
    TH1D hZDCM22_unpacker("hZDCM22_unpacker", "ZDC Minus 22", 1024, -0.5, 1023.5);
    TH1D hZDCP22_emu_trig("hZDCP22_emu_trig", "ZDC Plus 22", 1024, -0.5, 1023.5);
    TH1D hZDCM22_emu_trig("hZDCM22_emu_trig", "ZDC Minus 22", 1024, -0.5, 1023.5);
    TH1D hZDCP22_unpacker_trig("hZDCP22_unpacker_trig", "ZDC Plus 22", 1024, -0.5, 1023.5);
    TH1D hZDCM22_unpacker_trig("hZDCM22_unpacker_trig", "ZDC Minus 22", 1024, -0.5, 1023.5);

    TH1D hJet28_emu("hJet28_emu", "Single Jet 28", 1024, -0.5, 1023.5);
    TH1D hJet28_emu_trig("hJet28_emu_trig", "Single Jet 28", 1024, -0.5, 1023.5);
    TH1D hJet28_unpacker("hJet28_unpacker", "Single Jet 28", 1024, -0.5, 1023.5);
    TH1D hJet28_unpacker_trig("hJet28_unpacker_trig", "Single Jet 28", 1024, -0.5, 1023.5);

    Long64_t totalEvents = l1uGTReader.GetEntries(true);
    // read in information from TTrees 
    for (Long64_t i = 0; i < totalEvents; i++) {
        l1uGTReader.Next();l1uGTEmuReader.Next();unpackerReader.Next();emuReader.Next();
        if (i % 200000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }

	//==============================================================================
	//for ZDCP and ZDCM trigger fires
	hZDCP22_emu.Fill((*emuSum)[4]*2);
	hZDCM22_emu.Fill((*emuSum)[5]*2);
	hZDCP22_unpacker.Fill((*unpackerSum)[5]*2);
	hZDCM22_unpacker.Fill((*unpackerSum)[4]*2);

        l1uGTZDCP22_emu = m_algoDecisionInitial_Emu.At(SeedBit["L1_ZDCP22"]);
        l1uGTZDCM22_emu = m_algoDecisionInitial_Emu.At(SeedBit["L1_ZDCM22"]);
        l1uGTZDCP22_unpacker = m_algoDecisionInitial_unpacker.At(SeedBit["L1_ZDCP22"]);
        l1uGTZDCM22_unpacker = m_algoDecisionInitial_unpacker.At(SeedBit["L1_ZDCM22"]);
        if (l1uGTZDCP22_emu) hZDCP22_emu_trig.Fill((*emuSum)[4]*2);
        if (l1uGTZDCM22_emu) hZDCM22_emu_trig.Fill((*emuSum)[5]*2);
	if (l1uGTZDCP22_unpacker) hZDCP22_unpacker_trig.Fill((*unpackerSum)[5]*2);
        if (l1uGTZDCM22_unpacker) hZDCM22_unpacker_trig.Fill((*unpackerSum)[4]*2);
	//==============================================================================

	//==============================================================================
	//for Single Jet trigger fires
	float emuMaxJetPt = -999;
	/* iterate through emu jets and find matched and unmatched jets with max pT */
    	for (size_t j = 0; j < (*emuJetPt).size(); ++j) {
            if ((*emuJetPt)[j] > emuMaxJetPt) emuMaxJetPt = (*emuJetPt)[j];
        }
	hJet28_emu.Fill(emuMaxJetPt);
	hJet28_unpacker.Fill(emuMaxJetPt);

	l1uGTJet28_emu = m_algoDecisionInitial_Emu.At(SeedBit["L1_SingleJet28"]);
	l1uGTJet28_unpacker = m_algoDecisionInitial_unpacker.At(SeedBit["L1_SingleJet28"]);
	if (l1uGTJet28_emu) hJet28_emu_trig.Fill(emuMaxJetPt);
	if (l1uGTJet28_unpacker) hJet28_unpacker_trig.Fill(emuMaxJetPt);
	//==============================================================================

    }

    // save histograms to file so I can look at them 
    TFile* fout = new TFile("results/checkTrigFire_ZBv2.root", "recreate");
    hZDCP22_emu.Write(); 
    hZDCM22_emu.Write(); 
    hZDCP22_emu_trig.Write(); 
    hZDCM22_emu_trig.Write();
    hZDCP22_unpacker.Write(); 
    hZDCM22_unpacker.Write(); 
    hZDCP22_unpacker_trig.Write(); 
    hZDCM22_unpacker_trig.Write(); 

    hJet28_emu.Write();
    hJet28_unpacker.Write();
    hJet28_emu_trig.Write();
    hJet28_unpacker_trig.Write();
    fout->Close();
   
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        return checkTrigFire(argv[1]);
    else {
        cout << "ERROR" << endl;
        return -1;
    }
}
