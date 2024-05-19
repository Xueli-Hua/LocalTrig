//
/*
 * checkUnpacker.C: Macro to check that the unpacked sums and the emulated sums are equivalent.
 * Input: Folder of L1Ntuples
 * Output: A root file of histograms relevant to checking the unpacking. 
 * Authors: Hannah Bossi <hannah.bossi@cern.ch>, Gian Michele Innocenti, <gian.michele.innocenti@cern.ch> 
 * 9/10/23
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

#include <string>
#include <vector>
#include <iostream>

using namespace std;


void GetFiles(char const *input, vector<string> &files) {
  TSystemDirectory dir(input, input);
  TList *list = dir.GetListOfFiles();

  if (list) {
    TSystemFile *file;
    string fname;
    TIter next(list);
    while ((file = (TSystemFile *)next())) {
      fname = file->GetName();

      if (file->IsDirectory() && (fname.find(".") == string::npos)) {
        string newDir = string(input) + fname + "/";
        GetFiles(newDir.c_str(), files);
      } else if ((fname.find(".root") != string::npos)) {
        files.push_back(string(input) + fname);
        cout << files.back() << endl;
      }
    }
  }

  return;
}

void FillChain(TChain &chain, vector<string> &files) {
  for (auto file : files) {
    chain.Add(file.c_str());
  }
}

// params
// input: folder of l1ntuples to read from - must contain the branches specified above
// verbose: whether or not to print out the full output (error messgaes will always be printed)
int checkTriggerFires(
		  char const *input = "/eos/cms/store/group/phys_heavyions/hbossi/RunPrep/Run374521/", bool verbose = false, bool doOnline = false) {


  /* read in all files in the input folder */
  vector<string> files;
  GetFiles(input, files);

  /* read in emulated information */
  TChain emuChain("l1UpgradeEmuTree/L1UpgradeTree");
  FillChain(emuChain, files);
  TTreeReader emuReader(&emuChain);
  TTreeReaderValue<vector<float> > emuSum(emuReader, "sumZDCEt");
  TTreeReaderValue<vector<short>> emuType(emuReader, "sumZDCType"); 
  TTreeReaderValue<vector<float>>	emuBx(emuReader, "sumZDCBx");

  /* read in the unpacked information */
  TChain unpackerChain("l1UpgradeTree/L1UpgradeTree");
  FillChain(unpackerChain, files);
  TTreeReader unpackerReader(&unpackerChain);
  TTreeReaderValue<vector<float> > unpackerSum(unpackerReader, "sumZDCEt");
  TTreeReaderValue<vector<short>>	unpackerType(unpackerReader, "sumZDCType");
  TTreeReaderValue<vector<float>>	unpackerBx(unpackerReader, "sumZDCBx");

    /* read in trigger information */
  TChain trigChain("hltanalysis/HltTree");
  FillChain(trigChain, files);
  TTreeReader trigReader(&trigChain);
  TTreeReaderValue<int> zb(trigReader, "L1_ZeroBias");
  
  // zdc and
  // TTreeReaderValue<int> zdcAND22(trigReader, "L1_ZDC22_AND_BptxAND"); 
  // TTreeReaderValue<int> zdcAND80(trigReader, "L1_ZDC80_AND_BptxAND"); 
  // TTreeReaderValue<int> zdcAND133(trigReader, "L1_ZDC133_AND_BptxAND");

  // zdc or 
  // TTreeReaderValue<int> zdcOR22(trigReader, "L1_ZDC22_OR_BptxAND"); 
  // TTreeReaderValue<int> zdcOR80(trigReader, "L1_ZDC80_OR_BptxAND"); 
  // TTreeReaderValue<int> zdcOR133(trigReader, "L1_ZDC133_OR_BptxAND");

  // // zdc minus
  // TTreeReaderValue<int> zdcMinusAND22(trigReader, "L1_ZDCM22");
  // TTreeReaderValue<int> zdcMinusAND80(trigReader, "L1_ZDCM80");
  // TTreeReaderValue<int> zdcMinusAND133(trigReader, "L1_ZDCM133");

  // // zdc plus
  // TTreeReaderValue<int> zdcPlusAND22(trigReader, "L1_ZDCP22");
  // TTreeReaderValue<int> zdcPlusAND80(trigReader, "L1_ZDCP80");
  // TTreeReaderValue<int> zdcPlusAND133(trigReader, "L1_ZDCP133");

  // // not HF AND 
  // TTreeReaderValue<int> notHFAND_Jet8(trigReader, "L1_SingleJet8_NotMinimumBiasHF2_AND_BptxAND");
  // TTreeReaderValue<int> notHFAND_Jet16(trigReader, "L1_SingleJet16_NotMinimumBiasHF2_AND_BptxAND");
  // TTreeReaderValue<int> notHFAND_Jet20(trigReader, "L1_SingleJet20_NotMinimumBiasHF2_AND_BptxAND");
  // TTreeReaderValue<int> notHFAND_Jet24(trigReader, "L1_SingleJet24_NotMinimumBiasHF2_AND_BptxAND");
  // TTreeReaderValue<int> notHFAND_Jet28(trigReader, "L1_SingleJet28_NotMinimumBiasHF2_AND_BptxAND");

  /* read in the event information to get the global Bx*/
  TChain eventChain("l1EventTree/L1EventTree");
  FillChain(eventChain, files);
  TTreeReader eventReader(&eventChain);
  TTreeReaderValue<UInt_t> eventbx(eventReader, "bx");

  /* read in information on the emulated trigger fires*/
  TChain l1TrigChain("l1uGTTree/L1uGTTree");
  FillChain(l1TrigChain, files);
  TTreeReader emuTrigReader(&l1TrigChain);
  TTreeReaderValue<vector<bool>> initialTrigDecision(emuTrigReader, "m_algoDecisionInitial");
 
 
  /* create histograms */
  TH1D* hZDCP22 = new TH1D("hZDCP22", "ZDC Plus 22", 1024, -0.5, 1023.5);
  TH1D* hZDCM22 = new TH1D("hZDCM22", "ZDC Minus 22", 1024, -0.5, 1023.5);



  Long64_t totalEvents = unpackerReader.GetEntries(true);
  for (Long64_t i = 0; i < totalEvents; i++) {
    emuReader.Next(); unpackerReader.Next(); eventReader.Next(); emuTrigReader.Next(); trigReader.Next();

    if(!doOnline){
      if((*initialTrigDecision)[500] == 1 ){

          hZDCP22->Fill((*emuSum)[4]*2);            
      }
      else if ((*initialTrigDecision)[503] == 1){
          hZDCM22->Fill((*emuSum)[5]*2);
      }

    }
    
  




  } // end loop over  the  number of events


  TFile* outfile = new TFile("Run374521_ZDCTriggerCheck.root","RECREATE");
  outfile->cd();
  hZDCP22->Write();
  hZDCM22->Write();
  outfile->Close();

  return 1;
}

int main(int argc, char *argv[]) {
  if (argc == 2)
    return checkTriggerFires(argv[1]);
  else {
    cout << "ERROR" << endl;
    return -1;
  }
}
