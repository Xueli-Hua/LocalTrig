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
// input: folder of l1ntuples to read from - must contain the branches specified above (change this to your own personal directory)
// verbose: whether or not to print out the full output (error messgaes will always be printed)
int checkUnpacker(
		  char const *input = "/afs/cern.ch/user/x/xueli/L1rate/CMSSW_13_2_4/src/ZDCFiles/", bool verbose = false) {


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

  /* read in the event information to get the global Bx*/
  TChain eventChain("l1EventTree/L1EventTree");
  FillChain(eventChain, files);
  TTreeReader eventReader(&eventChain);
  TTreeReaderValue<UInt_t> eventbx(eventReader, "bx");
 
 
  /* create histograms */

  //histogram to give error rate in summary
  const Int_t nx = 2;
  const char *labels[nx] = {"Unpacked and Emulated Sums Match","Unpacked and Emulated Sums Don't Match"};
  TH1F *unpackerTest = new TH1F("unpackerTest","unpackerTest",nx,0,nx);
  unpackerTest->SetFillColor(38);
  unpackerTest->SetStats(0);
  for(int i=1;i<=nx;i++) unpackerTest->GetXaxis()->SetBinLabel(i,labels[i-1]);

  // correlation plots for exact agreement 
  TH2D *sumPlusCorr = new TH2D("sumPlusCorr", "sumPlusCorr", 5000, 0, 500,  5000, 0, 500 ); 
  TH2D *sumMinusCorr = new TH2D("sumMinusCorr", "sumMiniusCorr", 5000, 0, 500,  5000, 0, 500 ); 

  // to compare emulated and unpacked sums as a function of the global Bx
  TH1D* sumMinusBx = new TH1D("sumMinusBx", "sumMinusBx",  3500,0, 3500);  
  TH1D* sumPlusBx = new TH1D("sumPlusBx", "sumPlusBx", 3500,0, 3500); 
  TH1D* sumMinusEmuBx = new TH1D("sumMinusEmuBx", "sumMinusEmuBx",  3500,0, 3500);  
  TH1D* sumPlusEmuBx = new TH1D("sumPlusEmuBx", "sumPlusEmuBx", 3500,0, 3500);   

  // plot the Bx alone 
  TH1D* bx = new TH1D("globalBx", "globalBx", 3500,0, 3500); 

  // bx for different thresholds (only fill these for the unpacked sums)
  TH2D* bxThresholdPlus = new TH2D("thresholdBxPlus", "thresholdBxPlus", 3500,0, 3500, 50, 0, 100); 
  TH2D* bxThresholdMinus = new TH2D("thresholdBxMinus", "thresholdBxMinus", 3500,0, 3500, 50, 0, 100); 

  // plot the local bx for sum minus unpacked and emulated
  TH2D* localBxSumMinus  = new TH2D("localBxSumMinus", "localBxSumMinus", 3,-1.5, 1.5, 50, 0, 100);
  TH2D* localBxSumMinusEmu  = new	TH2D("localBxSumMinusEmu", "localBxSumMinusEmu", 3,-1.5, 1.5, 50, 0, 100);
  
  int unpackerErrors = 0;
  Long64_t totalEvents = unpackerReader.GetEntries(true);
  for (Long64_t i = 0; i < totalEvents; i++) {
    emuReader.Next(); unpackerReader.Next(); eventReader.Next();
    if(verbose)std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~ Event " << i << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    int eventBx      = *eventbx;

    // loop over the unpacked sums
    int max_threshold = 100; 
    int step = 2;
    // loop of the unpacker
    for(size_t j=0; j<(*unpackerSum).size(); j++){

      int unpackedBx   = (*unpackerBx)[j];
      int unpackedType = (*unpackerType)[j];
      int unpackedSum  = (*unpackerSum)[j]*2;
      int globalBx     = eventBx + unpackedBx; 
      
      if(verbose){
	std::cout << " Unpacked Bx: " << unpackedBx << " of Type " << unpackedType <<  " with sum = "<< unpackedSum << std::endl;
      }

      if(unpackedBx == -2) continue; 
     
      bx->Fill(globalBx);

      // fill the threshold histograms
      for(int threshold = 0; threshold < max_threshold; threshold += step) {
	if(unpackedSum > threshold){
	  // if we are looking at the sum minus
	  if (unpackedType == 28){
            bxThresholdMinus->Fill(globalBx, threshold, 1);
          }
	  // if we are looking at the sum plus
	  else if (unpackedType == 27){
            bxThresholdPlus->Fill(*eventbx + (*unpackerBx)[j], threshold, 1);
          } 
        }	
      }// end loop for threshold histgogram

      // now fill the other histograms
      if(unpackedType == 28){
        sumMinusBx->Fill(globalBx, unpackedSum);
	localBxSumMinus->Fill(unpackedBx, unpackedSum);
      }
      else{
        sumPlusBx->Fill(globalBx, unpackedSum);
      }
 
    } // end the unpacker loop
 
    // do the same as above, but for emulated sums
    for(size_t j = 0; j<(*emuSum).size(); j++){

      int emBx       = (*emuBx)[j];
      int emType     = (*emuType)[j];
      int emSum      = ((*emuSum)[j])*2;
      int globalBxEmu = eventBx + emBx;

      // account for the fact that the the indices are swapped in the case of minus and plus
      int unpackedSumToCompareMinus = (*unpackerSum)[j-1]*2;
      int unpackedSumToComparePlus  = (*unpackerSum)[j+1]*2;
    
      
      if(verbose){
	std::cout << " Emu Bx: " << emBx << " of Type " << emType <<  " with sum = "<< emSum <<  std::endl;
      }

      // don't look at emu bx = 3 as there is no unpacked equivalent
      // don't look at emu bx = -2 as there is not expected agreement (no possibility to apply unpacked sums)
      bool isInBxWindow = (emBx != 3) && (emBx != -2); 

      if(isInBxWindow){
	// check the zdc minus
        if(emType == 28){
	  sumMinusEmuBx->Fill(globalBxEmu, emSum);
	  localBxSumMinusEmu->Fill(emBx, emSum);
	  sumMinusCorr->Fill(emSum, unpackedSumToCompareMinus);
	  
	  if(emSum != unpackedSumToCompareMinus && unpackedSumToCompareMinus!=1023){
	    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	    std::cout << "[Error]: Unpacked sums and emulated sums do not match" << std::endl;
	    std::cout << "Emulated Sum Minus: " << emSum << " Emu Bx: " << emBx <<  " Unpacked Sum Minus: " << unpackedSumToCompareMinus << std::endl;
	    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	    unpackerTest->Fill(1.5);
	    unpackerErrors++; 
	  }
	  else{
	    unpackerTest->Fill(0.5);
	  }
	} // end check of the zdc minus

	// check the zdc plus
	else if(emType == 27){
	  
          sumPlusEmuBx->Fill(globalBxEmu, emSum);
          sumPlusCorr->Fill(emSum, unpackedSumToComparePlus);

	  if(emSum != unpackedSumToComparePlus && unpackedSumToComparePlus!=1023){
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << "[Error]: Unpacked sums and emulated sums do not match" << std::endl;
            std::cout << "Emulated Sum Plus: " << emSum << " Emu Bx: " << emBx <<  " Unpacked Sum Plus: " << unpackedSumToComparePlus << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	    unpackerTest->Fill(1.5);
            unpackerErrors++;
          }
	  else{
	    unpackerTest->Fill(0.5);
	  }
	} // end zdc plus loop
      } // end if is in Bx window
    }// end loop over the emualted sums
    if(verbose)std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;  


  } // end loop over  the  number of events
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "Finished Processing " << totalEvents << " events" << std::endl;
  std::cout << "Number of events with unpacker errors: " << unpackerErrors << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  TFile* outfile = new TFile("/eos/cms/store/group/phys_heavyions/xueli/HIZeroBias0/ZDCEmu_HIZeroBias_LUTv1_OOTPU_2/unpackerTest_Run374106_FullSample.root","RECREATE");
  unpackerTest->Write();
  sumPlusCorr->Write(); 
  sumMinusCorr->Write();
  localBxSumMinusEmu->Write();
  localBxSumMinus->Write(); 
  sumMinusBx->Write();
  sumPlusBx->Write(); 
  sumMinusEmuBx->Write();
  sumPlusEmuBx->Write(); 
  bxThresholdPlus->Write();
  bxThresholdMinus->Write();  
  bx->Write();
  outfile->Close();

  return 1;
}

int main(int argc, char *argv[]) {
  if (argc == 2)
    return checkUnpacker(argv[1]);
  else {
    cout << "ERROR" << endl;
    return -1;
  }
}
