#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// ********************** Tau3/2 SFs **************************
  TFile* f_Tau3toTau2   = new TFile("../../../../../../TTbar_Offline_SFs_2023/Tau3OverTau2_SF.root");
  TH1D * _SF_Tau3toTau2 = (TH1D*)f_Tau3toTau2->Get("SF");
// ********************** TXbb SFs **************************
  TFile* f_Xbb   = new TFile("../../../../../../TTbar_Offline_SFs_2023/Xbb_SF.root");
  TH1D * _SF_Xbb = (TH1D*)f_Xbb->Get("SF");
// ********************** PTJJ SFs **************************
  TFile* f_ptjj   = new TFile("../../../../../../TTbar_Offline_SFs_2023/PTJJ_SF.root");
  TH1D * _SF_PTJJ = (TH1D*)f_ptjj->Get("SF");
// ************************************************************

void Make_TTbar_Corr()
{
 TFile *f = new TFile("Histograms_TTbar_Corr.root","RECREATE");

 #include "Histograms.txt"

 TFile *f0 = TFile::Open("../../../../../NTuples/Tree_TTbar.root");
 TTree * T_ntuples    = (TTree*)f0->Get("tree");

 #include "Input_variables.txt"

 double pt_MIN=250.0;

 for(int i=0;i<T_ntuples->GetEntries();i++)
   {
    T_ntuples-> GetEntry(i);

    #include "Selection.txt"

    #include "No_Trigger_SF.txt"

    #include "Offlie_TTbar_Corr.txt"

    #include "Fill_Histo.txt"
   }


  f->Write();

}

