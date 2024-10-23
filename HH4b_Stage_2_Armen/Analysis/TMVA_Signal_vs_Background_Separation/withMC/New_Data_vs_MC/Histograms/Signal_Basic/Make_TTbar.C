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

void Make_TTbar()
{
 TFile *f = new TFile("Histograms_TTbar.root","RECREATE");

 #include "Histograms.txt"

 TFile *f0 = TFile::Open("../../../../../NTuples/Tree_TTbar.root");
 TTree * T_ntuples    = (TTree*)f0->Get("tree");

 #include "Input_variables.txt"

 for(int i=0;i<T_ntuples->GetEntries();i++)
   {
    T_ntuples-> GetEntry(i);

    #include "Selection.txt"

    #include "No_Trigger_SF.txt"

    #include "Fill_Histo.txt"
   }


  f->Write();
}

