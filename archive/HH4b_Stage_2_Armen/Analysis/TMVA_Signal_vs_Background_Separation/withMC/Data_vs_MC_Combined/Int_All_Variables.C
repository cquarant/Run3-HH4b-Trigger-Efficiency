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



void Int_All_Variables()
{
 TString Var_Name[]={"T_fatJet1_OppositeHemisphereHasBJet","T_fatJet2_OppositeHemisphereHasBJet", "T_NJets", "T_nBTaggedJets"};

   int Nbin[]         = { 2,  2,  10,  10 };
   double Xmin[]      = { 0,  0,   0,   0 };
   double Xmax[]      = { 2,  2,  10,  10 };

   int nVar = sizeof(Var_Name)/sizeof(Var_Name[0]);

for(int variable = 0; variable < nVar; variable++)
 {
// ***************************************************************************************************************
// Signal

   TFile *f0 = TFile::Open("../../../NTuples/Tree_Signal.root");
   TTree * ntps    = (TTree*)f0->Get("tree");

   Float_t MyS_weight;
//   Int_t MyS_Var;
   Bool_t MyS_Var;
   ntps->SetBranchAddress( "T_weight",         &MyS_weight);
   ntps->SetBranchAddress( Var_Name[variable], &MyS_Var);

// QCD
   TFile *f1 = TFile::Open("../../../NTuples/Tree_QCD.root");
//TTbar
   TFile *f2 = TFile::Open("../../../NTuples/Tree_TTbar_Corr.root");
//VJets
   TFile *f3 = TFile::Open("../../../NTuples/Tree_VJ.root");
//VV
   TFile *f4 = TFile::Open("../../../NTuples/Tree_VV.root");
// Full Data
   TFile *f5 = TFile::Open("../../../NTuples/Data_2022_PostEE.root");

   TTree * ntpb[] = { (TTree*)f1->Get("tree"),(TTree*)f2->Get("tree"),(TTree*)f3->Get("tree"),(TTree*)f4->Get("tree"),(TTree*)f5->Get("tree")};

   int nBProc = sizeof(ntpb)/sizeof(ntpb[0]);
   Float_t MyB_weight[nBProc];
   Bool_t MyB_Var[nBProc];
//   Int_t MyB_Var[nBProc];
   for(int k=0;k<nBProc;k++)
   {
    ntpb[k]->SetBranchAddress( "T_weight",         &MyB_weight[k] );
    ntpb[k]->SetBranchAddress( Var_Name[variable], &MyB_Var[k] );
   }

//***************************************************************************************
  TString Name = Var_Name[variable];

  TFile *f = new TFile(Name+"_Output.root","RECREATE");
  gStyle->SetOptFile(0);
  gStyle->SetOptStat("mre");
  gStyle->SetPaintTextFormat("1.2e");

  int nbin=Nbin[variable];
  double xmin=Xmin[variable];
  double xmax=Xmax[variable];

  TH1D * _Var_Signal  = new TH1D(Name+"_Signal",Name+"_Signal",nbin,xmin,xmax);
  TH1D * _Var_QCD     = new TH1D(Name+"_QCD",Name+"_QCD",nbin,xmin,xmax);
  TH1D * _Var_tt      = new TH1D(Name+"_tt",Name+"_tt",nbin,xmin,xmax);
  TH1D * _Var_VJets   = new TH1D(Name+"_VJets",Name+"_VJets",nbin,xmin,xmax);
  TH1D * _Var_VV      = new TH1D(Name+"_VV",Name+"_VV",nbin,xmin,xmax);
  TH1D * _Var_DATA    = new TH1D(Name+"_DATA",Name+"_DATA",nbin,xmin,xmax);
//***************************************************************************************

   for(int i=0;i<ntps->GetEntries();i++)
    {
      ntps -> GetEntry(i);
      _Var_Signal->Fill(MyS_Var,MyS_weight);
    }

  for(int k=0;k<nBProc;k++)
   for(int i=0;i<ntpb[k]->GetEntries();i++)
    {
     ntpb[k]    -> GetEntry(i);
     if(k==0)
      _Var_QCD->Fill(MyB_Var[k],MyB_weight[k]);
     else if(k==1)
      _Var_tt->Fill(MyB_Var[k],MyB_weight[k]);
     else if(k==2)
      _Var_VJets->Fill(MyB_Var[k],MyB_weight[k]);
     else if(k==3)
      _Var_VV->Fill(MyB_Var[k],MyB_weight[k]);
     else
      _Var_DATA->Fill(MyB_Var[k],MyB_weight[k]);
    }



  f->Write();

//  delete f00;
//  delete f01;

  delete f0;
  delete f1;
  delete f2;
  delete f3;
  delete f4;
  delete f5;

 }   // End of variable loop
}

