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

void All_Variables()
{
   TString Var_Name[] = {"PT1","PT2","MSD1","MSD2","ETA1","ETA2","Xbb1", "Xbb2"
                        };

   int Nbin[]         = { 100,  100,  100,  100,   100,   100,   100,   100};
   double Xmin[]      = { 200,  200,  40,  40, -3.0, -3.0, 0.0, 0.0};
   double Xmax[]      = { 1000, 1000, 400, 400,  3.0,  3.0,  1.0,  1.0};
   int nVar = sizeof(Var_Name)/sizeof(Var_Name[0]);

for(int variable = 0; variable < nVar; variable++)
 {

// ***************************************************************************************************************
// Signal

   TFile *f0 = TFile::Open("../../../NTuples/Tree_Signal.root");
   TTree * ntps    = (TTree*)f0->Get("tree");

   Float_t MyS_weight, MyS_Var;
   ntps->SetBranchAddress( "T_weight",     &MyS_weight);
   if(Var_Name[variable]=="PT1")
      ntps->SetBranchAddress( "T_fatJet1_pt",   &MyS_Var);
    else if(Var_Name[variable]=="PT2")
      ntps->SetBranchAddress( "T_fatJet2_pt",   &MyS_Var);
    else if(Var_Name[variable]=="MSD1")
      ntps->SetBranchAddress( "T_fatJet1_msoftdrop",   &MyS_Var );
    else if(Var_Name[variable]=="MSD2")
      ntps->SetBranchAddress( "T_fatJet2_msoftdrop",   &MyS_Var );
    else if(Var_Name[variable]=="ETA1")
      ntps->SetBranchAddress( "T_fatJet1_eta",   &MyS_Var);
    else if(Var_Name[variable]=="ETA2")
      ntps->SetBranchAddress( "T_fatJet2_eta",   &MyS_Var);
    else if(Var_Name[variable]=="Xbb1")
      ntps->SetBranchAddress( "T_fatJet1_particleNet_XbbVsQCD",   &MyS_Var);
    else if(Var_Name[variable]=="Xbb2")
      ntps->SetBranchAddress( "T_fatJet2_particleNet_XbbVsQCD",   &MyS_Var);



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
   Float_t MyB_weight[nBProc],MyB_Var[nBProc];
   for(int k=0;k<nBProc;k++)
   {
    ntpb[k]->SetBranchAddress( "T_weight",     &MyB_weight[k] );
    if(Var_Name[variable]=="PT1")
      ntpb[k]->SetBranchAddress( "T_fatJet1_pt",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PT2")
      ntpb[k]->SetBranchAddress( "T_fatJet2_pt",   &MyB_Var[k] );
    else if(Var_Name[variable]=="MSD1")
      ntpb[k]->SetBranchAddress( "T_fatJet1_msoftdrop",   &MyB_Var[k] );
    else if(Var_Name[variable]=="MSD2")
      ntpb[k]->SetBranchAddress( "T_fatJet2_msoftdrop",   &MyB_Var[k] );
    else if(Var_Name[variable]=="ETA1")
      ntpb[k]->SetBranchAddress( "T_fatJet1_eta",   &MyB_Var[k] );
    else if(Var_Name[variable]=="ETA2")
      ntpb[k]->SetBranchAddress( "T_fatJet2_eta",   &MyB_Var[k] );
    else if(Var_Name[variable]=="Xbb1")
      ntpb[k]->SetBranchAddress( "T_fatJet1_particleNet_XbbVsQCD",   &MyB_Var[k] );
    else if(Var_Name[variable]=="Xbb2")
      ntpb[k]->SetBranchAddress( "T_fatJet2_particleNet_XbbVsQCD",   &MyB_Var[k] );
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

