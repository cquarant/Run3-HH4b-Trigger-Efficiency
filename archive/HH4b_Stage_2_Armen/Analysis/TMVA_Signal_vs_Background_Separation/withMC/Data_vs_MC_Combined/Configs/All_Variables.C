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
   TString Var_Name[] = {"PT1","PT2","PT3","PT4","ETAb1","ETAb2","ETAq1","ETAq2","PTb1","PTb2","PTq1","PTq2","MINVQQ","DETAQQ",
                         "BTG1_DeepCSV","BTG1_DeepJet","BTG2_DeepCSV","BTG2_DeepJet","MINVBB"};
   int Nbin[]         = {50,50,50,50,94,94,94,94,50,50,50,50,100,100,100,100,100,100,100};
   double Xmin[]      = {80,70,60,20,-4.7,-4.7,-4.7,-4.7,30,30,30,30,400,3.5,0,0,0,0,0};
   double Xmax[]      = {480,370,310,220,4.7,4.7,4.7,4.7,330,330,330,330,3000,10,1,1,1,1,500};
   int nVar = sizeof(Var_Name)/sizeof(Var_Name[0]);

for(int variable = 0; variable < nVar; variable++)
 {

// ***************************************************************************************************************
// Signal

   TFile *f00 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_VBFHToBB_M_125_.root");
   TFile *f01 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_GluGluHToBB_M125_.root");
//   TFile *f00 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_VBFHToBB_M_125_dipoleRecoilOn_.root");
//   TFile *f01 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_GluGluHToBB_M125_POWHEG_.root");
   TTree * ntps[2]    = {(TTree*)f00->Get("Sig"), (TTree*)f01->Get("Sig")};
   int nSProc = sizeof(ntps)/sizeof(ntps[0]);

   Float_t MyS_weight[nSProc],MyS_Var[nSProc];
   for(int k=0;k<nSProc;k++)
   {
    ntps[k]->SetBranchAddress( "T_weight",     &MyS_weight[k] );
    if(Var_Name[variable]=="PT1")
      ntps[k]->SetBranchAddress( "T_pt0",   &MyS_Var[k] );
    else if(Var_Name[variable]=="PT2")
      ntps[k]->SetBranchAddress( "T_pt1",   &MyS_Var[k] );
    else if(Var_Name[variable]=="PT3")
      ntps[k]->SetBranchAddress( "T_pt2",   &MyS_Var[k] );
    else if(Var_Name[variable]=="PT4")
      ntps[k]->SetBranchAddress( "T_pt3",   &MyS_Var[k] );
    else if(Var_Name[variable]=="ETAb1")
      ntps[k]->SetBranchAddress( "T_etab1",   &MyS_Var[k] );
    else if(Var_Name[variable]=="ETAb2")
      ntps[k]->SetBranchAddress( "T_etab2",   &MyS_Var[k] );
    else if(Var_Name[variable]=="ETAq1")
      ntps[k]->SetBranchAddress( "T_etaq1",   &MyS_Var[k] );
    else if(Var_Name[variable]=="ETAq2")
      ntps[k]->SetBranchAddress( "T_etaq2",   &MyS_Var[k] );
    else if(Var_Name[variable]=="PTb1")
      ntps[k]->SetBranchAddress( "T_ptb1",   &MyS_Var[k] );
    else if(Var_Name[variable]=="PTb2")
      ntps[k]->SetBranchAddress( "T_ptb2",   &MyS_Var[k] );
    else if(Var_Name[variable]=="PTq1")
      ntps[k]->SetBranchAddress( "T_ptq1",   &MyS_Var[k] );
    else if(Var_Name[variable]=="PTq2")
      ntps[k]->SetBranchAddress( "T_ptq2",   &MyS_Var[k] );
    else if(Var_Name[variable]=="MINVQQ")
      ntps[k]->SetBranchAddress( "T_mqq",      &MyS_Var[k] );
    else if(Var_Name[variable]=="DETAQQ")
      ntps[k]->SetBranchAddress( "T_dETAqq",   &MyS_Var[k] );
    else if(Var_Name[variable]=="BTG1_DeepCSV")
      ntps[k]->SetBranchAddress( "T_btgb1_DeepCSV",  &MyS_Var[k] );
    else if(Var_Name[variable]=="BTG1_DeepJet")
      ntps[k]->SetBranchAddress( "T_btgb1_DeepJet",  &MyS_Var[k] );
    else if(Var_Name[variable]=="BTG2_DeepCSV")
      ntps[k]->SetBranchAddress( "T_btgb2_DeepCSV",  &MyS_Var[k] );
    else if(Var_Name[variable]=="BTG2_DeepJet")
      ntps[k]->SetBranchAddress( "T_btgb2_DeepJet",  &MyS_Var[k] );
    else if(Var_Name[variable]=="MINVBB")
      ntps[k]->SetBranchAddress( "T_mbb",   &MyS_Var[k] );
   }

// QCD
   TFile *f0 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT100to200_TuneCP5_.root");
   TFile *f1 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT200to300_TuneCP5_.root");
   TFile *f2 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT300to500_TuneCP5_.root");
   TFile *f3 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT500to700_TuneCP5_.root");
   TFile *f4 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT700to1000_TuneCP5_.root");
   TFile *f5 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT1000to1500_TuneCP5_.root");
   TFile *f6 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT1500to2000_TuneCP5_.root");
   TFile *f7 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT2000toInf_TuneCP5_.root");
//Single_Top
   TFile *f8  = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_t_channel_top_4f_InclusiveDecays_TuneCP5_.root");
   TFile *f9  = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_t_channel_antitop_4f_InclusiveDecays_TuneCP5_.root");
   TFile *f10  = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_tW_top_5f_inclusiveDecays_TuneCP5_.root");
   TFile *f11 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_.root");
//ttbar
   TFile *f12 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_TTTo2L2Nu_TuneCP5_.root");
   TFile *f13 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_TTToHadronic_TuneCP5_.root");
   TFile *f14 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_TTToSemiLeptonic_TuneCP5_.root");
//WJets
   TFile *f15 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT200to400_qc19_3j_TuneCP5_.root");
   TFile *f16 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT400to600_qc19_3j_TuneCP5_.root");
   TFile *f17 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT600to800_qc19_3j_TuneCP5_.root");
   TFile *f18 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT800toInf_qc19_3j_TuneCP5_.root");
//ZJets
   TFile *f19 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT200to400_qc19_4j_TuneCP5_.root");
   TFile *f20 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_.root");
   TFile *f21 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_.root");
   TFile *f22= TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT800toInf_qc19_4j_TuneCP5_.root");
// Full Data
   TFile *f23 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/Data_Derived/NTuples_Full_Data_2016_.root");

   TTree * ntpb[] = { (TTree*)f0->Get("Bkg"),(TTree*)f1->Get("Bkg"),(TTree*)f2->Get("Bkg"),(TTree*)f3->Get("Bkg"),(TTree*)f4->Get("Bkg"),(TTree*)f5->Get("Bkg"),(TTree*)f6->Get("Bkg"),(TTree*)f7->Get("Bkg"),
                      (TTree*)f8->Get("Bkg"),(TTree*)f9->Get("Bkg"),(TTree*)f10->Get("Bkg"),(TTree*)f11->Get("Bkg"),
                      (TTree*)f12->Get("Bkg"),(TTree*)f13->Get("Bkg"),(TTree*)f14->Get("Bkg"),
                      (TTree*)f15->Get("Bkg"),(TTree*)f16->Get("Bkg"),(TTree*)f17->Get("Bkg"),(TTree*)f18->Get("Bkg"),
                      (TTree*)f19->Get("Bkg"),(TTree*)f20->Get("Bkg"),(TTree*)f21->Get("Bkg"),(TTree*)f22->Get("Bkg"),
                      (TTree*)f23->Get("Bkg")
                    };


   int nBProc = sizeof(ntpb)/sizeof(ntpb[0]);

   Float_t MyB_weight[nBProc],MyB_Var[nBProc];
   for(int k=0;k<nBProc;k++)
   {
    ntpb[k]->SetBranchAddress( "T_weight",     &MyB_weight[k] );
    if(Var_Name[variable]=="PT1")
      ntpb[k]->SetBranchAddress( "T_pt0",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PT2")
      ntpb[k]->SetBranchAddress( "T_pt1",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PT3")
      ntpb[k]->SetBranchAddress( "T_pt2",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PT4")
      ntpb[k]->SetBranchAddress( "T_pt3",   &MyB_Var[k] );
    else if(Var_Name[variable]=="ETAb1")
      ntpb[k]->SetBranchAddress( "T_etab1",   &MyB_Var[k] );
    else if(Var_Name[variable]=="ETAb2")
      ntpb[k]->SetBranchAddress( "T_etab2",   &MyB_Var[k] );
    else if(Var_Name[variable]=="ETAq1")
      ntpb[k]->SetBranchAddress( "T_etaq1",   &MyB_Var[k] );
    else if(Var_Name[variable]=="ETAq2")
      ntpb[k]->SetBranchAddress( "T_etaq2",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PTb1")
      ntpb[k]->SetBranchAddress( "T_ptb1",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PTb2")
      ntpb[k]->SetBranchAddress( "T_ptb2",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PTq1")
      ntpb[k]->SetBranchAddress( "T_ptq1",   &MyB_Var[k] );
    else if(Var_Name[variable]=="PTq2")
      ntpb[k]->SetBranchAddress( "T_ptq2",   &MyB_Var[k] );
    else if(Var_Name[variable]=="MINVQQ")
      ntpb[k]->SetBranchAddress( "T_mqq",      &MyB_Var[k] );
    else if(Var_Name[variable]=="DETAQQ")
      ntpb[k]->SetBranchAddress( "T_dETAqq",   &MyB_Var[k] );
    else if(Var_Name[variable]=="BTG1_DeepCSV")
      ntpb[k]->SetBranchAddress( "T_btgb1_DeepCSV",  &MyB_Var[k] );
    else if(Var_Name[variable]=="BTG1_DeepJet")
      ntpb[k]->SetBranchAddress( "T_btgb1_DeepJet",  &MyB_Var[k] );
    else if(Var_Name[variable]=="BTG2_DeepCSV")
      ntpb[k]->SetBranchAddress( "T_btgb2_DeepCSV",  &MyB_Var[k] );
    else if(Var_Name[variable]=="BTG2_DeepJet")
      ntpb[k]->SetBranchAddress( "T_btgb2_DeepJet",  &MyB_Var[k] );
    else if(Var_Name[variable]=="MINVBB")
      ntpb[k]->SetBranchAddress( "T_mbb",   &MyB_Var[k] );
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

  TH1D * _Var_VBF_Hbb  = new TH1D(Name+"_VBF_Hbb",Name+"_VBF_Hbb",nbin,xmin,xmax);
  TH1D * _Var_ggF_Hbb  = new TH1D(Name+"_ggF_Hbb",Name+"_ggF_Hbb",nbin,xmin,xmax);
  TH1D * _Var_QCD  = new TH1D(Name+"_QCD",Name+"_QCD",nbin,xmin,xmax);
  TH1D * _Var_Single_Top  = new TH1D(Name+"_Single_Top",Name+"_Single_Top",nbin,xmin,xmax);
  TH1D * _Var_tt  = new TH1D(Name+"_tt",Name+"_tt",nbin,xmin,xmax);
  TH1D * _Var_WJets  = new TH1D(Name+"_WJets",Name+"_WJets",nbin,xmin,xmax);
  TH1D * _Var_ZJets  = new TH1D(Name+"_ZJets",Name+"_ZJets",nbin,xmin,xmax);
  TH1D * _Var_DATA   = new TH1D(Name+"_DATA",Name+"_DATA",nbin,xmin,xmax);
//***************************************************************************************

  for(int k=0;k<nSProc;k++)
   for(int i=0;i<ntps[k]->GetEntries();i++)
    {
     ntps[k]    -> GetEntry(i);
     if(k==0)
      _Var_VBF_Hbb->Fill(MyS_Var[k],MyS_weight[k]);
     if(k==1)
      _Var_ggF_Hbb->Fill(MyS_Var[k],MyS_weight[k]);
    }

  for(int k=0;k<nBProc;k++)
   for(int i=0;i<ntpb[k]->GetEntries();i++)
    {
     ntpb[k]    -> GetEntry(i);
     if(k<8)
      _Var_QCD->Fill(MyB_Var[k],MyB_weight[k]);
     else if(k>=8 && k<12)
      _Var_Single_Top->Fill(MyB_Var[k],MyB_weight[k]);
     else if(k>=12 && k<15)
      _Var_tt->Fill(MyB_Var[k],MyB_weight[k]);
     else if(k>=15 && k<19)
      _Var_WJets->Fill(MyB_Var[k],MyB_weight[k]);
     else if(k>=19 && k<23)
      _Var_ZJets->Fill(MyB_Var[k],MyB_weight[k]);
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
  delete f6;
  delete f7;
  delete f8;
  delete f9;
  delete f10;
  delete f11;
  delete f12;
  delete f13;
  delete f14;
  delete f15;
  delete f16;
  delete f17;
  delete f18;
  delete f19;
  delete f20;
  delete f21;
  delete f22;
  delete f23;

 }   // End of variable loop
}

