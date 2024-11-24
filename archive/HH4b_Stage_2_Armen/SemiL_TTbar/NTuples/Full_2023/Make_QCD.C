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


void Make_QCD()
{
 TFile *f = new TFile("Histograms_QCD.root","RECREATE");

 TH1D *_FatJet1_pt             = new TH1D("FatJet1_pt","FatJet1_pt",200,0,1200);
 TH1D *_FatJet1_eta            = new TH1D("FatJet1_eta","FatJet1_eta",100,-5,5);
 TH1D *_FatJet1_MassSD         = new TH1D("FatJet1_MassSD","FatJet1_MassSD",500,0,500);
 TH1D *_FatJet1_Xbb            = new TH1D("FatJet1_Xbb","FatJet1_Xbb",100,0.0,1.0);
 TH1D *_FatJet1_Xbb_Legacy     = new TH1D("FatJet1_Xbb_Legacy","FatJet1_Xbb_Legacy",100,0.0,1.0);
 TH1D *_FatJet1_Xbb_Legacy_AN  = new TH1D("FatJet1_Xbb_Legacy_AN","FatJet1_Xbb_Legacy_AN",100,0.0,1.0);
 TH1D *_FatJet1_Tau3OverTau2   = new TH1D("FatJet1_Tau3OverTau2","FatJet1_Tau3OverTau2",100,0.0,1.0);
 TH1D *_FatJet2_pt             = new TH1D("FatJet2_pt","FatJet2_pt",200,0,1200);
 TH1D *_FatJet2_eta            = new TH1D("FatJet2_eta","FatJet2_eta",100,-5,5);
 TH1D *_FatJet2_MassSD         = new TH1D("FatJet2_MassSD","FatJet2_MassSD",500,0,500);

 TH1D *_MET                    = new TH1D("MET","MET",100,0,500);
 TH1D *_Lep1_Pt                = new TH1D("Lep1_Pt","Lep1_Pt",300,0,300);
 TH1D *_Dr_LFJ                 = new TH1D("Dr_LFJ","Dr_LFJ",100,0,10);
 TH1D *_Dr_J1FJ                = new TH1D("Dr_J1FJ","Dr_J1FJ",100,0,10);
 TH1D *_Dr_J2FJ                = new TH1D("Dr_J2FJ","Dr_J2FJ",100,0,10);
 TH1D *_Dr_JmaxL               = new TH1D("Dr_JmaxL","Dr_JmaxL",100,0,10);

 TFile *f0 = TFile::Open("Tree_QCD.root");
 TTree * T_ntuples    = (TTree*)f0->Get("tree");

 Float_t weight;
 Float_t fatJet1_pt, fatJet1_eta, fatJet1_phi,fatJet1_msoftdrop;
 Float_t fatJet1_particleNet_XbbVsQCD, fatJet1_particleNet_Xbb_Legacy, fatJet1_particleNet_Xbb_Legacy_AN, fatJet1_Tau3OverTau2;
 Float_t fatJet2_pt, fatJet2_eta, fatJet2_msoftdrop;
 Float_t MET, Lep1_Pt, Dr_LFJ, Dr_J1FJ, Dr_J2FJ, Dr_JmaxL;

 T_ntuples->SetBranchAddress( "T_weight",                            &weight);
 T_ntuples->SetBranchAddress( "T_fatJet1_pt",                        &fatJet1_pt);
 T_ntuples->SetBranchAddress( "T_fatJet1_eta",                       &fatJet1_eta);
 T_ntuples->SetBranchAddress( "T_fatJet1_msoftdrop",                 &fatJet1_msoftdrop);
 T_ntuples->SetBranchAddress( "T_fatJet1_particleNet_XbbVsQCD",      &fatJet1_particleNet_XbbVsQCD);
 T_ntuples->SetBranchAddress( "T_fatJet1_particleNet_Xbb_Legacy",    &fatJet1_particleNet_Xbb_Legacy);
 T_ntuples->SetBranchAddress( "T_fatJet1_particleNet_Xbb_Legacy_AN", &fatJet1_particleNet_Xbb_Legacy_AN);
 T_ntuples->SetBranchAddress( "T_fatJet1_Tau3OverTau2",              &fatJet1_Tau3OverTau2);
 T_ntuples->SetBranchAddress( "T_fatJet2_pt",                        &fatJet2_pt);
 T_ntuples->SetBranchAddress( "T_fatJet2_eta",                       &fatJet2_eta);
 T_ntuples->SetBranchAddress( "T_fatJet2_msoftdrop",                 &fatJet2_msoftdrop);
 T_ntuples->SetBranchAddress( "T_MET",                               &MET);
 T_ntuples->SetBranchAddress( "T_Lep1_Pt",                           &Lep1_Pt);
 T_ntuples->SetBranchAddress( "T_Dr_LFJ",                            &Dr_LFJ);
 T_ntuples->SetBranchAddress( "T_Dr_J1FJ",                           &Dr_J1FJ);
 T_ntuples->SetBranchAddress( "T_Dr_J2FJ",                           &Dr_J2FJ);
 T_ntuples->SetBranchAddress( "T_Dr_JmaxL",                          &Dr_JmaxL);

 double pt_MIN=250.0;

 for(int i=0;i<T_ntuples->GetEntries();i++)
   {
    T_ntuples-> GetEntry(i);

    #include "Selection.txt"
     weight = weight*0.98; //normalizing to full 2022 (missed jobs)

    _FatJet1_pt->Fill(fatJet1_pt,weight);
    _FatJet1_eta->Fill(fatJet1_eta,weight);
    _FatJet1_MassSD->Fill(fatJet1_msoftdrop,weight);
    _FatJet1_Xbb->Fill(fatJet1_particleNet_XbbVsQCD,weight);
    _FatJet1_Xbb_Legacy->Fill(fatJet1_particleNet_Xbb_Legacy,weight);
    _FatJet1_Xbb_Legacy_AN->Fill(fatJet1_particleNet_Xbb_Legacy_AN,weight);
    _FatJet1_Tau3OverTau2->Fill(fatJet1_Tau3OverTau2,weight);
    _FatJet2_pt->Fill(fatJet2_pt,weight);
    _FatJet2_eta->Fill(fatJet2_eta,weight);
    _FatJet2_MassSD->Fill(fatJet2_msoftdrop,weight);
    _MET->Fill(MET,weight);
    _Lep1_Pt->Fill(Lep1_Pt,weight);
    _Dr_LFJ->Fill(Dr_LFJ,weight);
    _Dr_J1FJ->Fill(Dr_J1FJ,weight);
    _Dr_J2FJ->Fill(Dr_J2FJ,weight);
    _Dr_JmaxL->Fill(Dr_JmaxL,weight);

   }


  f->Write();

}

