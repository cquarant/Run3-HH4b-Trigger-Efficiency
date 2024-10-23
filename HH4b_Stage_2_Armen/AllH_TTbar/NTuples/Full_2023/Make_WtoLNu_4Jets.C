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


void Make_WtoLNu_4Jets()
{
 TFile *f = new TFile("Histograms_WtoLNu_4Jets.root","RECREATE");

 TH1D *_FatJet1_pt             = new TH1D("FatJet1_pt","FatJet1_pt",100,0,1000);
 TH1D *_FatJet1_eta            = new TH1D("FatJet1_eta","FatJet1_eta",100,-5,5);
 TH1D *_FatJet1_MassSD         = new TH1D("FatJet1_MassSD","FatJet1_MassSD",500,0,500);
 TH1D *_FatJet1_Xbb            = new TH1D("FatJet1_Xbb","FatJet1_Xbb",100,0.0,1.0);
 TH1D *_FatJet1_Xbb_Legacy     = new TH1D("FatJet1_Xbb_Legacy","FatJet1_Xbb_Legacy",100,0.0,1.0);
 TH1D *_FatJet1_Xbb_Legacy_AN  = new TH1D("FatJet1_Xbb_Legacy_AN","FatJet1_Xbb_Legacy_AN",100,0.0,1.0);
 TH1D *_FatJet1_Tau3OverTau2   = new TH1D("FatJet1_Tau3OverTau2","FatJet1_Tau3OverTau2",100,0.0,1.0);

 TH1D *_FatJet2_pt             = new TH1D("FatJet2_pt","FatJet2_pt",100,0,1000);
 TH1D *_FatJet2_eta            = new TH1D("FatJet2_eta","FatJet2_eta",100,-5,5);
 TH1D *_FatJet2_MassSD         = new TH1D("FatJet2_MassSD","FatJet2_MassSD",500,0,500);
 TH1D *_FatJet2_Xbb            = new TH1D("FatJet2_Xbb","FatJet2_Xbb",100,0.0,1.0);
 TH1D *_FatJet2_Xbb_Legacy     = new TH1D("FatJet2_Xbb_Legacy","FatJet2_Xbb_Legacy",100,0.0,1.0);
 TH1D *_FatJet2_Xbb_Legacy_AN  = new TH1D("FatJet2_Xbb_Legacy_AN","FatJet2_Xbb_Legacy_AN",100,0.0,1.0);
 TH1D *_FatJet2_Tau3OverTau2   = new TH1D("FatJet2_Tau3OverTau2","FatJet2_Tau3OverTau2",100,0.0,1.0);

 TH1D *_MET                    = new TH1D("MET","MET",100,0,500);
 TH1D *_PTjj                   = new TH1D("PTjj","PTjj",100,0,1000);

 TFile *f0 = TFile::Open("Tree_WtoLNu_4Jets.root");
 TTree * T_ntuples    = (TTree*)f0->Get("tree");

 Float_t weight;
 Float_t fatJet1_pt, fatJet1_eta, fatJet1_phi,fatJet1_msoftdrop;
 Float_t fatJet1_particleNet_XbbVsQCD, fatJet1_particleNet_Xbb_Legacy, fatJet1_particleNet_Xbb_Legacy_AN, fatJet1_Tau3OverTau2;
 Float_t fatJet2_pt, fatJet2_eta, fatJet2_phi,fatJet2_msoftdrop;
 Float_t fatJet2_particleNet_XbbVsQCD, fatJet2_particleNet_Xbb_Legacy, fatJet2_particleNet_Xbb_Legacy_AN, fatJet2_Tau3OverTau2;

 Float_t MET, PTjj;

 T_ntuples->SetBranchAddress( "T_weight",                            &weight);
 T_ntuples->SetBranchAddress( "T_fatJet1_pt",                        &fatJet1_pt);
 T_ntuples->SetBranchAddress( "T_fatJet1_eta",                       &fatJet1_eta);
 T_ntuples->SetBranchAddress( "T_fatJet1_phi",                       &fatJet1_phi);
 T_ntuples->SetBranchAddress( "T_fatJet1_msoftdrop",                 &fatJet1_msoftdrop);
 T_ntuples->SetBranchAddress( "T_fatJet1_particleNet_XbbVsQCD",      &fatJet1_particleNet_XbbVsQCD);
 T_ntuples->SetBranchAddress( "T_fatJet1_particleNet_Xbb_Legacy",    &fatJet1_particleNet_Xbb_Legacy);
 T_ntuples->SetBranchAddress( "T_fatJet1_particleNet_Xbb_Legacy_AN", &fatJet1_particleNet_Xbb_Legacy_AN);
 T_ntuples->SetBranchAddress( "T_fatJet1_Tau3OverTau2",              &fatJet1_Tau3OverTau2);

 T_ntuples->SetBranchAddress( "T_fatJet2_pt",                        &fatJet2_pt);
 T_ntuples->SetBranchAddress( "T_fatJet2_eta",                       &fatJet2_eta);
 T_ntuples->SetBranchAddress( "T_fatJet2_phi",                       &fatJet2_phi);
 T_ntuples->SetBranchAddress( "T_fatJet2_msoftdrop",                 &fatJet2_msoftdrop);
 T_ntuples->SetBranchAddress( "T_fatJet2_particleNet_XbbVsQCD",      &fatJet2_particleNet_XbbVsQCD);
 T_ntuples->SetBranchAddress( "T_fatJet2_particleNet_Xbb_Legacy",    &fatJet2_particleNet_Xbb_Legacy);
 T_ntuples->SetBranchAddress( "T_fatJet2_particleNet_Xbb_Legacy_AN", &fatJet2_particleNet_Xbb_Legacy_AN);
 T_ntuples->SetBranchAddress( "T_fatJet2_Tau3OverTau2",              &fatJet2_Tau3OverTau2);

 T_ntuples->SetBranchAddress( "T_MET",                               &MET);
 T_ntuples->SetBranchAddress( "T_PTjj",                              &PTjj);

 double pt_MIN=250.0;

 for(int i=0;i<T_ntuples->GetEntries();i++)
   {
    T_ntuples-> GetEntry(i);

     #include "Selection.txt"
     weight = weight*21.7/27.1; //normalizing to full 2022 (missed PreEE)

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
    _FatJet2_Xbb->Fill(fatJet2_particleNet_XbbVsQCD,weight);
    _FatJet2_Xbb_Legacy->Fill(fatJet2_particleNet_Xbb_Legacy,weight);
    _FatJet2_Xbb_Legacy_AN->Fill(fatJet2_particleNet_Xbb_Legacy_AN,weight);
    _FatJet2_Tau3OverTau2->Fill(fatJet2_Tau3OverTau2,weight);

    _MET->Fill(MET,weight);
    _PTjj->Fill(PTjj,weight);

   }


  f->Write();

}

