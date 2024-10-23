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

void bdt()
{
   Float_t ptFJ12, etaFJ12, massFJ12, MET, fatJet1_Tau3OverTau2, fatJet2_Tau3OverTau2, fatJet1_msoftdrop, fatJet1_pt, fatJet1_eta, fatJet2_pt;
   Float_t fatJet1_particleNet_QCD, fatJet1_particleNet_QCD0HF, fatJet1_particleNet_QCD1HF, fatJet1_particleNet_QCD2HF, FJ1_PtOverMFJ12, FJ2_PtOverMFJ12, ptFJ1OverptFJ2;

   Int_t NJets, nBTaggedJets;
   Float_t alpha_FJ12, fatJet1_particleNetWithMass_TvsQCD, fatJet2_particleNetWithMass_TvsQCD;

   Float_t fatJet2_msoftdrop;

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   reader->AddVariable( "T_ptFJ12",                              &ptFJ12);
   reader->AddVariable( "T_etaFJ12",                             &etaFJ12);
   reader->AddVariable( "T_massFJ12",                            &massFJ12);
   reader->AddVariable( "T_MET",                                 &MET);
   reader->AddVariable( "T_fatJet1_Tau3OverTau2",                &fatJet1_Tau3OverTau2);
   reader->AddVariable( "T_fatJet2_Tau3OverTau2",                &fatJet2_Tau3OverTau2);
   reader->AddVariable( "T_fatJet1_msoftdrop",                   &fatJet1_msoftdrop);
   reader->AddVariable( "T_fatJet1_pt",                          &fatJet1_pt);
   reader->AddVariable( "T_fatJet1_eta",                         &fatJet1_eta);
   reader->AddVariable( "T_fatJet1_particleNet_QCD",             &fatJet1_particleNet_QCD);
   reader->AddVariable( "T_fatJet1_particleNet_QCD0HF",          &fatJet1_particleNet_QCD0HF);
   reader->AddVariable( "T_fatJet1_particleNet_QCD1HF",          &fatJet1_particleNet_QCD1HF);
   reader->AddVariable( "T_fatJet1_particleNet_QCD2HF",          &fatJet1_particleNet_QCD2HF);
   reader->AddVariable( "T_fatJet2_pt",                          &fatJet2_pt);
   reader->AddVariable( "T_FJ1_PtOverMFJ12",                     &FJ1_PtOverMFJ12);
   reader->AddVariable( "T_FJ2_PtOverMFJ12",                     &FJ2_PtOverMFJ12);
   reader->AddVariable( "T_ptFJ1OverptFJ2",                      &ptFJ1OverptFJ2);

   reader->AddVariable( "T_NJets",                               &NJets);
   reader->AddVariable( "T_nBTaggedJets",                        &nBTaggedJets);
   reader->AddVariable( "T_alpha_FJ12",                          &alpha_FJ12);
   reader->AddVariable( "T_fatJet1_particleNetWithMass_TvsQCD",  &fatJet1_particleNetWithMass_TvsQCD);
   reader->AddVariable( "T_fatJet2_particleNetWithMass_TvsQCD",  &fatJet2_particleNetWithMass_TvsQCD);

   reader->AddSpectator( "T_fatJet2_msoftdrop",                  &fatJet2_msoftdrop);

   reader->BookMVA("TMVAClassification_BDTG", "dataset/weights/TMVAClassification_BDTG.weights.xml" );

// ***************************************************************************************************************
// Signal
   TFile *f0 = TFile::Open("../../../NTuples/Tree_Signal.root");
   TTree * ntps = (TTree*)f0->Get("tree");

   Float_t S_weight;
   Float_t S_ptFJ12, S_etaFJ12, S_massFJ12, S_MET, S_fatJet1_Tau3OverTau2, S_fatJet2_Tau3OverTau2, S_fatJet1_msoftdrop, S_fatJet1_pt, S_fatJet1_eta, S_fatJet2_pt;
   Float_t S_fatJet1_particleNet_QCD, S_fatJet1_particleNet_QCD0HF, S_fatJet1_particleNet_QCD1HF, S_fatJet1_particleNet_QCD2HF, S_FJ1_PtOverMFJ12, S_FJ2_PtOverMFJ12, S_ptFJ1OverptFJ2;

   Int_t S_NJets, S_nBTaggedJets;
   Float_t S_alpha_FJ12, S_fatJet1_particleNetWithMass_TvsQCD, S_fatJet2_particleNetWithMass_TvsQCD;
   Float_t S_fatJet2_msoftdrop;


   ntps->SetBranchAddress( "T_weight",                              &S_weight);
   ntps->SetBranchAddress( "T_ptFJ12",                              &S_ptFJ12);
   ntps->SetBranchAddress( "T_etaFJ12",                             &S_etaFJ12);
   ntps->SetBranchAddress( "T_massFJ12",                            &S_massFJ12);
   ntps->SetBranchAddress( "T_MET",                                 &S_MET);
   ntps->SetBranchAddress( "T_fatJet1_Tau3OverTau2",                &S_fatJet1_Tau3OverTau2);
   ntps->SetBranchAddress( "T_fatJet2_Tau3OverTau2",                &S_fatJet2_Tau3OverTau2);
   ntps->SetBranchAddress( "T_fatJet1_msoftdrop",                   &S_fatJet1_msoftdrop);
   ntps->SetBranchAddress( "T_fatJet1_pt",                          &S_fatJet1_pt);
   ntps->SetBranchAddress( "T_fatJet1_eta",                         &S_fatJet1_eta);
   ntps->SetBranchAddress( "T_fatJet2_pt",                          &S_fatJet2_pt);
   ntps->SetBranchAddress( "T_fatJet1_particleNet_QCD",             &S_fatJet1_particleNet_QCD);
   ntps->SetBranchAddress( "T_fatJet1_particleNet_QCD0HF",          &S_fatJet1_particleNet_QCD0HF);
   ntps->SetBranchAddress( "T_fatJet1_particleNet_QCD1HF",          &S_fatJet1_particleNet_QCD1HF);
   ntps->SetBranchAddress( "T_fatJet1_particleNet_QCD2HF",          &S_fatJet1_particleNet_QCD2HF);
   ntps->SetBranchAddress( "T_FJ1_PtOverMFJ12",                     &S_FJ1_PtOverMFJ12);
   ntps->SetBranchAddress( "T_FJ2_PtOverMFJ12",                     &S_FJ2_PtOverMFJ12);
   ntps->SetBranchAddress( "T_ptFJ1OverptFJ2",                      &S_ptFJ1OverptFJ2);

   ntps->SetBranchAddress( "T_NJets",                               &S_NJets);
   ntps->SetBranchAddress( "T_nBTaggedJets",                        &S_nBTaggedJets);
   ntps->SetBranchAddress( "T_alpha_FJ12",                          &S_alpha_FJ12);
   ntps->SetBranchAddress( "T_fatJet1_particleNetWithMass_TvsQCD",  &S_fatJet1_particleNetWithMass_TvsQCD);
   ntps->SetBranchAddress( "T_fatJet2_particleNetWithMass_TvsQCD",  &S_fatJet2_particleNetWithMass_TvsQCD);

   ntps->SetBranchAddress( "T_fatJet2_msoftdrop",                   &S_fatJet2_msoftdrop);


// QCD
   TFile *f1 = TFile::Open("../../../NTuples/Tree_QCD.root");
// TTbar
   TFile *f2 = TFile::Open("../../../NTuples/Tree_TTbar_Corr.root");
// VJ
   TFile *f3 = TFile::Open("../../../NTuples/Tree_VJ.root");
// VV
   TFile *f4 = TFile::Open("../../../NTuples/Tree_VV.root");
// Data
   TFile *f5 = TFile::Open("../../../NTuples/Data_2022_PostEE.root");

   TTree * ntpb[] =  {(TTree*)f1->Get("tree"), (TTree*)f2->Get("tree"), (TTree*)f3->Get("tree"), (TTree*)f4->Get("tree"), (TTree*)f5->Get("tree")};

   int nBProc = sizeof(ntpb)/sizeof(ntpb[0]);

   Float_t B_weight[nBProc];
   Float_t B_ptFJ12[nBProc], B_etaFJ12[nBProc], B_massFJ12[nBProc], B_MET[nBProc], B_fatJet1_Tau3OverTau2[nBProc], B_fatJet2_Tau3OverTau2[nBProc], B_fatJet1_msoftdrop[nBProc], B_fatJet1_pt[nBProc], B_fatJet1_eta[nBProc], B_fatJet2_pt[nBProc];
   Float_t B_fatJet1_particleNet_QCD[nBProc], B_fatJet1_particleNet_QCD0HF[nBProc], B_fatJet1_particleNet_QCD1HF[nBProc], B_fatJet1_particleNet_QCD2HF[nBProc], B_FJ1_PtOverMFJ12[nBProc], B_FJ2_PtOverMFJ12[nBProc], B_ptFJ1OverptFJ2[nBProc];

   Int_t B_NJets[nBProc], B_nBTaggedJets[nBProc];
   Float_t B_alpha_FJ12[nBProc], B_fatJet1_particleNetWithMasB_TvsQCD[nBProc], B_fatJet2_particleNetWithMasB_TvsQCD[nBProc];
   Float_t B_fatJet2_msoftdrop[nBProc];


   for(int k=0;k<nBProc;k++)
   {
    ntpb[k]->SetBranchAddress( "T_weight",                              &B_weight[k]);
    ntpb[k]->SetBranchAddress( "T_ptFJ12",                              &B_ptFJ12[k]);
    ntpb[k]->SetBranchAddress( "T_etaFJ12",                             &B_etaFJ12[k]);
    ntpb[k]->SetBranchAddress( "T_massFJ12",                            &B_massFJ12[k]);
    ntpb[k]->SetBranchAddress( "T_MET",                                 &B_MET[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_Tau3OverTau2",                &B_fatJet1_Tau3OverTau2[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet2_Tau3OverTau2",                &B_fatJet2_Tau3OverTau2[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_msoftdrop",                   &B_fatJet1_msoftdrop[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_pt",                          &B_fatJet1_pt[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_eta",                         &B_fatJet1_eta[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet2_pt",                          &B_fatJet2_pt[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_particleNet_QCD",             &B_fatJet1_particleNet_QCD[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_particleNet_QCD0HF",          &B_fatJet1_particleNet_QCD0HF[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_particleNet_QCD1HF",          &B_fatJet1_particleNet_QCD1HF[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet1_particleNet_QCD2HF",          &B_fatJet1_particleNet_QCD2HF[k]);
    ntpb[k]->SetBranchAddress( "T_FJ1_PtOverMFJ12",                     &B_FJ1_PtOverMFJ12);
    ntpb[k]->SetBranchAddress( "T_FJ2_PtOverMFJ12",                     &B_FJ2_PtOverMFJ12[k]);
    ntpb[k]->SetBranchAddress( "T_ptFJ1OverptFJ2",                      &B_ptFJ1OverptFJ2[k]);

    ntpb[k]->SetBranchAddress( "T_NJets",                               &B_NJets[k]);
    ntpb[k]->SetBranchAddress( "T_nBTaggedJets",                        &B_nBTaggedJets[k]);
    ntpb[k]->SetBranchAddress( "T_alpha_FJ12",                          &B_alpha_FJ12[k]);


                                                                         B_fatJet1_particleNetWithMasB_TvsQCD
    ntpb[k]->SetBranchAddress( "T_fatJet1_particleNetWithMasB_TvsQCD",  &B_fatJet1_particleNetWithMass_TvsQCD[k]);
    ntpb[k]->SetBranchAddress( "T_fatJet2_particleNetWithMass_TvsQCD",  &B_fatJet2_particleNetWithMass_TvsQCD[k]);




    ntpb[k]->SetBranchAddress( "T_fatJet2_msoftdrop",                   &B_fatJet2_msoftdrop[k]);
   }

//***************************************************************************************
  TFile *f = new TFile("BDT_Output.root","RECREATE");
  gStyle->SetOptFile(0);
  gStyle->SetOptStat("mre");
  gStyle->SetPaintTextFormat("1.2e");

  double xmin=-1;
  double xmax=1;


  TH1D * _BDT_Signal  = new TH1D("BDT_Signal","BDT_Signal",100,xmin,xmax);
  TH1D * _BDT_QCD     = new TH1D("BDT_QCD","BDT_QCD",100,xmin,xmax);
  TH1D * _BDT_tt      = new TH1D("BDT_tt","BDT_tt",100,xmin,xmax);
  TH1D * _BDT_VJ      = new TH1D("BDT_VJ","BDT_VJ",100,xmin,xmax);
  TH1D * _BDT_VV      = new TH1D("BDT_VV","BDT_VV",100,xmin,xmax);
  TH1D * _BDT_DATA    = new TH1D("BDT_DATA","BDT_DATA",100,xmin,xmax);

//***************************************************************************************

   for(int i=0;i<ntps->GetEntries();i++)
    {
     ntps->GetEntry(i);

     ptFJ12 = S_ptFJ12; etaFJ12 = S_etaFJ12; massFJ12 = S_massFJ12; MET = S_MET; fatJet1_Tau3OverTau2 = S_fatJet1_Tau3OverTau2; fatJet2_Tau3OverTau2 = S_fatJet2_Tau3OverTau2; 
     fatJet1_msoftdrop = S_fatJet1_msoftdrop; fatJet1_pt = S_fatJet1_pt; fatJet1_eta = S_fatJet1_eta; fatJet2_pt = S_fatJet2_pt;
     fatJet1_particleNet_QCD = S_fatJet1_particleNet_QCD; fatJet1_particleNet_QCD0HF = S_fatJet1_particleNet_QCD0HF; fatJet1_particleNet_QCD1HF = S_fatJet1_particleNet_QCD1HF;
     fatJet1_particleNet_QCD2HF = S_fatJet1_particleNet_QCD2HF; FJ1_PtOverMFJ12 = S_FJ1_PtOverMFJ12; FJ2_PtOverMFJ12 = S_FJ2_PtOverMFJ12; ptFJ1OverptFJ2 = S_ptFJ1OverptFJ2;

     NJets = S_NJets; nBTaggedJets = S_nBTaggedJets; alpha_FJ12 = S_alpha_FJ12; fatJet1_particleNetWithMass_TvsQCD = S_fatJet1_particleNetWithMass_TvsQCD; fatJet2_particleNetWithMass_TvsQCD = S_fatJet2_particleNetWithMass_TvsQCD;

     fatJet2_msoftdrop = S_fatJet2_msoftdrop;

     double bdt_out=reader->EvaluateMVA("TMVAClassification_BDTG");

    _BDT_Signal->Fill(reader->EvaluateMVA("TMVAClassification_BDTG"),S_weight);
    }

  for(int k=0;k<nBProc;k++)
   for(int i=0;i<ntpb[k]->GetEntries();i++)
    {

     ntpb[k]->GetEntry(i);

     ptFJ12 = B_ptFJ12[k]; etaFJ12 = B_etaFJ12[k]; massFJ12 = B_massFJ12[k]; MET = B_MET[k]; fatJet1_Tau3OverTau2 = B_fatJet1_Tau3OverTau2[k]; fatJet2_Tau3OverTau2 = B_fatJet2_Tau3OverTau2[k]; 
     fatJet1_msoftdrop = B_fatJet1_msoftdrop[k]; fatJet1_pt = B_fatJet1_pt[k]; fatJet1_eta = B_fatJet1_eta[k]; fatJet2_pt = B_fatJet2_pt[k];
     fatJet1_particleNet_QCD = B_fatJet1_particleNet_QCD[k]; fatJet1_particleNet_QCD0HF = B_fatJet1_particleNet_QCD0HF[k]; fatJet1_particleNet_QCD1HF = B_fatJet1_particleNet_QCD1HF[k];
     fatJet1_particleNet_QCD2HF = B_fatJet1_particleNet_QCD2HF[k]; FJ1_PtOverMFJ12 = B_FJ1_PtOverMFJ12[k]; FJ2_PtOverMFJ12 = B_FJ2_PtOverMFJ12[k]; ptFJ1OverptFJ2 = B_ptFJ1OverptFJ2[k];

     NJets = B_NJets[k]; nBTaggedJets = B_nBTaggedJets[k]; alpha_FJ12 = B_alpha_FJ12[k]; fatJet1_particleNetWithMass_TvsQCD = B_fatJet1_particleNetWithMass_TvsQCD[k]; fatJet2_particleNetWithMass_TvsQCD = B_fatJet2_particleNetWithMass_TvsQCD[k];

     fatJet2_msoftdrop = B_fatJet2_msoftdrop[k];

     if(k==0)  // exclude QCD_100to200
      _BDT_QCD->Fill(reader->EvaluateMVA("TMVAClassification_BDTG"),B_weight[k]);
     else if(k==1)
      _BDT_tt->Fill(reader->EvaluateMVA("TMVAClassification_BDTG"),B_weight[k]);
     else if(k==2)
      _BDT_VJ->Fill(reader->EvaluateMVA("TMVAClassification_BDTG"),B_weight[k]);
     else if(k==3)
      _BDT_VV->Fill(reader->EvaluateMVA("TMVAClassification_BDTG"),B_weight[k]);
     else if(k==4)
      _BDT_DATA->Fill(reader->EvaluateMVA("TMVAClassification_BDTG"),B_weight[k]);
    }

  f->Write();

  delete f0;
  delete f1;
  delete f2;
  delete f3;
  delete f4;
  delete f5;
}