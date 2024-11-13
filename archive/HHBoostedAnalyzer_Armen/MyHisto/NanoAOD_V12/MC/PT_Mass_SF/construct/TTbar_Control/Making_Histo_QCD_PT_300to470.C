#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include "TCanvas.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <vector>
#include <map>

// D-phi
 double phi_dist(double a, double b){
   if(fabs(a - b) > 3.14159265)
   {
    return 6.2831853 - fabs(a - b);
   }
   return fabs(a - b);
 }

bool inRange(int low, int high, int x)
{
    return (low <= x && x <= high);
}

double PU_Rew[100]={10.8497,4.42194,4.65723,20.3724,23.8022,13.7557,7.99435,4.19632,2.25956,1.1874,0.760477,0.484277,0.24048,0.172868,0.143258,0.064479,0.0499436,0.0266,0.0319528,0.0192531,0.0837497,0.0850733,0.15549,0.117273,0.130173,0.138538,0.131814,0.178734,0.249517,0.380582,0.413954,0.464502,0.464131,0.471535,0.53024,0.60752,0.571946,0.566939,0.480804,0.563873,0.512671,0.544993,0.539009,0.547109,0.576103,0.645572,0.719275,0.761689,0.874711,1.21862,1.51221,2.36156,3.95115,9.95135,16.477,18.3475,10.4246,3.23251,1.75302,2.70381,4.29273,6.57119,4.28892,0.986744,1.64368,5.21225,9.35711,9.75092,10.751,16.5981,30.7719,50.8975,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

// ******************************************

#include "/afs/cern.ch/work/t/tumasyan/HHTo4B/CMSSW_9_4_2/src/HHBoostedAnalyzer/MyHisto/MC/postEE/parameters.txt"

void Making_Histo_QCD_PT_300to470()
{
 TFile *f = new TFile("Histograms_QCD_PT_300to470.root","RECREATE");

TH1D *_FatJet1_pt  = new TH1D("FatJet1_pt","FatJet1_pt",200,0,1000);
TH1D *_FatJet1_eta  = new TH1D("FatJet1_eta","FatJet1_eta",100,-5,5);
TH1D *_FatJet1_phi  = new TH1D("FatJet1_phi","FatJet1_phi",100,-5,5);
TH2D *_FatJet1_eta_phi = new TH2D("FatJet1_eta_phi","FatJet1_eta_phi",100,-5,5,100,-5,5);
TH1D *_FatJet1_Mass  = new TH1D("FatJet1_Mass","FatJet1_Mass",500,0,500);
TH1D *_FatJet1_MassSD  = new TH1D("FatJet1_MassSD","FatJet1_MassSD",500,0,500);
TH1D *_FatJet1DDBTaggerV2  = new TH1D("FatJet1DDBTaggerV2","FatJet1DDBTaggerV2",100,0,1.0);
TH1D *_FatJet1PNetMD_QCD  = new TH1D("FatJet1PNetMD_QCD","FatJet1PNetMD_QCD",100,0,1.0);
TH1D *_FatJet1PNetMD_Xbb  = new TH1D("FatJet1PNetMD_Xbb","FatJet1PNetMD_Xbb",100,0,1.0);
TH1D *_FatJet1PNetMD_Xcc  = new TH1D("FatJet1PNetMD_Xcc","FatJet1PNetMD_Xcc",100,0,1.0);
TH1D *_FatJet1PNetMD_Xqq  = new TH1D("FatJet1PNetMD_Xqq","FatJet1PNetMD_Xqq",100,0,1.0);
TH1D *_FatJet1PNetHbbvsQCD  = new TH1D("FatJet1PNetHbbvsQCD","FatJet1PNetHbbvsQCD",100,0,1.0);
TH1D *_FatJet1PNet_QCD  = new TH1D("FatJet1PNet_QCD","FatJet1PNet_QCD",100,0,1.0);
TH1D *_FatJet1PNet_TvsQCD  = new TH1D("FatJet1PNet_TvsQCD","FatJet1PNet_TvsQCD",100,0,1.0);
TH1D *_FatJet1PNet_WvsQCD  = new TH1D("FatJet1PNet_WvsQCD","FatJet1PNet_WvsQCD",100,0,1.0);
TH1D *_FatJet1PNet_ZvsQCD  = new TH1D("FatJet1PNet_ZvsQCD","FatJet1PNet_ZvsQCD",100,0,1.0);
TH1D *_FatJet1PNetMD_Xbb_Modified_as_in_AN_Xbb  = new TH1D("FatJet1PNetMD_Xbb_Modified_as_in_AN_Xbb","FatJet1PNetMD_Xbb_Modified_as_in_AN_Xbb",100,0,1.0);
TH1D *_FatJet1PNet_mass  = new TH1D("FatJet1PNet_mass","FatJet1PNet_mass",500,0,500);
TH1D *_FatJet1_Tau3OverTau2  = new TH1D("FatJet1_Tau3OverTau2","FatJet1_Tau3OverTau2",100,0.0,1.0);

TH1D *_FatJet2_pt  = new TH1D("FatJet2_pt","FatJet2_pt",200,0,1000);
TH1D *_FatJet2_eta  = new TH1D("FatJet2_eta","FatJet2_eta",100,-5,5);
TH1D *_FatJet2_phi  = new TH1D("FatJet2_phi","FatJet2_phi",100,-5,5);
TH2D *_FatJet2_eta_phi = new TH2D("FatJet2_eta_phi","FatJet2_eta_phi",100,-5,5,100,-5,5);
TH1D *_FatJet2_Mass  = new TH1D("FatJet2_Mass","FatJet2_Mass",500,0,500);
TH1D *_FatJet2_MassSD  = new TH1D("FatJet2_MassSD","FatJet2_MassSD",500,0,500);
TH1D *_FatJet2DDBTaggerV2  = new TH1D("FatJet2DDBTaggerV2","FatJet2DDBTaggerV2",100,0,1.0);
TH1D *_FatJet2PNetMD_QCD  = new TH1D("FatJet2PNetMD_QCD","FatJet2PNetMD_QCD",100,0,1.0);
TH1D *_FatJet2PNetMD_Xbb  = new TH1D("FatJet2PNetMD_Xbb","FatJet2PNetMD_Xbb",100,0,1.0);
TH1D *_FatJet2PNetMD_Xcc  = new TH1D("FatJet2PNetMD_Xcc","FatJet2PNetMD_Xcc",100,0,1.0);
TH1D *_FatJet2PNetMD_Xqq  = new TH1D("FatJet2PNetMD_Xqq","FatJet2PNetMD_Xqq",100,0,1.0);
TH1D *_FatJet2PNetHbbvsQCD  = new TH1D("FatJet2PNetHbbvsQCD","FatJet2PNetHbbvsQCD",100,0,1.0);
TH1D *_FatJet2PNet_QCD  = new TH1D("FatJet2PNet_QCD","FatJet2PNet_QCD",100,0,1.0);
TH1D *_FatJet2PNet_TvsQCD  = new TH1D("FatJet2PNet_TvsQCD","FatJet2PNet_TvsQCD",100,0,1.0);
TH1D *_FatJet2PNet_WvsQCD  = new TH1D("FatJet2PNet_WvsQCD","FatJet2PNet_WvsQCD",100,0,1.0);
TH1D *_FatJet2PNet_ZvsQCD  = new TH1D("FatJet2PNet_ZvsQCD","FatJet2PNet_ZvsQCD",100,0,1.0);
TH1D *_FatJet2PNetMD_Xbb_Modified_as_in_AN_Xbb  = new TH1D("FatJet2PNetMD_Xbb_Modified_as_in_AN_Xbb","FatJet2PNetMD_Xbb_Modified_as_in_AN_Xbb",100,0,1.0);
TH1D *_FatJet2PNet_mass  = new TH1D("FatJet2PNet_mass","FatJet2PNet_mass",500,0,500);
TH1D *_FatJet2_Tau3OverTau2  = new TH1D("FatJet2_Tau3OverTau2","FatJet2_Tau3OverTau2",100,0.0,1.0);

TH1D *_MET = new TH1D("MET","MET",100,0,500);

/*
Float_t Lower_m[18]   = {0,50,100,150,200,250,300,350,400,450,500,600,700,800,1000,1200,1600,2000};
Float_t Lower_eta[25] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.2,4.6,5.0,5.5,7.0};
TH2D *_VBF_mINV_dETA = new TH2D("UnOrd_VBF_mINV_dETA","UnOrd_VBF_mINV_dETA",24,Lower_eta,17,Lower_m);

Float_t Lower_btg1[] ={0.0, 0.05, 0.10, 0.15, 0.20, 0.28, 0.36, 0.46, 0.60, 0.70, 0.80, 0.85, 0.90, 0.93, 0.97, 1.00};
Float_t Lower_btg2[] ={0.0, 0.05, 0.10, 0.15, 0.20, 0.28, 0.36, 0.46, 0.60, 0.70, 0.80, 0.85, 0.90, 0.93, 0.97, 1.00};
TH2D * _BTG1_BTG2 = new TH2D("UnOrd_BTG1_BTG2","UnOrd_BTG1_BTG2",15,Lower_btg1,15,Lower_btg2);
TH2D * _HBTG_btg_DeepJet_matchedtoL3BTag_2D = new TH2D("HBTG_btg_DeepJet_matchedtoL3BTag_2D","HBTG_btg_DeepJet_matchedtoL3BTag_2D",15,Lower_btg1,15,Lower_btg2);
*/

TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/MC/postEE_SFs/QCD_PT_300to470.root");

TH1F  *NEvents = (TH1F*)f1->Get("NEvents");
double SumGenWeights = NEvents->GetBinContent(1);
TTree *InputTree = (TTree*)f1->Get("tree");

Float_t    weight;
UInt_t     run;
UInt_t     lumi;
Float_t    npu;
Int_t      isVBFtag;

Bool_t     HLT_Ele32_WPTight_Gsf;
Bool_t     HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;

Bool_t     HLT_Mu50;
Bool_t     HLT_IsoMu27;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;

Bool_t     HLT_AK8PFJet230_SoftDropMass40;

Float_t    MET;

Float_t    lep1_Pt;
Float_t    lep1_Eta;
Float_t    lep1_Phi;
Int_t      lep1_Id;

Float_t    lep2_Pt;
Float_t    lep2_Eta;
Float_t    lep2_Phi;
Int_t      lep2_Id;

Float_t    FatJet1_pt;
Float_t    FatJet1_eta;
Float_t    FatJet1_phi;
Float_t    FatJet1_Mass;
Float_t    FatJet1_MassSD;
Float_t    FatJet1DDBTaggerV2;
Float_t    FatJet1PNetMD_QCD;
Float_t    FatJet1PNetMD_Xbb;
Float_t    FatJet1PNetMD_Xcc;
Float_t    FatJet1PNetMD_Xqq;
Float_t    FatJet1PNetHbbvsQCD;
Float_t    FatJet1PNet_QCD;
Float_t    FatJet1PNet_TvsQCD;
Float_t    FatJet1PNet_WvsQCD;
Float_t    FatJet1PNet_ZvsQCD;
Float_t    FatJet1PNet_mass;
Float_t    FatJet1Tau3OverTau2;

Float_t    FatJet2_pt;
Float_t    FatJet2_eta;
Float_t    FatJet2_phi;
Float_t    FatJet2_Mass;
Float_t    FatJet2_MassSD;
Float_t    FatJet2DDBTaggerV2;
Float_t    FatJet2PNetMD_QCD;
Float_t    FatJet2PNetMD_Xbb;
Float_t    FatJet2PNetMD_Xcc;
Float_t    FatJet2PNetMD_Xqq;
Float_t    FatJet2PNetHbbvsQCD;
Float_t    FatJet2PNet_QCD;
Float_t    FatJet2PNet_TvsQCD;
Float_t    FatJet2PNet_WvsQCD;
Float_t    FatJet2PNet_ZvsQCD;
Float_t    FatJet2PNet_mass;
Float_t    FatJet2Tau3OverTau2;

Float_t    FatJet3_pt;
Float_t    FatJet3_MassSD;


InputTree->SetBranchAddress("weight",&weight);
InputTree->SetBranchAddress("run",&run);
InputTree->SetBranchAddress("lumi",&lumi);
InputTree->SetBranchAddress("npu",&npu);
InputTree->SetBranchAddress("isVBFtag",&isVBFtag);

InputTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35);

InputTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
InputTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35);

InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);

InputTree->SetBranchAddress("lep1Pt",&lep1_Pt);
InputTree->SetBranchAddress("lep1Eta",&lep1_Eta);
InputTree->SetBranchAddress("lep1Phi",&lep1_Phi);
InputTree->SetBranchAddress("lep1Id",&lep1_Id);

InputTree->SetBranchAddress("lep2Pt",&lep2_Pt);
InputTree->SetBranchAddress("lep2Eta",&lep2_Eta);
InputTree->SetBranchAddress("lep2Phi",&lep2_Phi);
InputTree->SetBranchAddress("lep2Id",&lep2_Id);

InputTree->SetBranchAddress("MET",&MET);

InputTree->SetBranchAddress("fatJet1Pt",&FatJet1_pt);
InputTree->SetBranchAddress("fatJet1Eta",&FatJet1_eta);
InputTree->SetBranchAddress("fatJet1Phi",&FatJet1_phi);
InputTree->SetBranchAddress("fatJet1Mass",&FatJet1_Mass);
InputTree->SetBranchAddress("fatJet1MassSD",&FatJet1_MassSD);
InputTree->SetBranchAddress("fatJet1DDBTaggerV2",&FatJet1DDBTaggerV2);
InputTree->SetBranchAddress("fatJet1PNetMD_QCD",&FatJet1PNetMD_QCD);
InputTree->SetBranchAddress("fatJet1PNetMD_Xbb",&FatJet1PNetMD_Xbb);
InputTree->SetBranchAddress("fatJet1PNetMD_Xcc",&FatJet1PNetMD_Xcc);
InputTree->SetBranchAddress("fatJet1PNetMD_Xqq",&FatJet1PNetMD_Xqq);
InputTree->SetBranchAddress("fatJet1PNetHbbvsQCD",&FatJet1PNetHbbvsQCD);
InputTree->SetBranchAddress("fatJet1PNet_QCD",&FatJet1PNet_QCD);
InputTree->SetBranchAddress("fatJet1PNet_TvsQCD",&FatJet1PNet_TvsQCD);
InputTree->SetBranchAddress("fatJet1PNet_WvsQCD",&FatJet1PNet_WvsQCD);
InputTree->SetBranchAddress("fatJet1PNet_ZvsQCD",&FatJet1PNet_ZvsQCD);
InputTree->SetBranchAddress("fatJet1PNet_mass",&FatJet1PNet_mass);
InputTree->SetBranchAddress("fatJet1Tau3OverTau2",&FatJet1Tau3OverTau2);

InputTree->SetBranchAddress("fatJet2Pt",&FatJet2_pt);
InputTree->SetBranchAddress("fatJet2Eta",&FatJet2_eta);
InputTree->SetBranchAddress("fatJet2Phi",&FatJet2_phi);
InputTree->SetBranchAddress("fatJet2Mass",&FatJet2_Mass);
InputTree->SetBranchAddress("fatJet2MassSD",&FatJet2_MassSD);
InputTree->SetBranchAddress("fatJet2DDBTaggerV2",&FatJet2DDBTaggerV2);
InputTree->SetBranchAddress("fatJet2PNetMD_QCD",&FatJet2PNetMD_QCD);
InputTree->SetBranchAddress("fatJet2PNetMD_Xbb",&FatJet2PNetMD_Xbb);
InputTree->SetBranchAddress("fatJet2PNetMD_Xcc",&FatJet2PNetMD_Xcc);
InputTree->SetBranchAddress("fatJet2PNetMD_Xqq",&FatJet2PNetMD_Xqq);
InputTree->SetBranchAddress("fatJet2PNetHbbvsQCD",&FatJet2PNetHbbvsQCD);
InputTree->SetBranchAddress("fatJet2PNet_QCD",&FatJet2PNet_QCD);
InputTree->SetBranchAddress("fatJet2PNet_TvsQCD",&FatJet2PNet_TvsQCD);
InputTree->SetBranchAddress("fatJet2PNet_WvsQCD",&FatJet2PNet_WvsQCD);
InputTree->SetBranchAddress("fatJet2PNet_ZvsQCD",&FatJet2PNet_ZvsQCD);
InputTree->SetBranchAddress("fatJet2PNet_mass",&FatJet2PNet_mass);
InputTree->SetBranchAddress("fatJet2Tau3OverTau2",&FatJet2Tau3OverTau2);

InputTree->SetBranchAddress("fatJet3Pt",&FatJet3_pt);
InputTree->SetBranchAddress("fatJet3MassSD",&FatJet3_MassSD);

  // Events Loop
  for(int i=0;i<InputTree->GetEntries();i++)
   {
    InputTree->GetEntry(i);
// ********************************************************** HLT Selection
//   if(HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35==0) continue;
//   if(HLT_AK8PFJet230_SoftDropMass40==0) continue;

// ********************************************************** FatJets Selection
//    if(FatJet1_pt<300 || FatJet2_pt<300 || FatJet1_MassSD<40 || FatJet2_MassSD<40) continue;
//    if(FatJet1_pt<250 || FatJet2_pt<250 || FatJet1_MassSD<50 || FatJet2_MassSD<50 || FatJet1PNetMD_Xbb<0.8) continue;
//    if(FatJet1_pt<300 || FatJet2_pt<300 || FatJet1_MassSD<40 || FatJet2_MassSD<40 || FatJet1PNetMD_Xbb<0.8) continue;
//    if(FatJet1_pt<300 || FatJet2_pt<300 || FatJet1_MassSD<50 || FatJet2_MassSD<50 || FatJet1PNetMD_Xbb<0.8) continue;
    // Control region
//    if(FatJet1_MassSD<105 || FatJet1_MassSD>135) continue;
//    if(FatJet2_MassSD>95  && FatJet2_MassSD<135) continue;

   if(HLT_IsoMu27==0 || HLT_Mu50==0) continue;
//   if(HLT_Ele32_WPTight_Gsf==0) continue;

// ********************************************************** VBFTag veto
   if(isVBFtag) continue;

   if (fabs(lep1_Id) !=13 ) continue;
//   if (fabs(lep1_Id) !=11 ) continue;

   if (lep1_Pt < 55 ) continue;
//   if (lep1_Pt >  150 ) continue;
//   if (lep2_Pt >  0 ) continue;

   if(FatJet2_pt > 180 || FatJet3_pt > 180 ) continue;
   if(FatJet1_MassSD < 80) continue;
//   if((FatJet2_pt > 130 && FatJet2_MassSD > 30) || (FatJet3_pt > 130 && FatJet3_MassSD > 30)) continue;
// ********************************************************** weight
   weight = (weight/SumGenWeights)*XSec_QCD_PT_300to470*Lumi;
   double PU_weight=PU_Rew[(int)npu];
//   cout<<npu<<"           "<<PU_weight<<endl;
   weight = weight*PU_weight;

// ********************************************************** Fill Histograms
 if(sqrt ( pow(lep1_Eta - FatJet1_eta,2) + pow(phi_dist(lep1_Phi,FatJet1_phi),2) ) <=2.0) continue;
// if (MET < 100) continue;
// if (FatJet1Tau3OverTau2 > 0.46) continue;
// if (FatJet1_MassSD <140) continue;

 _FatJet1_pt  ->Fill(FatJet1_pt,weight);
 _FatJet1_eta ->Fill(FatJet1_eta,weight);
 _FatJet1_phi ->Fill(FatJet1_phi,weight);
 _FatJet1_eta_phi ->Fill(FatJet1_eta, FatJet1_phi, weight);
 _FatJet1_Mass ->Fill(FatJet1_Mass,weight);
 _FatJet1_MassSD ->Fill(FatJet1_MassSD,weight);
 _FatJet1DDBTaggerV2 ->Fill(FatJet1DDBTaggerV2,weight);
 _FatJet1PNetMD_QCD ->Fill(FatJet1PNetMD_QCD,weight);
 _FatJet1PNetMD_Xbb ->Fill(FatJet1PNetMD_Xbb,weight);
 _FatJet1PNetMD_Xcc ->Fill(FatJet1PNetMD_Xcc,weight);
 _FatJet1PNetMD_Xqq ->Fill(FatJet1PNetMD_Xqq,weight);
 _FatJet1PNetHbbvsQCD ->Fill(FatJet1PNetHbbvsQCD,weight);
 _FatJet1PNet_QCD ->Fill(FatJet1PNet_QCD,weight);
 _FatJet1PNet_TvsQCD ->Fill(FatJet1PNet_TvsQCD,weight);
 _FatJet1PNet_WvsQCD ->Fill(FatJet1PNet_WvsQCD,weight);
 _FatJet1PNet_ZvsQCD ->Fill(FatJet1PNet_ZvsQCD,weight);
  double FatJet1PNetMD_Xbb_Modified_as_in_AN_Xbb = FatJet1PNetMD_Xbb/(1 - FatJet1PNetMD_Xcc - FatJet1PNetMD_Xqq);
 _FatJet1PNetMD_Xbb_Modified_as_in_AN_Xbb->Fill(FatJet1PNetMD_Xbb_Modified_as_in_AN_Xbb,weight);
 _FatJet1PNet_mass ->Fill(FatJet1PNet_mass,weight);
 _FatJet1_Tau3OverTau2 ->Fill(FatJet1Tau3OverTau2,weight);

//  if(HLT_IsoMu50_AK8PFJet230_SoftDropMass40)
  if(HLT_AK8PFJet230_SoftDropMass40)
  {
     _FatJet2_pt  ->Fill(FatJet1_pt,weight);
     _FatJet2_eta ->Fill(FatJet1_eta,weight);
     _FatJet2_phi ->Fill(FatJet1_phi,weight);
     _FatJet2_eta_phi ->Fill(FatJet1_eta, FatJet1_phi, weight);
     _FatJet2_Mass ->Fill(FatJet1_Mass,weight);
     _FatJet2_MassSD ->Fill(FatJet1_MassSD,weight);
     _FatJet2DDBTaggerV2 ->Fill(FatJet1DDBTaggerV2,weight);
     _FatJet2PNetMD_QCD ->Fill(FatJet1PNetMD_QCD,weight);
     _FatJet2PNetMD_Xbb ->Fill(FatJet1PNetMD_Xbb,weight);
     _FatJet2PNetMD_Xcc ->Fill(FatJet1PNetMD_Xcc,weight);
     _FatJet2PNetMD_Xqq ->Fill(FatJet1PNetMD_Xqq,weight);
     _FatJet2PNetHbbvsQCD ->Fill(FatJet1PNetHbbvsQCD,weight);
     _FatJet2PNet_QCD ->Fill(FatJet1PNet_QCD,weight);
     _FatJet2PNet_TvsQCD ->Fill(FatJet1PNet_TvsQCD,weight);
     _FatJet2PNet_WvsQCD ->Fill(FatJet1PNet_WvsQCD,weight);
     _FatJet2PNet_ZvsQCD ->Fill(FatJet1PNet_ZvsQCD,weight);
      double FatJet2PNetMD_Xbb_Modified_as_in_AN_Xbb = FatJet2PNetMD_Xbb/(1 - FatJet2PNetMD_Xcc - FatJet2PNetMD_Xqq);
     _FatJet2PNetMD_Xbb_Modified_as_in_AN_Xbb->Fill(FatJet1PNetMD_Xbb_Modified_as_in_AN_Xbb,weight);
     _FatJet2PNet_mass ->Fill(FatJet1PNet_mass,weight);
     _FatJet2_Tau3OverTau2 ->Fill(FatJet1Tau3OverTau2,weight);
   }

 _MET -> Fill(MET,weight);

 } // end event loop

f->Write();

}
