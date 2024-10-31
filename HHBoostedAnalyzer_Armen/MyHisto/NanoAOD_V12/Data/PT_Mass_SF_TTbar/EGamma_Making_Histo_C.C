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
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

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


// ******************************************
TString Run = "C";

void EGamma_Making_Histo_C()
{
 map<int,vector<pair<int,int>>> Good_Lumis
  {
   #include "/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/Data/GoodLumiList_Map_eraC.txt"
  };

// JEC
gSystem->Load("libFWCoreFWLite.so");
vector<JetCorrectorParameters> vPar_AK8_Cv123;
vPar_AK8_Cv123.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2Relative_AK8PFPuppi.txt"));
vPar_AK8_Cv123.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK8_Cv123 = new FactorizedJetCorrector(vPar_AK8_Cv123);

vector<JetCorrectorParameters> vPar_AK8_Cv4;
vPar_AK8_Cv4.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2Relative_AK8PFPuppi.txt"));
vPar_AK8_Cv4.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK8_Cv4 = new FactorizedJetCorrector(vPar_AK8_Cv4);

TFile *f = new TFile("EGamma_Histograms_"+Run+".root","RECREATE");

    // Modification begin: New variables
Float_t Lower_m[16]   = {0,5,10,20,30,40,50,60,80,100,120,150,200,250,300,350};
Float_t Lower_pt[46]  = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,320,340,360,380,400,420,440,460,480,500,550,600,700,800,1000};

TH1D *_FatJet1_pt  = new TH1D("FatJet1_pt","FatJet1_pt",200,0,1000);
TH1D *_FatJet1_eta  = new TH1D("FatJet1_eta","FatJet1_eta",100,-5,5);
TH1D *_FatJet1_phi  = new TH1D("FatJet1_phi","FatJet1_phi",100,-5,5);
TH2D *_FatJet1_eta_phi = new TH2D("FatJet1_eta_phi","FatJet1_eta_phi",100,-5,5,100,-5,5);
TH1D *_FatJet1_Mass  = new TH1D("FatJet1_Mass","FatJet1_Mass",500,0,500);
TH1D *_FatJet1_MassSD  = new TH1D("FatJet1_MassSD","FatJet1_MassSD",500,0,500);
TH2D *_FatJet1_Pt_Mass    = new TH2D("FatJet1_Pt_Mass","FatJet1_Pt_Mass",45,Lower_pt,15,Lower_m);

TH2D *_FatJet2_Pt_Mass_M    = new TH2D("FatJet2_Pt_Mass_M","FatJet2_Pt_Mass_M",45,Lower_pt,15,Lower_m);

TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/PreBPix/Run2023"+Run+"_EGamma.root");

TTree *InputTree = (TTree*)f1->Get("tree");
TH1F  *NEvents = (TH1F*)f1->Get("nPU_True");
double Tot_Events = NEvents->GetEntries();

UInt_t     run;
UInt_t     lumi;
Int_t      isVBFtag;

Bool_t     HLT_Ele32_WPTight_Gsf;
Bool_t     HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

Bool_t     HLT_Mu50;
Bool_t     HLT_IsoMu27;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06;

Bool_t     HLT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

Float_t    MET;

Float_t    lep1_Pt;
Float_t    lep1_Eta;
Float_t    lep1_Phi;
Int_t      lep1_Id;
Float_t    lep2_Pt;
Float_t    lep2_Eta;
Float_t    lep2_Phi;
Int_t      lep2_Id;

Float_t    Jet1_Pt;
Float_t    Jet1_Eta;
Float_t    Jet1_Phi;

Float_t    Jet2_Pt;
Float_t    Jet2_Eta;
Float_t    Jet2_Phi;

Float_t    FatJet1_pt;
Float_t    FatJet1_eta;
Float_t    FatJet1_phi;
Float_t    FatJet1_Mass;
Float_t    FatJet1_MassSD;
Float_t    FatJet1_rawFactor;

Float_t    FatJet2_pt;
Float_t    FatJet2_eta;
Float_t    FatJet2_phi;
Float_t    FatJet2_Mass;
Float_t    FatJet2_MassSD;
Float_t    FatJet2_rawFactor;

Float_t    FatJet3_pt;
Float_t    FatJet3_rawFactor;

InputTree->SetBranchAddress("run",&run);
InputTree->SetBranchAddress("lumi",&lumi);
InputTree->SetBranchAddress("isVBFtag",&isVBFtag);

InputTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
InputTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("lep1Pt",&lep1_Pt);
InputTree->SetBranchAddress("lep1Eta",&lep1_Eta);
InputTree->SetBranchAddress("lep1Phi",&lep1_Phi);
InputTree->SetBranchAddress("lep1Id",&lep1_Id);
InputTree->SetBranchAddress("lep2Pt",&lep2_Pt);
InputTree->SetBranchAddress("lep2Eta",&lep2_Eta);
InputTree->SetBranchAddress("lep2Phi",&lep2_Phi);
InputTree->SetBranchAddress("lep2Id",&lep2_Id);
InputTree->SetBranchAddress("MET",&MET);

InputTree->SetBranchAddress("jet1Pt",&Jet1_Pt);
InputTree->SetBranchAddress("jet1Eta",&Jet1_Eta);
InputTree->SetBranchAddress("jet1Phi",&Jet1_Phi);
InputTree->SetBranchAddress("jet2Pt",&Jet2_Pt);
InputTree->SetBranchAddress("jet2Eta",&Jet2_Eta);
InputTree->SetBranchAddress("jet2Phi",&Jet2_Phi);

InputTree->SetBranchAddress("fatJet1_pt",&FatJet1_pt);
InputTree->SetBranchAddress("fatJet1_eta",&FatJet1_eta);
InputTree->SetBranchAddress("fatJet1_phi",&FatJet1_phi);
InputTree->SetBranchAddress("fatJet1_mass",&FatJet1_Mass);
InputTree->SetBranchAddress("fatJet1_msoftdrop",&FatJet1_MassSD);
InputTree->SetBranchAddress("fatJet1_rawFactor",&FatJet1_rawFactor);

InputTree->SetBranchAddress("fatJet2_pt",&FatJet2_pt);
InputTree->SetBranchAddress("fatJet2_eta",&FatJet2_eta);
InputTree->SetBranchAddress("fatJet2_phi",&FatJet2_phi);
InputTree->SetBranchAddress("fatJet2_mass",&FatJet2_Mass);
InputTree->SetBranchAddress("fatJet2_msoftdrop",&FatJet2_MassSD);
InputTree->SetBranchAddress("fatJet2_rawFactor",&FatJet2_rawFactor);

InputTree->SetBranchAddress("fatJet3_pt",&FatJet3_pt);
InputTree->SetBranchAddress("fatJet3_rawFactor",&FatJet3_rawFactor);

// Trigger Objects
TTree *InputTree_TrgObj = (TTree*)f1->Get("tree_TrgObj");
Int_t NTrigger_Objects;
Float_t Trigger_Object_pt[20];
Float_t Trigger_Object_eta[20];
Float_t Trigger_Object_phi[20];
Int_t   Trigger_Object_bit[20];
InputTree_TrgObj->SetBranchAddress("NTrigger_Objects",   &NTrigger_Objects);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_pt",   Trigger_Object_pt);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_eta",  Trigger_Object_eta);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_phi",  Trigger_Object_phi);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_bit",  Trigger_Object_bit);

  // Events Loop
  for(int i=0;i<InputTree_TrgObj->GetEntries();i++)
   {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);
// ********************************************************** HLT Selection
   if ( ! (HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) ==11) ) continue;
   if (lep1_Pt < 50) continue;

// ********************************************************** FatJets correction and selection
   if(run <367765)
    {
    if(FatJet1_pt > 0)
     {
      double Raw_FatJet1_pt = FatJet1_pt*(1.0 - FatJet1_rawFactor);
      corrector_AK8_Cv123->setJetPt(Raw_FatJet1_pt);
      corrector_AK8_Cv123->setJetEta(FatJet1_eta);
      corrector_AK8_Cv123->setJetPhi(FatJet1_phi);
      double This_correction = corrector_AK8_Cv123->getCorrection();
      FatJet1_pt     = Raw_FatJet1_pt * This_correction;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
     }
    if(FatJet2_pt > 0)
     {
      double Raw_FatJet2_pt = FatJet2_pt*(1.0 - FatJet2_rawFactor);
      corrector_AK8_Cv123->setJetPt(Raw_FatJet2_pt);
      corrector_AK8_Cv123->setJetEta(FatJet2_eta);
      corrector_AK8_Cv123->setJetPhi(FatJet2_phi);
      double This_correction = corrector_AK8_Cv123->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt* This_correction;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
     }
    }

   if(run >=367765)
    {
    if(FatJet1_pt > 0)
     {
      double Raw_FatJet1_pt = FatJet1_pt*(1.0 - FatJet1_rawFactor);
      corrector_AK8_Cv4->setJetPt(Raw_FatJet1_pt);
      corrector_AK8_Cv4->setJetEta(FatJet1_eta);
      corrector_AK8_Cv4->setJetPhi(FatJet1_phi);
      double This_correction = corrector_AK8_Cv4->getCorrection();
      FatJet1_pt     = Raw_FatJet1_pt * This_correction;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
     }
    if(FatJet2_pt > 0)
     {
      double Raw_FatJet2_pt = FatJet2_pt*(1.0 - FatJet2_rawFactor);
      corrector_AK8_Cv4->setJetPt(Raw_FatJet2_pt);
      corrector_AK8_Cv4->setJetEta(FatJet2_eta);
      corrector_AK8_Cv4->setJetPhi(FatJet2_phi);
      double This_correction = corrector_AK8_Cv4->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt* This_correction;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
     }
    }
// ********************************************************** FatJets Selection
    if(FatJet2_pt > 180) continue;
    if(FatJet1_pt < 160) continue;
    if (lep1_Pt < 55) continue;
    if (lep2_Pt > 30) continue;
    if(phi_dist(FatJet1_phi,lep1_Phi) < 2.0) continue;
//    if(fabs(FatJet1_eta) > 1.4) continue;
//    if(fabs(FatJet1_eta) > 2.5 || fabs(FatJet1_eta) < 1.4) continue;
   if (MET < 50) continue;

   double Dr_J1FJ =-1;
   double Dr_J1L =10;
   if(Jet1_Pt > 40)
     {
      Dr_J1FJ = sqrt ( pow(Jet1_Eta - FatJet1_eta,2) + pow(phi_dist(Jet1_Phi,FatJet1_phi),2) );
      Dr_J1L  = sqrt ( pow(Jet1_Eta - lep1_Eta,2) + pow(phi_dist(Jet1_Phi,lep1_Phi),2) );
     }
   double Dr_J2FJ =-1;
   double Dr_J2L =10;
   if(Jet2_Pt > 40)
    {
      Dr_J2FJ = sqrt ( pow(Jet2_Eta - FatJet1_eta,2) + pow(phi_dist(Jet2_Phi,FatJet1_phi),2) );
      Dr_J2L  = sqrt ( pow(Jet2_Eta - lep1_Eta,2) + pow(phi_dist(Jet2_Phi,lep1_Phi),2) );
    }

   double Dr_JFJ_Max = Dr_J1FJ;
   double Dr_JFJ_Min = Dr_J2FJ;
   double Dr_JmaxL   = Dr_J1L;
   if(Dr_J2FJ > Dr_J1FJ)
     {
      Dr_JFJ_Max = Dr_J2FJ;
      Dr_JFJ_Min = Dr_J1FJ;
      Dr_JmaxL   = Dr_J2L;
     }

  if(Dr_JFJ_Max < 0) continue;
  if( Dr_J1L <= 0.4 || Dr_J2L <= 0.4 ) continue;
  if(Dr_JmaxL > 3.5 ) continue;

// ********************************************************** VBFTag veto
//   if(isVBFtag) continue;
// ********************************************************** JSON Certification
    bool Certified=false;
    for (pair<int, vector<pair<int,int>>> Run_Lumi : Good_Lumis)
      {
       if(Run_Lumi.first>run) break;
       if(Run_Lumi.first==run)
       for(auto LumiBlok : Run_Lumi.second)
         if ( inRange(LumiBlok.first,LumiBlok.second,lumi) ) { Certified=true; break; }
      }
    if(!Certified) continue;

// ********************************************************** Lepton Selection or Veto
//   if (fabs(lep1_Id) !=11 ) continue;

// ********************************************************** Trigger Objects and Matchings
  // Probe Matched
  bool Probe_Matched = false;

  bool matched_to_AK8PFJet230_SoftDropMass40 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if((Trigger_Object_bit[itrg] & 4) == 4)
       if(sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet1_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 100)
         {matched_to_AK8PFJet230_SoftDropMass40 = true; break;}

  if (matched_to_AK8PFJet230_SoftDropMass40)
     Probe_Matched=true;

// ********************************************************** Fill Histograms

 _FatJet1_pt  ->Fill(FatJet1_pt);
 _FatJet1_eta ->Fill(FatJet1_eta);
 _FatJet1_phi ->Fill(FatJet1_phi);
 _FatJet1_eta_phi ->Fill(FatJet1_eta, FatJet1_phi);
 _FatJet1_Mass ->Fill(FatJet1_Mass);
 _FatJet1_MassSD ->Fill(FatJet1_MassSD);
 _FatJet1_Pt_Mass  -> Fill(FatJet1_pt,FatJet1_MassSD);

  if(Probe_Matched)
  {
   _FatJet2_Pt_Mass_M  -> Fill(FatJet1_pt,FatJet1_MassSD);
  }


 } // end event loop

f->Write();

}
