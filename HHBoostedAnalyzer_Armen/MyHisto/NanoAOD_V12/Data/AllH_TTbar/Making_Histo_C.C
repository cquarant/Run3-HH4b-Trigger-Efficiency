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

void Making_Histo_C()
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

TFile *f = new TFile("JetMET_Histograms_"+Run+".root","RECREATE");
TTree *outputTree = new TTree("tree", "");

Float_t T_weight;
Float_t T_fatJet1_pt, T_fatJet1_eta, T_fatJet1_phi, T_fatJet1_msoftdrop, T_fatJet1_particleNet_XbbVsQCD, T_fatJet1_Tau3OverTau2, T_fatJet1_particleNet_Xbb_Legacy, T_fatJet1_particleNet_Xbb_Legacy_AN;
Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_phi, T_fatJet2_msoftdrop, T_fatJet2_particleNet_XbbVsQCD, T_fatJet2_Tau3OverTau2, T_fatJet2_particleNet_Xbb_Legacy, T_fatJet2_particleNet_Xbb_Legacy_AN;
Float_t T_PTjj;

outputTree->Branch("T_weight",                           &T_weight,      "T_weight/F");

outputTree->Branch("T_fatJet1_pt",                        &T_fatJet1_pt,                        "T_fatJet1_pt/F");
outputTree->Branch("T_fatJet1_eta",                       &T_fatJet1_eta,                       "T_fatJet1_eta/F");
outputTree->Branch("T_fatJet1_phi",                       &T_fatJet1_phi,                       "T_fatJet1_phi/F");
outputTree->Branch("T_fatJet1_msoftdrop",                 &T_fatJet1_msoftdrop,                 "T_fatJet1_msoftdrop/F");
outputTree->Branch("T_fatJet1_particleNet_XbbVsQCD",      &T_fatJet1_particleNet_XbbVsQCD,      "T_fatJet1_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet1_Tau3OverTau2",              &T_fatJet1_Tau3OverTau2,              "T_fatJet1_Tau3OverTau2/F");
outputTree->Branch("T_fatJet1_particleNet_Xbb_Legacy",    &T_fatJet1_particleNet_Xbb_Legacy,    "T_fatJet1_particleNet_Xbb_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_Xbb_Legacy_AN", &T_fatJet1_particleNet_Xbb_Legacy_AN, "T_fatJet1_particleNet_Xbb_Legacy_AN/F");

outputTree->Branch("T_fatJet2_pt",                        &T_fatJet2_pt,                        "T_fatJet2_pt/F");
outputTree->Branch("T_fatJet2_eta",                       &T_fatJet2_eta,                       "T_fatJet2_eta/F");
outputTree->Branch("T_fatJet2_phi",                       &T_fatJet2_phi,                       "T_fatJet2_phi/F");
outputTree->Branch("T_fatJet2_msoftdrop",                 &T_fatJet2_msoftdrop,                 "T_fatJet2_msoftdrop/F");
outputTree->Branch("T_fatJet2_particleNet_XbbVsQCD",      &T_fatJet2_particleNet_XbbVsQCD,      "T_fatJet2_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet2_Tau3OverTau2",              &T_fatJet2_Tau3OverTau2,              "T_fatJet2_Tau3OverTau2/F");
outputTree->Branch("T_fatJet2_particleNet_Xbb_Legacy",    &T_fatJet2_particleNet_Xbb_Legacy,    "T_fatJet2_particleNet_Xbb_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_Xbb_Legacy_AN", &T_fatJet2_particleNet_Xbb_Legacy_AN, "T_fatJet2_particleNet_Xbb_Legacy_AN/F");

outputTree->Branch("T_PTjj",                              &T_PTjj,         "T_PTjj/F");

TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/Legacy/PreBPix/Run2023"+Run+"_JetMET.root");

TTree *InputTree = (TTree*)f1->Get("tree");
TH1F  *NEvents = (TH1F*)f1->Get("nPU_True");
double Tot_Events = NEvents->GetEntries();

UInt_t     run;
UInt_t     lumi;
Int_t      isVBFtag;

Bool_t     HLT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_AK8PFJet425_SoftDropMass40;

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
Float_t    FatJet1PNetMD_Xbb;
Float_t    FatJet1_HbbvsQCD;
Float_t    FatJet1_Tau3OverTau2;
Float_t    FatJet1_rawFactor;
Float_t    FatJet1PNetMD_Xbb_Legacy;
Float_t    FatJet1PNetMD_QCD_Legacy;

Float_t    FatJet2_pt;
Float_t    FatJet2_eta;
Float_t    FatJet2_phi;
Float_t    FatJet2_Mass;
Float_t    FatJet2_MassSD;
Float_t    FatJet2PNetMD_Xbb;
Float_t    FatJet2_Tau3OverTau2;
Float_t    FatJet2_rawFactor;
Float_t    FatJet2PNetMD_Xbb_Legacy;
Float_t    FatJet2PNetMD_QCD_Legacy;

Float_t    FatJet3_pt;
Float_t    FatJet3_rawFactor;

InputTree->SetBranchAddress("run",&run);
InputTree->SetBranchAddress("lumi",&lumi);
InputTree->SetBranchAddress("isVBFtag",&isVBFtag);

InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass40", &HLT_AK8PFJet425_SoftDropMass40);

InputTree->SetBranchAddress("fatJet1_pt",&FatJet1_pt);
InputTree->SetBranchAddress("fatJet1_eta",&FatJet1_eta);
InputTree->SetBranchAddress("fatJet1_phi",&FatJet1_phi);
InputTree->SetBranchAddress("fatJet1_mass",&FatJet1_Mass);
InputTree->SetBranchAddress("fatJet1_msoftdrop",&FatJet1_MassSD);
InputTree->SetBranchAddress("fatJet1_rawFactor",&FatJet1_rawFactor);
InputTree->SetBranchAddress("fatJet1_Tau3OverTau2",&FatJet1_Tau3OverTau2);
InputTree->SetBranchAddress("fatJet1_particleNet_XbbVsQCD",&FatJet1PNetMD_Xbb);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xbb",&FatJet1PNetMD_Xbb_Legacy);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCD",&FatJet1PNetMD_QCD_Legacy);

InputTree->SetBranchAddress("fatJet2_pt",&FatJet2_pt);
InputTree->SetBranchAddress("fatJet2_eta",&FatJet2_eta);
InputTree->SetBranchAddress("fatJet2_phi",&FatJet2_phi);
InputTree->SetBranchAddress("fatJet2_mass",&FatJet2_Mass);
InputTree->SetBranchAddress("fatJet2_msoftdrop",&FatJet2_MassSD);
InputTree->SetBranchAddress("fatJet2_rawFactor",&FatJet2_rawFactor);
InputTree->SetBranchAddress("fatJet2_Tau3OverTau2",&FatJet2_Tau3OverTau2);
InputTree->SetBranchAddress("fatJet2_particleNet_XbbVsQCD",&FatJet2PNetMD_Xbb);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xbb",&FatJet2PNetMD_Xbb_Legacy);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCD",&FatJet2PNetMD_QCD_Legacy);

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
   if(HLT_AK8PFJet425_SoftDropMass40==0) continue;

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

// ********************************************************** FatJets correction and selection
   if (run<367661) continue;
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

      if(FatJet1_pt < 450 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50 ) continue;
      if(FatJet2_pt < 450 || fabs(FatJet2_eta) > 2.4 || FatJet2_MassSD < 50 ) continue;

// ********************************************************** Lepton Selection or Veto
//   if (fabs(lep1_Id) !=11 ) continue;

// ********************************************************** Trigger Objects and Matchings
  // Matching 1st
  bool matched_TRG_1=false;

  bool matched_to_AK8PFJet230_SoftDropMass40_1 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if((Trigger_Object_bit[itrg] & 4) == 4)
       if(sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet1_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 100)
         {matched_to_AK8PFJet230_SoftDropMass40_1 = true; break;}

  bool matched_to_AK8PFJet250_1 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if(Trigger_Object_bit[itrg] == 1)
       if(sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet1_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 425)
         {matched_to_AK8PFJet250_1 = true; break;}

  if (matched_to_AK8PFJet230_SoftDropMass40_1 && matched_to_AK8PFJet250_1) matched_TRG_1=true;

  // Matching 2nd
  bool matched_TRG_2=false;

  bool matched_to_AK8PFJet230_SoftDropMass40_2 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if((Trigger_Object_bit[itrg] & 4) == 4)
       if(sqrt(pow((FatJet2_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet2_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 100)
         {matched_to_AK8PFJet230_SoftDropMass40_2 = true; break;}

  bool matched_to_AK8PFJet250_2 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if(Trigger_Object_bit[itrg] == 1)
       if(sqrt(pow((FatJet2_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet2_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 425)
         {matched_to_AK8PFJet250_2 = true; break;}

  if (matched_to_AK8PFJet230_SoftDropMass40_2 && matched_to_AK8PFJet250_2) matched_TRG_2=true;

  if (!(matched_TRG_1 || matched_TRG_2)) continue;

// ********************************************** PTJJ

   double FatJet1_px = FatJet1_pt*cos(FatJet1_phi), FatJet1_py = FatJet1_pt*sin(FatJet1_phi);
   double FatJet2_px = FatJet2_pt*cos(FatJet2_phi), FatJet2_py = FatJet2_pt*sin(FatJet2_phi);
   double PTjj = sqrt(pow((FatJet1_px +FatJet2_px),2) + pow((FatJet1_py +FatJet2_py),2));

// ********************************************************** Fill Histograms
  T_weight   = 1.0;

  T_fatJet1_pt                         = FatJet1_pt;
  T_fatJet1_eta                        = FatJet1_eta;
  T_fatJet1_phi                        = FatJet1_phi;
  T_fatJet1_msoftdrop                  = FatJet1_MassSD;
  T_fatJet1_particleNet_XbbVsQCD       = FatJet1PNetMD_Xbb;
  T_fatJet1_Tau3OverTau2               = FatJet1_Tau3OverTau2;
  T_fatJet1_particleNet_Xbb_Legacy     = FatJet1PNetMD_Xbb_Legacy;
  double FatJet1PNetMD_Xbb_Legacy_AN   = FatJet1PNetMD_Xbb_Legacy/(FatJet1PNetMD_Xbb_Legacy + FatJet1PNetMD_QCD_Legacy);
  T_fatJet1_particleNet_Xbb_Legacy_AN  = FatJet1PNetMD_Xbb_Legacy_AN;

  T_fatJet2_pt                         = FatJet2_pt;
  T_fatJet2_eta                        = FatJet2_eta;
  T_fatJet2_phi                        = FatJet2_phi;
  T_fatJet2_msoftdrop                  = FatJet2_MassSD;
  T_fatJet2_particleNet_XbbVsQCD       = FatJet2PNetMD_Xbb;
  T_fatJet2_Tau3OverTau2               = FatJet2_Tau3OverTau2;
  T_fatJet2_particleNet_Xbb_Legacy     = FatJet2PNetMD_Xbb_Legacy;
  double FatJet2PNetMD_Xbb_Legacy_AN   = FatJet2PNetMD_Xbb_Legacy/(FatJet2PNetMD_Xbb_Legacy + FatJet2PNetMD_QCD_Legacy);
  T_fatJet2_particleNet_Xbb_Legacy_AN  = FatJet2PNetMD_Xbb_Legacy_AN;

  T_PTjj = PTjj;

  outputTree->Fill();

 } // end event loop

f->Write();

}


