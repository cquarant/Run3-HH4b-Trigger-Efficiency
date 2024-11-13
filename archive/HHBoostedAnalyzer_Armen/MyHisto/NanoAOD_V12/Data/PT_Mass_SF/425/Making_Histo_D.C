#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include <iostream>
#include <map>
#include <math.h>
#include <vector>

// D-phi
double phi_dist(double a, double b) {
  if (fabs(a - b) > 3.14159265) {
    return 6.2831853 - fabs(a - b);
  }
  return fabs(a - b);
}

bool inRange(int low, int high, int x) { return (low <= x && x <= high); }

// ******************************************
TString Run = "D";

void Making_Histo_D() {
  map<int, vector<pair<int, int>>> Good_Lumis{
#include "/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/Data/GoodLumiList_Map_eraD.txt"
  };

  // JEC
  gSystem->Load("libFWCoreFWLite.so");
  vector<JetCorrectorParameters> vPar_AK8;
  vPar_AK8.push_back(JetCorrectorParameters(
      "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      "Summer23BPixPrompt23_RunD_V1_DATA/"
      "Summer23BPixPrompt23_RunD_V1_DATA_L2Relative_AK8PFPuppi.txt"));
  vPar_AK8.push_back(JetCorrectorParameters(
      "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      "Summer23BPixPrompt23_RunD_V1_DATA/"
      "Summer23BPixPrompt23_RunD_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
  FactorizedJetCorrector *corrector_AK8 = new FactorizedJetCorrector(vPar_AK8);

  TFile *f = new TFile("Histograms_" + Run + ".root", "RECREATE");

  // Modification begin: New variables
  Float_t Lower_m[13] = {20,  30,  40,  50,  60,  80, 100,
                         120, 150, 200, 250, 300, 350};
  Float_t Lower_pt[13] = {300, 320, 350, 370, 400,  450, 500,
                          550, 600, 700, 800, 1000, 1200};

  TH2D *_FatJet1_Pt_Mass =
      new TH2D("FatJet1_Pt_Mass", "FatJet1_Pt_Mass", 12, Lower_pt, 12, Lower_m);
  TH2D *_FatJet1_Pt_Mass_M = new TH2D("FatJet1_Pt_Mass_M", "FatJet1_Pt_Mass_M",
                                      12, Lower_pt, 12, Lower_m);

  TH2D *_FatJet2_Pt_Mass =
      new TH2D("FatJet2_Pt_Mass", "FatJet2_Pt_Mass", 12, Lower_pt, 12, Lower_m);
  TH2D *_FatJet2_Pt_Mass_M = new TH2D("FatJet2_Pt_Mass_M", "FatJet2_Pt_Mass_M",
                                      12, Lower_pt, 12, Lower_m);

  // TFile *f1 = new
  // TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/PostBPix/Run2023"+Run+"_JetMET.root");
  TFile *f1 = new TFile(
      "/eos/cms/store/group/phys_higgs/nonresonant_HH/bbbb/sixie/Run3Analysis/"
      "HH/HHTo4BNtupler/ArmenVersion_ICHEP2024/Data_2023/PostBPix/Run2023" +
      Run + "_JetMET.root");

  TTree *InputTree = (TTree *)f1->Get("tree");
  TH1F *NEvents = (TH1F *)f1->Get("nPU_True");
  double Tot_Events = NEvents->GetEntries();

  UInt_t run;
  UInt_t lumi;
  Int_t isVBFtag;

  Bool_t HLT_AK8PFJet230_SoftDropMass40;
  Bool_t HLT_AK8PFJet425_SoftDropMass40;

  Float_t FatJet1_pt;
  Float_t FatJet1_eta;
  Float_t FatJet1_phi;
  Float_t FatJet1_MassSD;
  Float_t FatJet1_rawFactor;

  Float_t FatJet2_pt;
  Float_t FatJet2_eta;
  Float_t FatJet2_phi;
  Float_t FatJet2_MassSD;
  Float_t FatJet2_rawFactor;

  Float_t FatJet3_pt;
  Float_t FatJet3_rawFactor;

  InputTree->SetBranchAddress("run", &run);
  InputTree->SetBranchAddress("lumi", &lumi);
  InputTree->SetBranchAddress("isVBFtag", &isVBFtag);

  InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40",
                              &HLT_AK8PFJet230_SoftDropMass40);
  InputTree->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass40",
                              &HLT_AK8PFJet425_SoftDropMass40);

  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);

  InputTree->SetBranchAddress("fatJet2_pt", &FatJet2_pt);
  InputTree->SetBranchAddress("fatJet2_eta", &FatJet2_eta);
  InputTree->SetBranchAddress("fatJet2_phi", &FatJet2_phi);
  InputTree->SetBranchAddress("fatJet2_msoftdrop", &FatJet2_MassSD);
  InputTree->SetBranchAddress("fatJet2_rawFactor", &FatJet2_rawFactor);

  InputTree->SetBranchAddress("fatJet3_pt", &FatJet3_pt);
  InputTree->SetBranchAddress("fatJet3_rawFactor", &FatJet3_rawFactor);

  // Trigger Objects
  TTree *InputTree_TrgObj = (TTree *)f1->Get("tree_TrgObj");
  Int_t NTrigger_Objects;
  Float_t Trigger_Object_pt[20];
  Float_t Trigger_Object_eta[20];
  Float_t Trigger_Object_phi[20];
  Int_t Trigger_Object_bit[20];
  InputTree_TrgObj->SetBranchAddress("NTrigger_Objects", &NTrigger_Objects);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_pt", Trigger_Object_pt);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_eta", Trigger_Object_eta);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_phi", Trigger_Object_phi);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_bit", Trigger_Object_bit);

  // Events Loop
  for (int i = 0; i < InputTree_TrgObj->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);
    // ********************************************************** HLT Selection
    //   if(HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35==0) continue;
    //   if(HLT_Ele32_WPTight_Gsf==0) continue;
    if (HLT_AK8PFJet230_SoftDropMass40 == 0)
      continue;

    // ********************************************************** FatJets
    // correction and selection

    if (FatJet1_pt > 0) {
      double Raw_FatJet1_pt = FatJet1_pt * (1.0 - FatJet1_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_pt);
      corrector_AK8->setJetEta(FatJet1_eta);
      corrector_AK8->setJetPhi(FatJet1_phi);
      double This_correction = corrector_AK8->getCorrection();
      FatJet1_pt = Raw_FatJet1_pt * This_correction;
      FatJet1_MassSD =
          FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
    }
    if (FatJet2_pt > 0) {
      double Raw_FatJet2_pt = FatJet2_pt * (1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet2_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      corrector_AK8->setJetPhi(FatJet2_phi);
      double This_correction = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt * This_correction;
      FatJet2_MassSD =
          FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
    }

    // ********************************************************** FatJets
    // Selection
    if (FatJet3_pt > 150)
      continue;
    if (FatJet1_pt < 300 || fabs(FatJet1_eta) > 2.5)
      continue;
    if (FatJet2_pt < 160)
      continue;
    if (phi_dist(FatJet1_phi, FatJet2_phi) < 2.5)
      continue;

    // ********************************************************** VBFTag veto
    if (isVBFtag)
      continue;
    // ********************************************************** JSON
    // Certification
    bool Certified = false;
    for (pair<int, vector<pair<int, int>>> Run_Lumi : Good_Lumis) {
      if (Run_Lumi.first > run)
        break;
      if (Run_Lumi.first == run)
        for (auto LumiBlok : Run_Lumi.second)
          if (inRange(LumiBlok.first, LumiBlok.second, lumi)) {
            Certified = true;
            break;
          }
    }
    if (!Certified)
      continue;

    // ********************************************************** Lepton
    // Selection or Veto
    //   if (fabs(lep1_Id) !=11 ) continue;

    // ********************************************************** Trigger
    // Objects and Matchings
    bool Tag_Matched = false;

    bool Tag_matched_to_AK8PFJet230_SoftDropMass40 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if ((Trigger_Object_bit[itrg] & 4) == 4) {
        if (sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet1_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 100) {
          Tag_matched_to_AK8PFJet230_SoftDropMass40 = true;
          break;
        }
      }

    bool Tag_matched_to_AK8PFJet425 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if (Trigger_Object_bit[itrg] == 1)
        if (sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet1_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 425) {
          Tag_matched_to_AK8PFJet425 = true;
          break;
        }

    if (Tag_matched_to_AK8PFJet230_SoftDropMass40 &&
        Tag_matched_to_AK8PFJet425 && HLT_AK8PFJet425_SoftDropMass40)
      Tag_Matched = true;

    // Probe Matched
    bool Probe_Matched = false;

    bool matched_to_AK8PFJet230_SoftDropMass40 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if ((Trigger_Object_bit[itrg] & 4) == 4)
        if (sqrt(pow((FatJet2_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet2_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 100) {
          matched_to_AK8PFJet230_SoftDropMass40 = true;
          break;
        }

    bool matched_to_AK8PFJet250 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if (Trigger_Object_bit[itrg] == 1)
        if (sqrt(pow((FatJet2_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet2_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 425) {
          matched_to_AK8PFJet250 = true;
          break;
        }

    if (matched_to_AK8PFJet230_SoftDropMass40 && matched_to_AK8PFJet250 &&
        HLT_AK8PFJet425_SoftDropMass40)
      Probe_Matched = true;

    /*
      if (Tag_matched_to_AK8PFJet425 && HLT_AK8PFJet425_SoftDropMass40)
         cout<<"good"<<endl;
      if (!Tag_matched_to_AK8PFJet425 && HLT_AK8PFJet425_SoftDropMass40)
         {
           cout<<"bad      "<<Tag_Matched<<"       "<<FatJet1_pt<<"
      "<<FatJet1_MassSD<<"        "<<FatJet1_eta<<endl; cout<<"bad
      "<<Probe_Matched<<"       "<<FatJet2_pt<<"        "<<FatJet2_MassSD<<"
      "<<FatJet2_eta<<endl;
         }
    */
    // ********************************************************** Fill
    // Histograms
    _FatJet1_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);
    if (Tag_Matched)
      _FatJet1_Pt_Mass_M->Fill(FatJet1_pt, FatJet1_MassSD);

    _FatJet2_Pt_Mass->Fill(FatJet2_pt, FatJet2_MassSD);
    if (Probe_Matched)
      _FatJet2_Pt_Mass_M->Fill(FatJet2_pt, FatJet2_MassSD);

  } // end event loop

  f->Write();
}
