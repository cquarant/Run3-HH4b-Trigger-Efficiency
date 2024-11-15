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
#include <cstdlib> // for std::getenv
#include <iostream>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <vector>

std::string getCMSSWBase() {
  const char *cmssw_base = std::getenv("CMSSW_BASE");
  if (!cmssw_base) {
    throw std::runtime_error("CMSSW_BASE environment variable not set! Did you "
                             "forget to run 'cmsenv'?");
  }
  return std::string(cmssw_base);
}

// D-phi
double phi_dist(double a, double b) {
  if (fabs(a - b) > 3.14159265) {
    return 6.2831853 - fabs(a - b);
  }
  return fabs(a - b);
}

bool inRange(int low, int high, int x) { return (low <= x && x <= high); }

// ******************************************
TString Run = "C";

void Making_Histo_C() {
  map<int, vector<pair<int, int>>> Good_Lumis{
#include "HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/Data/GoodLumiList_Map_eraC.txt"
  };

  // JEC
  std::string jec_path = getCMSSWBase() + "/src/JECs/";
  gSystem->Load("libFWCoreFWLite.so");
  vector<JetCorrectorParameters> vPar_AK8_Cv123;
  vPar_AK8_Cv123.push_back(JetCorrectorParameters(
      jec_path + "Summer23Prompt23_RunCv123_V1_DATA/"
               "Summer23Prompt23_RunCv123_V1_DATA_L2Relative_AK8PFPuppi.txt"));
  vPar_AK8_Cv123.push_back(JetCorrectorParameters(
      jec_path + 
      "Summer23Prompt23_RunCv123_V1_DATA/"
      "Summer23Prompt23_RunCv123_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
  FactorizedJetCorrector *corrector_AK8_Cv123 =
      new FactorizedJetCorrector(vPar_AK8_Cv123);

  vector<JetCorrectorParameters> vPar_AK8_Cv4;
  vPar_AK8_Cv4.push_back(JetCorrectorParameters(
      // "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      jec_path + "Summer23Prompt23_RunCv4_V1_DATA/"
               "Summer23Prompt23_RunCv4_V1_DATA_L2Relative_AK8PFPuppi.txt"));
  vPar_AK8_Cv4.push_back(JetCorrectorParameters(
      // "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      jec_path + "Summer23Prompt23_RunCv4_V1_DATA/"
               "Summer23Prompt23_RunCv4_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
  FactorizedJetCorrector *corrector_AK8_Cv4 =
      new FactorizedJetCorrector(vPar_AK8_Cv4);

  TFile *f = new TFile("JetMET_Histograms_" + Run + ".root", "RECREATE");

  // Modification begin: New variables
  Float_t Lower_m[16] = {0,  5,   10,  20,  30,  40,  50,  60,
                         80, 100, 120, 150, 200, 250, 300, 350};
  Float_t Lower_pt[46] = {0,   10,  20,  30,  40,  50,  60,  70,  80,  90,
                          100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
                          200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
                          300, 320, 340, 360, 380, 400, 420, 440, 460, 480,
                          500, 550, 600, 700, 800, 1000};

  TH1D *_FatJet1_pt = new TH1D("FatJet1_pt", "FatJet1_pt", 200, 0, 1000);
  TH1D *_FatJet1_eta = new TH1D("FatJet1_eta", "FatJet1_eta", 100, -5, 5);
  TH1D *_FatJet1_phi = new TH1D("FatJet1_phi", "FatJet1_phi", 100, -5, 5);
  TH2D *_FatJet1_eta_phi =
      new TH2D("FatJet1_eta_phi", "FatJet1_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet1_Mass = new TH1D("FatJet1_Mass", "FatJet1_Mass", 500, 0, 500);
  TH1D *_FatJet1_MassSD =
      new TH1D("FatJet1_MassSD", "FatJet1_MassSD", 500, 0, 500);
  TH2D *_FatJet1_Pt_Mass =
      new TH2D("FatJet1_Pt_Mass", "FatJet1_Pt_Mass", 45, Lower_pt, 15, Lower_m);

  TH1D *_FatJet2_pt = new TH1D("FatJet2_pt", "FatJet2_pt", 200, 0, 1000);
  TH1D *_FatJet2_eta = new TH1D("FatJet2_eta", "FatJet2_eta", 100, -5, 5);
  TH1D *_FatJet2_phi = new TH1D("FatJet2_phi", "FatJet2_phi", 100, -5, 5);
  TH2D *_FatJet2_eta_phi =
      new TH2D("FatJet2_eta_phi", "FatJet2_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet2_Mass = new TH1D("FatJet2_Mass", "FatJet2_Mass", 500, 0, 500);
  TH1D *_FatJet2_MassSD =
      new TH1D("FatJet2_MassSD", "FatJet2_MassSD", 500, 0, 500);
  TH2D *_FatJet2_Pt_Mass =
      new TH2D("FatJet2_Pt_Mass", "FatJet2_Pt_Mass", 45, Lower_pt, 15, Lower_m);
  TH2D *_FatJet2_Pt_Mass_M = new TH2D("FatJet2_Pt_Mass_M", "FatJet2_Pt_Mass_M",
                                      45, Lower_pt, 15, Lower_m);

  Float_t Lower_pt_N[9] = {230, 240, 250, 270, 300, 350, 500, 700, 1000};
  TH2D *_FatJet2_Pt_MassN = new TH2D("FatJet2_Pt_MassN", "FatJet2_Pt_MassN", 8,
                                     Lower_pt_N, 15, Lower_m);
  TH2D *_FatJet2_Pt_MassN_M = new TH2D(
      "FatJet2_Pt_MassN_M", "FatJet2_Pt_MassN_M", 8, Lower_pt_N, 15, Lower_m);

  // TODO: Modify the path if needed
  TFile *f1 = new TFile(
      "/eos/cms/store/group/phys_higgs/nonresonant_HH/bbbb/sixie/Run3Analysis/"
      "HH/HHTo4BNtupler/ArmenVersion_ICHEP2024/Data_2023/PreBPix/Run2023" +
      Run + "_JetMET.root");

  TTree *InputTree = (TTree *)f1->Get("tree");
  TH1F *NEvents = (TH1F *)f1->Get("nPU_True");
  double Tot_Events = NEvents->GetEntries();

  UInt_t run;
  UInt_t lumi;
  Int_t isVBFtag;

  Bool_t HLT_Ele32_WPTight_Gsf;
  Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;
  Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

  Bool_t HLT_Mu50;
  Bool_t HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
  Bool_t HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06;

  Bool_t HLT_AK8PFJet230_SoftDropMass40;
  Bool_t HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

  Float_t MET;

  Float_t lep1_Pt;
  Float_t lep1_Eta;
  Float_t lep1_Phi;
  Int_t lep1_Id;

  Float_t FatJet1_pt;
  Float_t FatJet1_eta;
  Float_t FatJet1_phi;
  Float_t FatJet1_Mass;
  Float_t FatJet1_MassSD;
  Float_t FatJet1DDBTaggerV2;
  Float_t FatJet1PNetMD_QCD;
  Float_t FatJet1PNetMD_Xbb;
  Float_t FatJet1PNetMD_Xcc;
  Float_t FatJet1PNetMD_Xqq;
  Float_t FatJet1PNetHbbvsQCD;
  Float_t FatJet1PNet_QCD;
  Float_t FatJet1PNet_TvsQCD;
  Float_t FatJet1PNet_WvsQCD;
  Float_t FatJet1PNet_ZvsQCD;
  Float_t FatJet1PNet_mass;
  Float_t FatJet1Tau3OverTau2;
  Float_t FatJet1_rawFactor;

  Float_t FatJet2_pt;
  Float_t FatJet2_eta;
  Float_t FatJet2_phi;
  Float_t FatJet2_Mass;
  Float_t FatJet2_MassSD;
  Float_t FatJet2DDBTaggerV2;
  Float_t FatJet2PNetMD_QCD;
  Float_t FatJet2PNetMD_Xbb;
  Float_t FatJet2PNetMD_Xcc;
  Float_t FatJet2PNetMD_Xqq;
  Float_t FatJet2PNetHbbvsQCD;
  Float_t FatJet2PNet_QCD;
  Float_t FatJet2PNet_TvsQCD;
  Float_t FatJet2PNet_WvsQCD;
  Float_t FatJet2PNet_ZvsQCD;
  Float_t FatJet2PNet_mass;
  Float_t FatJet2Tau3OverTau2;
  Float_t FatJet2_rawFactor;

  Float_t FatJet3_pt;
  Float_t FatJet3_rawFactor;

  InputTree->SetBranchAddress("run", &run);
  InputTree->SetBranchAddress("lumi", &lumi);
  InputTree->SetBranchAddress("isVBFtag", &isVBFtag);

  InputTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
  InputTree->SetBranchAddress(
      "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40",
      &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
  InputTree->SetBranchAddress(
      "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06",
      &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

  InputTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
  InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40",
                              &HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
  InputTree->SetBranchAddress(
      "HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06",
      &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06);

  InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40",
                              &HLT_AK8PFJet230_SoftDropMass40);
  InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06",
                              &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

  InputTree->SetBranchAddress("lep1Pt", &lep1_Pt);
  InputTree->SetBranchAddress("lep1Eta", &lep1_Eta);
  InputTree->SetBranchAddress("lep1Phi", &lep1_Phi);
  InputTree->SetBranchAddress("lep1Id", &lep1_Id);
  InputTree->SetBranchAddress("MET", &MET);

  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_mass", &FatJet1_Mass);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);

  InputTree->SetBranchAddress("fatJet2_pt", &FatJet2_pt);
  InputTree->SetBranchAddress("fatJet2_eta", &FatJet2_eta);
  InputTree->SetBranchAddress("fatJet2_phi", &FatJet2_phi);
  InputTree->SetBranchAddress("fatJet2_mass", &FatJet2_Mass);
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
  std::cout << "Total Events: " << InputTree_TrgObj->GetEntries() << std::endl;
  for (int i = 0; i < InputTree_TrgObj->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);
    // ********************************************************** HLT Selection
    //   if(HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06==0) continue;
    //   if(HLT_Ele32_WPTight_Gsf==0) continue;
    // TODO: this is the difference
    if (HLT_AK8PFJet230_SoftDropMass40 == 0)
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

    // ********************************************************** FatJets
    // correction and selection
    if (run < 367765) {
      if (FatJet1_pt > 0) {
        double Raw_FatJet1_pt = FatJet1_pt * (1.0 - FatJet1_rawFactor);
        corrector_AK8_Cv123->setJetPt(Raw_FatJet1_pt);
        corrector_AK8_Cv123->setJetEta(FatJet1_eta);
        corrector_AK8_Cv123->setJetPhi(FatJet1_phi);
        double This_correction = corrector_AK8_Cv123->getCorrection();
        FatJet1_pt = Raw_FatJet1_pt * This_correction;
        FatJet1_MassSD =
            FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
      }
      if (FatJet2_pt > 0) {
        double Raw_FatJet2_pt = FatJet2_pt * (1.0 - FatJet2_rawFactor);
        corrector_AK8_Cv123->setJetPt(Raw_FatJet2_pt);
        corrector_AK8_Cv123->setJetEta(FatJet2_eta);
        corrector_AK8_Cv123->setJetPhi(FatJet2_phi);
        double This_correction = corrector_AK8_Cv123->getCorrection();
        FatJet2_pt = Raw_FatJet2_pt * This_correction;
        FatJet2_MassSD =
            FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
      }
    }

    if (run >= 367765) {
      if (FatJet1_pt > 0) {
        double Raw_FatJet1_pt = FatJet1_pt * (1.0 - FatJet1_rawFactor);
        corrector_AK8_Cv4->setJetPt(Raw_FatJet1_pt);
        corrector_AK8_Cv4->setJetEta(FatJet1_eta);
        corrector_AK8_Cv4->setJetPhi(FatJet1_phi);
        double This_correction = corrector_AK8_Cv4->getCorrection();
        FatJet1_pt = Raw_FatJet1_pt * This_correction;
        FatJet1_MassSD =
            FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
      }
      if (FatJet2_pt > 0) {
        double Raw_FatJet2_pt = FatJet2_pt * (1.0 - FatJet2_rawFactor);
        corrector_AK8_Cv4->setJetPt(Raw_FatJet2_pt);
        corrector_AK8_Cv4->setJetEta(FatJet2_eta);
        corrector_AK8_Cv4->setJetPhi(FatJet2_phi);
        double This_correction = corrector_AK8_Cv4->getCorrection();
        FatJet2_pt = Raw_FatJet2_pt * This_correction;
        FatJet2_MassSD =
            FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
      }
    }

    // ********************************************************** FatJets
    // Selection
    if (FatJet3_pt > 150)
      continue;
    if (FatJet1_pt < 300 || fabs(FatJet1_eta) > 2.5 || FatJet1_MassSD < 80)
      continue;
    if (FatJet2_pt < 160)
      continue;
    if (phi_dist(FatJet1_phi, FatJet2_phi) < 2.5)
      continue;
    //    if(fabs(FatJet2_eta) > 1.4) continue;
    //    if(fabs(FatJet2_eta) > 2.5 || fabs(FatJet2_eta) < 1.4) continue;

    // ********************************************************** VBFTag veto
    //   if(isVBFtag) continue;

    // ********************************************************** Lepton
    // Selection or Veto
    //   if (fabs(lep1_Id) !=11 ) continue;

    // ********************************************************** Trigger
    // Objects and Matchings
    bool Tag_Matched = false;

    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if ((Trigger_Object_bit[itrg] & 4) == 4) {
        if (sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet1_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 100) {
          Tag_Matched = true;
          break;
        } else
          break;
      }
    if (!Tag_Matched)
      continue;

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

    if (matched_to_AK8PFJet230_SoftDropMass40)
      Probe_Matched = true;

    // ********************************************************** Fill
    // Histograms

    _FatJet1_pt->Fill(FatJet1_pt);
    _FatJet1_eta->Fill(FatJet1_eta);
    _FatJet1_phi->Fill(FatJet1_phi);
    _FatJet1_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
    _FatJet1_Mass->Fill(FatJet1_Mass);
    _FatJet1_MassSD->Fill(FatJet1_MassSD);
    _FatJet1_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);

    _FatJet2_pt->Fill(FatJet2_pt);
    _FatJet2_eta->Fill(FatJet2_eta);
    _FatJet2_phi->Fill(FatJet2_phi);
    _FatJet2_eta_phi->Fill(FatJet2_eta, FatJet2_phi);
    _FatJet2_Mass->Fill(FatJet2_Mass);
    _FatJet2_MassSD->Fill(FatJet2_MassSD);
    _FatJet2_Pt_Mass->Fill(FatJet2_pt, FatJet2_MassSD);
    _FatJet2_Pt_MassN->Fill(FatJet2_pt, FatJet2_MassSD);

    if (Probe_Matched) {
      _FatJet2_Pt_Mass_M->Fill(FatJet2_pt, FatJet2_MassSD);
      _FatJet2_Pt_MassN_M->Fill(FatJet2_pt, FatJet2_MassSD);
    }

  } // end event loop

  f->Write();
  std::cout << "Done. Writing to file Histograms_" + Run + ".root" << std::endl;
}
