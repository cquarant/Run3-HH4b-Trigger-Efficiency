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
TString Run = "D";

void Muon_Making_Histo_D() {
  map<int, vector<pair<int, int>>> Good_Lumis{
#include "HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/Data/GoodLumiList_Map_eraD.txt"
  };

  // JEC
  std::string jec_path = getCMSSWBase() + "/src/JECs/";
  gSystem->Load("libFWCoreFWLite.so");
  vector<JetCorrectorParameters> vPar_AK8;
  vPar_AK8.push_back(JetCorrectorParameters(
      // "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      jec_path + 
      "Summer23BPixPrompt23_RunD_V1_DATA/"
      "Summer23BPixPrompt23_RunD_V1_DATA_L2Relative_AK8PFPuppi.txt"));
  vPar_AK8.push_back(JetCorrectorParameters(
      // "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      jec_path + 
      "Summer23BPixPrompt23_RunD_V1_DATA/"
      "Summer23BPixPrompt23_RunD_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
  FactorizedJetCorrector *corrector_AK8 = new FactorizedJetCorrector(vPar_AK8);

  TFile *f = new TFile("Muon_Histograms_" + Run + ".root", "RECREATE");

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
  TH1D *_FatJet1PNetMD_Xbb =
      new TH1D("FatJet1PNetMD_Xbb", "FatJet1PNetMD_Xbb", 100, 0, 1.0);
  TH1D *_FatJet1PNetMD_Xbb_Legacy = new TH1D(
      "FatJet1PNetMD_Xbb_Legacy", "FatJet1PNetMD_Xbb_Legacy", 100, 0, 1.0);
  TH1D *_FatJet1PNetMD_Xbb_Legacy_AN =
      new TH1D("FatJet1PNetMD_Xbb_Legacy_AN", "FatJet1PNetMD_Xbb_Legacy_AN",
               100, 0, 1.0);

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
  TH1D *_FatJet2PNetMD_Xbb =
      new TH1D("FatJet2PNetMD_Xbb", "FatJet2PNetMD_Xbb", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_Legacy = new TH1D(
      "FatJet2PNetMD_Xbb_Legacy", "FatJet2PNetMD_Xbb_Legacy", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_Legacy_AN =
      new TH1D("FatJet2PNetMD_Xbb_Legacy_AN", "FatJet2PNetMD_Xbb_Legacy_AN",
               100, 0, 1.0);

  TH1D *_FatJet3_pt = new TH1D("FatJet3_pt", "FatJet3_pt", 200, 0, 1000);
  TH1D *_FatJet3_eta = new TH1D("FatJet3_eta", "FatJet3_eta", 100, -5, 5);
  TH1D *_FatJet3_phi = new TH1D("FatJet3_phi", "FatJet3_phi", 100, -5, 5);
  TH2D *_FatJet3_eta_phi =
      new TH2D("FatJet3_eta_phi", "FatJet3_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet3_Mass = new TH1D("FatJet3_Mass", "FatJet3_Mass", 500, 0, 500);
  TH1D *_FatJet3_MassSD =
      new TH1D("FatJet3_MassSD", "FatJet3_MassSD", 500, 0, 500);
  TH2D *_FatJet3_Pt_Mass =
      new TH2D("FatJet3_Pt_Mass", "FatJet3_Pt_Mass", 45, Lower_pt, 15, Lower_m);
  TH1D *_FatJet3PNetMD_Xbb =
      new TH1D("FatJet3PNetMD_Xbb", "FatJet3PNetMD_Xbb", 100, 0, 1.0);
  TH1D *_FatJet3PNetMD_Xbb_Legacy = new TH1D(
      "FatJet3PNetMD_Xbb_Legacy", "FatJet3PNetMD_Xbb_Legacy", 100, 0, 1.0);
  TH1D *_FatJet3PNetMD_Xbb_Legacy_AN =
      new TH1D("FatJet3PNetMD_Xbb_Legacy_AN", "FatJet3PNetMD_Xbb_Legacy_AN",
               100, 0, 1.0);

  TH1D *_MET = new TH1D("MET", "MET", 100, 0, 500);
  TH1D *_Lep1_Pt = new TH1D("Lep1_Pt", "Lep1_Pt", 150, 0, 300);
  TH1D *_Dr_LFJ = new TH1D("Dr_LFJ", "Dr_LFJ", 100, -1.0, 9.0);
  TH1D *_Dr_J1FJ = new TH1D("Dr_J1FJ", "Dr_J1FJ", 100, -1.0, 9.0);
  TH1D *_Dr_J2FJ = new TH1D("Dr_J2FJ", "Dr_J2FJ", 100, -1.0, 9.0);
  TH1D *_Dr_JmaxL = new TH1D("Dr_JmaxL", "Dr_JmaxL", 100, -1.0, 9.0);

  // TODO: Update this if necessary
  // TFile *f1 = new TFile(
  //     "/eos/home-t/tumasyan/HHTo4B/Data_2023/Legacy/PostBPix/Run2023" + Run +
  //     "_Muon.root");
  TFile *f1 = new TFile(
      "/eos/cms/store/group/phys_higgs/nonresonant_HH/bbbb/sixie/Run3Analysis/"
      "HH/HHTo4BNtupler/ArmenVersion_ICHEP2024/Data_2023/PostBPix/Run2023" +
      Run + "_Muon.root");

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
  Bool_t HLT_IsoMu27;
  Bool_t HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
  Bool_t HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06;

  Bool_t HLT_AK8PFJet230_SoftDropMass40;
  Bool_t HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

  Float_t MET;

  Float_t lep1_Pt;
  Float_t lep1_Eta;
  Float_t lep1_Phi;
  Int_t lep1_Id;
  Float_t lep2_Pt;
  Float_t lep2_Eta;
  Float_t lep2_Phi;
  Int_t lep2_Id;

  Float_t Jet1_Pt;
  Float_t Jet1_Eta;
  Float_t Jet1_Phi;

  Float_t Jet2_Pt;
  Float_t Jet2_Eta;
  Float_t Jet2_Phi;

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
  Int_t FatJet1_hadronFlavour;
  Int_t FatJet1_nBHadrons;
  Int_t FatJet1_nCHadrons;
  Float_t FatJet1PNetMD_Xbb_Legacy;
  Float_t FatJet1PNetMD_QCD_Legacy;

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
  Int_t FatJet2_hadronFlavour;
  Int_t FatJet2_nBHadrons;
  Int_t FatJet2_nCHadrons;
  Float_t FatJet2PNetMD_Xbb_Legacy;
  Float_t FatJet2PNetMD_QCD_Legacy;

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
  InputTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
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
  InputTree->SetBranchAddress("lep2Pt", &lep2_Pt);
  InputTree->SetBranchAddress("lep2Eta", &lep2_Eta);
  InputTree->SetBranchAddress("lep2Phi", &lep2_Phi);
  InputTree->SetBranchAddress("lep2Id", &lep2_Id);
  InputTree->SetBranchAddress("MET", &MET);

  InputTree->SetBranchAddress("jet1Pt", &Jet1_Pt);
  InputTree->SetBranchAddress("jet1Eta", &Jet1_Eta);
  InputTree->SetBranchAddress("jet1Phi", &Jet1_Phi);
  InputTree->SetBranchAddress("jet2Pt", &Jet2_Pt);
  InputTree->SetBranchAddress("jet2Eta", &Jet2_Eta);
  InputTree->SetBranchAddress("jet2Phi", &Jet2_Phi);

  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_mass", &FatJet1_Mass);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);
  InputTree->SetBranchAddress("fatJet1_particleNet_XbbVsQCD",
                              &FatJet1PNetMD_Xbb);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xbb",
                              &FatJet1PNetMD_Xbb_Legacy);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCD",
                              &FatJet1PNetMD_QCD_Legacy);

  InputTree->SetBranchAddress("fatJet2_pt", &FatJet2_pt);
  InputTree->SetBranchAddress("fatJet2_eta", &FatJet2_eta);
  InputTree->SetBranchAddress("fatJet2_phi", &FatJet2_phi);
  InputTree->SetBranchAddress("fatJet2_mass", &FatJet2_Mass);
  InputTree->SetBranchAddress("fatJet2_msoftdrop", &FatJet2_MassSD);
  InputTree->SetBranchAddress("fatJet2_rawFactor", &FatJet2_rawFactor);
  InputTree->SetBranchAddress("fatJet2_particleNet_XbbVsQCD",
                              &FatJet2PNetMD_Xbb);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xbb",
                              &FatJet2PNetMD_Xbb_Legacy);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCD",
                              &FatJet2PNetMD_QCD_Legacy);

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
  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);
    // ********************************************************** HLT Selection
    if (!(HLT_IsoMu27 && fabs(lep1_Id) == 13))
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

    if (FatJet1_pt < 250 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50)
      continue;

    // ********************************************************** VBFTag veto
    //   if(isVBFtag) continue;

    // ********************************************************** Lepton
    // Selection or Veto
    if (lep1_Pt < 50)
      continue;
    if (lep2_Pt > 30)
      continue;
    if (FatJet2_pt > 200 && FatJet2_MassSD > 50)
      continue;
    double Dr_LFJ = sqrt(pow(lep1_Eta - FatJet1_eta, 2) +
                         pow(phi_dist(lep1_Phi, FatJet1_phi), 2));
    if (Dr_LFJ < 1.5)
      continue;
    if (MET < 50)
      continue;
    //   if (FatJet1Tau3OverTau2 > 0.5) continue;
    if (FatJet1_MassSD < 50)
      continue;
    double Dr_J1FJ = -1;
    double Dr_J1L = 10;
    if (Jet1_Pt > 40) {
      Dr_J1FJ = sqrt(pow(Jet1_Eta - FatJet1_eta, 2) +
                     pow(phi_dist(Jet1_Phi, FatJet1_phi), 2));
      Dr_J1L = sqrt(pow(Jet1_Eta - lep1_Eta, 2) +
                    pow(phi_dist(Jet1_Phi, lep1_Phi), 2));
    }
    double Dr_J2FJ = -1;
    double Dr_J2L = 10;
    if (Jet2_Pt > 40) {
      Dr_J2FJ = sqrt(pow(Jet2_Eta - FatJet1_eta, 2) +
                     pow(phi_dist(Jet2_Phi, FatJet1_phi), 2));
      Dr_J2L = sqrt(pow(Jet2_Eta - lep1_Eta, 2) +
                    pow(phi_dist(Jet2_Phi, lep1_Phi), 2));
    }

    double Dr_JFJ_Max = Dr_J1FJ;
    double Dr_JFJ_Min = Dr_J2FJ;
    double Dr_JmaxL = Dr_J1L;
    if (Dr_J2FJ > Dr_J1FJ) {
      Dr_JFJ_Max = Dr_J2FJ;
      Dr_JFJ_Min = Dr_J1FJ;
      Dr_JmaxL = Dr_J2L;
    }

    if (Dr_JFJ_Max < 1.5)
      continue;
    if (Dr_J1L <= 0.4 || Dr_J2L <= 0.4)
      continue;
    if (Dr_JmaxL > 3.5)
      continue;

    // ********************************************************** Trigger
    // Objects and Matchings Matching 1st
    bool matched_TRG_1 = false;

    bool matched_to_AK8PFJet230_SoftDropMass40 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if ((Trigger_Object_bit[itrg] & 4) == 4)
        if (sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet1_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 100) {
          matched_to_AK8PFJet230_SoftDropMass40 = true;
          break;
        }

    if (matched_to_AK8PFJet230_SoftDropMass40)
      matched_TRG_1 = true;

    // Matching 2nd
    bool matched_TRG_2 = false;
    if (HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06)
      matched_TRG_2 = true;

    // ********************************************************** Fill
    // Histograms

    _FatJet1_pt->Fill(FatJet1_pt);
    _FatJet1_eta->Fill(FatJet1_eta);
    _FatJet1_phi->Fill(FatJet1_phi);
    _FatJet1_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
    _FatJet1_Mass->Fill(FatJet1_Mass);
    _FatJet1_MassSD->Fill(FatJet1_MassSD);
    _FatJet1_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);
    _FatJet1PNetMD_Xbb->Fill(FatJet1PNetMD_Xbb);
    _FatJet1PNetMD_Xbb_Legacy->Fill(FatJet1PNetMD_Xbb_Legacy);
    double FatJet1PNetMD_Xbb_Legacy_AN =
        FatJet1PNetMD_Xbb_Legacy /
        (FatJet1PNetMD_Xbb_Legacy + FatJet1PNetMD_QCD_Legacy);
    _FatJet1PNetMD_Xbb_Legacy_AN->Fill(FatJet1PNetMD_Xbb_Legacy_AN);

    if (matched_TRG_1) {
      _FatJet2_pt->Fill(FatJet1_pt);
      _FatJet2_eta->Fill(FatJet1_eta);
      _FatJet2_phi->Fill(FatJet1_phi);
      _FatJet2_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
      _FatJet2_Mass->Fill(FatJet1_Mass);
      _FatJet2_MassSD->Fill(FatJet1_MassSD);
      _FatJet2_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);
      _FatJet2PNetMD_Xbb->Fill(FatJet1PNetMD_Xbb);
      _FatJet2PNetMD_Xbb_Legacy->Fill(FatJet1PNetMD_Xbb_Legacy);
      _FatJet2PNetMD_Xbb_Legacy_AN->Fill(FatJet1PNetMD_Xbb_Legacy_AN);
    }

    if (matched_TRG_1 && matched_TRG_2) {
      _FatJet3_pt->Fill(FatJet1_pt);
      _FatJet3_eta->Fill(FatJet1_eta);
      _FatJet3_phi->Fill(FatJet1_phi);
      _FatJet3_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
      _FatJet3_Mass->Fill(FatJet1_Mass);
      _FatJet3_MassSD->Fill(FatJet1_MassSD);
      _FatJet3_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);
      _FatJet3PNetMD_Xbb->Fill(FatJet1PNetMD_Xbb);
      _FatJet3PNetMD_Xbb_Legacy->Fill(FatJet1PNetMD_Xbb_Legacy);
      _FatJet3PNetMD_Xbb_Legacy_AN->Fill(FatJet1PNetMD_Xbb_Legacy_AN);
    }

    _MET->Fill(MET);
    _Lep1_Pt->Fill(lep1_Pt);
    _Dr_LFJ->Fill(Dr_LFJ);
    _Dr_J1FJ->Fill(Dr_JFJ_Max);
    _Dr_J2FJ->Fill(Dr_JFJ_Min);
    _Dr_JmaxL->Fill(Dr_JmaxL);

  } // end event loop

  f->Write();
  std::cout << "Done. Written to Muon_Histograms_" + Run + ".root" << std::endl;
}
