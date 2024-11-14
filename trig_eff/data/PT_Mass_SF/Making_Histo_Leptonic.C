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
#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <vector>

#define ARR_SIZE 10000

// D-phi
double phi_dist(double a, double b) {
  if (fabs(a - b) > 3.14159265) {
    return 6.2831853 - fabs(a - b);
  }
  return fabs(a - b);
}

double get_dR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = phi_dist(phi1, phi2);
  return sqrt(pow(deta, 2) + pow(dphi, 2));
}

bool inRange(int low, int high, int x) { return (low <= x && x <= high); }

void Making_Histo_Leptonic(
    const std::string &run_tag,              // e.g. 2023C, 2023D
    const std::string &channel,              // e.g. Muon, EGamma
    const std::string &sample_path,          // path to the root file
    const std::string &output_path,          // path to the output root file
    const std::string &jec_path_L2Relative,  // path to the AK4 JEC txt file
    const std::string &jec_path_L2L3Residual // path to the AK8 JEC txt file
) {
  // Load appropriate Good Lumi List based on era
  // TODO: add more eras if needed
  map<int, vector<pair<int, int>>> Good_Lumis;
  if (run_tag == "2023C") {
    Good_Lumis = {
#include "GoodLumiList/GoodLumiList_Map_2023C.txt"
    };
  } else if (run_tag == "2023D") {
    Good_Lumis = {
#include "GoodLumiList/GoodLumiList_Map_2023D.txt"
    };
  } else {
    throw std::runtime_error("Unsupported era for lumi dict: " + run_tag);
  }

  gSystem->Load("libFWCoreFWLite.so");

  // JEC
  vector<JetCorrectorParameters> vPar_AK8;
  vPar_AK8.push_back(JetCorrectorParameters(jec_path_L2Relative.c_str()));
  vPar_AK8.push_back(JetCorrectorParameters(jec_path_L2L3Residual.c_str()));
  FactorizedJetCorrector *corrector_AK8 = new FactorizedJetCorrector(vPar_AK8);

  TFile *f = new TFile(output_path.c_str(), "RECREATE");

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

  TH2D *_FatJet2_Pt_Mass_M = new TH2D("FatJet2_Pt_Mass_M", "FatJet2_Pt_Mass_M",
                                      45, Lower_pt, 15, Lower_m);

  TFile *f1 = new TFile(sample_path.c_str());

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
  Float_t FatJet1_rawFactor;

  Float_t FatJet2_pt;
  Float_t FatJet2_eta;
  Float_t FatJet2_phi;
  Float_t FatJet2_Mass;
  Float_t FatJet2_MassSD;
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
  for (int i = 0; i < InputTree_TrgObj->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);
    
    // HLT Selection
    if (channel == "EGamma") {
      if (!(HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) == 11)) {
            continue;
      }
    } else if (channel == "Muon") {
      if (!(HLT_IsoMu27 && fabs(lep1_Id) == 13)) {
        continue;
      }
    } else {
      throw std::invalid_argument("Invalid channel");
    }
    
    if (lep1_Pt < 50)
      continue;

    // FatJets correction and selection
    if (FatJet1_pt > 0) {
      double Raw_FatJet1_pt = FatJet1_pt*(1.0 - FatJet1_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_pt);
      corrector_AK8->setJetEta(FatJet1_eta);
      corrector_AK8->setJetPhi(FatJet1_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet1_pt     = Raw_FatJet1_pt * corr;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * corr;
     }
    if (FatJet2_pt > 0) {
      double Raw_FatJet2_pt = FatJet2_pt*(1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet2_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      corrector_AK8->setJetPhi(FatJet2_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt* corr;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * corr;
     }

    // FatJets selection
    if (FatJet2_pt > 180)
      continue;
    if (FatJet1_pt < 160)
      continue;
    if (lep1_Pt < 55)
      continue;
    if (lep2_Pt > 30)
      continue;
    if (phi_dist(FatJet1_phi, lep1_Phi) < 2.0)
      continue;
    // if (fabs(FatJet1_eta) > 1.4) continue;
    // if (fabs(FatJet1_eta) > 2.5 || fabs(FatJet1_eta) < 1.4) continue;
    if (MET < 50)
      continue;

    double dR_J1FJ = -1;
    double dR_J1L = 10;
    if (Jet1_Pt > 40) {
      dR_J1FJ = get_dR(Jet1_Eta, Jet1_Phi, FatJet1_eta, FatJet1_phi);
      dR_J1L = get_dR(Jet1_Eta, Jet1_Phi, lep1_Eta, lep1_Phi);
    }
    double dR_J2FJ = -1;
    double dR_J2L = 10;
    if (Jet2_Pt > 40) {
      dR_J2FJ = get_dR(Jet2_Eta, Jet2_Phi, FatJet1_eta, FatJet1_phi);
      dR_J2L = get_dR(Jet2_Eta, Jet2_Phi, lep1_Eta, lep1_Phi);
    }

    double dR_JFJ_Max = dR_J1FJ;
    double dR_JFJ_Min = dR_J2FJ;
    double dR_JmaxL = dR_J1L;
    if (dR_J2FJ > dR_J1FJ) {
      dR_JFJ_Max = dR_J2FJ;
      dR_JFJ_Min = dR_J1FJ;
      dR_JmaxL = dR_J2L;
    }

    if (dR_JFJ_Max < 0)
      continue;
    if (dR_J1L <= 0.4 || dR_J2L <= 0.4)
      continue;
    if (dR_JmaxL > 3.5)
      continue;

    // VBFTag veto
    // if(isVBFtag) continue;
    // JSON certification
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

    // Lepton selection or veto
    //   if (fabs(lep1_Id) !=11 ) continue;

    // Trigger objects and Matchings Probe Matched
    bool Probe_Matched = false;

    bool matched_to_AK8PFJet230_SoftDropMass40 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if ((Trigger_Object_bit[itrg] & 4) == 4) {
        double dR = get_dR(FatJet1_eta, FatJet1_phi, Trigger_Object_eta[itrg],
                           Trigger_Object_phi[itrg]);
        if (dR < 0.4 && Trigger_Object_pt[itrg] > 100) {
          matched_to_AK8PFJet230_SoftDropMass40 = true;
          break;
        }
      }

    if (matched_to_AK8PFJet230_SoftDropMass40) {
      Probe_Matched = true;
    }

    // Fill histograms
    _FatJet1_pt->Fill(FatJet1_pt);
    _FatJet1_eta->Fill(FatJet1_eta);
    _FatJet1_phi->Fill(FatJet1_phi);
    _FatJet1_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
    _FatJet1_Mass->Fill(FatJet1_Mass);
    _FatJet1_MassSD->Fill(FatJet1_MassSD);
    _FatJet1_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);

    if (Probe_Matched) {
      _FatJet2_Pt_Mass_M->Fill(FatJet1_pt, FatJet1_MassSD);
    }

  } // end event loop

  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
