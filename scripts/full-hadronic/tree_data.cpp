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

#define ARR_SIZE 1000

std::string to_lower(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  return str;
}

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

struct JetCorrectionResult {
  double corrected_pt;
  double corrected_massSD;
};

JetCorrectionResult applyJEC(Double_t jet_pt, Double_t jet_eta, Double_t jet_phi,
                             Double_t jet_rawFactor, Double_t jet_massSD,
                             FactorizedJetCorrector *corrector) {
  JetCorrectionResult result{jet_pt, jet_massSD};

  if (jet_pt > 0) {
    double raw_pt = jet_pt * (1.0 - jet_rawFactor);
    corrector->setJetPt(raw_pt);
    corrector->setJetEta(jet_eta);
    corrector->setJetPhi(jet_phi);
    double correction_factor = corrector->getCorrection();

    result.corrected_pt = raw_pt * correction_factor;
    result.corrected_massSD =
        jet_massSD * (1.0 - jet_rawFactor) * correction_factor;
  }

  return result;
}

bool trigger_matching(int trigger_bit,
                      Int_t NTrigger_Objects, 
                      Int_t *Trigger_Object_bit,
                      Float_t *Trigger_Object_pt,
                      Float_t *Trigger_Object_eta,
                      Float_t *Trigger_Object_phi,
                      Float_t FatJet_eta,
                      Float_t FatJet_phi,
                      Float_t trigger_pt_threshold = 100.0,
                      Float_t trigger_dR_threshold = 0.4) {
  bool matched = false;
  for (int itrg=0; itrg < NTrigger_Objects; itrg++) {
    if ((Trigger_Object_bit[itrg] & trigger_bit) == trigger_bit) {
      double dR = get_dR(FatJet_eta, FatJet_phi, Trigger_Object_eta[itrg], Trigger_Object_phi[itrg]);
      if (dR < trigger_dR_threshold && Trigger_Object_pt[itrg] > trigger_pt_threshold) {
        matched = true; break;
      }
    }
  }
  return matched;
}


void tree_data(
    const std::string &year,                 // 2022, 2023
    const std::string &run_tag,              // e.g. 2023C, 2023D
    const std::string &channel,              // e.g. Muon, EGamma, JetMET
    const std::string &sample_path,          // path to the root file
    const std::string &output_path,          // path to the output root file
    const std::string &jec_path_L2Relative,  // path to the AK4 JEC txt file
    const std::string &jec_path_L2L3Residual // path to the AK8 JEC txt file
) {
  // Load appropriate Good Lumi List based on era
  // TODO: add more eras if needed
  map<int, vector<pair<int, int>>> Good_Lumis;
  if (run_tag == "2022C") {
    Good_Lumis = {
#include "GoodLumiList/GoodLumiList_Map_2022C.txt"
    };
  } else if (run_tag == "2022D") {
    Good_Lumis = {
#include "GoodLumiList/GoodLumiList_Map_2022D.txt"
    };
  } else if (run_tag == "2022E") {
    Good_Lumis = {
#include "GoodLumiList/GoodLumiList_Map_2022E.txt"
    };
  } else if (run_tag == "2022F") {
    Good_Lumis = {
#include "GoodLumiList/GoodLumiList_Map_2022F.txt"
    };
  } else if (run_tag == "2022G") {
    Good_Lumis = {
#include "GoodLumiList/GoodLumiList_Map_2022G.txt"
    };
  } else if (run_tag == "2023C") {
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

  std::string channel_lower = to_lower(channel);
  std::cout << "Run tag: " << run_tag << std::endl;
  std::cout << "Channel: " << channel << std::endl;
  std::cout << "Sample Path: " << sample_path << std::endl;
  std::cout << "Output Path: " << output_path << std::endl;
  std::cout << "L2Relative JEC Path: " << jec_path_L2Relative << std::endl;
  std::cout << "L2L3Residual JEC Path: " << jec_path_L2L3Residual << std::endl;

  // JEC
  vector<JetCorrectorParameters> vPar_AK8;
  vPar_AK8.push_back(JetCorrectorParameters(jec_path_L2Relative.c_str()));
  vPar_AK8.push_back(JetCorrectorParameters(jec_path_L2L3Residual.c_str()));
  FactorizedJetCorrector *corrector = new FactorizedJetCorrector(vPar_AK8);

  TFile *f = new TFile(output_path.c_str(), "RECREATE");
  TTree *outputTree = new TTree("tree", "");

  Float_t T_weight;
  Float_t T_fatJet1_pt, T_fatJet1_eta, T_fatJet1_phi, T_fatJet1_msoftdrop,
      T_fatJet1_GloParT_massVis, T_fatJet1_GloParT_massRes, T_fatJet1_Tau3OverTau2;
  Float_t T_fatJet1_ParticleNetLegacy_Xbb, T_fatJet1_ParticleNetLegacy_XbbVsQCD,
      T_fatJet1_globalParT_XbbVsQCD;
  Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_msoftdrop, T_fatJet2_globalParT_XbbVsQCD, 
      T_fatJet2_Tau3OverTau2, T_fatJet2_GloParT_massVis, T_fatJet2_GloParT_massRes, 
      T_fatJet2_GloParT_XbbVsQCD;
  Float_t T_pTjj;

  outputTree->Branch("weight", &T_weight, "weight/F");
  outputTree->Branch("fatJet1_pt", &T_fatJet1_pt, "fatJet1_pt/F");
  outputTree->Branch("fatJet1_eta", &T_fatJet1_eta, "fatJet1_eta/F");
  outputTree->Branch("fatJet1_phi", &T_fatJet1_phi, "fatJet1_phi/F");
  outputTree->Branch("fatJet1_msoftdrop", &T_fatJet1_msoftdrop,
                     "fatJet1_msoftdrop/F");
  outputTree->Branch("fatJet1_GloParT_massVis", &T_fatJet1_GloParT_massVis,
                     "fatJet1_GloParT_massVis/F");
  outputTree->Branch("fatJet1_GloParT_massRes", &T_fatJet1_GloParT_massRes,
                     "fatJet1_GloParT_massRes/F");
  outputTree->Branch("fatJet1_Tau3OverTau2", &T_fatJet1_Tau3OverTau2,
                     "fatJet1_Tau3OverTau2/F");
  outputTree->Branch("fatJet1_ParticleNetLegacy_Xbb",
                     &T_fatJet1_ParticleNetLegacy_Xbb,
                     "fatJet1_ParticleNetLegacy_Xbb/F");
  outputTree->Branch("fatJet1_ParticleNetLegacy_XbbVsQCD",
                     &T_fatJet1_ParticleNetLegacy_XbbVsQCD,
                     "T_fatJet1_ParticleNetLegacy_XbbVsQCD/F");
  outputTree->Branch("fatJet1_GloParT_XbbVsQCD", &T_fatJet1_globalParT_XbbVsQCD,
                     "fatJet1_GloParT_XbbVsQCD/F");
  outputTree->Branch("fatJet2_pt", &T_fatJet2_pt, "fatJet2_pt/F");
  outputTree->Branch("fatJet2_eta", &T_fatJet2_eta, "fatJet2_eta/F");
  outputTree->Branch("fatJet2_msoftdrop", &T_fatJet2_msoftdrop,
                     "fatJet2_msoftdrop/F");
  outputTree->Branch("fatJet2_GloParT_massVis", &T_fatJet2_GloParT_massVis,
                     "fatJet2_GloParT_massVis/F");
  outputTree->Branch("fatJet2_GloParT_massRes", &T_fatJet2_GloParT_massRes,
                     "fatJet2_GloParT_massRes/F");
  outputTree->Branch("fatJet2_Tau3OverTau2", &T_fatJet2_Tau3OverTau2,
                     "fatJet2_Tau3OverTau2/F");
  outputTree->Branch("fatJet2_GloParT_XbbVsQCD", &T_fatJet2_globalParT_XbbVsQCD,
                     "fatJet2_GloParT_XbbVsQCD/F");
  outputTree->Branch("pTjj", &T_pTjj, "pTjj/F");

  TFile *f1 = new TFile(sample_path.c_str());

  TTree *InputTree = (TTree *)f1->Get("tree");
  TH1F *NEvents = (TH1F *)f1->Get("nPU_True");

  UInt_t run;
  UInt_t lumi;
  Int_t isVBFtag;

  Bool_t HLT_AK8PFJet230_SoftDropMass40;
  Bool_t HLT_AK8PFJet425_SoftDropMass40;

  Float_t Jet1_Pt;
  Float_t Jet1_Eta;
  Float_t Jet1_Phi;

  Float_t Jet2_Pt;
  Float_t Jet2_Eta;
  Float_t Jet2_Phi;

  // FatJet 1 kinematics
  Float_t FatJet1_pt;
  Float_t FatJet1_eta;
  Float_t FatJet1_phi;
  Float_t FatJet1_Mass;
  Float_t FatJet1_MassSD;
  Float_t FatJet1_rawFactor;
  Float_t FatJet1_Tau3OverTau2;
  // Fatjet 1 ParticleNetLegacy scores
  Float_t FatJet1_PNetLegacy_Xbb;
  Float_t FatJet1_PNetLegacy_Xcc;
  Float_t FatJet1_PNetLegacy_Xqq;
  Float_t FatJet1_PNetLegacy_QCD;
  Float_t FatJet1_PNetLegacy_QCDb;
  Float_t FatJet1_PNetLegacy_QCDbb;
  Float_t FatJet1_PNetLegacy_QCDothers;
  // Fatjet 1 GloParT scores
  Float_t FatJet1_GloParT_QCD0HF;
  Float_t FatJet1_GloParT_QCD1HF;
  Float_t FatJet1_GloParT_QCD2HF;
  Float_t FatJet1_GloParT_Xbb;
  Float_t FatJet1_GloParT_Xcc;
  Float_t FatJet1_GloParT_Xqq;
  Float_t FatJet1_GloParT_XbbVsQCD;
  Float_t FatJet1_GloParT_massRes;
  Float_t FatJet1_GloParT_massVis;

  // FatJet 2 kinematics
  Float_t FatJet2_pt;
  Float_t FatJet2_eta;
  Float_t FatJet2_phi;
  Float_t FatJet2_Mass;
  Float_t FatJet2_MassSD;
  Float_t FatJet2_GloParT_XbbVsQCD;
  Float_t FatJet2_GloParT_massRes;
  Float_t FatJet2_GloParT_massVis;
  Float_t FatJet2_rawFactor;
  Float_t FatJet2_Tau3OverTau2;

  InputTree->SetBranchAddress("run", &run);
  InputTree->SetBranchAddress("lumi", &lumi);
  InputTree->SetBranchAddress("isVBFtag", &isVBFtag);

  InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);
  InputTree->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass40", &HLT_AK8PFJet425_SoftDropMass40);


  // FatJet 1 kinematics
  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_mass", &FatJet1_Mass);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);
  InputTree->SetBranchAddress("fatJet1_Tau3OverTau2", &FatJet1_Tau3OverTau2);
  // FatJet 1 ParticleNetLegacy scores
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xbb",
                              &FatJet1_PNetLegacy_Xbb);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xcc",
                              &FatJet1_PNetLegacy_Xcc);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xqq",
                              &FatJet1_PNetLegacy_Xqq);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCD",
                              &FatJet1_PNetLegacy_QCD);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDb",
                              &FatJet1_PNetLegacy_QCDb);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDbb",
                              &FatJet1_PNetLegacy_QCDbb);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDothers",
                              &FatJet1_PNetLegacy_QCDothers);
  // FatJet 1 GloParT scores
  InputTree->SetBranchAddress("fatJet1_globalParT_QCD0HF",
                              &FatJet1_GloParT_QCD0HF);
  InputTree->SetBranchAddress("fatJet1_globalParT_QCD1HF",
                              &FatJet1_GloParT_QCD1HF);
  InputTree->SetBranchAddress("fatJet1_globalParT_QCD2HF",
                              &FatJet1_GloParT_QCD2HF);
  InputTree->SetBranchAddress("fatJet1_globalParT_Xbb", &FatJet1_GloParT_Xbb);
  InputTree->SetBranchAddress("fatJet1_globalParT_Xcc", &FatJet1_GloParT_Xcc);
  InputTree->SetBranchAddress("fatJet1_globalParT_Xqq", &FatJet1_GloParT_Xqq);
  InputTree->SetBranchAddress("fatJet1_globalParT_XbbVsQCD",
                              &FatJet1_GloParT_XbbVsQCD);
  InputTree->SetBranchAddress("fatJet1_globalParT_massRes",
                              &FatJet1_GloParT_massRes);
  InputTree->SetBranchAddress("fatJet1_globalParT_massVis",
                              &FatJet1_GloParT_massVis);

  // FatJet 2 kinematics
  InputTree->SetBranchAddress("fatJet2_pt", &FatJet2_pt);
  InputTree->SetBranchAddress("fatJet2_eta", &FatJet2_eta);
  InputTree->SetBranchAddress("fatJet2_phi", &FatJet2_phi);
  InputTree->SetBranchAddress("fatJet2_mass", &FatJet2_Mass);
  InputTree->SetBranchAddress("fatJet2_msoftdrop", &FatJet2_MassSD);
  InputTree->SetBranchAddress("fatJet2_rawFactor", &FatJet2_rawFactor);
  InputTree->SetBranchAddress("fatJet2_Tau3OverTau2", &FatJet2_Tau3OverTau2);
  InputTree->SetBranchAddress("fatJet2_globalParT_XbbVsQCD",
                              &FatJet2_GloParT_XbbVsQCD);
  InputTree->SetBranchAddress("fatJet2_globalParT_massRes", &FatJet2_GloParT_massRes);
  InputTree->SetBranchAddress("fatJet2_globalParT_massVis", &FatJet2_GloParT_massVis);

  // Trigger Objects
  TTree *InputTree_TrgObj = (TTree *)f1->Get("tree_TrgObj");
  Int_t NTrigger_Objects;
  Float_t Trigger_Object_pt[ARR_SIZE];
  Float_t Trigger_Object_eta[ARR_SIZE];
  Float_t Trigger_Object_phi[ARR_SIZE];
  Int_t Trigger_Object_bit[ARR_SIZE];
  InputTree_TrgObj->SetBranchAddress("NTrigger_Objects", &NTrigger_Objects);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_pt", Trigger_Object_pt);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_eta", Trigger_Object_eta);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_phi", Trigger_Object_phi);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_bit", Trigger_Object_bit);

  // Events Loop
  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);

    // Trigger
    if (HLT_AK8PFJet425_SoftDropMass40 == 0) continue;

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
    if (!Certified) {
      continue;
    }

    // FatJets correction (JEC)
    JetCorrectionResult jec1 =
        applyJEC(FatJet1_pt, FatJet1_eta, FatJet1_phi, FatJet1_rawFactor,
                 FatJet1_MassSD, corrector);
    FatJet1_pt = jec1.corrected_pt;
    FatJet1_MassSD = jec1.corrected_massSD;

    JetCorrectionResult jec2 =
        applyJEC(FatJet2_pt, FatJet2_eta, FatJet2_phi, FatJet2_rawFactor,
                 FatJet2_MassSD, corrector);
    FatJet2_pt = jec2.corrected_pt;
    FatJet2_MassSD = jec2.corrected_massSD;

    // Selections
    // // pT > 450, |eta| < 2.5, massSD > 50, TXbb > 0.1, tau3/tau2 < 0.46
    if (FatJet1_pt <= 450 || fabs(FatJet1_eta) >= 2.5 || FatJet1_MassSD <= 50 || FatJet1_GloParT_XbbVsQCD <= 0.1 || FatJet1_Tau3OverTau2 >= 0.46) continue;
    if (FatJet2_pt <= 450 || fabs(FatJet2_eta) >= 2.5 || FatJet2_MassSD <= 50 || FatJet2_GloParT_XbbVsQCD <= 0.1 || FatJet2_Tau3OverTau2 >= 0.46) continue;
    // // pT > 450, |eta| < 2.5, massSD > 50, TXbb > 0.1, tau3/tau2 < 0.46
    // if (FatJet1_pt < 450 || fabs(FatJet1_eta) > 2.5 || FatJet1_MassSD < 50) continue;
    // if (FatJet2_pt < 450 || fabs(FatJet2_eta) > 2.5 || FatJet2_MassSD < 50) continue;

    int trigger_bit = 4;
    bool matched1 = trigger_matching(trigger_bit, NTrigger_Objects,
                                     Trigger_Object_bit, Trigger_Object_pt, Trigger_Object_eta, Trigger_Object_phi,
                                     FatJet1_eta, FatJet1_phi);
    if (!matched1) continue;
    bool matched2 = trigger_matching(trigger_bit, NTrigger_Objects,
                                     Trigger_Object_bit, Trigger_Object_pt, Trigger_Object_eta, Trigger_Object_phi,
                                     FatJet2_eta, FatJet2_phi);
    if (!matched2) continue;

    // Fill the tree
    T_weight = 1.0;
    T_fatJet1_pt = FatJet1_pt;
    T_fatJet1_eta = FatJet1_eta;
    T_fatJet1_phi = FatJet1_phi;
    T_fatJet1_msoftdrop = FatJet1_MassSD;
    T_fatJet1_GloParT_massVis = FatJet1_GloParT_massVis * FatJet1_Mass * (1.0 - FatJet1_rawFactor);
    T_fatJet1_GloParT_massRes = FatJet1_GloParT_massRes * FatJet1_Mass * (1.0 - FatJet1_rawFactor);
    T_fatJet1_Tau3OverTau2 = FatJet1_Tau3OverTau2;
    T_fatJet1_ParticleNetLegacy_Xbb = FatJet1_PNetLegacy_Xbb;
    double FatJet1_PNetLegacy_XbbVsQCD =
        FatJet1_PNetLegacy_Xbb / (FatJet1_PNetLegacy_Xbb + FatJet1_PNetLegacy_QCD);
    T_fatJet1_ParticleNetLegacy_XbbVsQCD = FatJet1_PNetLegacy_XbbVsQCD;
    T_fatJet1_globalParT_XbbVsQCD = FatJet1_GloParT_XbbVsQCD;
    T_fatJet2_pt = FatJet2_pt;
    T_fatJet2_eta = FatJet2_eta;
    T_fatJet2_msoftdrop = FatJet2_MassSD;
    T_fatJet2_GloParT_massVis = FatJet2_GloParT_massVis * FatJet2_Mass * (1.0 - FatJet2_rawFactor);
    T_fatJet2_GloParT_massRes = FatJet2_GloParT_massRes * FatJet2_Mass * (1.0 - FatJet2_rawFactor);
    T_fatJet2_Tau3OverTau2 = FatJet2_Tau3OverTau2;
    T_fatJet2_globalParT_XbbVsQCD = FatJet2_GloParT_XbbVsQCD;

    Float_t fatjet1_px = FatJet1_pt * cos(FatJet1_phi);
    Float_t fatjet1_py = FatJet1_pt * sin(FatJet1_phi);
    Float_t fatjet2_px = FatJet2_pt * cos(FatJet2_phi);
    Float_t fatjet2_py = FatJet2_pt * sin(FatJet2_phi);
    T_pTjj = sqrt(pow(fatjet1_px + fatjet2_px, 2) + pow(fatjet1_py + fatjet2_py, 2));

    outputTree->Fill();

  } // end event loop

  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
