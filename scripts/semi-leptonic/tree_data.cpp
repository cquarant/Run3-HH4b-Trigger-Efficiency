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
  Float_t corrected_pt;
  Float_t corrected_massSD;
};

JetCorrectionResult applyJEC(Float_t jet_pt, Float_t jet_eta, Float_t jet_phi,
                             Float_t jet_rawFactor, Float_t jet_massSD,
                             FactorizedJetCorrector *corrector) {
  JetCorrectionResult result{jet_pt, jet_massSD};

  if (jet_pt > 0) {
    Float_t raw_pt = jet_pt * (1.0 - jet_rawFactor);
    corrector->setJetPt(raw_pt);
    corrector->setJetEta(jet_eta);
    corrector->setJetPhi(jet_phi);
    Float_t correction_factor = corrector->getCorrection();

    result.corrected_pt = raw_pt * correction_factor;
    result.corrected_massSD =
        jet_massSD * (1.0 - jet_rawFactor) * correction_factor;
  }

  return result;
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
      T_fatJet1_Tau3OverTau2;
  Float_t T_fatJet1_ParticleNetLegacy_Xbb, T_fatJet1_ParticleNetLegacy_XbbVsQCD,
      T_fatJet1_globalParT_XbbVsQCD;
  Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_msoftdrop;
  Float_t T_MET, T_lep1_pt, T_lep1_eta, T_lep1_phi, T_dR_LFJ, T_dR_JFJ_max,
      T_dR_JFJ_min, T_dR_JmaxL;

  outputTree->Branch("weight", &T_weight, "weight/F");
  outputTree->Branch("fatJet1_pt", &T_fatJet1_pt, "fatJet1_pt/F");
  outputTree->Branch("fatJet1_eta", &T_fatJet1_eta, "fatJet1_eta/F");
  outputTree->Branch("fatJet1_phi", &T_fatJet1_phi, "fatJet1_phi/F");
  outputTree->Branch("fatJet1_msoftdrop", &T_fatJet1_msoftdrop,
                     "fatJet1_msoftdrop/F");
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
  outputTree->Branch("MET", &T_MET, "MET/F");
  outputTree->Branch("lep1_pt", &T_lep1_pt, "lep1_pt/F");
  outputTree->Branch("lep1_eta", &T_lep1_eta, "lep1_eta/F");
  outputTree->Branch("lep1_phi", &T_lep1_phi, "lep1_phi/F");
  outputTree->Branch("dR_LFJ", &T_dR_LFJ, "dR_LFJ/F");
  outputTree->Branch("dR_J1FJ", &T_dR_JFJ_max, "dR_J1FJ/F");
  outputTree->Branch("dR_J2FJ", &T_dR_JFJ_min, "dR_J2FJ/F");
  outputTree->Branch("dR_JmaxL", &T_dR_JmaxL, "dR_JmaxL/F");

  TFile *f1 = new TFile(sample_path.c_str());

  TTree *InputTree = (TTree *)f1->Get("tree");
  TH1F *NEvents = (TH1F *)f1->Get("nPU_True");

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

  Float_t lep1_pt;
  Float_t lep1_eta;
  Float_t lep1_phi;
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

  // FatJet 1 kinematics
  Float_t FatJet1_pt;
  Float_t FatJet1_eta;
  Float_t FatJet1_phi;
  Float_t FatJet1_Mass;
  Float_t FatJet1_MassSD;
  Float_t FatJet1_rawFactor;
  Float_t FatJet1_Tau3OverTau2;
  // Fatjet 1 ParticleNet scores
  Float_t FatJet1PNet_QCD;
  Float_t FatJet1PNet_QCD0HF;
  Float_t FatJet1PNet_QCD1HF;
  Float_t FatJet1PNet_QCD2HF;
  Float_t FatJet1PNet_XbbVsQCD;
  Float_t FatJet1PNet_XccVsQCD;
  Float_t FatJet1PNet_XggVsQCD;
  Float_t FatJet1PNet_XqqVsQCD;
  // Fatjet 1 ParticleNetLegacy scores
  Float_t FatJet1PNetLegacy_Xbb;
  Float_t FatJet1PNetLegacy_Xcc;
  Float_t FatJet1PNetLegacy_Xqq;
  Float_t FatJet1PNetLegacy_QCD;
  Float_t FatJet1PNetLegacy_QCDb;
  Float_t FatJet1PNetLegacy_QCDbb;
  Float_t FatJet1PNetLegacy_QCDothers;
  // Fatjet 1 GloParT scores
  Float_t FatJet1GloParT_QCD0HF;
  Float_t FatJet1GloParT_QCD1HF;
  Float_t FatJet1GloParT_QCD2HF;
  Float_t FatJet1GloParT_Xbb;
  Float_t FatJet1GloParT_Xcc;
  Float_t FatJet1GloParT_Xqq;
  Float_t FatJet1GloParT_XbbVsQCD;
  Float_t FatJet1GloParT_massRes;
  Float_t FatJet1GloParT_massVis;

  // FatJet 2 kinematics
  Float_t FatJet2_pt;
  Float_t FatJet2_eta;
  Float_t FatJet2_phi;
  Float_t FatJet2_Mass;
  Float_t FatJet2_MassSD;
  Float_t FatJet2_rawFactor;
  Float_t FatJet2_Tau3OverTau2;
  // Fatjet 2 ParticleNet scores
  Float_t FatJet2PNet_QCD;
  Float_t FatJet2PNet_QCD0HF;
  Float_t FatJet2PNet_QCD1HF;
  Float_t FatJet2PNet_QCD2HF;
  Float_t FatJet2PNet_XbbVsQCD;
  Float_t FatJet2PNet_XccVsQCD;
  Float_t FatJet2PNet_XggVsQCD;
  Float_t FatJet2PNet_XqqVsQCD;
  // Fatjet 2 ParticleNetLegacy scores
  Float_t FatJet2PNetLegacy_Xbb;
  Float_t FatJet2PNetLegacy_Xcc;
  Float_t FatJet2PNetLegacy_Xqq;
  Float_t FatJet2PNetLegacy_QCD;
  Float_t FatJet2PNetLegacy_QCDb;
  Float_t FatJet2PNetLegacy_QCDbb;
  Float_t FatJet2PNetLegacy_QCDothers;
  // Fatjet 2 GloParT scores
  Float_t FatJet2GloParT_QCD0HF;
  Float_t FatJet2GloParT_QCD1HF;
  Float_t FatJet2GloParT_QCD2HF;
  Float_t FatJet2GloParT_Xbb;
  Float_t FatJet2GloParT_Xcc;
  Float_t FatJet2GloParT_Xqq;
  Float_t FatJet2GloParT_XbbVsQCD;
  Float_t FatJet2GloParT_massRes;
  Float_t FatJet2GloParT_massVis;

  // FatJet 3 kinematics
  Float_t FatJet3_pt;
  Float_t FatJet3_eta;
  Float_t FatJet3_phi;
  Float_t FatJet3_Mass;
  Float_t FatJet3_MassSD;
  Float_t FatJet3_rawFactor;
  Float_t FatJet3_Tau3OverTau2;

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

  InputTree->SetBranchAddress("lep1Pt", &lep1_pt);
  InputTree->SetBranchAddress("lep1Eta", &lep1_eta);
  InputTree->SetBranchAddress("lep1Phi", &lep1_phi);
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

  // FatJet 1 kinematics
  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_mass", &FatJet1_Mass);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);
  InputTree->SetBranchAddress("fatJet1_Tau3OverTau2", &FatJet1_Tau3OverTau2);
  // FatJet 1 ParticleNet scores
  InputTree->SetBranchAddress("fatJet1_particleNet_QCD", &FatJet1PNet_QCD);
  InputTree->SetBranchAddress("fatJet1_particleNet_QCD0HF",
                              &FatJet1PNet_QCD0HF);
  InputTree->SetBranchAddress("fatJet1_particleNet_QCD1HF",
                              &FatJet1PNet_QCD1HF);
  InputTree->SetBranchAddress("fatJet1_particleNet_QCD2HF",
                              &FatJet1PNet_QCD2HF);
  InputTree->SetBranchAddress("fatJet1_particleNet_XbbVsQCD",
                              &FatJet1PNet_XbbVsQCD);
  InputTree->SetBranchAddress("fatJet1_particleNet_XccVsQCD",
                              &FatJet1PNet_XccVsQCD);
  InputTree->SetBranchAddress("fatJet1_particleNet_XggVsQCD",
                              &FatJet1PNet_XggVsQCD);
  InputTree->SetBranchAddress("fatJet1_particleNet_XqqVsQCD",
                              &FatJet1PNet_XqqVsQCD);
  // FatJet 1 ParticleNetLegacy scores
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xbb",
                              &FatJet1PNetLegacy_Xbb);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xcc",
                              &FatJet1PNetLegacy_Xcc);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xqq",
                              &FatJet1PNetLegacy_Xqq);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCD",
                              &FatJet1PNetLegacy_QCD);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDb",
                              &FatJet1PNetLegacy_QCDb);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDbb",
                              &FatJet1PNetLegacy_QCDbb);
  InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDothers",
                              &FatJet1PNetLegacy_QCDothers);
  // FatJet 1 GloParT scores
  InputTree->SetBranchAddress("fatJet1_globalParT_QCD0HF",
                              &FatJet1GloParT_QCD0HF);
  InputTree->SetBranchAddress("fatJet1_globalParT_QCD1HF",
                              &FatJet1GloParT_QCD1HF);
  InputTree->SetBranchAddress("fatJet1_globalParT_QCD2HF",
                              &FatJet1GloParT_QCD2HF);
  InputTree->SetBranchAddress("fatJet1_globalParT_Xbb", &FatJet1GloParT_Xbb);
  InputTree->SetBranchAddress("fatJet1_globalParT_Xcc", &FatJet1GloParT_Xcc);
  InputTree->SetBranchAddress("fatJet1_globalParT_Xqq", &FatJet1GloParT_Xqq);
  InputTree->SetBranchAddress("fatJet1_globalParT_XbbVsQCD",
                              &FatJet1GloParT_XbbVsQCD);
  InputTree->SetBranchAddress("fatJet1_globalParT_massRes",
                              &FatJet1GloParT_massRes);
  InputTree->SetBranchAddress("fatJet1_globalParT_massVis",
                              &FatJet1GloParT_massVis);

  // FatJet 2 kinematics
  InputTree->SetBranchAddress("fatJet2_pt", &FatJet2_pt);
  InputTree->SetBranchAddress("fatJet2_eta", &FatJet2_eta);
  InputTree->SetBranchAddress("fatJet2_phi", &FatJet2_phi);
  InputTree->SetBranchAddress("fatJet2_mass", &FatJet2_Mass);
  InputTree->SetBranchAddress("fatJet2_msoftdrop", &FatJet2_MassSD);
  InputTree->SetBranchAddress("fatJet2_rawFactor", &FatJet2_rawFactor);
  InputTree->SetBranchAddress("fatJet2_Tau3OverTau2", &FatJet2_Tau3OverTau2);
  // FatJet 2 ParticleNet scores
  InputTree->SetBranchAddress("fatJet2_particleNet_QCD", &FatJet2PNet_QCD);
  InputTree->SetBranchAddress("fatJet2_particleNet_QCD0HF",
                              &FatJet2PNet_QCD0HF);
  InputTree->SetBranchAddress("fatJet2_particleNet_QCD1HF",
                              &FatJet2PNet_QCD1HF);
  InputTree->SetBranchAddress("fatJet2_particleNet_QCD2HF",
                              &FatJet2PNet_QCD2HF);
  InputTree->SetBranchAddress("fatJet2_particleNet_XbbVsQCD",
                              &FatJet2PNet_XbbVsQCD);
  InputTree->SetBranchAddress("fatJet2_particleNet_XccVsQCD",
                              &FatJet2PNet_XccVsQCD);
  InputTree->SetBranchAddress("fatJet2_particleNet_XggVsQCD",
                              &FatJet2PNet_XggVsQCD);
  InputTree->SetBranchAddress("fatJet2_particleNet_XqqVsQCD",
                              &FatJet2PNet_XqqVsQCD);
  // FatJet 2 ParticleNetLegacy scores
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xbb",
                              &FatJet2PNetLegacy_Xbb);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xcc",
                              &FatJet2PNetLegacy_Xcc);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xqq",
                              &FatJet2PNetLegacy_Xqq);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCD",
                              &FatJet2PNetLegacy_QCD);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCDb",
                              &FatJet2PNetLegacy_QCDb);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCDbb",
                              &FatJet2PNetLegacy_QCDbb);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCDothers",
                              &FatJet2PNetLegacy_QCDothers);
  // FatJet 2 GloParT scores
  InputTree->SetBranchAddress("fatJet2_globalParT_QCD0HF",
                              &FatJet2GloParT_QCD0HF);
  InputTree->SetBranchAddress("fatJet2_globalParT_QCD1HF",
                              &FatJet2GloParT_QCD1HF);
  InputTree->SetBranchAddress("fatJet2_globalParT_QCD2HF",
                              &FatJet2GloParT_QCD2HF);
  InputTree->SetBranchAddress("fatJet2_globalParT_Xbb", &FatJet2GloParT_Xbb);
  InputTree->SetBranchAddress("fatJet2_globalParT_Xcc", &FatJet2GloParT_Xcc);
  InputTree->SetBranchAddress("fatJet2_globalParT_Xqq", &FatJet2GloParT_Xqq);
  InputTree->SetBranchAddress("fatJet2_globalParT_XbbVsQCD",
                              &FatJet2GloParT_XbbVsQCD);
  InputTree->SetBranchAddress("fatJet2_globalParT_massRes",
                              &FatJet2GloParT_massRes);
  InputTree->SetBranchAddress("fatJet2_globalParT_massVis",
                              &FatJet2GloParT_massVis);

  // FatJet 3 kinematics
  InputTree->SetBranchAddress("fatJet3_pt", &FatJet3_pt);
  InputTree->SetBranchAddress("fatJet3_eta", &FatJet3_eta);
  InputTree->SetBranchAddress("fatJet3_phi", &FatJet3_phi);
  InputTree->SetBranchAddress("fatJet3_mass", &FatJet3_Mass);
  InputTree->SetBranchAddress("fatJet3_msoftdrop", &FatJet3_MassSD);
  InputTree->SetBranchAddress("fatJet3_rawFactor", &FatJet3_rawFactor);
  InputTree->SetBranchAddress("fatJet3_Tau3OverTau2", &FatJet3_Tau3OverTau2);

  // Events Loop
  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);

    if (year == "2023") {
      // for 2023 HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06 started operations
      // from 367661
      if (run < 367661) {
        continue;
      }
    }

    // HLT Selection
    bool HLT_ele = (HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) == 11);
    bool HLT_mu = (HLT_IsoMu27 && fabs(lep1_Id) == 13);
    if (channel_lower == "egamma" or channel_lower == "electron") {
      if (!HLT_ele) {
        continue;
      }
    } else if (channel_lower == "muon") {
      if (!HLT_mu) {
        continue;
      }
    } else if (channel_lower == "lepton" or channel_lower == "leptonic") {
      if (!HLT_ele && !HLT_mu) {
        continue;
      }
    } else {
      throw std::invalid_argument("Invalid channel: " + channel);
    }

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

    JetCorrectionResult jec3 =
        applyJEC(FatJet3_pt, FatJet3_eta, FatJet3_phi, FatJet3_rawFactor,
                 FatJet3_MassSD, corrector);
    FatJet3_pt = jec3.corrected_pt;
    FatJet3_MassSD = jec3.corrected_massSD;

    // Selection
    // EGamma, Muon, Lepton (EGamma + Muon)
    if (FatJet1_pt < 250 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50) {
      continue;
    }
    if (lep1_pt < 50 || lep2_Pt > 30) {
      continue;
    }
    if (FatJet2_pt > 200 && FatJet2_MassSD > 50) {
      continue;
    }
    if (MET < 50) {
      continue;
    }

    double dR_LFJ = get_dR(lep1_eta, lep1_phi, FatJet1_eta, FatJet1_phi);
    if (dR_LFJ < 1.5) {
      continue;
    }

    double dR_J1FJ = -1;
    double dR_J1L = 10;
    if (Jet1_Pt > 40) {
      dR_J1FJ = get_dR(Jet1_Eta, Jet1_Phi, FatJet1_eta, FatJet1_phi);
      dR_J1L = get_dR(Jet1_Eta, Jet1_Phi, lep1_eta, lep1_phi);
    }
    double dR_J2FJ = -1;
    double dR_J2L = 10;
    if (Jet2_Pt > 40) {
      dR_J2FJ = get_dR(Jet2_Eta, Jet2_Phi, FatJet1_eta, FatJet1_phi);
      dR_J2L = get_dR(Jet2_Eta, Jet2_Phi, lep1_eta, lep1_phi);
    }

    double dR_JFJ_max = dR_J1FJ;
    double dR_JFJ_min = dR_J2FJ;
    double dR_JmaxL = dR_J1L;
    if (dR_J2FJ > dR_J1FJ) {
      dR_JFJ_max = dR_J2FJ;
      dR_JFJ_min = dR_J1FJ;
      dR_JmaxL = dR_J2L;
    }

    if (dR_J1L <= 0.4 || dR_J2L <= 0.4) {
      continue;
    }

    T_weight = 1.0;
    T_fatJet1_pt = FatJet1_pt;
    T_fatJet1_eta = FatJet1_eta;
    T_fatJet1_phi = FatJet1_phi;
    T_fatJet1_msoftdrop = FatJet1_MassSD;
    T_fatJet1_Tau3OverTau2 = FatJet1_Tau3OverTau2;
    T_fatJet1_ParticleNetLegacy_Xbb = FatJet1PNetLegacy_Xbb;
    double FatJet1PNetLegacy_XbbVsQCD =
        FatJet1PNetLegacy_Xbb / (FatJet1PNetLegacy_Xbb + FatJet1PNetLegacy_QCD);
    T_fatJet1_ParticleNetLegacy_XbbVsQCD = FatJet1PNetLegacy_XbbVsQCD;
    T_fatJet1_globalParT_XbbVsQCD = FatJet1GloParT_XbbVsQCD;
    T_fatJet2_pt = FatJet2_pt;
    T_fatJet2_eta = FatJet2_eta;
    T_fatJet2_msoftdrop = FatJet2_MassSD;

    T_MET = MET;
    T_lep1_pt = lep1_pt;
    T_lep1_eta = lep1_eta;
    T_lep1_phi = lep1_phi;
    T_dR_LFJ = dR_LFJ;
    T_dR_JFJ_max = dR_JFJ_max;
    T_dR_JFJ_min = dR_JFJ_min;
    T_dR_JmaxL = dR_JmaxL;

    outputTree->Fill();

  } // end event loop

  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
