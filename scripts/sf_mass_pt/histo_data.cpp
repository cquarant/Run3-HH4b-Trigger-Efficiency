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

bool checkTriggerMatching(Int_t NTrigger_Objects, 
                          const Float_t* Trigger_Object_pt, const Float_t* Trigger_Object_eta,
                          const Float_t* Trigger_Object_phi, const Int_t* Trigger_Object_bit,
                          float jet_eta, float jet_phi, int required_bits, float pt_threshold) {
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++) {
        if ((Trigger_Object_bit[itrg] & required_bits) == required_bits) {
            double dR = get_dR(jet_eta, jet_phi, Trigger_Object_eta[itrg], Trigger_Object_phi[itrg]);
            if (dR < 0.4 && Trigger_Object_pt[itrg] > pt_threshold) {
                return true;
            }
        }
    }
    return false;
}

bool checkTriggerMatching(Int_t NTrigger_Objects, 
                          const Float_t* Trigger_Object_pt, const Float_t* Trigger_Object_eta,
                          const Float_t* Trigger_Object_phi, const Int_t* Trigger_Object_bit,
                          float jet_eta, float jet_phi, float pt_threshold) {
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++) {
        if (Trigger_Object_bit[itrg] == 1) {
            double dR = get_dR(jet_eta, jet_phi, Trigger_Object_eta[itrg], Trigger_Object_phi[itrg]);
            if (dR < 0.4 && Trigger_Object_pt[itrg] > pt_threshold) {
                return true;
            }
        }
    }
    return false;
}

bool checkAK4JetRequirements(Float_t j1_pt, Float_t j1_eta, Float_t j1_phi,
                             Float_t j2_pt, Float_t j2_eta, Float_t j2_phi,
                             Float_t fj_eta, Float_t fj_phi, Float_t lep_eta,
                             Float_t lep_phi) {
  double dR_J1FJ = -1;
  double dR_J1L = 10;
  if (j1_pt > 40) {
    dR_J1FJ = get_dR(j1_eta, j1_phi, fj_eta, fj_phi);
    dR_J1L = get_dR(j1_eta, j1_phi, lep_eta, lep_phi);
  }

  double dR_J2FJ = -1;
  double dR_J2L = 10;
  if (j2_pt > 40) {
    dR_J2FJ = get_dR(j2_eta, j2_phi, fj_eta, fj_phi);
    dR_J2L = get_dR(j2_eta, j2_phi, lep_eta, lep_phi);
  }

  // Find max and min distances
  double dR_JFJ_Max = dR_J1FJ;
  double dR_JFJ_Min = dR_J2FJ;
  double dR_JmaxL = dR_J1L;
  if (dR_J2FJ > dR_J1FJ) {
    dR_JFJ_Max = dR_J2FJ;
    dR_JFJ_Min = dR_J1FJ;
    dR_JmaxL = dR_J2L;
  }

  // Apply cuts
  if (dR_JFJ_Max < 1.5)
    return false;
  if (dR_J1L <= 0.4 || dR_J2L <= 0.4)
    return false;
  if (dR_JmaxL > 3.5)
    return false;

  return true;
}

struct JetCorrectionResult {
  Float_t corrected_pt;
  Float_t corrected_massSD;
};

JetCorrectionResult applyJEC(Float_t jet_pt, Float_t jet_eta, Float_t jet_phi,
                             Float_t jet_rawFactor, Float_t jet_massSD,
                             FactorizedJetCorrector *corrector,
                             bool correct_mass = true) {
  JetCorrectionResult result{jet_pt, jet_massSD};

  if (jet_pt > 0) {
    Float_t raw_pt = jet_pt * (1.0 - jet_rawFactor);
    corrector->setJetPt(raw_pt);
    corrector->setJetEta(jet_eta);
    corrector->setJetPhi(jet_phi);
    Float_t correction_factor = corrector->getCorrection();

    result.corrected_pt = raw_pt * correction_factor;
    if (correct_mass) {
      result.corrected_massSD =
          jet_massSD * (1.0 - jet_rawFactor) * correction_factor;
    }
  }

  return result;
}

void histo_data(
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

  // Float_t bins_pt[9] = {300, 350, 400, 450, 500, 600, 700, 850, 1000};
  // int num_pt_bins = 8;

  // Float_t bins_m[8] = {60, 90, 120, 150, 180, 210, 240, 300};
  // int num_m_bins = 7;
  // Float_t bins_pt[46] =
  // {0.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,
  //                        80.0,  90.0,  100.0, 110.0, 120.0, 130.0, 140.0,
  //                        150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0,
  //                        220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0,
  //                        290.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0,
  //                        420.0, 440.0, 460.0, 480.0, 500.0, 550.0, 600.0,
  //                        700.0, 800.0, 1000.0};
  // int num_pt_bins = 45;

  // Float_t bins_m[16] = {0.0,  5.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,
  //                       80.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0,
  //                       350.0};
  // int num_m_bins = 15;

  Float_t bins_pt[9] = {250, 275, 300, 350, 400, 450, 500, 600, 1000};
  int num_pt_bins = 8;

  Float_t bins_m[10] = {50, 60, 80, 100, 120, 150, 200, 250, 300, 350};
  int num_m_bins = 9;

  // probe FatJet 1 kinematics
  TH1D *_ProbeJet_pass_pt =
      new TH1D("ProbeJet_pass_pt", "ProbeJet_pass_pt", 200, 0, 1000);
  TH1D *_ProbeJet_pass_eta =
      new TH1D("ProbeJet_pass_eta", "ProbeJet_pass_eta", 100, -5, 5);
  TH1D *_ProbeJet_pass_phi =
      new TH1D("ProbeJet_pass_phi", "ProbeJet_pass_phi", 100, -5, 5);
  TH2D *_ProbeJet_pass_eta_phi = new TH2D(
      "ProbeJet_pass_eta_phi", "ProbeJet_pass_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_ProbeJet_pass_Mass =
      new TH1D("ProbeJet_pass_Mass", "ProbeJet_pass_Mass", 500, 0, 500);
  TH1D *_ProbeJet_pass_MassSD =
      new TH1D("ProbeJet_pass_MassSD", "ProbeJet_pass_MassSD", 500, 0, 500);
  TH2D *_ProbeJet_pass_Mass_Pt =
      new TH2D("ProbeJet_pass_Mass_Pt", "ProbeJet_pass_Mass_Pt", num_m_bins,
               bins_m, num_pt_bins, bins_pt);
  // probe FatJet 1 ParticleNet scores
  TH1D *_ProbeJet_pass_PNet_QCD =
      new TH1D("ProbeJet_pass_PNet_QCD", "ProbeJet_pass_PNet_QCD", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_QCD0HF = new TH1D(
      "ProbeJet_pass_PNet_QCD0HF", "ProbeJet_pass_PNet_QCD0HF", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_QCD1HF = new TH1D(
      "ProbeJet_pass_PNet_QCD1HF", "ProbeJet_pass_PNet_QCD1HF", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_QCD2HF = new TH1D(
      "ProbeJet_pass_PNet_QCD2HF", "ProbeJet_pass_PNet_QCD2HF", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XbbVsQCD =
      new TH1D("ProbeJet_pass_PNet_XbbVsQCD", "ProbeJet_pass_PNet_XbbVsQCD",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XccVsQCD =
      new TH1D("ProbeJet_pass_PNet_XccVsQCD", "ProbeJet_pass_PNet_XccVsQCD",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XggVsQCD =
      new TH1D("ProbeJet_pass_PNet_XggVsQCD", "ProbeJet_pass_PNet_XggVsQCD",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XqqVsQCD =
      new TH1D("ProbeJet_pass_PNet_XqqVsQCD", "ProbeJet_pass_PNet_XqqVsQCD",
               100, 0, 1.0);
  // probe FatJet 1 ParticleNetLegacy scores
  TH1D *_ProbeJet_pass_PNetLegacy_Xbb =
      new TH1D("ProbeJet_pass_PNetLegacy_Xbb", "ProbeJet_pass_PNetLegacy_Xbb",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_XbbVsQCD =
      new TH1D("ProbeJet_pass_PNetLegacy_XbbVsQCD",
               "ProbeJet_pass_PNetLegacy_XbbVsQCD", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_Xcc =
      new TH1D("ProbeJet_pass_PNetLegacy_Xcc", "ProbeJet_pass_PNetLegacy_Xcc",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_Xqq =
      new TH1D("ProbeJet_pass_PNetLegacy_Xqq", "ProbeJet_pass_PNetLegacy_Xqq",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCD =
      new TH1D("ProbeJet_pass_PNetLegacy_QCD", "ProbeJet_pass_PNetLegacy_QCD",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCDb =
      new TH1D("ProbeJet_pass_PNetLegacy_QCDb", "ProbeJet_pass_PNetLegacy_QCDb",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCDbb =
      new TH1D("ProbeJet_pass_PNetLegacy_QCDbb",
               "ProbeJet_pass_PNetLegacy_QCDbb", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCDothers =
      new TH1D("ProbeJet_pass_PNetLegacy_QCDothers",
               "ProbeJet_pass_PNetLegacy_QCDothers", 100, 0, 1.0);
  // probe FatJet 1 GloParT scores
  TH1D *_ProbeJet_pass_GloParT_QCD0HF =
      new TH1D("ProbeJet_pass_GloParT_QCD0HF", "ProbeJet_pass_GloParT_QCD0HF",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_QCD1HF =
      new TH1D("ProbeJet_pass_GloParT_QCD1HF", "ProbeJet_pass_GloParT_QCD1HF",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_QCD2HF =
      new TH1D("ProbeJet_pass_GloParT_QCD2HF", "ProbeJet_pass_GloParT_QCD2HF",
               100, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_Xbb = new TH1D(
      "ProbeJet_pass_GloParT_Xbb", "ProbeJet_pass_GloParT_Xbb", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_Xcc = new TH1D(
      "ProbeJet_pass_GloParT_Xcc", "ProbeJet_pass_GloParT_Xcc", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_Xqq = new TH1D(
      "ProbeJet_pass_GloParT_Xqq", "ProbeJet_pass_GloParT_Xqq", 100, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_XbbVsQCD =
      new TH1D("ProbeJet_pass_GloParT_XbbVsQCD",
               "ProbeJet_pass_GloParT_XbbVsQCD", 100, 0, 1.0);

  // tag FatJet 1 kinematics
  TH1D *_ProbeJet_all_pt =
      new TH1D("ProbeJet_all_pt", "ProbeJet_all_pt", 200, 0, 1000);
  TH1D *_ProbeJet_all_eta =
      new TH1D("ProbeJet_all_eta", "ProbeJet_all_eta", 100, -5, 5);
  TH1D *_ProbeJet_all_phi =
      new TH1D("ProbeJet_all_phi", "ProbeJet_all_phi", 100, -5, 5);
  TH2D *_ProbeJet_all_eta_phi = new TH2D(
      "ProbeJet_all_eta_phi", "ProbeJet_all_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_ProbeJet_all_Mass =
      new TH1D("ProbeJet_all_Mass", "ProbeJet_all_Mass", 500, 0, 500);
  TH1D *_ProbeJet_all_MassSD =
      new TH1D("ProbeJet_all_MassSD", "ProbeJet_all_MassSD", 500, 0, 500);
  TH2D *_ProbeJet_all_Mass_Pt =
      new TH2D("ProbeJet_all_Mass_Pt", "ProbeJet_all_Mass_Pt", num_m_bins,
               bins_m, num_pt_bins, bins_pt);
  // tag FatJet 1 ParticleNet scores
  TH1D *_ProbeJet_all_PNet_QCD =
      new TH1D("ProbeJet_all_PNet_QCD", "ProbeJet_all_PNet_QCD", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_QCD0HF = new TH1D(
      "ProbeJet_all_PNet_QCD0HF", "ProbeJet_all_PNet_QCD0HF", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_QCD1HF = new TH1D(
      "ProbeJet_all_PNet_QCD1HF", "ProbeJet_all_PNet_QCD1HF", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_QCD2HF = new TH1D(
      "ProbeJet_all_PNet_QCD2HF", "ProbeJet_all_PNet_QCD2HF", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XbbVsQCD = new TH1D(
      "ProbeJet_all_PNet_XbbVsQCD", "ProbeJet_all_PNet_XbbVsQCD", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XccVsQCD = new TH1D(
      "ProbeJet_all_PNet_XccVsQCD", "ProbeJet_all_PNet_XccVsQCD", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XggVsQCD = new TH1D(
      "ProbeJet_all_PNet_XggVsQCD", "ProbeJet_all_PNet_XggVsQCD", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XqqVsQCD = new TH1D(
      "ProbeJet_all_PNet_XqqVsQCD", "ProbeJet_all_PNet_XqqVsQCD", 100, 0, 1.0);
  // tag FatJet 1 ParticleNetLegacy scores
  TH1D *_ProbeJet_all_PNetLegacy_Xbb =
      new TH1D("ProbeJet_all_PNetLegacy_Xbb", "ProbeJet_all_PNetLegacy_Xbb",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_XbbVsQCD =
      new TH1D("ProbeJet_all_PNetLegacy_XbbVsQCD",
               "ProbeJet_all_PNetLegacy_XbbVsQCD", 100, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_Xcc =
      new TH1D("ProbeJet_all_PNetLegacy_Xcc", "ProbeJet_all_PNetLegacy_Xcc",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_Xqq =
      new TH1D("ProbeJet_all_PNetLegacy_Xqq", "ProbeJet_all_PNetLegacy_Xqq",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCD =
      new TH1D("ProbeJet_all_PNetLegacy_QCD", "ProbeJet_all_PNetLegacy_QCD",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCDb =
      new TH1D("ProbeJet_all_PNetLegacy_QCDb", "ProbeJet_all_PNetLegacy_QCDb",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCDbb =
      new TH1D("ProbeJet_all_PNetLegacy_QCDbb", "ProbeJet_all_PNetLegacy_QCDbb",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCDothers =
      new TH1D("ProbeJet_all_PNetLegacy_QCDothers",
               "ProbeJet_all_PNetLegacy_QCDothers", 100, 0, 1.0);
  // tag FatJet 1 GloParT scores
  TH1D *_ProbeJet_all_GloParT_QCD0HF =
      new TH1D("ProbeJet_all_GloParT_QCD0HF", "ProbeJet_all_GloParT_QCD0HF",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_QCD1HF =
      new TH1D("ProbeJet_all_GloParT_QCD1HF", "ProbeJet_all_GloParT_QCD1HF",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_QCD2HF =
      new TH1D("ProbeJet_all_GloParT_QCD2HF", "ProbeJet_all_GloParT_QCD2HF",
               100, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_Xbb = new TH1D(
      "ProbeJet_all_GloParT_Xbb", "ProbeJet_all_GloParT_Xbb", 100, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_Xcc = new TH1D(
      "ProbeJet_all_GloParT_Xcc", "ProbeJet_all_GloParT_Xcc", 100, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_Xqq = new TH1D(
      "ProbeJet_all_GloParT_Xqq", "ProbeJet_all_GloParT_Xqq", 100, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_XbbVsQCD =
      new TH1D("ProbeJet_all_GloParT_XbbVsQCD", "ProbeJet_all_GloParT_XbbVsQCD",
               100, 0, 1.0);

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

  // probe FatJet kinematics
  Float_t ProbeJet_pt;
  Float_t ProbeJet_eta;
  Float_t ProbeJet_phi;
  Float_t ProbeJet_Mass;
  Float_t ProbeJet_MassSD;
  Float_t ProbeJet_rawFactor;
  // probe Fatjet ParticleNet scores
  Float_t ProbeJetPNet_QCD;
  Float_t ProbeJetPNet_QCD0HF;
  Float_t ProbeJetPNet_QCD1HF;
  Float_t ProbeJetPNet_QCD2HF;
  Float_t ProbeJetPNet_XbbVsQCD;
  Float_t ProbeJetPNet_XccVsQCD;
  Float_t ProbeJetPNet_XggVsQCD;
  Float_t ProbeJetPNet_XqqVsQCD;
  // probe Fatjet ParticleNetLegacy scores
  Float_t ProbeJetPNetLegacy_Xbb;
  Float_t ProbeJetPNetLegacy_Xcc;
  Float_t ProbeJetPNetLegacy_Xqq;
  Float_t ProbeJetPNetLegacy_QCD;
  Float_t ProbeJetPNetLegacy_QCDb;
  Float_t ProbeJetPNetLegacy_QCDbb;
  Float_t ProbeJetPNetLegacy_QCDothers;
  // probe Fatjet GloParT scores
  Float_t ProbeJetGloParT_QCD0HF;
  Float_t ProbeJetGloParT_QCD1HF;
  Float_t ProbeJetGloParT_QCD2HF;
  Float_t ProbeJetGloParT_Xbb;
  Float_t ProbeJetGloParT_Xcc;
  Float_t ProbeJetGloParT_Xqq;
  Float_t ProbeJetGloParT_XbbVsQCD;
  Float_t ProbeJetGloParT_massRes;
  Float_t ProbeJetGloParT_massVis;

  // FatJet 1 kinematics
  Float_t FatJet1_pt;
  Float_t FatJet1_eta;
  Float_t FatJet1_phi;
  Float_t FatJet1_Mass;
  Float_t FatJet1_MassSD;
  Float_t FatJet1_rawFactor;
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

  // FatJet 1 kinematics
  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_mass", &FatJet1_Mass);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);
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
  for (int i = 0; i < InputTree_TrgObj->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);

    // HLT Selection
    bool HLT_QCD = HLT_AK8PFJet230_SoftDropMass40;
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
    } else if (channel_lower == "jetmet" or channel_lower == "qcd") {
      if (!HLT_QCD) {
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

    // Tag and Probe
    if (channel_lower == "jetmet" or channel_lower == "qcd") {
      // Tag jet requirements
      if (FatJet1_pt < 300 || fabs(FatJet1_eta) > 2.5 ||
          FatJet1_MassSD < 80) {
        continue;
      }
      // lepton requirements
      if (lep1_Pt > 20.0) {
        continue;
      }

      // Back-to-back requirement
      if (phi_dist(FatJet1_phi, FatJet2_phi) < 2.5) {
        continue;
      }

      // No additional jets
      if (FatJet3_pt > 150) {
        continue;
      }

      bool tag_matched = checkTriggerMatching(
          NTrigger_Objects, Trigger_Object_pt, Trigger_Object_eta, Trigger_Object_phi,
          Trigger_Object_bit, FatJet1_eta, FatJet1_phi, 4, 100);
      if (!tag_matched) {
        continue;
      }

      // Assign probe jet variables (using FatJet2 for QCD)
      ProbeJet_pt = FatJet2_pt;
      ProbeJet_eta = FatJet2_eta;
      ProbeJet_phi = FatJet2_phi;
      ProbeJet_Mass = FatJet2_Mass;
      ProbeJet_MassSD = FatJet2_MassSD;
      ProbeJet_rawFactor = FatJet2_rawFactor;
      // probe Fatjet ParticleNet scores
      ProbeJetPNet_QCD = FatJet2PNet_QCD;
      ProbeJetPNet_QCD0HF = FatJet2PNet_QCD0HF;
      ProbeJetPNet_QCD1HF = FatJet2PNet_QCD1HF;
      ProbeJetPNet_QCD2HF = FatJet2PNet_QCD2HF;
      ProbeJetPNet_XbbVsQCD = FatJet2PNet_XbbVsQCD;
      ProbeJetPNet_XccVsQCD = FatJet2PNet_XccVsQCD;
      ProbeJetPNet_XggVsQCD = FatJet2PNet_XggVsQCD;
      ProbeJetPNet_XqqVsQCD = FatJet2PNet_XqqVsQCD;
      // probe Fatjet ParticleNetLegacy scores
      ProbeJetPNetLegacy_Xbb = FatJet2PNetLegacy_Xbb;
      ProbeJetPNetLegacy_Xcc = FatJet2PNetLegacy_Xcc;
      ProbeJetPNetLegacy_Xqq = FatJet2PNetLegacy_Xqq;
      ProbeJetPNetLegacy_QCD = FatJet2PNetLegacy_QCD;
      ProbeJetPNetLegacy_QCDb = FatJet2PNetLegacy_QCDb;
      ProbeJetPNetLegacy_QCDbb = FatJet2PNetLegacy_QCDbb;
      ProbeJetPNetLegacy_QCDothers = FatJet2PNetLegacy_QCDothers;
      // probe Fatjet GloParT scores
      ProbeJetGloParT_QCD0HF = FatJet2GloParT_QCD0HF;
      ProbeJetGloParT_QCD1HF = FatJet2GloParT_QCD1HF;
      ProbeJetGloParT_QCD2HF = FatJet2GloParT_QCD2HF;
      ProbeJetGloParT_Xbb = FatJet2GloParT_Xbb;
      ProbeJetGloParT_Xcc = FatJet2GloParT_Xcc;
      ProbeJetGloParT_Xqq = FatJet2GloParT_Xqq;
      ProbeJetGloParT_XbbVsQCD = FatJet2GloParT_XbbVsQCD;
      ProbeJetGloParT_massRes = FatJet2GloParT_massRes;
      ProbeJetGloParT_massVis = FatJet2GloParT_massVis;

      if (ProbeJet_pt < 160 || fabs(ProbeJet_eta) > 2.5) {
        continue;
      }

    } else {
      // Leptonic channel
      // EGamma, Muon, Lepton (EGamma + Muon)
      if (lep1_Pt < 55 || lep2_Pt > 30 || fabs(lep1_Eta) > 2.4) {
        continue;
      }

      // Back-to-back requirement
      if (phi_dist(FatJet1_phi, lep1_Phi) < 2.0) {
        continue;
      }

      // No additional jets
      if (FatJet2_pt > 180) {
        continue;
      }
      if (MET < 50) {
        continue;
      }

      if (!checkAK4JetRequirements(Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet2_Pt,
                                   Jet2_Eta, Jet2_Phi, FatJet1_eta, FatJet1_phi,
                                   lep1_Eta, lep1_Phi)) {
        continue;
      }

      // Assign probe jet variables (using FatJet2 for QCD)
      ProbeJet_pt = FatJet1_pt;
      ProbeJet_eta = FatJet1_eta;
      ProbeJet_phi = FatJet1_phi;
      ProbeJet_Mass = FatJet1_Mass;
      ProbeJet_MassSD = FatJet1_MassSD;
      ProbeJet_rawFactor = FatJet1_rawFactor;
      // probe Fatjet ParticleNet scores
      ProbeJetPNet_QCD = FatJet1PNet_QCD;
      ProbeJetPNet_QCD0HF = FatJet1PNet_QCD0HF;
      ProbeJetPNet_QCD1HF = FatJet1PNet_QCD1HF;
      ProbeJetPNet_QCD2HF = FatJet1PNet_QCD2HF;
      ProbeJetPNet_XbbVsQCD = FatJet1PNet_XbbVsQCD;
      ProbeJetPNet_XccVsQCD = FatJet1PNet_XccVsQCD;
      ProbeJetPNet_XggVsQCD = FatJet1PNet_XggVsQCD;
      ProbeJetPNet_XqqVsQCD = FatJet1PNet_XqqVsQCD;
      // probe Fatjet ParticleNetLegacy scores
      ProbeJetPNetLegacy_Xbb = FatJet1PNetLegacy_Xbb;
      ProbeJetPNetLegacy_Xcc = FatJet1PNetLegacy_Xcc;
      ProbeJetPNetLegacy_Xqq = FatJet1PNetLegacy_Xqq;
      ProbeJetPNetLegacy_QCD = FatJet1PNetLegacy_QCD;
      ProbeJetPNetLegacy_QCDb = FatJet1PNetLegacy_QCDb;
      ProbeJetPNetLegacy_QCDbb = FatJet1PNetLegacy_QCDbb;
      ProbeJetPNetLegacy_QCDothers = FatJet1PNetLegacy_QCDothers;
      // probe Fatjet GloParT scores
      ProbeJetGloParT_QCD0HF = FatJet1GloParT_QCD0HF;
      ProbeJetGloParT_QCD1HF = FatJet1GloParT_QCD1HF;
      ProbeJetGloParT_QCD2HF = FatJet1GloParT_QCD2HF;
      ProbeJetGloParT_Xbb = FatJet1GloParT_Xbb;
      ProbeJetGloParT_Xcc = FatJet1GloParT_Xcc;
      ProbeJetGloParT_Xqq = FatJet1GloParT_Xqq;
      ProbeJetGloParT_XbbVsQCD = FatJet1GloParT_XbbVsQCD;
      ProbeJetGloParT_massRes = FatJet1GloParT_massRes;
      ProbeJetGloParT_massVis = FatJet1GloParT_massVis;

      if (ProbeJet_pt < 160 || fabs(ProbeJet_eta) > 2.5) {
        continue;
      }

    } // end if

    // Probe jet matching
    bool probe_pass = false;
    bool matched_to_AK8PFJet230_SoftDropMass40 = checkTriggerMatching(
        NTrigger_Objects, Trigger_Object_pt, Trigger_Object_eta, Trigger_Object_phi,
        Trigger_Object_bit, ProbeJet_eta, ProbeJet_phi, 4, 100);

    if (year == "2022") {
      // $CMSSW_RELEASE_BASE/src/PhysicsTools/NanoAOD/python/triggerObjects_cff.py
      if (matched_to_AK8PFJet230_SoftDropMass40) {
        bool matched_to_AK8PFJet250 = checkTriggerMatching(
            NTrigger_Objects, Trigger_Object_pt, Trigger_Object_eta, Trigger_Object_phi,
            Trigger_Object_bit, ProbeJet_eta, ProbeJet_phi, 250);
        probe_pass =
            matched_to_AK8PFJet230_SoftDropMass40 & matched_to_AK8PFJet250;
      } else {
        probe_pass = false;
      }

    } else if (year == "2023") {
      probe_pass = matched_to_AK8PFJet230_SoftDropMass40;
    } else {
      throw std::invalid_argument("Invalid year: " + year);
    }

    // Fill histograms
    // FatJet 1
    _ProbeJet_all_pt->Fill(ProbeJet_pt);
    _ProbeJet_all_eta->Fill(ProbeJet_eta);
    _ProbeJet_all_phi->Fill(ProbeJet_phi);
    _ProbeJet_all_eta_phi->Fill(ProbeJet_eta, ProbeJet_phi);
    _ProbeJet_all_Mass->Fill(ProbeJet_Mass);
    _ProbeJet_all_MassSD->Fill(ProbeJet_MassSD);
    _ProbeJet_all_Mass_Pt->Fill(ProbeJet_MassSD, ProbeJet_pt);

    // ParticleNet
    _ProbeJet_all_PNet_QCD->Fill(ProbeJetPNet_QCD);
    _ProbeJet_all_PNet_QCD0HF->Fill(ProbeJetPNet_QCD0HF);
    _ProbeJet_all_PNet_QCD1HF->Fill(ProbeJetPNet_QCD1HF);
    _ProbeJet_all_PNet_QCD2HF->Fill(ProbeJetPNet_QCD2HF);
    _ProbeJet_all_PNet_XbbVsQCD->Fill(ProbeJetPNet_XbbVsQCD);
    _ProbeJet_all_PNet_XccVsQCD->Fill(ProbeJetPNet_XccVsQCD);
    _ProbeJet_all_PNet_XggVsQCD->Fill(ProbeJetPNet_XggVsQCD);
    _ProbeJet_all_PNet_XqqVsQCD->Fill(ProbeJetPNet_XqqVsQCD);
    // ParticleNetLegacy
    Float_t ProbeJetPNetLegacy_XbbVsQCD =
        ProbeJetPNetLegacy_Xbb /
        (ProbeJetPNetLegacy_Xbb + ProbeJetPNetLegacy_QCD + 1e-12);
    _ProbeJet_all_PNetLegacy_Xbb->Fill(ProbeJetPNetLegacy_Xbb);
    _ProbeJet_all_PNetLegacy_XbbVsQCD->Fill(ProbeJetPNetLegacy_XbbVsQCD);
    _ProbeJet_all_PNetLegacy_Xcc->Fill(ProbeJetPNetLegacy_Xcc);
    _ProbeJet_all_PNetLegacy_Xqq->Fill(ProbeJetPNetLegacy_Xqq);
    _ProbeJet_all_PNetLegacy_QCD->Fill(ProbeJetPNetLegacy_QCD);
    _ProbeJet_all_PNetLegacy_QCDb->Fill(ProbeJetPNetLegacy_QCDb);
    _ProbeJet_all_PNetLegacy_QCDbb->Fill(ProbeJetPNetLegacy_QCDbb);
    _ProbeJet_all_PNetLegacy_QCDothers->Fill(ProbeJetPNetLegacy_QCDothers);
    // GloParT
    _ProbeJet_all_GloParT_QCD0HF->Fill(ProbeJetGloParT_QCD0HF);
    _ProbeJet_all_GloParT_QCD1HF->Fill(ProbeJetGloParT_QCD1HF);
    _ProbeJet_all_GloParT_QCD2HF->Fill(ProbeJetGloParT_QCD2HF);
    _ProbeJet_all_GloParT_Xbb->Fill(ProbeJetGloParT_Xbb);
    _ProbeJet_all_GloParT_Xcc->Fill(ProbeJetGloParT_Xcc);
    _ProbeJet_all_GloParT_Xqq->Fill(ProbeJetGloParT_Xqq);
    _ProbeJet_all_GloParT_XbbVsQCD->Fill(ProbeJetGloParT_XbbVsQCD);

    if (probe_pass) {
      _ProbeJet_pass_pt->Fill(ProbeJet_pt);
      _ProbeJet_pass_eta->Fill(ProbeJet_eta);
      _ProbeJet_pass_phi->Fill(ProbeJet_phi);
      _ProbeJet_pass_eta_phi->Fill(ProbeJet_eta, ProbeJet_phi);
      _ProbeJet_pass_Mass->Fill(ProbeJet_Mass);
      _ProbeJet_pass_MassSD->Fill(ProbeJet_MassSD);
      _ProbeJet_pass_Mass_Pt->Fill(ProbeJet_MassSD, ProbeJet_pt);

      // ParticleNet
      _ProbeJet_pass_PNet_QCD->Fill(ProbeJetPNet_QCD);
      _ProbeJet_pass_PNet_QCD0HF->Fill(ProbeJetPNet_QCD0HF);
      _ProbeJet_pass_PNet_QCD1HF->Fill(ProbeJetPNet_QCD1HF);
      _ProbeJet_pass_PNet_QCD2HF->Fill(ProbeJetPNet_QCD2HF);
      _ProbeJet_pass_PNet_XbbVsQCD->Fill(ProbeJetPNet_XbbVsQCD);
      _ProbeJet_pass_PNet_XccVsQCD->Fill(ProbeJetPNet_XccVsQCD);
      _ProbeJet_pass_PNet_XggVsQCD->Fill(ProbeJetPNet_XggVsQCD);
      _ProbeJet_pass_PNet_XqqVsQCD->Fill(ProbeJetPNet_XqqVsQCD);
      // ParticleNetLegacy
      _ProbeJet_pass_PNetLegacy_Xbb->Fill(ProbeJetPNetLegacy_Xbb);
      _ProbeJet_pass_PNetLegacy_XbbVsQCD->Fill(ProbeJetPNetLegacy_XbbVsQCD);
      _ProbeJet_pass_PNetLegacy_Xcc->Fill(ProbeJetPNetLegacy_Xcc);
      _ProbeJet_pass_PNetLegacy_Xqq->Fill(ProbeJetPNetLegacy_Xqq);
      _ProbeJet_pass_PNetLegacy_QCD->Fill(ProbeJetPNetLegacy_QCD);
      _ProbeJet_pass_PNetLegacy_QCDb->Fill(ProbeJetPNetLegacy_QCDb);
      _ProbeJet_pass_PNetLegacy_QCDbb->Fill(ProbeJetPNetLegacy_QCDbb);
      _ProbeJet_pass_PNetLegacy_QCDothers->Fill(ProbeJetPNetLegacy_QCDothers);
      // GloParT
      _ProbeJet_pass_GloParT_QCD0HF->Fill(ProbeJetGloParT_QCD0HF);
      _ProbeJet_pass_GloParT_QCD1HF->Fill(ProbeJetGloParT_QCD1HF);
      _ProbeJet_pass_GloParT_QCD2HF->Fill(ProbeJetGloParT_QCD2HF);
      _ProbeJet_pass_GloParT_Xbb->Fill(ProbeJetGloParT_Xbb);
      _ProbeJet_pass_GloParT_Xcc->Fill(ProbeJetGloParT_Xcc);
      _ProbeJet_pass_GloParT_Xqq->Fill(ProbeJetGloParT_Xqq);
      _ProbeJet_pass_GloParT_XbbVsQCD->Fill(ProbeJetGloParT_XbbVsQCD);
    }

  } // end event loop

  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
