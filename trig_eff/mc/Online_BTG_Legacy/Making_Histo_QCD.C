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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#define ARR_SIZE 10000

struct ParamDict {
  double Lumi;
  double XSec_QCD_HT_100to200;
  double XSec_QCD_HT_200to400;
  double XSec_QCD_HT_400to600;
  double XSec_QCD_HT_600to800;
  double XSec_QCD_HT_800to1000;
  double XSec_QCD_HT_1000to1200;
  double XSec_QCD_HT_1200to1500;
  double XSec_QCD_HT_1500to2000;
  double XSec_QCD_HT_2000toInf;
};

void loadParamDict(ParamDict *param_dict, const std::string &param_path) {
  std::ifstream file(param_path);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open parameters file: " + param_path);
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '/')
      continue;

    // Look for double declarations
    if (line.find("double") != std::string::npos) {
      std::istringstream iss(line);
      std::string type, name, equals;
      double value;

      // Parse line of format: double XSec_QCD_HT_100to200    =  25220000.00;
      if (!(iss >> type >> name >> equals >> value)) {
        continue;
      }

      // Remove semicolon if present
      if (name.back() == ';')
        name = name.substr(0, name.size() - 1);

      if (name == "Lumi") {
        param_dict->Lumi = value;
      } else if (name == "XSec_QCD_HT_100to200") {
        param_dict->XSec_QCD_HT_100to200 = value;
      } else if (name == "XSec_QCD_HT_200to400") {
        param_dict->XSec_QCD_HT_200to400 = value;
      } else if (name == "XSec_QCD_HT_400to600") {
        param_dict->XSec_QCD_HT_400to600 = value;
      } else if (name == "XSec_QCD_HT_600to800") {
        param_dict->XSec_QCD_HT_600to800 = value;
      } else if (name == "XSec_QCD_HT_800to1000") {
        param_dict->XSec_QCD_HT_800to1000 = value;
      } else if (name == "XSec_QCD_HT_1000to1200") {
        param_dict->XSec_QCD_HT_1000to1200 = value;
      } else if (name == "XSec_QCD_HT_1200to1500") {
        param_dict->XSec_QCD_HT_1200to1500 = value;
      } else if (name == "XSec_QCD_HT_1500to2000") {
        param_dict->XSec_QCD_HT_1500to2000 = value;
      } else if (name == "XSec_QCD_HT_2000toInf") {
        param_dict->XSec_QCD_HT_2000toInf = value;
      }
    }
  }
}

double getXSec(ParamDict *param_dict, std::string ht_bin) {
  if (ht_bin == "100to200") {
    return param_dict->XSec_QCD_HT_100to200;
  } else if (ht_bin == "200to400") {
    return param_dict->XSec_QCD_HT_200to400;
  } else if (ht_bin == "400to600") {
    return param_dict->XSec_QCD_HT_400to600;
  } else if (ht_bin == "600to800") {
    return param_dict->XSec_QCD_HT_600to800;
  } else if (ht_bin == "800to1000") {
    return param_dict->XSec_QCD_HT_800to1000;
  } else if (ht_bin == "1000to1200") {
    return param_dict->XSec_QCD_HT_1000to1200;
  } else if (ht_bin == "1200to1500") {
    return param_dict->XSec_QCD_HT_1200to1500;
  } else if (ht_bin == "1500to2000") {
    return param_dict->XSec_QCD_HT_1500to2000;
  } else if (ht_bin == "2000toInf") {
    return param_dict->XSec_QCD_HT_2000toInf;
  } else {
    throw std::runtime_error("Invalid HT bin: " + ht_bin);
  }
}

std::vector<double> loadPUReweighting(const std::string &pu_file) {
  std::vector<double> values;
  std::ifstream file(pu_file);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open PU reweighting file: " + pu_file);
  }

  double value;
  while (file >> value) {
    values.push_back(value);
  }

  if (values.size() != 100) { // Typically expecting 100 bins
    throw std::runtime_error("Expected 100 PU reweighting values, got " +
                             std::to_string(values.size()));
  }

  return values;
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

void Making_Histo_QCD(
    const std::string &ht_bin,       // QCD HT bin, e.g. "100to200"
    const std::string &sample_path,  // path to the root file
    const std::string &output_path,  // path to the output root file
    const std::string &pu_path,      // path to the pileup reweighting file
    const std::string &sf_path,      // path to the PT-Mass-SFs root file
    const std::string &param_path,   // path to the parameters file
    const std::string &jec_path_ak4, // path to the AK4 JEC txt file
    const std::string &jec_path_ak8  // path to the AK8 JEC txt file
) {
  gSystem->Load("libFWCoreFWLite.so");

  // PT-Mass-SFs
  TFile *f_PT_Mass_SF = new TFile(sf_path.c_str());
  TH2D *_Eff_Data = (TH2D *)f_PT_Mass_SF->Get("Eff_Data_ETA0");
  TH2D *_Eff_MC = (TH2D *)f_PT_Mass_SF->Get("Eff_MC_ETA0");

  // parse ParamDict
  ParamDict param_dict;
  loadParamDict(&param_dict, param_path);
  double xsec = getXSec(&param_dict, ht_bin);

  // pu weight
  // std::vector<double> PU_Rew = loadPUReweighting(pu_path);
  std::vector<double> PU_Rew_vec = loadPUReweighting(pu_path);
  double PU_Rew[100];
  std::copy(PU_Rew_vec.begin(), PU_Rew_vec.end(), PU_Rew);

  // JEC
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(JetCorrectorParameters(jec_path_ak4.c_str()));
  FactorizedJetCorrector *corrector = new FactorizedJetCorrector(vPar);

  vector<JetCorrectorParameters> vParAK8;
  vParAK8.push_back(JetCorrectorParameters(jec_path_ak8.c_str()));
  FactorizedJetCorrector *corrector_AK8 = new FactorizedJetCorrector(vParAK8);

  /* (should be added when available)
  // JER
  std::string resptstr =
  "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_PtResolution_AK4PFPuppi.txt";
  std::string resptstr_sf =
  "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_SF_AK4PFPuppi.txt";
  JME::JetResolution resolution_pt = JME::JetResolution(resptstr.c_str());
  JME::JetResolutionScaleFactor resolution_pt_sf =
  JME::JetResolutionScaleFactor(resptstr_sf.c_str());

  std::string resptstr_AK8 =
  "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_PtResolution_AK8PFPuppi.txt";
  std::string resptstr_sf_AK8 =
  "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_SF_AK8PFPuppi.txt";
  JME::JetResolution resolution_pt_AK8 =
  JME::JetResolution(resptstr_AK8.c_str()); JME::JetResolutionScaleFactor
  resolution_pt_sf_AK8 = JME::JetResolutionScaleFactor(resptstr_sf_AK8.c_str());
  */

  TFile *f = new TFile(output_path.c_str(), "RECREATE");

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
  TH1D *_FatJet1PNetMD_Xbb_UDSG =
      new TH1D("FatJet1PNetMD_Xbb_UDSG", "FatJet1PNetMD_Xbb_UDSG", 100, 0, 1.0);
  TH1D *_FatJet1PNetMD_Xbb_C =
      new TH1D("FatJet1PNetMD_Xbb_C", "FatJet1PNetMD_Xbb_C", 100, 0, 1.0);
  TH1D *_FatJet1PNetMD_Xbb_B =
      new TH1D("FatJet1PNetMD_Xbb_B", "FatJet1PNetMD_Xbb_B", 100, 0, 1.0);
  TH1D *_FatJet1PNetMD_Xbb_B_1 =
      new TH1D("FatJet1PNetMD_Xbb_B_1", "FatJet1PNetMD_Xbb_B_1", 100, 0, 1.0);
  TH1D *_FatJet1PNetMD_Xbb_B_2 =
      new TH1D("FatJet1PNetMD_Xbb_B_2", "FatJet1PNetMD_Xbb_B_2", 100, 0, 1.0);
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
  TH1D *_FatJet2PNetMD_Xbb_UDSG =
      new TH1D("FatJet2PNetMD_Xbb_UDSG", "FatJet2PNetMD_Xbb_UDSG", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_C =
      new TH1D("FatJet2PNetMD_Xbb_C", "FatJet2PNetMD_Xbb_C", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_B =
      new TH1D("FatJet2PNetMD_Xbb_B", "FatJet2PNetMD_Xbb_B", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_B_1 =
      new TH1D("FatJet2PNetMD_Xbb_B_1", "FatJet2PNetMD_Xbb_B_1", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_B_2 =
      new TH1D("FatJet2PNetMD_Xbb_B_2", "FatJet2PNetMD_Xbb_B_2", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_Legacy = new TH1D(
      "FatJet2PNetMD_Xbb_Legacy", "FatJet2PNetMD_Xbb_Legacy", 100, 0, 1.0);
  TH1D *_FatJet2PNetMD_Xbb_Legacy_AN =
      new TH1D("FatJet2PNetMD_Xbb_Legacy_AN", "FatJet2PNetMD_Xbb_Legacy_AN",
               100, 0, 1.0);

  TFile *f1 = new TFile(sample_path.c_str());

  TH1F *NEvents = (TH1F *)f1->Get("NEvents");
  double SumGenWeights = NEvents->GetBinContent(1);
  TTree *InputTree = (TTree *)f1->Get("tree");

  Float_t weight;
  UInt_t run;
  UInt_t lumi;
  Float_t npu;
  Int_t isVBFtag;
  Float_t rho;

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

  Int_t nGenJet;
  Float_t GenJet_eta[ARR_SIZE];
  Float_t GenJet_phi[ARR_SIZE];
  Float_t GenJet_pt[ARR_SIZE];
  Int_t nGenJetAK8;
  Float_t GenJetAK8_eta[ARR_SIZE];
  Float_t GenJetAK8_phi[ARR_SIZE];
  Float_t GenJetAK8_pt[ARR_SIZE];

  InputTree->SetBranchAddress("weight", &weight);
  InputTree->SetBranchAddress("run", &run);
  InputTree->SetBranchAddress("lumi", &lumi);
  InputTree->SetBranchAddress("npu", &npu);
  InputTree->SetBranchAddress("rho", &rho);
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
  InputTree->SetBranchAddress("fatJet1_particleNet_XbbVsQCD",
                              &FatJet1PNetMD_Xbb);
  InputTree->SetBranchAddress("fatJet1_hadronFlavour", &FatJet1_hadronFlavour);
  InputTree->SetBranchAddress("fatJet1_nBHadrons", &FatJet1_nBHadrons);
  InputTree->SetBranchAddress("fatJet1_nCHadrons", &FatJet1_nCHadrons);
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
  InputTree->SetBranchAddress("fatJet2_hadronFlavour", &FatJet2_hadronFlavour);
  InputTree->SetBranchAddress("fatJet2_nBHadrons", &FatJet2_nBHadrons);
  InputTree->SetBranchAddress("fatJet2_nCHadrons", &FatJet2_nCHadrons);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xbb",
                              &FatJet2PNetMD_Xbb_Legacy);
  InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCD",
                              &FatJet2PNetMD_QCD_Legacy);

  InputTree->SetBranchAddress("fatJet3_pt", &FatJet3_pt);
  InputTree->SetBranchAddress("fatJet3_rawFactor", &FatJet3_rawFactor);

  InputTree->SetBranchAddress("nGenJet", &nGenJet);
  InputTree->SetBranchAddress("GenJet_eta", GenJet_eta);
  InputTree->SetBranchAddress("GenJet_phi", GenJet_phi);
  InputTree->SetBranchAddress("GenJet_pt", GenJet_pt);
  InputTree->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
  InputTree->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta);
  InputTree->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi);
  InputTree->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt);

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

  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);

    // HLT Selection
    if (HLT_AK8PFJet230_SoftDropMass40 == 0) {
      continue;
    }

    // FatJets correction and selection
    if (FatJet1_pt > 0) {
      double Raw_FatJet1_pt = FatJet1_pt * (1.0 - FatJet1_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_pt);
      corrector_AK8->setJetEta(FatJet1_eta);
      corrector_AK8->setJetPhi(FatJet1_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet1_pt = Raw_FatJet1_pt * corr;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * corr;
    }
    if (FatJet2_pt > 0) {
      double Raw_FatJet2_pt = FatJet2_pt * (1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet2_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      corrector_AK8->setJetPhi(FatJet2_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt * corr;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * corr;
    }

    /*
    // Jet Smearing
    double res_pt_1;
    double res_pt_sf_1;
    JME::JetParameters JerPARAM_1 = {{JME::Binning::JetPt, FatJet1_pt},
                                     {JME::Binning::JetEta, FatJet1_eta},
                                     {JME::Binning::Rho, rho}};
    JME::JetParameters JerSFPARAM_1;
    JerSFPARAM_1.set(JME::Binning::JetPt, FatJet1_pt);
    JerSFPARAM_1.set(JME::Binning::JetEta, FatJet1_eta);
    JerSFPARAM_1.set(JME::Binning::Rho, rho);
    res_pt_1 = resolution_pt_AK8.getResolution(JerPARAM_1);
    res_pt_sf_1 = resolution_pt_sf_AK8.getScaleFactor(JerSFPARAM_1);

    double res_pt_2;
    double res_pt_sf_2;
    JME::JetParameters JerPARAM_2 = {{JME::Binning::JetPt, FatJet2_pt},
                                     {JME::Binning::JetEta, FatJet2_eta},
                                     {JME::Binning::Rho, rho}};
    JME::JetParameters JerSFPARAM_2;
    JerSFPARAM_2.set(JME::Binning::JetPt, FatJet2_pt);
    JerSFPARAM_2.set(JME::Binning::JetEta, FatJet2_eta);
    JerSFPARAM_2.set(JME::Binning::Rho, rho);
    res_pt_2 = resolution_pt_AK8.getResolution(JerPARAM_2);
    res_pt_sf_2 = resolution_pt_sf_AK8.getScaleFactor(JerSFPARAM_2);

    double SmearFactor_1 = 1;
    bool GenJetMatched_1 = false;
    double SmearFactor_2 = 1;
    bool GenJetMatched_2 = false;

    for (int nGJAK8 = 0; nGJAK8 < nGenJetAK8; nGJAK8++) {
      if (!GenJetMatched_1 &&
          sqrt(pow(FatJet1_eta - GenJetAK8_eta[nGJAK8], 2) +
               pow(phi_dist(FatJet1_phi, GenJetAK8_phi[nGJAK8]), 2)) < 0.2 &&
          (fabs(FatJet1_pt - GenJetAK8_pt[nGJAK8]) / FatJet1_pt <
           3 * res_pt_1)) {
        SmearFactor_1 = 1.0 + (res_pt_sf_1 - 1.0) *
                                  (FatJet1_pt - GenJetAK8_pt[nGJAK8]) /
                                  FatJet1_pt;
        GenJetMatched_1 = true;
      }
      if (!GenJetMatched_2 &&
          sqrt(pow(FatJet2_eta - GenJetAK8_eta[nGJAK8], 2) +
               pow(phi_dist(FatJet2_phi, GenJetAK8_phi[nGJAK8]), 2)) < 0.2 &&
          (fabs(FatJet2_pt - GenJetAK8_pt[nGJAK8]) / FatJet2_pt <
           3 * res_pt_2)) {
        SmearFactor_2 = 1.0 + (res_pt_sf_2 - 1.0) *
                                  (FatJet2_pt - GenJetAK8_pt[nGJAK8]) /
                                  FatJet2_pt;
        GenJetMatched_2 = true;
      }
    }

    // if(!GenJetMatched && res_pt_sf[nJ] > 1.0)
    //   {
    //    double sigma = res_pt[nJ] * sqrt(res_pt_sf[nJ]*res_pt_sf[nJ] - 1);
    //    normal_distribution<> d(0, sigma);
    //    SmearFactor = 1.0 + d(m_random_generator);
    //   }

    // Smear
    FatJet1_pt = FatJet1_pt * SmearFactor_1;
    FatJet1_MassSD = FatJet1_MassSD * SmearFactor_1;
    FatJet2_pt = FatJet2_pt * SmearFactor_2;
    FatJet2_MassSD = FatJet2_MassSD * SmearFactor_2;
    */

    // FatJets selection
    if (FatJet1_pt < 250 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50) {
      continue;
    }

    // VBFTag veto
    //   if(isVBFtag) continue;

    // Lepton Selection or Veto
    if (lep1_Pt > 20.0) {
      continue;
    }

    // Trigger Objects and Matchings
    // Matching 1st
    bool matched_TRG_1 = false;

    bool matched_to_AK8PFJet230_SoftDropMass40 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++) {
      if ((Trigger_Object_bit[itrg] & 4) == 4) {
        double dR = get_dR(FatJet1_eta, FatJet1_phi, Trigger_Object_eta[itrg],
                           Trigger_Object_phi[itrg]);
        if (dR < 0.4 && Trigger_Object_pt[itrg] > 100) {
          matched_to_AK8PFJet230_SoftDropMass40 = true;
          break;
        }
      }
    }

    if (matched_to_AK8PFJet230_SoftDropMass40)
      matched_TRG_1 = true;
    else
      continue;

    bool matched_to_AK8PFJet230_SoftDropMass40_fJ2 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++) {
      if ((Trigger_Object_bit[itrg] & 4) == 4) {
        double dR = get_dR(FatJet2_eta, FatJet2_phi, Trigger_Object_eta[itrg],
                           Trigger_Object_phi[itrg]);
        if (dR < 0.4 && Trigger_Object_pt[itrg] > 100) {
          matched_to_AK8PFJet230_SoftDropMass40_fJ2 = true;
          break;
        }
      }
    }

    if (matched_to_AK8PFJet230_SoftDropMass40_fJ2)
      continue;
    if (FatJet3_pt > 200)
      continue;

    // Matching 2nd
    bool matched_TRG_2 = false;
    if (HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06)
      matched_TRG_2 = true;

    // weight
    weight = (weight / SumGenWeights) * xsec * param_dict.Lumi;
    double PU_weight = PU_Rew[(int)npu];
    if (PU_weight < 20.0) {
      weight = weight * PU_weight;
    }

    // Add Trigger PT_Mass and PNet Scale Factors
    double Eff_Data_1 = 0;
    double Eff_MC_1 = 0;

    if (matched_TRG_1) {
      Int_t bin_PT_1 = _Eff_Data->GetXaxis()->FindBin(FatJet1_pt);
      Int_t bin_Mass_1 = _Eff_Data->GetYaxis()->FindBin(FatJet1_MassSD);
      Eff_Data_1 = 1.0;
      if (_Eff_Data->GetBinContent(bin_PT_1, bin_Mass_1) > 0)
        Eff_Data_1 = _Eff_Data->GetBinContent(bin_PT_1, bin_Mass_1);
      Eff_MC_1 = 1.0;
      if (_Eff_MC->GetBinContent(bin_PT_1, bin_Mass_1) > 0)
        Eff_MC_1 = _Eff_MC->GetBinContent(bin_PT_1, bin_Mass_1);
    }

    double Tot_Data = 1 - (1 - Eff_Data_1);
    double Tot_MC = 1 - (1 - Eff_MC_1);
    double PT_Mass_BTG_SF = Tot_Data / Tot_MC;

    if (PT_Mass_BTG_SF > 0) {
      weight = weight * PT_Mass_BTG_SF;
    }

    // Fill histograms
    _FatJet1_pt->Fill(FatJet1_pt, weight);
    _FatJet1_eta->Fill(FatJet1_eta, weight);
    _FatJet1_phi->Fill(FatJet1_phi, weight);
    _FatJet1_eta_phi->Fill(FatJet1_eta, FatJet1_phi, weight);
    _FatJet1_Mass->Fill(FatJet1_Mass, weight);
    _FatJet1_MassSD->Fill(FatJet1_MassSD, weight);
    _FatJet1_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD, weight);
    _FatJet1PNetMD_Xbb->Fill(FatJet1PNetMD_Xbb, weight);
    _FatJet1PNetMD_Xbb_Legacy->Fill(FatJet1PNetMD_Xbb_Legacy, weight);
    double FatJet1PNetMD_Xbb_Legacy_AN =
        FatJet1PNetMD_Xbb_Legacy /
        (FatJet1PNetMD_Xbb_Legacy + FatJet1PNetMD_QCD_Legacy);
    _FatJet1PNetMD_Xbb_Legacy_AN->Fill(FatJet1PNetMD_Xbb_Legacy_AN, weight);

    if (FatJet1_hadronFlavour == 0)
      _FatJet1PNetMD_Xbb_UDSG->Fill(FatJet1PNetMD_Xbb, weight);
    if (FatJet1_hadronFlavour == 4)
      _FatJet1PNetMD_Xbb_C->Fill(FatJet1PNetMD_Xbb, weight);
    if (FatJet1_hadronFlavour == 5) {
      _FatJet1PNetMD_Xbb_B->Fill(FatJet1PNetMD_Xbb, weight);
      if (FatJet1_nBHadrons <= 1)
        _FatJet1PNetMD_Xbb_B_1->Fill(FatJet1PNetMD_Xbb, weight);
      if (FatJet1_nBHadrons >= 2)
        _FatJet1PNetMD_Xbb_B_2->Fill(FatJet1PNetMD_Xbb, weight);
    }

    if (matched_TRG_2) {
      _FatJet2_pt->Fill(FatJet1_pt, weight);
      _FatJet2_eta->Fill(FatJet1_eta, weight);
      _FatJet2_phi->Fill(FatJet1_phi, weight);
      _FatJet2_eta_phi->Fill(FatJet1_eta, FatJet1_phi, weight);
      _FatJet2_Mass->Fill(FatJet1_Mass, weight);
      _FatJet2_MassSD->Fill(FatJet1_MassSD, weight);
      _FatJet2_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD, weight);
      _FatJet2PNetMD_Xbb->Fill(FatJet1PNetMD_Xbb, weight);
      _FatJet2PNetMD_Xbb_Legacy->Fill(FatJet1PNetMD_Xbb_Legacy, weight);
      _FatJet2PNetMD_Xbb_Legacy_AN->Fill(FatJet1PNetMD_Xbb_Legacy_AN, weight);

      if (FatJet1_hadronFlavour == 0)
        _FatJet2PNetMD_Xbb_UDSG->Fill(FatJet1PNetMD_Xbb, weight);
      if (FatJet1_hadronFlavour == 4)
        _FatJet2PNetMD_Xbb_C->Fill(FatJet1PNetMD_Xbb, weight);
      if (FatJet1_hadronFlavour == 5) {
        _FatJet2PNetMD_Xbb_B->Fill(FatJet1PNetMD_Xbb, weight);
        if (FatJet1_nBHadrons <= 1)
          _FatJet2PNetMD_Xbb_B_1->Fill(FatJet1PNetMD_Xbb, weight);
        if (FatJet1_nBHadrons >= 2)
          _FatJet2PNetMD_Xbb_B_2->Fill(FatJet1PNetMD_Xbb, weight);
      }
    }

  } // end event loop

  f->Write();

  std::cout << "Done with QCD HT " << ht_bin << std::endl;
}