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
#include <bits/stdc++.h>
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

std::string to_lower(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  return str;
}

struct ParamDict {
  double Lumi;
  double XSec_TTto4Q;
  double XSec_TTtoLNu2Q;
  double XSec_TTto2L2Nu;
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
    // if (line.find("double") != std::string::npos) {
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
    } else if (name == "XSec_TTtoLNu2Q") {
      param_dict->XSec_TTtoLNu2Q = value;
    } else if (name == "XSec_TTto4Q") {
      param_dict->XSec_TTto4Q = value;
    } else if (name == "XSec_TTto2L2Nu") {
      param_dict->XSec_TTto2L2Nu = value;
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
    } else {
      // warn if unknown parameter
      std::cerr << "Unknown parameter: " << name << std::endl;
    }
  }
}

double getXSec(ParamDict *param_dict, std::string data_type) {
  if (data_type == "TTtoLNu2Q") {
    return param_dict->XSec_TTtoLNu2Q;
  } else if (data_type == "TTto4Q") {
    return param_dict->XSec_TTto4Q;
  } else if (data_type == "TTto2L2Nu") {
    return param_dict->XSec_TTto2L2Nu;
  } else if (data_type == "QCD_HT_100to200") {
    return param_dict->XSec_QCD_HT_100to200;
  } else if (data_type == "QCD_HT_200to400") {
    return param_dict->XSec_QCD_HT_200to400;
  } else if (data_type == "QCD_HT_400to600") {
    return param_dict->XSec_QCD_HT_400to600;
  } else if (data_type == "QCD_HT_600to800") {
    return param_dict->XSec_QCD_HT_600to800;
  } else if (data_type == "QCD_HT_800to1000") {
    return param_dict->XSec_QCD_HT_800to1000;
  } else if (data_type == "QCD_HT_1000to1200") {
    return param_dict->XSec_QCD_HT_1000to1200;
  } else if (data_type == "QCD_HT_1200to1500") {
    return param_dict->XSec_QCD_HT_1200to1500;
  } else if (data_type == "QCD_HT_1500to2000") {
    return param_dict->XSec_QCD_HT_1500to2000;
  } else if (data_type == "QCD_HT_2000toInf") {
    return param_dict->XSec_QCD_HT_2000toInf;
  } else {
    throw std::invalid_argument("Invalid data_type: " + data_type);
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

void histo_mc(
    const std::string &data_type,    // Example: TTtoLNu2Q, QCD_HT100to200
    const std::string &channel,      // Muon, EGamma, leptonic
    const std::string &sample_path,  // path to the root file
    const std::string &output_path,  // path to the output root file
    const std::string &pu_path,      // path to the pileup reweighting file
    const std::string &param_path,   // path to the parameters file
    const std::string &jec_path_ak4, // path to th  e AK4 JEC txt file
    const std::string &jec_path_ak8  // path to the AK8 JEC txt file
) {
  gSystem->Load("libFWCoreFWLite.so");

  std::string channel_lower = to_lower(channel);
  std::cout << "Data Type: " << data_type << std::endl;
  std::cout << "Channel: " << channel << std::endl;
  std::cout << "Sample Path: " << sample_path << std::endl;
  std::cout << "Output Path: " << output_path << std::endl;
  std::cout << "PU Path: " << pu_path << std::endl;
  std::cout << "Param Path: " << param_path << std::endl;
  std::cout << "AK4 JEC Path: " << jec_path_ak4 << std::endl;
  std::cout << "AK8 JEC Path: " << jec_path_ak8 << std::endl;

  // parse ParamDict
  ParamDict param_dict;
  loadParamDict(&param_dict, param_path);
  double xsec = getXSec(&param_dict, data_type);

  // pu weight
  std::vector<double> PU_Rew = loadPUReweighting(pu_path);

  // JEC
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(JetCorrectorParameters(jec_path_ak4.c_str()));
  FactorizedJetCorrector *corrector = new FactorizedJetCorrector(vPar);

  vector<JetCorrectorParameters> vParAK8;
  vParAK8.push_back(JetCorrectorParameters(jec_path_ak8.c_str()));
  FactorizedJetCorrector *corrector_AK8 = new FactorizedJetCorrector(vParAK8);

  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  // Float_t bins_pt[9] = {300, 350, 400, 450, 500, 600, 700, 850, 1000};
  // int num_pt_bins = 8;

  // Float_t bins_m[8] = {60, 90, 120, 150, 180, 210, 240, 300};
  // int num_m_bins = 7;
  // Float_t bins_pt[46] = {0.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,
  //                        80.0,  90.0,  100.0, 110.0, 120.0, 130.0, 140.0, 150.0,
  //                        160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0,
  //                        240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 320.0,
  //                        340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0,
  //                        500.0, 550.0, 600.0, 700.0, 800.0, 1000.0};
  // int num_pt_bins = 45;

  // Float_t bins_m[16] = {0.0,  5.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,
  //                       80.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 350.0};
  // int num_m_bins = 15;

  Float_t bins_pt[9] = {250, 275, 300, 350, 400, 450, 500, 600, 1000};
  int num_pt_bins = 8;

  Float_t bins_m[10] = {50, 60, 80, 100, 120, 150, 200, 250, 300, 350};
  int num_m_bins = 9;

  // tag FatJet 1 kinematics
  TH1D *_FatJet1_tag_pt =
      new TH1D("FatJet1_tag_pt", "FatJet1_tag_pt", 200, 0, 1000);
  TH1D *_FatJet1_tag_eta =
      new TH1D("FatJet1_tag_eta", "FatJet1_tag_eta", 100, -5, 5);
  TH1D *_FatJet1_tag_phi =
      new TH1D("FatJet1_tag_phi", "FatJet1_tag_phi", 100, -5, 5);
  TH2D *_FatJet1_tag_eta_phi = new TH2D(
      "FatJet1_tag_eta_phi", "FatJet1_tag_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet1_tag_Mass =
      new TH1D("FatJet1_tag_Mass", "FatJet1_tag_Mass", 500, 0, 500);
  TH1D *_FatJet1_tag_MassSD =
      new TH1D("FatJet1_tag_MassSD", "FatJet1_tag_MassSD", 500, 0, 500);
  TH2D *_FatJet1_tag_Mass_Pt =
      new TH2D("FatJet1_tag_Mass_Pt", "FatJet1_tag_Mass_Pt", num_m_bins, bins_m,
               num_pt_bins, bins_pt);
  // tag FatJet 1 ParticleNet scores
  TH1D *_FatJet1_tag_PNet_QCD =
      new TH1D("FatJet1_tag_PNet_QCD", "FatJet1_tag_PNet_QCD", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNet_QCD0HF = new TH1D(
      "FatJet1_tag_PNet_QCD0HF", "FatJet1_tag_PNet_QCD0HF", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNet_QCD1HF = new TH1D(
      "FatJet1_tag_PNet_QCD1HF", "FatJet1_tag_PNet_QCD1HF", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNet_QCD2HF = new TH1D(
      "FatJet1_tag_PNet_QCD2HF", "FatJet1_tag_PNet_QCD2HF", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNet_XbbVsQCD = new TH1D(
      "FatJet1_tag_PNet_XbbVsQCD", "FatJet1_tag_PNet_XbbVsQCD", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNet_XccVsQCD = new TH1D(
      "FatJet1_tag_PNet_XccVsQCD", "FatJet1_tag_PNet_XccVsQCD", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNet_XggVsQCD = new TH1D(
      "FatJet1_tag_PNet_XggVsQCD", "FatJet1_tag_PNet_XggVsQCD", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNet_XqqVsQCD = new TH1D(
      "FatJet1_tag_PNet_XqqVsQCD", "FatJet1_tag_PNet_XqqVsQCD", 100, 0, 1.0);
  // tag FatJet 1 ParticleNetLegacy scores
  TH1D *_FatJet1_tag_PNetLegacy_Xbb = new TH1D(
      "FatJet1_tag_PNetLegacy_Xbb", "FatJet1_tag_PNetLegacy_Xbb", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNetLegacy_XbbVsQCD = new TH1D(
      "FatJet1_tag_PNetLegacy_XbbVsQCD", "FatJet1_tag_PNetLegacy_XbbVsQCD", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNetLegacy_Xcc = new TH1D(
      "FatJet1_tag_PNetLegacy_Xcc", "FatJet1_tag_PNetLegacy_Xcc", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNetLegacy_Xqq = new TH1D(
      "FatJet1_tag_PNetLegacy_Xqq", "FatJet1_tag_PNetLegacy_Xqq", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNetLegacy_QCD = new TH1D(
      "FatJet1_tag_PNetLegacy_QCD", "FatJet1_tag_PNetLegacy_QCD", 100, 0, 1.0);
  TH1D *_FatJet1_tag_PNetLegacy_QCDb =
      new TH1D("FatJet1_tag_PNetLegacy_QCDb", "FatJet1_tag_PNetLegacy_QCDb",
               100, 0, 1.0);
  TH1D *_FatJet1_tag_PNetLegacy_QCDbb =
      new TH1D("FatJet1_tag_PNetLegacy_QCDbb", "FatJet1_tag_PNetLegacy_QCDbb",
               100, 0, 1.0);
  TH1D *_FatJet1_tag_PNetLegacy_QCDothers =
      new TH1D("FatJet1_tag_PNetLegacy_QCDothers",
               "FatJet1_tag_PNetLegacy_QCDothers", 100, 0, 1.0);
  // tag FatJet 1 GloParT scores
  TH1D *_FatJet1_tag_GloParT_QCD0HF = new TH1D(
      "FatJet1_tag_GloParT_QCD0HF", "FatJet1_tag_GloParT_QCD0HF", 100, 0, 1.0);
  TH1D *_FatJet1_tag_GloParT_QCD1HF = new TH1D(
      "FatJet1_tag_GloParT_QCD1HF", "FatJet1_tag_GloParT_QCD1HF", 100, 0, 1.0);
  TH1D *_FatJet1_tag_GloParT_QCD2HF = new TH1D(
      "FatJet1_tag_GloParT_QCD2HF", "FatJet1_tag_GloParT_QCD2HF", 100, 0, 1.0);
  TH1D *_FatJet1_tag_GloParT_Xbb = new TH1D(
      "FatJet1_tag_GloParT_Xbb", "FatJet1_tag_GloParT_Xbb", 100, 0, 1.0);
  TH1D *_FatJet1_tag_GloParT_Xcc = new TH1D(
      "FatJet1_tag_GloParT_Xcc", "FatJet1_tag_GloParT_Xcc", 100, 0, 1.0);
  TH1D *_FatJet1_tag_GloParT_Xqq = new TH1D(
      "FatJet1_tag_GloParT_Xqq", "FatJet1_tag_GloParT_Xqq", 100, 0, 1.0);
  TH1D *_FatJet1_tag_GloParT_XbbVsQCD =
      new TH1D("FatJet1_tag_GloParT_XbbVsQCD", "FatJet1_tag_GloParT_XbbVsQCD",
               100, 0, 1.0);

  // probe FatJet 1 kinematics
  TH1D *_FatJet1_probe_pt =
      new TH1D("FatJet1_probe_pt", "FatJet1_probe_pt", 200, 0, 1000);
  TH1D *_FatJet1_probe_eta =
      new TH1D("FatJet1_probe_eta", "FatJet1_probe_eta", 100, -5, 5);
  TH1D *_FatJet1_probe_phi =
      new TH1D("FatJet1_probe_phi", "FatJet1_probe_phi", 100, -5, 5);
  TH2D *_FatJet1_probe_eta_phi = new TH2D(
      "FatJet1_probe_eta_phi", "FatJet1_probe_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet1_probe_Mass =
      new TH1D("FatJet1_probe_Mass", "FatJet1_probe_Mass", 500, 0, 500);
  TH1D *_FatJet1_probe_MassSD =
      new TH1D("FatJet1_probe_MassSD", "FatJet1_probe_MassSD", 500, 0, 500);
  TH2D *_FatJet1_probe_Mass_Pt =
      new TH2D("FatJet1_probe_Mass_Pt", "FatJet1_probe_Mass_Pt", num_m_bins,
               bins_m, num_pt_bins, bins_pt);
  // probe FatJet 1 ParticleNet scores
  TH1D *_FatJet1_probe_PNet_QCD =
      new TH1D("FatJet1_probe_PNet_QCD", "FatJet1_probePNet_QCD", 100, 0, 1.0);
  TH1D *_FatJet1_probe_PNet_QCD0HF = new TH1D(
      "FatJet1_probe_PNet_QCD0HF", "FatJet1_probe_PNet_QCD0HF", 100, 0, 1.0);
  TH1D *_FatJet1_probe_PNet_QCD1HF = new TH1D(
      "FatJet1_probe_PNet_QCD1HF", "FatJet1_probe_PNet_QCD1HF", 100, 0, 1.0);
  TH1D *_FatJet1_probe_PNet_QCD2HF = new TH1D(
      "FatJet1_probe_PNet_QCD2HF", "FatJet1_probe_PNet_QCD2HF", 100, 0, 1.0);
  TH1D *_FatJet1_probe_PNet_XbbVsQCD =
      new TH1D("FatJet1_probe_PNet_XbbVsQCD", "FatJet1_probe_PNet_XbbVsQCD",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNet_XccVsQCD =
      new TH1D("FatJet1_probe_PNet_XccVsQCD", "FatJet1_probe_PNet_XccVsQCD",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNet_XggVsQCD =
      new TH1D("FatJet1_probe_PNet_XggVsQCD", "FatJet1_probe_PNet_XggVsQCD",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNet_XqqVsQCD =
      new TH1D("FatJet1_probe_PNet_XqqVsQCD", "FatJet1_probe_PNet_XqqVsQCD",
               100, 0, 1.0);
  // probe FatJet 1 ParticleNetLegacy scores
  TH1D *_FatJet1_probe_PNetLegacy_Xbb =
      new TH1D("FatJet1_probe_PNetLegacy_Xbb", "FatJet1_probe_PNetLegacy_Xbb",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNetLegacy_XbbVsQCD =
      new TH1D("FatJet1_probe_PNetLegacy_XbbVsQCD", "FatJet1_probe_PNetLegacy_XbbVsQCD",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNetLegacy_Xcc =
      new TH1D("FatJet1_probe_PNetLegacy_Xcc", "FatJet1_probe_PNetLegacy_Xcc",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNetLegacy_Xqq =
      new TH1D("FatJet1_probe_PNetLegacy_Xqq", "FatJet1_probe_PNetLegacy_Xqq",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNetLegacy_QCD =
      new TH1D("FatJet1_probe_PNetLegacy_QCD", "FatJet1_probe_PNetLegacy_QCD",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNetLegacy_QCDb =
      new TH1D("FatJet1_probe_PNetLegacy_QCDb", "FatJet1_probe_PNetLegacy_QCDb",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_PNetLegacy_QCDbb =
      new TH1D("FatJet1_probe_PNetLegacy_QCDbb",
               "FatJet1_probe_PNetLegacy_QCDbb", 100, 0, 1.0);
  TH1D *_FatJet1_probe_PNetLegacy_QCDothers =
      new TH1D("FatJet1_probe_PNetLegacy_QCDothers",
               "FatJet1_probe_PNetLegacy_QCDothers", 100, 0, 1.0);
  // probe FatJet 1 GloParT scores
  TH1D *_FatJet1_probe_GloParT_QCD0HF =
      new TH1D("FatJet1_probe_GloParT_QCD0HF", "FatJet1_probe_GloParT_QCD0HF",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_GloParT_QCD1HF =
      new TH1D("FatJet1_probe_GloParT_QCD1HF", "FatJet1_probe_GloParT_QCD1HF",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_GloParT_QCD2HF =
      new TH1D("FatJet1_probe_GloParT_QCD2HF", "FatJet1_probe_GloParT_QCD2HF",
               100, 0, 1.0);
  TH1D *_FatJet1_probe_GloParT_Xbb = new TH1D(
      "FatJet1_probe_GloParT_Xbb", "FatJet1_probe_GloParT_Xbb", 100, 0, 1.0);
  TH1D *_FatJet1_probe_GloParT_Xcc = new TH1D(
      "FatJet1_probe_GloParT_Xcc", "FatJet1_probe_GloParT_Xcc", 100, 0, 1.0);
  TH1D *_FatJet1_probe_GloParT_Xqq = new TH1D(
      "FatJet1_probe_GloParT_Xqq", "FatJet1_probe_GloParT_Xqq", 100, 0, 1.0);
  TH1D *_FatJet1_probe_GloParT_XbbVsQCD =
      new TH1D("FatJet1_probe_GloParT_XbbVsQCD",
               "FatJet1_probe_GloParT_XbbVsQCD", 100, 0, 1.0);

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

  // Events Loop
  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);

    // HLT Selection
    bool EGamma = (HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) == 11);
    bool Muon = (HLT_IsoMu27 && fabs(lep1_Id) == 13);
    if (channel_lower == "egamma" or channel_lower == "electron") {
      if (!EGamma) {
        continue;
      }
    } else if (channel_lower == "muon") {
      if (!Muon) {
        continue;
      }
    } else if (channel_lower == "lepton" or channel_lower == "leptonic") {
      if (!EGamma && !Muon) {
        continue;
      }
    } else if (channel_lower == "jetmet" or channel_lower == "qcd") {
      if (HLT_AK8PFJet230_SoftDropMass40 == 0) {
        continue;
      }
    } else {
      throw std::invalid_argument("Invalid channel: " + channel);
    }

    // FatJets correction and selection
    if (FatJet1_pt > 0) {
      double Raw_FatJet1_tag_pt = FatJet1_pt * (1.0 - FatJet1_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_tag_pt);
      corrector_AK8->setJetEta(FatJet1_eta);
      corrector_AK8->setJetPhi(FatJet1_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet1_pt = Raw_FatJet1_tag_pt * corr;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * corr;
    }
    if (FatJet2_pt > 0) {
      double Raw_FatJet1_probe_pt = FatJet2_pt * (1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_probe_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      corrector_AK8->setJetPhi(FatJet2_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet1_probe_pt * corr;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * corr;
    }

    // FatJets selection
    bool probe_match = false;
    if (channel_lower == "jetmet" or channel_lower == "qcd") {
      if (FatJet1_pt <= 300 || fabs(FatJet1_eta) >= 2.5 ||
          FatJet1_MassSD <= 80) {
        continue;
      }
      if (FatJet2_pt <= 160 || fabs(FatJet2_eta) >= 2.5) {
        continue;
      }
      if (phi_dist(FatJet1_phi, FatJet2_phi) <= 2.5) {
        continue;
      }
      if (FatJet3_pt > 160) {
        continue;
      }

      // Trigger objects and matchings
      bool tag_match = false;
      for (int itrg = 0; itrg < NTrigger_Objects; itrg++) {
        if ((Trigger_Object_bit[itrg] & 4) == 4) {
          double dR = get_dR(FatJet1_eta, FatJet1_phi, Trigger_Object_eta[itrg],
                             Trigger_Object_phi[itrg]);
          if (dR < 0.4 && Trigger_Object_pt[itrg] > 100) {
            tag_match = true;
            break;
          }
        }
      }
      if (!tag_match) {
        continue;
      }

      // Probe matching
      bool matched_to_AK8PFJet230_SoftDropMass40 = false;
      for (int itrg = 0; itrg < NTrigger_Objects; itrg++) {
        if ((Trigger_Object_bit[itrg] & 4) == 4) {
          double dR = get_dR(FatJet2_eta, FatJet2_phi, Trigger_Object_eta[itrg],
                             Trigger_Object_phi[itrg]);
          if (dR < 0.4 && Trigger_Object_pt[itrg] > 100) {
            matched_to_AK8PFJet230_SoftDropMass40 = true;
            break;
          }
        }
      }
      if (matched_to_AK8PFJet230_SoftDropMass40) {
        probe_match = true;
      }
    } else {
      // EGamma, Muon, Lepton (EGamma + Muon)
      if (lep1_Pt <= 50 || fabs(lep1_Eta) >= 2.4) {
        continue;
      }
      if (FatJet1_pt <= 160 || fabs(FatJet1_eta) >= 2.5) {
        continue;
      }
      if (phi_dist(FatJet1_phi, lep1_Phi) <= 2.0) {
        continue;
      }
      if (FatJet2_pt > 160) {
        continue;
      }
      if (lep2_Pt > 30) {
        continue;
      }
      if (MET <= 50) {
        continue;
      }

      bool has_valid_tagged_AK4Jet = false;
      if (Jet1_Pt > 40 && fabs(Jet1_Eta) < 2.4) {
        double dR_J1_FJ = get_dR(Jet1_Eta, Jet1_Phi, FatJet1_eta, FatJet1_phi);
        double dR_J1_L = get_dR(Jet1_Eta, Jet1_Phi, lep1_Eta, lep1_Phi);

        if (dR_J1_FJ > 1.5 && dR_J1_L > 0.5) {
          has_valid_tagged_AK4Jet = true;
        }
      }

      if (!has_valid_tagged_AK4Jet && Jet2_Pt > 40 && fabs(Jet2_Eta) < 2.4) {
        double dR_J2_FJ = get_dR(Jet2_Eta, Jet2_Phi, FatJet1_eta, FatJet1_phi);
        double dR_J2_L = get_dR(Jet2_Eta, Jet2_Phi, lep1_Eta, lep1_Phi);

        if (dR_J2_FJ > 1.5 && dR_J2_L > 0.5) {
          has_valid_tagged_AK4Jet = true;
        }
      }

      if (!has_valid_tagged_AK4Jet) {
        continue;
      }

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
      if (matched_to_AK8PFJet230_SoftDropMass40) {
        probe_match = true;
      }
    }

    // weight
    weight = (weight / SumGenWeights) * xsec * param_dict.Lumi;
    double PU_weight = PU_Rew[(int)npu];
    if (PU_weight < 20.0) {
      weight = weight * PU_weight;
    }

    // Fill histograms
    // Kinematics
    _FatJet1_tag_pt->Fill(FatJet1_pt, weight);
    _FatJet1_tag_eta->Fill(FatJet1_eta, weight);
    _FatJet1_tag_phi->Fill(FatJet1_phi, weight);
    _FatJet1_tag_eta_phi->Fill(FatJet1_eta, FatJet1_phi, weight);
    _FatJet1_tag_Mass->Fill(FatJet1_Mass, weight);
    _FatJet1_tag_MassSD->Fill(FatJet1_MassSD, weight);
    _FatJet1_tag_Mass_Pt->Fill(FatJet1_MassSD, FatJet1_pt, weight);

    // ParticleNet
    _FatJet1_tag_PNet_QCD->Fill(FatJet1PNet_QCD, weight);
    _FatJet1_tag_PNet_QCD0HF->Fill(FatJet1PNet_QCD0HF, weight);
    _FatJet1_tag_PNet_QCD1HF->Fill(FatJet1PNet_QCD1HF, weight);
    _FatJet1_tag_PNet_QCD2HF->Fill(FatJet1PNet_QCD2HF, weight);
    _FatJet1_tag_PNet_XbbVsQCD->Fill(FatJet1PNet_XbbVsQCD, weight);
    _FatJet1_tag_PNet_XccVsQCD->Fill(FatJet1PNet_XccVsQCD, weight);
    _FatJet1_tag_PNet_XggVsQCD->Fill(FatJet1PNet_XggVsQCD, weight);
    _FatJet1_tag_PNet_XqqVsQCD->Fill(FatJet1PNet_XqqVsQCD, weight);
    // ParticleNetLegacy
    Float_t FatJet1PNetLegacy_XbbVsQCD =
        FatJet1PNetLegacy_Xbb / (FatJet1PNetLegacy_Xbb + FatJet1PNetLegacy_QCD + 1e-12);
    _FatJet1_tag_PNetLegacy_Xbb->Fill(FatJet1PNetLegacy_Xbb, weight);
    _FatJet1_tag_PNetLegacy_XbbVsQCD->Fill(FatJet1PNetLegacy_XbbVsQCD, weight);
    _FatJet1_tag_PNetLegacy_Xcc->Fill(FatJet1PNetLegacy_Xcc, weight);
    _FatJet1_tag_PNetLegacy_Xqq->Fill(FatJet1PNetLegacy_Xqq, weight);
    _FatJet1_tag_PNetLegacy_QCD->Fill(FatJet1PNetLegacy_QCD, weight);
    _FatJet1_tag_PNetLegacy_QCDb->Fill(FatJet1PNetLegacy_QCDb, weight);
    _FatJet1_tag_PNetLegacy_QCDbb->Fill(FatJet1PNetLegacy_QCDbb, weight);
    _FatJet1_tag_PNetLegacy_QCDothers->Fill(FatJet1PNetLegacy_QCDothers,
                                            weight);
    // GloParT
    _FatJet1_tag_GloParT_QCD0HF->Fill(FatJet1GloParT_QCD0HF, weight);
    _FatJet1_tag_GloParT_QCD1HF->Fill(FatJet1GloParT_QCD1HF, weight);
    _FatJet1_tag_GloParT_QCD2HF->Fill(FatJet1GloParT_QCD2HF, weight);
    _FatJet1_tag_GloParT_Xbb->Fill(FatJet1GloParT_Xbb, weight);
    _FatJet1_tag_GloParT_Xcc->Fill(FatJet1GloParT_Xcc, weight);
    _FatJet1_tag_GloParT_Xqq->Fill(FatJet1GloParT_Xqq, weight);
    _FatJet1_tag_GloParT_XbbVsQCD->Fill(FatJet1GloParT_XbbVsQCD, weight);

    if (probe_match) {
      // Kinematics
      _FatJet1_probe_pt->Fill(FatJet1_pt, weight);
      _FatJet1_probe_eta->Fill(FatJet1_eta, weight);
      _FatJet1_probe_phi->Fill(FatJet1_phi, weight);
      _FatJet1_probe_eta_phi->Fill(FatJet1_eta, FatJet1_phi, weight);
      _FatJet1_probe_Mass->Fill(FatJet1_Mass, weight);
      _FatJet1_probe_MassSD->Fill(FatJet1_MassSD, weight);
      _FatJet1_probe_Mass_Pt->Fill(FatJet1_MassSD, FatJet1_pt, weight);
      // ParticleNet
      _FatJet1_probe_PNet_QCD->Fill(FatJet1PNet_QCD, weight);
      _FatJet1_probe_PNet_QCD0HF->Fill(FatJet1PNet_QCD0HF, weight);
      _FatJet1_probe_PNet_QCD1HF->Fill(FatJet1PNet_QCD1HF, weight);
      _FatJet1_probe_PNet_QCD2HF->Fill(FatJet1PNet_QCD2HF, weight);
      _FatJet1_probe_PNet_XbbVsQCD->Fill(FatJet1PNet_XbbVsQCD, weight);
      _FatJet1_probe_PNet_XccVsQCD->Fill(FatJet1PNet_XccVsQCD, weight);
      _FatJet1_probe_PNet_XggVsQCD->Fill(FatJet1PNet_XggVsQCD, weight);
      _FatJet1_probe_PNet_XqqVsQCD->Fill(FatJet1PNet_XqqVsQCD, weight);
      // ParticleNetLegacy
      _FatJet1_probe_PNetLegacy_Xbb->Fill(FatJet1PNetLegacy_Xbb, weight);
      _FatJet1_probe_PNetLegacy_XbbVsQCD->Fill(FatJet1PNetLegacy_XbbVsQCD, weight);
      _FatJet1_probe_PNetLegacy_Xcc->Fill(FatJet1PNetLegacy_Xcc, weight);
      _FatJet1_probe_PNetLegacy_Xqq->Fill(FatJet1PNetLegacy_Xqq, weight);
      _FatJet1_probe_PNetLegacy_QCD->Fill(FatJet1PNetLegacy_QCD, weight);
      _FatJet1_probe_PNetLegacy_QCDb->Fill(FatJet1PNetLegacy_QCDb, weight);
      _FatJet1_probe_PNetLegacy_QCDbb->Fill(FatJet1PNetLegacy_QCDbb, weight);
      _FatJet1_probe_PNetLegacy_QCDothers->Fill(FatJet1PNetLegacy_QCDothers,
                                                weight);
      // GloParT
      _FatJet1_probe_GloParT_QCD0HF->Fill(FatJet1GloParT_QCD0HF, weight);
      _FatJet1_probe_GloParT_QCD1HF->Fill(FatJet1GloParT_QCD1HF, weight);
      _FatJet1_probe_GloParT_QCD2HF->Fill(FatJet1GloParT_QCD2HF, weight);
      _FatJet1_probe_GloParT_Xbb->Fill(FatJet1GloParT_Xbb, weight);
      _FatJet1_probe_GloParT_Xcc->Fill(FatJet1GloParT_Xcc, weight);
      _FatJet1_probe_GloParT_Xqq->Fill(FatJet1GloParT_Xqq, weight);
      _FatJet1_probe_GloParT_XbbVsQCD->Fill(FatJet1GloParT_XbbVsQCD, weight);
    }

  } // end event loop

  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
