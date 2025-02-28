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

struct ParamDict {
  double Lumi;
  // QCD HT-binned samples
  double XSec_QCD_HT_100to200;
  double XSec_QCD_HT_200to400;
  double XSec_QCD_HT_400to600;
  double XSec_QCD_HT_600to800;
  double XSec_QCD_HT_800to1000;
  double XSec_QCD_HT_1000to1200;
  double XSec_QCD_HT_1200to1500;
  double XSec_QCD_HT_1500to2000;
  double XSec_QCD_HT_2000toInf;

  // TT samples
  double XSec_TTto4Q;
  double XSec_TTtoLNu2Q;
  double XSec_TTto2L2Nu;

  // ttHto2B sample
  double XSec_ttHto2B_M_125;

  // Diboson samples
  double XSec_WW;
  double XSec_WZ;
  double XSec_ZZ;

  // Z+jets HT-binned samples
  double XSec_Zto2Q_4Jets_HT_200to400;
  double XSec_Zto2Q_4Jets_HT_400to600;
  double XSec_Zto2Q_4Jets_HT_600to800;
  double XSec_Zto2Q_4Jets_HT_800toInf;

  // W+jets HT-binned samples
  double XSec_Wto2Q_3Jets_HT_200to400;
  double XSec_Wto2Q_3Jets_HT_400to600;
  double XSec_Wto2Q_3Jets_HT_600to800;
  double XSec_Wto2Q_3Jets_HT_800toInf;

  // W+jets PTQQ-binned samples
  double XSec_Wto2Q_2Jets_PTQQ_100to200_1J;
  double XSec_Wto2Q_2Jets_PTQQ_100to200_2J;
  double XSec_Wto2Q_2Jets_PTQQ_200to400_1J;
  double XSec_Wto2Q_2Jets_PTQQ_200to400_2J;
  double XSec_Wto2Q_2Jets_PTQQ_400to600_1J;
  double XSec_Wto2Q_2Jets_PTQQ_400to600_2J;
  double XSec_Wto2Q_2Jets_PTQQ_600_1J;
  double XSec_Wto2Q_2Jets_PTQQ_600_2J;

  // Z+jets PTQQ-binned samples
  double XSec_Zto2Q_2Jets_PTQQ_100to200_1J;
  double XSec_Zto2Q_2Jets_PTQQ_100to200_2J;
  double XSec_Zto2Q_2Jets_PTQQ_200to400_1J;
  double XSec_Zto2Q_2Jets_PTQQ_200to400_2J;
  double XSec_Zto2Q_2Jets_PTQQ_400to600_1J;
  double XSec_Zto2Q_2Jets_PTQQ_400to600_2J;
  double XSec_Zto2Q_2Jets_PTQQ_600_1J;
  double XSec_Zto2Q_2Jets_PTQQ_600_2J;

  // DY+jets samples
  double XSec_DYto2L_2Jets_MLL_50_0J;
  double XSec_DYto2L_2Jets_MLL_50_1J;
  double XSec_DYto2L_2Jets_MLL_50_2J;

  // W+leptons samples
  double XSec_WtoLNu_2Jets_0J;
  double XSec_WtoLNu_2Jets_1J;
  double XSec_WtoLNu_2Jets_2J;
  double XSec_WtoLNu_4Jets_1J;
  double XSec_WtoLNu_4Jets_2J;
  double XSec_WtoLNu_4Jets_3J;

  // HH sample
  double XSec_GluGlutoHHto4B;
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

      // Parse line of format: double XSec_QCD_HT_100to200 = 25220000.00;
      if (!(iss >> type >> name >> equals >> value)) {
        continue;
      }

      // Remove semicolon if present
      if (name.back() == ';')
        name = name.substr(0, name.size() - 1);

      // Use pointer to member to make assignment more concise
      using MemberPtr = double ParamDict::*;
      std::map<std::string, MemberPtr> memberMap = {
          {"Lumi", &ParamDict::Lumi},
          // QCD
          {"XSec_QCD_HT_100to200", &ParamDict::XSec_QCD_HT_100to200},
          {"XSec_QCD_HT_200to400", &ParamDict::XSec_QCD_HT_200to400},
          {"XSec_QCD_HT_400to600", &ParamDict::XSec_QCD_HT_400to600},
          {"XSec_QCD_HT_600to800", &ParamDict::XSec_QCD_HT_600to800},
          {"XSec_QCD_HT_800to1000", &ParamDict::XSec_QCD_HT_800to1000},
          {"XSec_QCD_HT_1000to1200", &ParamDict::XSec_QCD_HT_1000to1200},
          {"XSec_QCD_HT_1200to1500", &ParamDict::XSec_QCD_HT_1200to1500},
          {"XSec_QCD_HT_1500to2000", &ParamDict::XSec_QCD_HT_1500to2000},
          {"XSec_QCD_HT_2000toInf", &ParamDict::XSec_QCD_HT_2000toInf},
          // TT
          {"XSec_TTto4Q", &ParamDict::XSec_TTto4Q},
          {"XSec_TTtoLNu2Q", &ParamDict::XSec_TTtoLNu2Q},
          {"XSec_TTto2L2Nu", &ParamDict::XSec_TTto2L2Nu},
          // ttH
          {"XSec_ttHto2B_M_125", &ParamDict::XSec_ttHto2B_M_125},
          // Diboson
          {"XSec_WW", &ParamDict::XSec_WW},
          {"XSec_WZ", &ParamDict::XSec_WZ},
          {"XSec_ZZ", &ParamDict::XSec_ZZ},
          // Z+jets HT
          {"XSec_Zto2Q_4Jets_HT_200to400",
           &ParamDict::XSec_Zto2Q_4Jets_HT_200to400},
          {"XSec_Zto2Q_4Jets_HT_400to600",
           &ParamDict::XSec_Zto2Q_4Jets_HT_400to600},
          {"XSec_Zto2Q_4Jets_HT_600to800",
           &ParamDict::XSec_Zto2Q_4Jets_HT_600to800},
          {"XSec_Zto2Q_4Jets_HT_800toInf",
           &ParamDict::XSec_Zto2Q_4Jets_HT_800toInf},
          // W+jets HT
          {"XSec_Wto2Q_3Jets_HT_200to400",
           &ParamDict::XSec_Wto2Q_3Jets_HT_200to400},
          {"XSec_Wto2Q_3Jets_HT_400to600",
           &ParamDict::XSec_Wto2Q_3Jets_HT_400to600},
          {"XSec_Wto2Q_3Jets_HT_600to800",
           &ParamDict::XSec_Wto2Q_3Jets_HT_600to800},
          {"XSec_Wto2Q_3Jets_HT_800toInf",
           &ParamDict::XSec_Wto2Q_3Jets_HT_800toInf},
          // W+jets PTQQ
          {"XSec_Wto2Q_2Jets_PTQQ_100to200_1J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_1J},
          {"XSec_Wto2Q_2Jets_PTQQ_100to200_2J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_2J},
          {"XSec_Wto2Q_2Jets_PTQQ_200to400_1J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_1J},
          {"XSec_Wto2Q_2Jets_PTQQ_200to400_2J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_2J},
          {"XSec_Wto2Q_2Jets_PTQQ_400to600_1J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_1J},
          {"XSec_Wto2Q_2Jets_PTQQ_400to600_2J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_2J},
          {"XSec_Wto2Q_2Jets_PTQQ_600_1J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_1J},
          {"XSec_Wto2Q_2Jets_PTQQ_600_2J",
           &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_2J},
          // Z+jets PTQQ
          {"XSec_Zto2Q_2Jets_PTQQ_100to200_1J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_1J},
          {"XSec_Zto2Q_2Jets_PTQQ_100to200_2J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_2J},
          {"XSec_Zto2Q_2Jets_PTQQ_200to400_1J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_1J},
          {"XSec_Zto2Q_2Jets_PTQQ_200to400_2J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_2J},
          {"XSec_Zto2Q_2Jets_PTQQ_400to600_1J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_1J},
          {"XSec_Zto2Q_2Jets_PTQQ_400to600_2J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_2J},
          {"XSec_Zto2Q_2Jets_PTQQ_600_1J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_600_1J},
          {"XSec_Zto2Q_2Jets_PTQQ_600_2J",
           &ParamDict::XSec_Zto2Q_2Jets_PTQQ_600_2J},
          // DY+jets
          {"XSec_DYto2L_2Jets_MLL_50_0J",
           &ParamDict::XSec_DYto2L_2Jets_MLL_50_0J},
          {"XSec_DYto2L_2Jets_MLL_50_1J",
           &ParamDict::XSec_DYto2L_2Jets_MLL_50_1J},
          {"XSec_DYto2L_2Jets_MLL_50_2J",
           &ParamDict::XSec_DYto2L_2Jets_MLL_50_2J},
          // W+leptons
          {"XSec_WtoLNu_2Jets_0J", &ParamDict::XSec_WtoLNu_2Jets_0J},
          {"XSec_WtoLNu_2Jets_1J", &ParamDict::XSec_WtoLNu_2Jets_1J},
          {"XSec_WtoLNu_2Jets_2J", &ParamDict::XSec_WtoLNu_2Jets_2J},
          {"XSec_WtoLNu_4Jets_1J", &ParamDict::XSec_WtoLNu_4Jets_1J},
          {"XSec_WtoLNu_4Jets_2J", &ParamDict::XSec_WtoLNu_4Jets_2J},
          {"XSec_WtoLNu_4Jets_3J", &ParamDict::XSec_WtoLNu_4Jets_3J},
          // HH
          {"XSec_GluGlutoHHto4B", &ParamDict::XSec_GluGlutoHHto4B}};

      auto it = memberMap.find(name);
      if (it != memberMap.end()) {
        param_dict->*(it->second) = value;
      } else {
        std::cerr << "Unknown parameter: " << name << std::endl;
      }
    }
  }
}

double getXSec(ParamDict *param_dict, std::string data_type) {
  std::map<std::string, double ParamDict::*> xsecMap = {
      // TT samples
      {"TTto4Q", &ParamDict::XSec_TTto4Q},
      {"TTtoLNu2Q", &ParamDict::XSec_TTtoLNu2Q},
      {"TTto2L2Nu", &ParamDict::XSec_TTto2L2Nu},
      // ttHto2B sample
      {"ttHto2B_M_125", &ParamDict::XSec_ttHto2B_M_125},
      // QCD samples
      {"QCD_HT_100to200", &ParamDict::XSec_QCD_HT_100to200},
      {"QCD_HT_200to400", &ParamDict::XSec_QCD_HT_200to400},
      {"QCD_HT_400to600", &ParamDict::XSec_QCD_HT_400to600},
      {"QCD_HT_600to800", &ParamDict::XSec_QCD_HT_600to800},
      {"QCD_HT_800to1000", &ParamDict::XSec_QCD_HT_800to1000},
      {"QCD_HT_1000to1200", &ParamDict::XSec_QCD_HT_1000to1200},
      {"QCD_HT_1200to1500", &ParamDict::XSec_QCD_HT_1200to1500},
      {"QCD_HT_1500to2000", &ParamDict::XSec_QCD_HT_1500to2000},
      {"QCD_HT_2000toInf", &ParamDict::XSec_QCD_HT_2000toInf},
      // Diboson samples
      {"WW", &ParamDict::XSec_WW},
      {"WZ", &ParamDict::XSec_WZ},
      {"ZZ", &ParamDict::XSec_ZZ},
      // Z+jets HT samples
      {"Zto2Q_4Jets_HT_200to400", &ParamDict::XSec_Zto2Q_4Jets_HT_200to400},
      {"Zto2Q_4Jets_HT_400to600", &ParamDict::XSec_Zto2Q_4Jets_HT_400to600},
      {"Zto2Q_4Jets_HT_600to800", &ParamDict::XSec_Zto2Q_4Jets_HT_600to800},
      {"Zto2Q_4Jets_HT_800toInf", &ParamDict::XSec_Zto2Q_4Jets_HT_800toInf},
      // W+jets HT samples
      {"Wto2Q_3Jets_HT_200to400", &ParamDict::XSec_Wto2Q_3Jets_HT_200to400},
      {"Wto2Q_3Jets_HT_400to600", &ParamDict::XSec_Wto2Q_3Jets_HT_400to600},
      {"Wto2Q_3Jets_HT_600to800", &ParamDict::XSec_Wto2Q_3Jets_HT_600to800},
      {"Wto2Q_3Jets_HT_800toInf", &ParamDict::XSec_Wto2Q_3Jets_HT_800toInf},
      // W+jets PTQQ samples
      {"Wto2Q_2Jets_PTQQ_100to200_1J",
       &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_1J},
      {"Wto2Q_2Jets_PTQQ_100to200_2J",
       &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_2J},
      {"Wto2Q_2Jets_PTQQ_200to400_1J",
       &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_1J},
      {"Wto2Q_2Jets_PTQQ_200to400_2J",
       &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_2J},
      {"Wto2Q_2Jets_PTQQ_400to600_1J",
       &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_1J},
      {"Wto2Q_2Jets_PTQQ_400to600_2J",
       &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_2J},
      {"Wto2Q_2Jets_PTQQ_600_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_1J},
      {"Wto2Q_2Jets_PTQQ_600_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_2J},
      // Z+jets PTQQ samples
      {"Zto2Q_2Jets_PTQQ_100to200_1J",
       &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_1J},
      {"Zto2Q_2Jets_PTQQ_100to200_2J",
       &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_2J},
      {"Zto2Q_2Jets_PTQQ_200to400_1J",
       &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_1J},
      {"Zto2Q_2Jets_PTQQ_200to400_2J",
       &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_2J},
      {"Zto2Q_2Jets_PTQQ_400to600_1J",
       &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_1J},
      {"Zto2Q_2Jets_PTQQ_400to600_2J",
       &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_2J},
      {"Zto2Q_2Jets_PTQQ_600_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_600_1J},
      {"Zto2Q_2Jets_PTQQ_600_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_600_2J},
      // DY+jets samples
      {"DYto2L_2Jets_MLL_50_0J", &ParamDict::XSec_DYto2L_2Jets_MLL_50_0J},
      {"DYto2L_2Jets_MLL_50_1J", &ParamDict::XSec_DYto2L_2Jets_MLL_50_1J},
      {"DYto2L_2Jets_MLL_50_2J", &ParamDict::XSec_DYto2L_2Jets_MLL_50_2J},
      // W+leptons samples
      {"WtoLNu_2Jets_0J", &ParamDict::XSec_WtoLNu_2Jets_0J},
      {"WtoLNu_2Jets_1J", &ParamDict::XSec_WtoLNu_2Jets_1J},
      {"WtoLNu_2Jets_2J", &ParamDict::XSec_WtoLNu_2Jets_2J},
      {"WtoLNu_4Jets_1J", &ParamDict::XSec_WtoLNu_4Jets_1J},
      {"WtoLNu_4Jets_2J", &ParamDict::XSec_WtoLNu_4Jets_2J},
      {"WtoLNu_4Jets_3J", &ParamDict::XSec_WtoLNu_4Jets_3J},
      // HH sample
      {"GluGlutoHHto4B", &ParamDict::XSec_GluGlutoHHto4B}};

  auto it = xsecMap.find(data_type);
  if (it != xsecMap.end()) {
    return param_dict->*(it->second);
  }
  throw std::invalid_argument("Invalid data_type: " + data_type);
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

// JEC and JER
struct JetCorrectionResult {
  double corrected_pt;
  double corrected_massSD;
};

JetCorrectionResult applyJEC(double jet_pt, double jet_eta, double jet_phi,
                             double jet_rawFactor, double jet_massSD,
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

JetCorrectionResult
applyJER(double jet_pt, double jet_eta, double jet_phi, double jet_massSD,
         const JME::JetResolution &resolution_pt,
         const JME::JetResolutionScaleFactor &resolution_pt_sf, double rho,
         const Float_t *genJetPts, const Float_t *genJetEtas,
         const Float_t *genJetPhis, int nGenJets) {

  JME::JetParameters JerPARAM = {{JME::Binning::JetPt, jet_pt},
                                 {JME::Binning::JetEta, jet_eta},
                                 {JME::Binning::Rho, rho}};

  JME::JetParameters JerSFPARAM;
  JerSFPARAM.set(JME::Binning::JetPt, jet_pt);
  JerSFPARAM.set(JME::Binning::JetEta, jet_eta);
  JerSFPARAM.set(JME::Binning::Rho, rho);

  double resolution = resolution_pt.getResolution(JerPARAM);
  double resolution_sf = resolution_pt_sf.getScaleFactor(JerSFPARAM);

  const double MAX_DELTA_R = 0.2;
  double smearFactor = 1.0;
  bool gen_matched = false;

  for (int nGJ = 0; nGJ < nGenJets; nGJ++) {
    double delta_R =
        get_dR(jet_eta, jet_phi, genJetEtas[nGJ], genJetPhis[nGJ]);
    double pt_diff_ratio = fabs(jet_pt - genJetPts[nGJ]) / jet_pt;
    double resolution_threshold = 3 * resolution;

    if (!gen_matched && delta_R < MAX_DELTA_R &&
        pt_diff_ratio < resolution_threshold) {
      smearFactor =
          1.0 + (resolution_sf - 1.0) * (jet_pt - genJetPts[nGJ]) / jet_pt;
      gen_matched = true;
      break;
    }
  }

  if (!gen_matched && resolution_sf > 1.0) {
    double sigma = resolution * sqrt(resolution_sf * resolution_sf - 1);
    smearFactor = 1.0 + gRandom->Gaus(0, sigma);
  }

  return JetCorrectionResult{jet_pt * smearFactor, jet_massSD * smearFactor};
}

void tree_mc(const std::string &year,      // 2022, 2023
             const std::string &data_type, // Example: TTtoLNu2Q, QCD_HT100to200
             const std::string &sample_path, // path to the root file
             const std::string &output_path, // path to the output root file
             const std::string &pu_path, // path to the pileup reweighting file
             const std::string &param_path, // path to the parameters file
             const std::string &jec_path,   // path to the AK8 JEC txt file
             const std::string &jer_path,   // path to the AK8 JER txt file
             const std::string &jer_path_sf // path to the AK8 JER SF txt file
) {
  gSystem->Load("libFWCoreFWLite.so");

  if (year != "2022" && year != "2023") {
    throw std::invalid_argument("Invalid year: " + year);
  }

  std::cout << "Data Type: " << data_type << std::endl;
  std::cout << "Sample Path: " << sample_path << std::endl;
  std::cout << "Output Path: " << output_path << std::endl;
  std::cout << "PU Path: " << pu_path << std::endl;
  std::cout << "Param Path: " << param_path << std::endl;
  std::cout << "JEC Path: " << jec_path << std::endl;
  std::cout << "JER Path: " << jer_path << std::endl;
  std::cout << "JER SF Path: " << jer_path_sf << std::endl;

  // parse ParamDict
  ParamDict param_dict;
  loadParamDict(&param_dict, param_path);
  std::cout << "Lumi: " << param_dict.Lumi << std::endl;
  double xsec = getXSec(&param_dict, data_type);
  std::cout << "XSec: " << xsec << std::endl;

  // pu weight
  std::vector<double> PU_Rew = loadPUReweighting(pu_path);

  // JEC
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(JetCorrectorParameters(jec_path.c_str()));
  FactorizedJetCorrector *corrector = new FactorizedJetCorrector(vPar);

  // JER
  JME::JetResolution resolution_pt = JME::JetResolution(jer_path.c_str());
  JME::JetResolutionScaleFactor resolution_pt_sf =
      JME::JetResolutionScaleFactor(jer_path_sf.c_str());

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

  TH1F *NEvents = (TH1F *)f1->Get("NEvents");
  double SumGenWeights = NEvents->GetBinContent(1);
  std::cout << "SumGenWeights: " << SumGenWeights << std::endl;
  TTree *InputTree = (TTree *)f1->Get("tree");

  Float_t weight;
  UInt_t run;
  UInt_t lumi;
  Float_t npu;
  Int_t isVBFtag;
  Float_t rho;

  Bool_t HLT_Ele32_WPTight_Gsf;
  Bool_t HLT_IsoMu27;

  Float_t MET;

  Float_t lep1_pt;
  Float_t lep1_eta;
  Float_t lep1_phi;
  Int_t lep1_Id;
  Float_t lep2_pt;
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

  // gen-level
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
  InputTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);

  InputTree->SetBranchAddress("lep1Pt", &lep1_pt);
  InputTree->SetBranchAddress("lep1Eta", &lep1_eta);
  InputTree->SetBranchAddress("lep1Phi", &lep1_phi);
  InputTree->SetBranchAddress("lep1Id", &lep1_Id);
  InputTree->SetBranchAddress("lep2Pt", &lep2_pt);
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

  // gen-level jets
  InputTree->SetBranchAddress("nGenJet", &nGenJet);
  InputTree->SetBranchAddress("GenJet_eta", GenJet_eta);
  InputTree->SetBranchAddress("GenJet_phi", GenJet_phi);
  InputTree->SetBranchAddress("GenJet_pt", GenJet_pt);
  InputTree->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
  InputTree->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta);
  InputTree->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi);
  InputTree->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt);

  // Events Loop
  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);

    // HLT Selection
    bool HLT_ele = (HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) == 11);
    bool HLT_mu = (HLT_IsoMu27 && fabs(lep1_Id) == 13);
    if (!HLT_ele && !HLT_mu) {
      continue;
    }

    // FatJets correction (JEC and JER)
    JetCorrectionResult jec1 =
        applyJEC(FatJet1_pt, FatJet1_eta, FatJet1_phi, FatJet1_rawFactor,
                 FatJet1_MassSD, corrector);
    FatJet1_pt = jec1.corrected_pt;
    FatJet1_MassSD = jec1.corrected_massSD;

    JetCorrectionResult jer1 =
        applyJER(FatJet1_pt, FatJet1_eta, FatJet1_phi, FatJet1_MassSD,
                 resolution_pt, resolution_pt_sf, rho, GenJetAK8_pt,
                 GenJetAK8_eta, GenJetAK8_phi, nGenJetAK8);
    FatJet1_pt = jer1.corrected_pt;
    FatJet1_MassSD = jer1.corrected_massSD;

    JetCorrectionResult jec2 =
        applyJEC(FatJet2_pt, FatJet2_eta, FatJet2_phi, FatJet2_rawFactor,
                 FatJet2_MassSD, corrector);
    FatJet2_pt = jec2.corrected_pt;
    FatJet2_MassSD = jec2.corrected_massSD;

    JetCorrectionResult jer2 =
        applyJER(FatJet2_pt, FatJet2_eta, FatJet2_phi, FatJet2_MassSD,
                 resolution_pt, resolution_pt_sf, rho, GenJetAK8_pt,
                 GenJetAK8_eta, GenJetAK8_phi, nGenJetAK8);
    FatJet2_pt = jer2.corrected_pt;
    FatJet2_MassSD = jer2.corrected_massSD;

    // Selection
    if (year == "2022") {
      if (FatJet1_pt < 270 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50) {
          continue;
      }
      if (FatJet2_pt > 250 && FatJet2_MassSD > 50) {
        continue;
      }
      if (lep1_pt < 50 || lep2_pt > 30) {
        continue;
      }
      if (MET < 50) {
        continue;
      }
    } else if (year == "2023") {
      if (FatJet1_pt < 250 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50) {
      continue;
      }
      if (FatJet2_pt > 200 && FatJet2_MassSD > 50) {
        continue;
      }
      if (lep1_pt < 50 || lep2_pt > 30) {
        continue;
      }
      if (MET < 50) {
        continue;
      }
    } else {
      throw std::invalid_argument("Invalid year: " + year);
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

    // weight
    weight = (weight / SumGenWeights) * xsec * param_dict.Lumi;
    double PU_weight = PU_Rew[(int)npu];
    if (PU_weight < 20.0) {
      weight = weight * PU_weight;
    }

    T_weight = weight;
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
