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
#include "TXbb.h"

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
                // ttHto2B
                {"XSec_ttHto2B_M_125", &ParamDict::XSec_ttHto2B_M_125},
                // Diboson
                {"XSec_WW", &ParamDict::XSec_WW},
                {"XSec_WZ", &ParamDict::XSec_WZ},
                {"XSec_ZZ", &ParamDict::XSec_ZZ},
                // Z+jets HT
                {"XSec_Zto2Q_4Jets_HT_200to400", &ParamDict::XSec_Zto2Q_4Jets_HT_200to400},
                {"XSec_Zto2Q_4Jets_HT_400to600", &ParamDict::XSec_Zto2Q_4Jets_HT_400to600},
                {"XSec_Zto2Q_4Jets_HT_600to800", &ParamDict::XSec_Zto2Q_4Jets_HT_600to800},
                {"XSec_Zto2Q_4Jets_HT_800toInf", &ParamDict::XSec_Zto2Q_4Jets_HT_800toInf},
                // W+jets HT
                {"XSec_Wto2Q_3Jets_HT_200to400", &ParamDict::XSec_Wto2Q_3Jets_HT_200to400},
                {"XSec_Wto2Q_3Jets_HT_400to600", &ParamDict::XSec_Wto2Q_3Jets_HT_400to600},
                {"XSec_Wto2Q_3Jets_HT_600to800", &ParamDict::XSec_Wto2Q_3Jets_HT_600to800},
                {"XSec_Wto2Q_3Jets_HT_800toInf", &ParamDict::XSec_Wto2Q_3Jets_HT_800toInf},
                // W+jets PTQQ
                {"XSec_Wto2Q_2Jets_PTQQ_100to200_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_1J},
                {"XSec_Wto2Q_2Jets_PTQQ_100to200_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_2J},
                {"XSec_Wto2Q_2Jets_PTQQ_200to400_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_1J},
                {"XSec_Wto2Q_2Jets_PTQQ_200to400_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_2J},
                {"XSec_Wto2Q_2Jets_PTQQ_400to600_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_1J},
                {"XSec_Wto2Q_2Jets_PTQQ_400to600_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_2J},
                {"XSec_Wto2Q_2Jets_PTQQ_600_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_1J},
                {"XSec_Wto2Q_2Jets_PTQQ_600_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_2J},
                // Z+jets PTQQ
                {"XSec_Zto2Q_2Jets_PTQQ_100to200_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_1J},
                {"XSec_Zto2Q_2Jets_PTQQ_100to200_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_2J},
                {"XSec_Zto2Q_2Jets_PTQQ_200to400_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_1J},
                {"XSec_Zto2Q_2Jets_PTQQ_200to400_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_2J},
                {"XSec_Zto2Q_2Jets_PTQQ_400to600_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_1J},
                {"XSec_Zto2Q_2Jets_PTQQ_400to600_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_2J},
                {"XSec_Zto2Q_2Jets_PTQQ_600_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_600_1J},
                {"XSec_Zto2Q_2Jets_PTQQ_600_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_600_2J},
                // DY+jets
                {"XSec_DYto2L_2Jets_MLL_50_0J", &ParamDict::XSec_DYto2L_2Jets_MLL_50_0J},
                {"XSec_DYto2L_2Jets_MLL_50_1J", &ParamDict::XSec_DYto2L_2Jets_MLL_50_1J},
                {"XSec_DYto2L_2Jets_MLL_50_2J", &ParamDict::XSec_DYto2L_2Jets_MLL_50_2J},
                // W+leptons
                {"XSec_WtoLNu_2Jets_0J", &ParamDict::XSec_WtoLNu_2Jets_0J},
                {"XSec_WtoLNu_2Jets_1J", &ParamDict::XSec_WtoLNu_2Jets_1J},
                {"XSec_WtoLNu_2Jets_2J", &ParamDict::XSec_WtoLNu_2Jets_2J},
                {"XSec_WtoLNu_4Jets_1J", &ParamDict::XSec_WtoLNu_4Jets_1J},
                {"XSec_WtoLNu_4Jets_2J", &ParamDict::XSec_WtoLNu_4Jets_2J},
                {"XSec_WtoLNu_4Jets_3J", &ParamDict::XSec_WtoLNu_4Jets_3J},
                // HH
                {"XSec_GluGlutoHHto4B", &ParamDict::XSec_GluGlutoHHto4B}
            };

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
        {"Wto2Q_2Jets_PTQQ_100to200_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_1J},
        {"Wto2Q_2Jets_PTQQ_100to200_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_100to200_2J},
        {"Wto2Q_2Jets_PTQQ_200to400_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_1J},
        {"Wto2Q_2Jets_PTQQ_200to400_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_200to400_2J},
        {"Wto2Q_2Jets_PTQQ_400to600_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_1J},
        {"Wto2Q_2Jets_PTQQ_400to600_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_400to600_2J},
        {"Wto2Q_2Jets_PTQQ_600_1J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_1J},
        {"Wto2Q_2Jets_PTQQ_600_2J", &ParamDict::XSec_Wto2Q_2Jets_PTQQ_600_2J},
        // Z+jets PTQQ samples
        {"Zto2Q_2Jets_PTQQ_100to200_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_1J},
        {"Zto2Q_2Jets_PTQQ_100to200_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_100to200_2J},
        {"Zto2Q_2Jets_PTQQ_200to400_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_1J},
        {"Zto2Q_2Jets_PTQQ_200to400_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_200to400_2J},
        {"Zto2Q_2Jets_PTQQ_400to600_1J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_1J},
        {"Zto2Q_2Jets_PTQQ_400to600_2J", &ParamDict::XSec_Zto2Q_2Jets_PTQQ_400to600_2J},
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
        {"GluGlutoHHto4B", &ParamDict::XSec_GluGlutoHHto4B}
    };

    auto it = xsecMap.find(data_type);
    if (it != xsecMap.end()) {
        return param_dict->*(it->second);
    }
    throw std::invalid_argument("Invalid data_type: " + data_type);
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
  double dR_LFJ = get_dR(fj_eta, fj_phi, lep_eta, lep_phi);
  if (dR_LFJ < 1.5)
    return false;

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

// JEC and JER
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

JetCorrectionResult
applyJER(Float_t jet_pt, Float_t jet_eta, Float_t jet_phi, Float_t jet_massSD,
         const JME::JetResolution &resolution_pt,
         const JME::JetResolutionScaleFactor &resolution_pt_sf, Float_t rho,
         const Float_t *genJetPts, const Float_t *genJetEtas,
         const Float_t *genJetPhis, int nGenJets) {

  JME::JetParameters JerPARAM = {
    {JME::Binning::JetPt, jet_pt},
    {JME::Binning::JetEta, jet_eta},
    {JME::Binning::Rho, rho}
  };
  
  JME::JetParameters JerSFPARAM;
  JerSFPARAM.set(JME::Binning::JetPt, jet_pt);
  JerSFPARAM.set(JME::Binning::JetEta, jet_eta);
  JerSFPARAM.set(JME::Binning::Rho, rho);

  Float_t resolution = resolution_pt.getResolution(JerPARAM);
  Float_t resolution_sf = resolution_pt_sf.getScaleFactor(JerSFPARAM); 


  const Float_t MAX_DELTA_R = 0.2;
  Float_t smearFactor = 1.0;
  bool gen_matched = false;

  for (int nGJ = 0; nGJ < nGenJets; nGJ++) {
    Float_t delta_R = get_dR(jet_eta, jet_phi, genJetEtas[nGJ], genJetPhis[nGJ]);
    Float_t pt_diff_ratio = fabs(jet_pt - genJetPts[nGJ]) / jet_pt;
    Float_t resolution_threshold = 3 * resolution;

    if (!gen_matched && delta_R < MAX_DELTA_R &&
        pt_diff_ratio < resolution_threshold) {
      smearFactor = 1.0 + (resolution_sf - 1.0) * (jet_pt - genJetPts[nGJ]) / jet_pt;
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

void histo_mc(
    const std::string &year,        // 2022, 2023
    const std::string &data_type,   // Example: TTtoLNu2Q, QCD_HT100to200
    const std::string &channel,     // Muon, EGamma, leptonic
    const std::string &sample_path, // path to the root file
    const std::string &output_path, // path to the output root file
    const std::string &pu_path,     // path to the pileup reweighting file
    const std::string &param_path,  // path to the parameters file
    const std::string &jec_path,    // path to the AK8 JEC txt file
    const std::string &jer_path,    // path to the AK8 JER txt file
    const std::string &jer_path_sf  // path to the AK8 JER SF txt file
) {
  gSystem->Load("libFWCoreFWLite.so");

  std::string channel_lower = to_lower(channel);
  std::cout << "Data Type: " << data_type << std::endl;
  std::cout << "Channel: " << channel << std::endl;
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
  double xsec = getXSec(&param_dict, data_type);

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

  // probe jet kinematics
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
  // probe jet ParticleNet scores
  TH1D *_ProbeJet_all_PNet_QCD =
      new TH1D("ProbeJet_all_PNet_QCD", "ProbeJet_all_PNet_QCD", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_QCD0HF = new TH1D(
      "ProbeJet_all_PNet_QCD0HF", "ProbeJet_all_PNet_QCD0HF", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_QCD1HF = new TH1D(
      "ProbeJet_all_PNet_QCD1HF", "ProbeJet_all_PNet_QCD1HF", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_QCD2HF = new TH1D(
      "ProbeJet_all_PNet_QCD2HF", "ProbeJet_all_PNet_QCD2HF", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XbbVsQCD = new TH1D(
      "ProbeJet_all_PNet_XbbVsQCD", "ProbeJet_all_PNet_XbbVsQCD", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XccVsQCD = new TH1D(
      "ProbeJet_all_PNet_XccVsQCD", "ProbeJet_all_PNet_XccVsQCD", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XggVsQCD = new TH1D(
      "ProbeJet_all_PNet_XggVsQCD", "ProbeJet_all_PNet_XggVsQCD", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNet_XqqVsQCD = new TH1D(
      "ProbeJet_all_PNet_XqqVsQCD", "ProbeJet_all_PNet_XqqVsQCD", N_TXbb, 0, 1.0);
  // probe jet ParticleNetLegacy scores
  TH1D *_ProbeJet_all_PNetLegacy_Xbb =
      new TH1D("ProbeJet_all_PNetLegacy_Xbb", "ProbeJet_all_PNetLegacy_Xbb",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_XbbVsQCD =
      new TH1D("ProbeJet_all_PNetLegacy_XbbVsQCD",
               "ProbeJet_all_PNetLegacy_XbbVsQCD", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_Xcc =
      new TH1D("ProbeJet_all_PNetLegacy_Xcc", "ProbeJet_all_PNetLegacy_Xcc",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_Xqq =
      new TH1D("ProbeJet_all_PNetLegacy_Xqq", "ProbeJet_all_PNetLegacy_Xqq",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCD =
      new TH1D("ProbeJet_all_PNetLegacy_QCD", "ProbeJet_all_PNetLegacy_QCD",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCDb =
      new TH1D("ProbeJet_all_PNetLegacy_QCDb", "ProbeJet_all_PNetLegacy_QCDb",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCDbb =
      new TH1D("ProbeJet_all_PNetLegacy_QCDbb", "ProbeJet_all_PNetLegacy_QCDbb",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_PNetLegacy_QCDothers =
      new TH1D("ProbeJet_all_PNetLegacy_QCDothers",
               "ProbeJet_all_PNetLegacy_QCDothers", N_TXbb, 0, 1.0);
  // probe jet GloParT scores
  TH1D *_ProbeJet_all_GloParT_QCD0HF =
      new TH1D("ProbeJet_all_GloParT_QCD0HF", "ProbeJet_all_GloParT_QCD0HF",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_QCD1HF =
      new TH1D("ProbeJet_all_GloParT_QCD1HF", "ProbeJet_all_GloParT_QCD1HF",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_QCD2HF =
      new TH1D("ProbeJet_all_GloParT_QCD2HF", "ProbeJet_all_GloParT_QCD2HF",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_Xbb = new TH1D(
      "ProbeJet_all_GloParT_Xbb", "ProbeJet_all_GloParT_Xbb", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_Xcc = new TH1D(
      "ProbeJet_all_GloParT_Xcc", "ProbeJet_all_GloParT_Xcc", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_Xqq = new TH1D(
      "ProbeJet_all_GloParT_Xqq", "ProbeJet_all_GloParT_Xqq", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_all_GloParT_XbbVsQCD =
      new TH1D("ProbeJet_all_GloParT_XbbVsQCD", "ProbeJet_all_GloParT_XbbVsQCD",
               N_TXbb, 0, 1.0);

  // passing prob jet kinematics
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
  // passing prob jet ParticleNet scores
  TH1D *_ProbeJet_pass_PNet_QCD =
      new TH1D("ProbeJet_pass_PNet_QCD", "ProbeJet_pass_PNet_QCD", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_QCD0HF = new TH1D(
      "ProbeJet_pass_PNet_QCD0HF", "ProbeJet_pass_PNet_QCD0HF", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_QCD1HF = new TH1D(
      "ProbeJet_pass_PNet_QCD1HF", "ProbeJet_pass_PNet_QCD1HF", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_QCD2HF = new TH1D(
      "ProbeJet_pass_PNet_QCD2HF", "ProbeJet_pass_PNet_QCD2HF", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XbbVsQCD =
      new TH1D("ProbeJet_pass_PNet_XbbVsQCD", "ProbeJet_pass_PNet_XbbVsQCD",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XccVsQCD =
      new TH1D("ProbeJet_pass_PNet_XccVsQCD", "ProbeJet_pass_PNet_XccVsQCD",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XggVsQCD =
      new TH1D("ProbeJet_pass_PNet_XggVsQCD", "ProbeJet_pass_PNet_XggVsQCD",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNet_XqqVsQCD =
      new TH1D("ProbeJet_pass_PNet_XqqVsQCD", "ProbeJet_pass_PNet_XqqVsQCD",
               N_TXbb, 0, 1.0);
  // passing prob jet ParticleNetLegacy scores
  TH1D *_ProbeJet_pass_PNetLegacy_Xbb =
      new TH1D("ProbeJet_pass_PNetLegacy_Xbb", "ProbeJet_pass_PNetLegacy_Xbb",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_XbbVsQCD =
      new TH1D("ProbeJet_pass_PNetLegacy_XbbVsQCD",
               "ProbeJet_pass_PNetLegacy_XbbVsQCD", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_Xcc =
      new TH1D("ProbeJet_pass_PNetLegacy_Xcc", "ProbeJet_pass_PNetLegacy_Xcc",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_Xqq =
      new TH1D("ProbeJet_pass_PNetLegacy_Xqq", "ProbeJet_pass_PNetLegacy_Xqq",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCD =
      new TH1D("ProbeJet_pass_PNetLegacy_QCD", "ProbeJet_pass_PNetLegacy_QCD",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCDb =
      new TH1D("ProbeJet_pass_PNetLegacy_QCDb", "ProbeJet_pass_PNetLegacy_QCDb",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCDbb =
      new TH1D("ProbeJet_pass_PNetLegacy_QCDbb",
               "ProbeJet_pass_PNetLegacy_QCDbb", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_PNetLegacy_QCDothers =
      new TH1D("ProbeJet_pass_PNetLegacy_QCDothers",
               "ProbeJet_pass_PNetLegacy_QCDothers", N_TXbb, 0, 1.0);
  // passing prob jet GloParT scores
  TH1D *_ProbeJet_pass_GloParT_QCD0HF =
      new TH1D("ProbeJet_pass_GloParT_QCD0HF", "ProbeJet_pass_GloParT_QCD0HF",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_QCD1HF =
      new TH1D("ProbeJet_pass_GloParT_QCD1HF", "ProbeJet_pass_GloParT_QCD1HF",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_QCD2HF =
      new TH1D("ProbeJet_pass_GloParT_QCD2HF", "ProbeJet_pass_GloParT_QCD2HF",
               N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_Xbb = new TH1D(
      "ProbeJet_pass_GloParT_Xbb", "ProbeJet_pass_GloParT_Xbb", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_Xcc = new TH1D(
      "ProbeJet_pass_GloParT_Xcc", "ProbeJet_pass_GloParT_Xcc", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_Xqq = new TH1D(
      "ProbeJet_pass_GloParT_Xqq", "ProbeJet_pass_GloParT_Xqq", N_TXbb, 0, 1.0);
  TH1D *_ProbeJet_pass_GloParT_XbbVsQCD =
      new TH1D("ProbeJet_pass_GloParT_XbbVsQCD",
               "ProbeJet_pass_GloParT_XbbVsQCD", N_TXbb, 0, 1.0);

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

  // gen-level jets
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

    JetCorrectionResult jec3 =
        applyJEC(FatJet3_pt, FatJet3_eta, FatJet3_phi, FatJet3_rawFactor,
                 FatJet3_MassSD, corrector);
    FatJet3_pt = jec3.corrected_pt;
    FatJet3_MassSD = jec3.corrected_massSD;

    JetCorrectionResult jer3 =
        applyJER(FatJet3_pt, FatJet3_eta, FatJet3_phi, FatJet3_MassSD,
                 resolution_pt, resolution_pt_sf, rho, GenJetAK8_pt,
                 GenJetAK8_eta, GenJetAK8_phi, nGenJetAK8);
    FatJet3_pt = jer3.corrected_pt;
    FatJet3_MassSD = jer3.corrected_massSD;

    // Tag and Probe
    if (channel_lower == "jetmet" or channel_lower == "qcd") {
      // Tag jet requirements
      if (FatJet1_pt < 300 || fabs(FatJet1_eta) > 2.5 ||
          FatJet1_MassSD < 80) {
        continue;
      }
      if (FatJet2_pt < 160) {
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

    } else {
      // Leptonic channel
      if (FatJet1_pt < 160) {
        continue;
      }
      // EGamma, Muon, Lepton (EGamma + Muon)
      if (year == "2022") {
        if (lep1_Pt < 55 || lep2_Pt > 30) {
          continue;
        }
      } else if (year == "2023") {
        if (lep1_Pt < 50 || lep2_Pt > 30) {
          continue;
        }
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

    } // end if

    // Probe jet matching
    bool probe_pass = false;
    bool matched_to_AK8PFJet230_SoftDropMass40 = checkTriggerMatching(
        NTrigger_Objects, Trigger_Object_pt, Trigger_Object_eta, Trigger_Object_phi,
        Trigger_Object_bit, ProbeJet_eta, ProbeJet_phi, 4, 100);

    if (year == "2022") {
      // $CMSSW_RELEASE_BASE/src/PhysicsTools/NanoAOD/python/triggerObjects_cff.py
      // matched_to_AK8PFJet230_SoftDropMass40 & matched_to_AK8PFJet250
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

    // weight
    weight = (weight / SumGenWeights) * xsec * param_dict.Lumi;
    double PU_weight = PU_Rew[(int)npu];
    if (PU_weight < 20.0) {
      weight = weight * PU_weight;
    }

    // Fill histograms
    // Kinematics
    _ProbeJet_all_pt->Fill(ProbeJet_pt, weight);
    _ProbeJet_all_eta->Fill(ProbeJet_eta, weight);
    _ProbeJet_all_phi->Fill(ProbeJet_phi, weight);
    _ProbeJet_all_eta_phi->Fill(ProbeJet_eta, ProbeJet_phi, weight);
    _ProbeJet_all_Mass->Fill(ProbeJet_Mass, weight);
    _ProbeJet_all_MassSD->Fill(ProbeJet_MassSD, weight);
    _ProbeJet_all_Mass_Pt->Fill(ProbeJet_MassSD, ProbeJet_pt, weight);

    // ParticleNet
    _ProbeJet_all_PNet_QCD->Fill(ProbeJetPNet_QCD, weight);
    _ProbeJet_all_PNet_QCD0HF->Fill(ProbeJetPNet_QCD0HF, weight);
    _ProbeJet_all_PNet_QCD1HF->Fill(ProbeJetPNet_QCD1HF, weight);
    _ProbeJet_all_PNet_QCD2HF->Fill(ProbeJetPNet_QCD2HF, weight);
    _ProbeJet_all_PNet_XbbVsQCD->Fill(ProbeJetPNet_XbbVsQCD, weight);
    _ProbeJet_all_PNet_XccVsQCD->Fill(ProbeJetPNet_XccVsQCD, weight);
    _ProbeJet_all_PNet_XggVsQCD->Fill(ProbeJetPNet_XggVsQCD, weight);
    _ProbeJet_all_PNet_XqqVsQCD->Fill(ProbeJetPNet_XqqVsQCD, weight);
    // ParticleNetLegacy
    Float_t ProbeJetPNetLegacy_XbbVsQCD =
        ProbeJetPNetLegacy_Xbb /
        (ProbeJetPNetLegacy_Xbb + ProbeJetPNetLegacy_QCD + 1e-12);
    _ProbeJet_all_PNetLegacy_Xbb->Fill(ProbeJetPNetLegacy_Xbb, weight);
    _ProbeJet_all_PNetLegacy_XbbVsQCD->Fill(ProbeJetPNetLegacy_XbbVsQCD,
                                            weight);
    _ProbeJet_all_PNetLegacy_Xcc->Fill(ProbeJetPNetLegacy_Xcc, weight);
    _ProbeJet_all_PNetLegacy_Xqq->Fill(ProbeJetPNetLegacy_Xqq, weight);
    _ProbeJet_all_PNetLegacy_QCD->Fill(ProbeJetPNetLegacy_QCD, weight);
    _ProbeJet_all_PNetLegacy_QCDb->Fill(ProbeJetPNetLegacy_QCDb, weight);
    _ProbeJet_all_PNetLegacy_QCDbb->Fill(ProbeJetPNetLegacy_QCDbb, weight);
    _ProbeJet_all_PNetLegacy_QCDothers->Fill(ProbeJetPNetLegacy_QCDothers,
                                             weight);
    // GloParT
    _ProbeJet_all_GloParT_QCD0HF->Fill(ProbeJetGloParT_QCD0HF, weight);
    _ProbeJet_all_GloParT_QCD1HF->Fill(ProbeJetGloParT_QCD1HF, weight);
    _ProbeJet_all_GloParT_QCD2HF->Fill(ProbeJetGloParT_QCD2HF, weight);
    _ProbeJet_all_GloParT_Xbb->Fill(ProbeJetGloParT_Xbb, weight);
    _ProbeJet_all_GloParT_Xcc->Fill(ProbeJetGloParT_Xcc, weight);
    _ProbeJet_all_GloParT_Xqq->Fill(ProbeJetGloParT_Xqq, weight);
    _ProbeJet_all_GloParT_XbbVsQCD->Fill(ProbeJetGloParT_XbbVsQCD, weight);

    if (probe_pass) {
      // Kinematics
      _ProbeJet_pass_pt->Fill(ProbeJet_pt, weight);
      _ProbeJet_pass_eta->Fill(ProbeJet_eta, weight);
      _ProbeJet_pass_phi->Fill(ProbeJet_phi, weight);
      _ProbeJet_pass_eta_phi->Fill(ProbeJet_eta, ProbeJet_phi, weight);
      _ProbeJet_pass_Mass->Fill(ProbeJet_Mass, weight);
      _ProbeJet_pass_MassSD->Fill(ProbeJet_MassSD, weight);
      _ProbeJet_pass_Mass_Pt->Fill(ProbeJet_MassSD, ProbeJet_pt, weight);
      // ParticleNet
      _ProbeJet_pass_PNet_QCD->Fill(ProbeJetPNet_QCD, weight);
      _ProbeJet_pass_PNet_QCD0HF->Fill(ProbeJetPNet_QCD0HF, weight);
      _ProbeJet_pass_PNet_QCD1HF->Fill(ProbeJetPNet_QCD1HF, weight);
      _ProbeJet_pass_PNet_QCD2HF->Fill(ProbeJetPNet_QCD2HF, weight);
      _ProbeJet_pass_PNet_XbbVsQCD->Fill(ProbeJetPNet_XbbVsQCD, weight);
      _ProbeJet_pass_PNet_XccVsQCD->Fill(ProbeJetPNet_XccVsQCD, weight);
      _ProbeJet_pass_PNet_XggVsQCD->Fill(ProbeJetPNet_XggVsQCD, weight);
      _ProbeJet_pass_PNet_XqqVsQCD->Fill(ProbeJetPNet_XqqVsQCD, weight);
      // ParticleNetLegacy
      _ProbeJet_pass_PNetLegacy_Xbb->Fill(ProbeJetPNetLegacy_Xbb, weight);
      _ProbeJet_pass_PNetLegacy_XbbVsQCD->Fill(ProbeJetPNetLegacy_XbbVsQCD,
                                               weight);
      _ProbeJet_pass_PNetLegacy_Xcc->Fill(ProbeJetPNetLegacy_Xcc, weight);
      _ProbeJet_pass_PNetLegacy_Xqq->Fill(ProbeJetPNetLegacy_Xqq, weight);
      _ProbeJet_pass_PNetLegacy_QCD->Fill(ProbeJetPNetLegacy_QCD, weight);
      _ProbeJet_pass_PNetLegacy_QCDb->Fill(ProbeJetPNetLegacy_QCDb, weight);
      _ProbeJet_pass_PNetLegacy_QCDbb->Fill(ProbeJetPNetLegacy_QCDbb, weight);
      _ProbeJet_pass_PNetLegacy_QCDothers->Fill(ProbeJetPNetLegacy_QCDothers,
                                                weight);
      // GloParT
      _ProbeJet_pass_GloParT_QCD0HF->Fill(ProbeJetGloParT_QCD0HF, weight);
      _ProbeJet_pass_GloParT_QCD1HF->Fill(ProbeJetGloParT_QCD1HF, weight);
      _ProbeJet_pass_GloParT_QCD2HF->Fill(ProbeJetGloParT_QCD2HF, weight);
      _ProbeJet_pass_GloParT_Xbb->Fill(ProbeJetGloParT_Xbb, weight);
      _ProbeJet_pass_GloParT_Xcc->Fill(ProbeJetGloParT_Xcc, weight);
      _ProbeJet_pass_GloParT_Xqq->Fill(ProbeJetGloParT_Xqq, weight);
      _ProbeJet_pass_GloParT_XbbVsQCD->Fill(ProbeJetGloParT_XbbVsQCD, weight);
    }

  } // end event loop

  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
