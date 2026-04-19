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
#include "TChain.h"
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

  // VH samples
  double XSec_WplusH_Hto2B_Wto2Q_M_125;
  double XSec_WplusH_Hto2B_WtoLNu_M_125;
  double XSec_WminusH_Hto2B_Wto2Q_M_125;
  double XSec_WminusH_Hto2B_WtoLNu_M_125;
  double XSec_ZH_Hto2B_Zto2L_M_125;
  double XSec_ZH_Hto2B_Zto2Nu_M_125;
  double XSec_ZH_Hto2B_Zto2Q_M_125;
  double XSec_ZH_Hto2C_Zto2Q_M_125;

  // ggZH samples
  double XSec_ggZH_Hto2B_Zto2L_M_125;
  double XSec_ggZH_Hto2B_Zto2Nu_M_125;
  double XSec_ggZH_Hto2B_Zto2Q_M_125;
  double XSec_ggZH_Hto2C_Zto2Q_M_125;

  // singleT samples
  double XSec_TBbarQ_t_channel_4FS;
  double XSec_TbarBQ_t_channel_4FS;
  double XSec_TWminusto4Q;
  double XSec_TWminustoLNu2Q;
  double XSec_TbarWplusto4Q;
  double XSec_TbarWplustoLNu2Q;

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
          // VH
          {"XSec_WplusH_Hto2B_Wto2Q_M_125",
           &ParamDict::XSec_WplusH_Hto2B_Wto2Q_M_125},
          {"XSec_WplusH_Hto2B_WtoLNu_M_125",
           &ParamDict::XSec_WplusH_Hto2B_WtoLNu_M_125},
          {"XSec_WminusH_Hto2B_Wto2Q_M_125",
           &ParamDict::XSec_WminusH_Hto2B_Wto2Q_M_125},
          {"XSec_WminusH_Hto2B_WtoLNu_M_125",
           &ParamDict::XSec_WminusH_Hto2B_WtoLNu_M_125},
          {"XSec_ZH_Hto2B_Zto2L_M_125", &ParamDict::XSec_ZH_Hto2B_Zto2L_M_125},
          {"XSec_ZH_Hto2B_Zto2Nu_M_125", &ParamDict::XSec_ZH_Hto2B_Zto2Nu_M_125},
          {"XSec_ZH_Hto2B_Zto2Q_M_125", &ParamDict::XSec_ZH_Hto2B_Zto2Q_M_125},
          {"XSec_ZH_Hto2C_Zto2Q_M_125", &ParamDict::XSec_ZH_Hto2C_Zto2Q_M_125},
          // ggZH
          {"XSec_ggZH_Hto2B_Zto2L_M_125", &ParamDict::XSec_ggZH_Hto2B_Zto2L_M_125},
          {"XSec_ggZH_Hto2B_Zto2Nu_M_125", &ParamDict::XSec_ggZH_Hto2B_Zto2Nu_M_125},
          {"XSec_ggZH_Hto2B_Zto2Q_M_125", &ParamDict::XSec_ggZH_Hto2B_Zto2Q_M_125},
          {"XSec_ggZH_Hto2C_Zto2Q_M_125", &ParamDict::XSec_ggZH_Hto2C_Zto2Q_M_125},
          // singleT
          {"XSec_TBbarQ_t_channel_4FS", &ParamDict::XSec_TBbarQ_t_channel_4FS},
          {"XSec_TbarBQ_t_channel_4FS", &ParamDict::XSec_TbarBQ_t_channel_4FS},
          {"XSec_TWminusto4Q", &ParamDict::XSec_TWminusto4Q},
          {"XSec_TWminustoLNu2Q", &ParamDict::XSec_TWminustoLNu2Q},
          {"XSec_TbarWplusto4Q", &ParamDict::XSec_TbarWplusto4Q},
          {"XSec_TbarWplustoLNu2Q", &ParamDict::XSec_TbarWplustoLNu2Q},
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
      // VH
      {"WplusH_Hto2B_Wto2Q_M_125",
        &ParamDict::XSec_WplusH_Hto2B_Wto2Q_M_125},
      {"WplusH_Hto2B_WtoLNu_M_125",
        &ParamDict::XSec_WplusH_Hto2B_WtoLNu_M_125},
      {"WminusH_Hto2B_Wto2Q_M_125",
        &ParamDict::XSec_WminusH_Hto2B_Wto2Q_M_125},
      {"WminusH_Hto2B_WtoLNu_M_125",
        &ParamDict::XSec_WminusH_Hto2B_WtoLNu_M_125},
      {"ZH_Hto2B_Zto2L_M_125", &ParamDict::XSec_ZH_Hto2B_Zto2L_M_125},
      {"ZH_Hto2B_Zto2Nu_M_125", &ParamDict::XSec_ZH_Hto2B_Zto2Nu_M_125},
      {"ZH_Hto2B_Zto2Q_M_125", &ParamDict::XSec_ZH_Hto2B_Zto2Q_M_125},
      {"ZH_Hto2C_Zto2Q_M_125", &ParamDict::XSec_ZH_Hto2C_Zto2Q_M_125},
      // ggZH
      {"ggZH_Hto2B_Zto2L_M_125", &ParamDict::XSec_ggZH_Hto2B_Zto2L_M_125},
      {"ggZH_Hto2B_Zto2Nu_M_125", &ParamDict::XSec_ggZH_Hto2B_Zto2Nu_M_125},
      {"ggZH_Hto2B_Zto2Q_M_125", &ParamDict::XSec_ggZH_Hto2B_Zto2Q_M_125},
      {"ggZH_Hto2C_Zto2Q_M_125", &ParamDict::XSec_ggZH_Hto2C_Zto2Q_M_125},
      // singleT
      {"TBbarQ_t_channel_4FS", &ParamDict::XSec_TBbarQ_t_channel_4FS},
      {"TbarBQ_t_channel_4FS", &ParamDict::XSec_TbarBQ_t_channel_4FS},
      {"TWminusto4Q", &ParamDict::XSec_TWminusto4Q},
      {"TWminustoLNu2Q", &ParamDict::XSec_TWminustoLNu2Q},
      {"TbarWplusto4Q", &ParamDict::XSec_TbarWplusto4Q},
      {"TbarWplustoLNu2Q", &ParamDict::XSec_TbarWplustoLNu2Q},
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

void tree_mc_hl(const std::string &lepton,      // ele, muon
                const std::string &year,      // 2022, 2023
                const std::string &data_type, // Example: TTtoLNu2Q, QCD_HT100to200
                const std::string &sample_path, // path to the root file
                const std::string &output_path, // path to the output root file
                const std::string &pu_path, // path to the pileup reweighting file
                const std::string &param_path, // path to the parameters file
                const std::string &jec_path,   // path to the AK8 JEC txt file
                const std::string &jer_path,   // path to the AK8 JER txt file
                const std::string &jer_path_sf, // path to the AK8 JER SF txt file
                const int begin_file_index=0,
                const int end_file_index=999999999
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


  std::string new_data_type = "QCD_HT_2000toInf";
  double xsec = 0.0;
  if (data_type == "QCD_HT_2000") {
    xsec = getXSec(&param_dict, new_data_type);
  } else {
    xsec = getXSec(&param_dict, data_type);
  }
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

  
  // Outputs
  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  TH1F *h_cutflow = new TH1F("cutflow", "cutflow", 10, 0, 10);
  h_cutflow->GetXaxis()->SetBinLabel(1, "Total Events");
  h_cutflow->GetXaxis()->SetBinLabel(2, "HLT trigger selection");
  h_cutflow->GetXaxis()->SetBinLabel(3, "fatjet0 kinematic selection");   
  h_cutflow->GetXaxis()->SetBinLabel(4, "fatjet1 kinematic selection");   
  h_cutflow->GetXaxis()->SetBinLabel(5, "lepton kinematic selection");
  h_cutflow->GetXaxis()->SetBinLabel(6, "MET selection");
  h_cutflow->GetXaxis()->SetBinLabel(7, "lepton-fatjet dR selection");

  TTree *outputTree = new TTree("tree", "");

  Float_t T_weight;
  Float_t T_fatJet1_pt, T_fatJet1_eta, T_fatJet1_phi, T_fatJet1_mass, T_fatJet1_msd,
      T_fatJet1_CAmass, T_fatJet1_CAmsoftdrop, T_fatJet1_CAmass_et,
      T_fatJet1_CAmsoftdrop_et, T_fatJet1_CAmass_fatjet_et,
      T_fatJet1_GloParT_massVis, T_fatJet1_GloParT_massRes, T_fatJet1_Tau3OverTau2;
  Float_t T_fatJet1_ParticleNetLegacy_Xbb, T_fatJet1_ParticleNetLegacy_XbbVsQCD;
  Float_t T_fatJet1_GloParT_Xbb, T_fatJet1_GloParT_XbbVsQCD, 
      T_fatJet1_GloParT_Xtauhtaue;
  Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_phi, T_fatJet2_mass, T_fatJet2_msd,
      T_fatJet2_CAmass, T_fatJet2_CAmsoftdrop, T_fatJet2_CAmass_et,
      T_fatJet2_CAmsoftdrop_et, T_fatJet2_CAmass_fatjet_et,
      T_fatJet2_GloParT_massVis, T_fatJet2_GloParT_massRes, T_fatJet2_Tau3OverTau2;
  Float_t T_MET, T_ele1_pt, T_ele1_eta, T_ele1_phi, T_ele2_pt, T_ele2_eta, T_ele2_phi,
      T_muon1_pt, T_muon1_eta, T_muon1_phi, T_muon2_pt, T_muon2_eta, T_muon2_phi,
      T_dR_LFJ, T_dR_JFJ, T_dR_JL;
  Float_t T_jet1_pt, T_jet1_eta, T_jet1_phi, T_jet1_mass, T_jet1_btag;
  Float_t T_jet2_pt, T_jet2_eta, T_jet2_phi, T_jet2_mass, T_jet2_btag;



  outputTree->Branch("weight", &T_weight, "weight/F");
  outputTree->Branch("fatJet1_pt", &T_fatJet1_pt, "fatJet1_pt/F");
  outputTree->Branch("fatJet1_eta", &T_fatJet1_eta, "fatJet1_eta/F");
  outputTree->Branch("fatJet1_phi", &T_fatJet1_phi, "fatJet1_phi/F");
  outputTree->Branch("fatJet1_mass", &T_fatJet1_mass, "fatJet1_mass/F");
  outputTree->Branch("fatJet1_msd", &T_fatJet1_msd, "fatJet1_msd/F");
  outputTree->Branch("fatJet1_CAmass", &T_fatJet1_CAmass, "fatJet1_CAmass/F");
  outputTree->Branch("fatJet1_CAmsoftdrop", &T_fatJet1_CAmsoftdrop,
                     "fatJet1_CAmsoftdrop/F");
  outputTree->Branch("fatJet1_CAmass_et", &T_fatJet1_CAmass_et,
                     "fatJet1_CAmass_et/F");
  outputTree->Branch("fatJet1_CAmsoftdrop_et", &T_fatJet1_CAmsoftdrop_et,
                     "fatJet1_CAmsoftdrop_et/F");
  outputTree->Branch("fatJet1_CAmass_fatjet_et", &T_fatJet1_CAmass_fatjet_et,
                     "fatJet1_CAmass_fatjet_et/F");
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

  outputTree->Branch("fatJet1_GloParT_Xbb", &T_fatJet1_GloParT_Xbb,
                     "fatJet1_GloParT_Xbb/F");
  outputTree->Branch("fatJet1_GloParT_XbbVsQCD", &T_fatJet1_GloParT_XbbVsQCD,
                     "fatJet1_GloParT_XbbVsQCD/F");
  outputTree->Branch("fatJet1_GloParT_Xtauhtaue",
                     &T_fatJet1_GloParT_Xtauhtaue,
                     "fatJet1_GloParT_Xtauhtaue/F");

  outputTree->Branch("fatJet2_pt", &T_fatJet2_pt, "fatJet2_pt/F");
  outputTree->Branch("fatJet2_eta", &T_fatJet2_eta, "fatJet2_eta/F");
  outputTree->Branch("fatJet2_phi", &T_fatJet2_phi, "fatJet2_phi/F");
  outputTree->Branch("fatJet2_mass", &T_fatJet2_mass, "fatJet2_mass/F");
  outputTree->Branch("fatJet2_msd", &T_fatJet2_msd, "fatJet2_msd/F");
  outputTree->Branch("fatJet2_CAmass", &T_fatJet2_CAmass, "fatJet2_CAmass/F");
  outputTree->Branch("fatJet2_CAmsoftdrop", &T_fatJet2_CAmsoftdrop,
                     "fatJet2_CAmsoftdrop/F");
  outputTree->Branch("fatJet2_CAmass_et", &T_fatJet2_CAmass_et,
                     "fatJet2_CAmass_et/F");
  outputTree->Branch("fatJet2_CAmsoftdrop_et", &T_fatJet2_CAmsoftdrop_et,
                     "fatJet2_CAmsoftdrop_et/F");
  outputTree->Branch("fatJet2_CAmass_fatjet_et", &T_fatJet2_CAmass_fatjet_et,
                     "fatJet2_CAmass_fatjet_et/F");
  outputTree->Branch("fatJet2_GloParT_massVis", &T_fatJet2_GloParT_massVis,
                     "fatJet2_GloParT_massVis/F");
  outputTree->Branch("fatJet2_GloParT_massRes", &T_fatJet2_GloParT_massRes,
                     "fatJet2_GloParT_massRes/F");
  outputTree->Branch("fatJet2_Tau3OverTau2", &T_fatJet2_Tau3OverTau2,
                     "fatJet2_Tau3OverTau2/F");
  outputTree->Branch("MET", &T_MET, "MET/F");
  outputTree->Branch("ele1_pt", &T_ele1_pt, "ele1_pt/F");
  outputTree->Branch("ele1_eta", &T_ele1_eta, "ele1_eta/F");
  outputTree->Branch("ele1_phi", &T_ele1_phi, "ele1_phi/F");
  outputTree->Branch("ele2_pt", &T_ele2_pt, "ele2_pt/F");
  outputTree->Branch("ele2_eta", &T_ele2_eta, "ele2_eta/F");
  outputTree->Branch("ele2_phi", &T_ele2_phi, "ele2_phi/F");
  outputTree->Branch("muon1_pt", &T_muon1_pt, "muon1_pt/F");
  outputTree->Branch("muon1_eta", &T_muon1_eta, "muon1_eta/F");
  outputTree->Branch("muon1_phi", &T_muon1_phi, "muon1_phi/F");
  outputTree->Branch("muon2_pt", &T_muon2_pt, "muon2_pt/F");
  outputTree->Branch("muon2_eta", &T_muon2_eta, "muon2_eta/F");
  outputTree->Branch("muon2_phi", &T_muon2_phi, "muon2_phi/F");
  outputTree->Branch("dR_LFJ", &T_dR_LFJ, "dR_LFJ/F");
  outputTree->Branch("dR_JFJ", &T_dR_JFJ, "dR_JFJ/F");
  outputTree->Branch("dR_JL", &T_dR_JL, "dR_JL/F");

  outputTree->Branch("jet1_pt", &T_jet1_pt, "jet1_pt/F");
  outputTree->Branch("jet1_eta", &T_jet1_eta, "jet1_eta/F");
  outputTree->Branch("jet1_phi", &T_jet1_phi, "jet1_phi/F");
  outputTree->Branch("jet1_mass", &T_jet1_mass, "jet1_mass/F");
  outputTree->Branch("jet1_btag", &T_jet1_btag, "jet1_btag/F");
  outputTree->Branch("jet2_pt", &T_jet2_pt, "jet2_pt/F");
  outputTree->Branch("jet2_eta", &T_jet2_eta, "jet2_eta/F");
  outputTree->Branch("jet2_phi", &T_jet2_phi, "jet2_phi/F");
  outputTree->Branch("jet2_mass", &T_jet2_mass, "jet2_mass/F");
  outputTree->Branch("jet2_btag", &T_jet2_btag, "jet2_btag/F");

 #include <filesystem>
  namespace fs = std::filesystem;

  TChain *InputTree = new TChain("Events");
  std::string dir_path = sample_path.c_str();
  
  // scan dir to get root files and add to TChain only selected ones
  std::vector<fs::path> files;
  for (const auto &e : fs::directory_iterator(dir_path))
  if (e.is_regular_file()) files.push_back(e.path());
  std::cout << "Total files: " << files.size() << std::endl;
  
  std::cout << "Adding files from " << begin_file_index << " to " << end_file_index << std::endl;
  int i=0;
  for (const auto &p : files){
    if (i < begin_file_index) {
        i++;
        continue;
    }
    else if (i > end_file_index) {
        break;
    } 
    std::cout << "Adding file: " << p << std::endl;
    InputTree->Add(p.c_str());
    i++;
  } 
  std::cout << "Total files added to TChain: " << i - begin_file_index << std::endl;
  
  
  Double_t weight;
  UInt_t run;
  UInt_t lumi;
  Int_t npu;
  Int_t isVBFtag;
  Float_t rho;

  Long64_t HLT_Ele30_WPTight_Gsf;
  Long64_t HLT_IsoMu24;


  Float_t METPt;

  Double_t ElectronPt0;
  Double_t ElectronEta0;
  Double_t ElectronPhi0;
  Double_t ElectronPt1;
  Double_t ElectronEta1;
  Double_t ElectronPhi1;

  Double_t MuonPt0;
  Double_t MuonEta0;
  Double_t MuonPhi0;
  Double_t MuonPt1;
  Double_t MuonEta1;
  Double_t MuonPhi1;

  Double_t ak4JetPt[4];
  Double_t ak4JetEta[4];
  Double_t ak4JetPhi[4];
  Double_t ak4JetMass[4];
  Double_t ak4JetbtagPNetB[4];
  
  // FatJet 0 kinematics
  Double_t ak8FatJetPt0;
  Double_t ak8FatJetEta0;
  Double_t ak8FatJetPhi0;
  Double_t ak8FatJetMass0;
  Double_t ak8FatJetMsd0;
  Double_t ak8FatJetCAmass0;
  Double_t ak8FatJetCAmsoftdrop0;
  Double_t ak8FatJetCAmass_et0;
  Double_t ak8FatJetCAmsoftdrop_et0;
  Double_t ak8FatJetCAmass_fatjet_et0;
  Double_t ak8FatJetrawFactor0;
  Double_t ak8FatJetTau3OverTau20;
  // Fatjet 0 ParticleNetLegacy scores
  Double_t ak8FatJetPNetXbbLegacy0;
  Double_t ak8FatJetPNetQCDLegacy0;
  Double_t ak8FatJetPNetQCDbLegacy0;
  Double_t ak8FatJetPNetQCDbbLegacy0;
  Double_t ak8FatJetPNetQCDothersLegacy0;
  Double_t ak8FatJetPNetXbbvsQCDLegacy0;
  // Fatjet 0 GloParT scores
  Double_t ak8FatJetParTQCD0HF0;
  Double_t ak8FatJetParTQCD1HF0;
  Double_t ak8FatJetParTQCD2HF0;
  Double_t ak8FatJetParTXbb0;
  Double_t ak8FatJetParTXcc0;
  Double_t ak8FatJetParTXqq0;
  Double_t ak8FatJetParTXbbvsQCD0;
  Double_t ak8FatJetCAglobalParT_massResApplied0;
  Double_t ak8FatJetCAglobalParT_massVisApplied0;
  Double_t ak8FatJetParTXtauhtaue0;

  // FatJet 1 kinematics
  Double_t ak8FatJetPt1;
  Double_t ak8FatJetEta1;
  Double_t ak8FatJetPhi1;
  Double_t ak8FatJetMass1;
  Double_t ak8FatJetMsd1;
  Double_t ak8FatJetCAmass1;
  Double_t ak8FatJetCAmsoftdrop1;
  Double_t ak8FatJetCAmass_et1;
  Double_t ak8FatJetCAmsoftdrop_et1;
  Double_t ak8FatJetCAmass_fatjet_et1;
  Double_t ak8FatJetrawFactor1;
  Double_t ak8FatJetTau3OverTau21;

  // Fatjet 1 GloParT scores
  Double_t ak8FatJetCAglobalParT_massResApplied1;
  Double_t ak8FatJetCAglobalParT_massVisApplied1;
  // // gen-level
  // Int_t nGenJet;
  // Float_t GenJet_eta[ARR_SIZE];
  // Float_t GenJet_phi[ARR_SIZE];
  // Float_t GenJet_pt[ARR_SIZE];
  // Int_t nGenJetAK8;
  // Float_t GenJetAK8_eta[ARR_SIZE];
  // Float_t GenJetAK8_phi[ARR_SIZE];
  // Float_t GenJetAK8_pt[ARR_SIZE];

  InputTree->SetBranchAddress("weight", &weight);
  InputTree->SetBranchAddress("run", &run);
  InputTree->SetBranchAddress("luminosityBlock", &lumi);
  InputTree->SetBranchAddress("nPU", &npu);

  // Lepton triggers
  InputTree->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf);
  InputTree->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
  
  //Lepton kinematics
  InputTree->SetBranchAddress("ElectronPt0", &ElectronPt0);
  InputTree->SetBranchAddress("ElectronEta0", &ElectronEta0);
  InputTree->SetBranchAddress("ElectronPhi0", &ElectronPhi0);
  InputTree->SetBranchAddress("ElectronPt1", &ElectronPt1);
  InputTree->SetBranchAddress("ElectronEta1", &ElectronEta1);
  InputTree->SetBranchAddress("ElectronPhi1", &ElectronPhi1);

  InputTree->SetBranchAddress("MuonPt0", &MuonPt0);
  InputTree->SetBranchAddress("MuonEta0", &MuonEta0);
  InputTree->SetBranchAddress("MuonPhi0", &MuonPhi0);
  InputTree->SetBranchAddress("MuonPt1", &MuonPt1);
  InputTree->SetBranchAddress("MuonEta1", &MuonEta1);
  InputTree->SetBranchAddress("MuonPhi1", &MuonPhi1);

  InputTree->SetBranchAddress("METPt", &METPt);

  for (int i = 0; i < 4; i++) {
    InputTree->SetBranchAddress(Form("ak4JetPt%d", i), &ak4JetPt[i]);
    InputTree->SetBranchAddress(Form("ak4JetEta%d", i), &ak4JetEta[i]);
    InputTree->SetBranchAddress(Form("ak4JetPhi%d", i), &ak4JetPhi[i]);
    InputTree->SetBranchAddress(Form("ak4JetMass%d", i), &ak4JetMass[i]);
    InputTree->SetBranchAddress(Form("ak4JetbtagPNetB%d", i), &ak4JetbtagPNetB[i]);
  } 

  // FatJet 0 kinematics
  InputTree->SetBranchAddress("ak8FatJetPt0", &ak8FatJetPt0);
  InputTree->SetBranchAddress("ak8FatJetEta0", &ak8FatJetEta0);
  InputTree->SetBranchAddress("ak8FatJetPhi0", &ak8FatJetPhi0);
  InputTree->SetBranchAddress("ak8FatJetMass0", &ak8FatJetMass0);  
  InputTree->SetBranchAddress("ak8FatJetMsd0", &ak8FatJetMsd0);  
  InputTree->SetBranchAddress("ak8FatJetCAmass0", &ak8FatJetCAmass0);
  InputTree->SetBranchAddress("ak8FatJetCAmsoftdrop0", &ak8FatJetCAmsoftdrop0);
  InputTree->SetBranchAddress("ak8FatJetCAmass_et0", &ak8FatJetCAmass_et0);
  InputTree->SetBranchAddress("ak8FatJetCAmsoftdrop_et0", &ak8FatJetCAmsoftdrop_et0);
  InputTree->SetBranchAddress("ak8FatJetCAmass_fatjet_et0", &ak8FatJetCAmass_fatjet_et0);
  InputTree->SetBranchAddress("ak8FatJetrawFactor0", &ak8FatJetrawFactor0);
  InputTree->SetBranchAddress("ak8FatJetTau3OverTau20", &ak8FatJetTau3OverTau20);
  
  // FatJet 0 ParticleNetLegacy scores
  InputTree->SetBranchAddress("ak8FatJetPNetXbbLegacy0",
                              &ak8FatJetPNetXbbLegacy0);
  InputTree->SetBranchAddress("ak8FatJetPNetQCDLegacy0",
                              &ak8FatJetPNetQCDLegacy0);
  InputTree->SetBranchAddress("ak8FatJetPNetQCDbLegacy0",
                              &ak8FatJetPNetQCDbLegacy0);
  InputTree->SetBranchAddress("ak8FatJetPNetQCDbbLegacy0",
                              &ak8FatJetPNetQCDbbLegacy0);
  InputTree->SetBranchAddress("ak8FatJetPNetQCDothersLegacy0",
                              &ak8FatJetPNetQCDothersLegacy0);
  InputTree->SetBranchAddress("ak8FatJetPNetXbbvsQCDLegacy0",
                              &ak8FatJetPNetXbbvsQCDLegacy0);
  // FatJet 0 GloParT scores
  InputTree->SetBranchAddress("ak8FatJetParTQCD0HF0",
                              &ak8FatJetParTQCD0HF0);
  InputTree->SetBranchAddress("ak8FatJetParTQCD1HF0",
                              &ak8FatJetParTQCD1HF0);
  InputTree->SetBranchAddress("ak8FatJetParTQCD2HF0",
                              &ak8FatJetParTQCD2HF0);
  InputTree->SetBranchAddress("ak8FatJetParTXbb0", &ak8FatJetParTXbb0);
  InputTree->SetBranchAddress("ak8FatJetParTXcc0", &ak8FatJetParTXcc0);
  InputTree->SetBranchAddress("ak8FatJetParTXqq0", &ak8FatJetParTXqq0);
  InputTree->SetBranchAddress("ak8FatJetParTXbbvsQCD0",
                              &ak8FatJetParTXbbvsQCD0);
  InputTree->SetBranchAddress("ak8FatJetCAglobalParT_massResApplied0",
                              &ak8FatJetCAglobalParT_massResApplied0);
  InputTree->SetBranchAddress("ak8FatJetCAglobalParT_massVisApplied0",
                              &ak8FatJetCAglobalParT_massVisApplied0);
  InputTree->SetBranchAddress("ak8FatJetParTXtauhtaue0", &ak8FatJetParTXtauhtaue0);

  // FatJet 1 kinematics
  InputTree->SetBranchAddress("ak8FatJetPt1", &ak8FatJetPt1);
  InputTree->SetBranchAddress("ak8FatJetEta1", &ak8FatJetEta1);
  InputTree->SetBranchAddress("ak8FatJetPhi1", &ak8FatJetPhi1);
  InputTree->SetBranchAddress("ak8FatJetMass1", &ak8FatJetMass1);  
  InputTree->SetBranchAddress("ak8FatJetCAmass1", &ak8FatJetCAmass1);
  InputTree->SetBranchAddress("ak8FatJetMsd1", &ak8FatJetMsd1);
  InputTree->SetBranchAddress("ak8FatJetCAmsoftdrop1", &ak8FatJetCAmsoftdrop1);
  InputTree->SetBranchAddress("ak8FatJetCAmass_et1", &ak8FatJetCAmass_et1);
  InputTree->SetBranchAddress("ak8FatJetCAmsoftdrop_et1", &ak8FatJetCAmsoftdrop_et1);
  InputTree->SetBranchAddress("ak8FatJetCAmass_fatjet_et1", &ak8FatJetCAmass_fatjet_et1);
  InputTree->SetBranchAddress("ak8FatJetrawFactor1", &ak8FatJetrawFactor1);
  InputTree->SetBranchAddress("ak8FatJetTau3OverTau21", &ak8FatJetTau3OverTau21);

  // FatJet 1 GloParT scores
  InputTree->SetBranchAddress("ak8FatJetCAglobalParT_massResApplied1",
                              &ak8FatJetCAglobalParT_massResApplied1);
  InputTree->SetBranchAddress("ak8FatJetCAglobalParT_massVisApplied1",
                              &ak8FatJetCAglobalParT_massVisApplied1);
  // // gen-level jets
  // InputTree->SetBranchAddress("nGenJet", &nGenJet);
  // InputTree->SetBranchAddress("GenJet_eta", GenJet_eta);
  // InputTree->SetBranchAddress("GenJet_phi", GenJet_phi);
  // InputTree->SetBranchAddress("GenJet_pt", GenJet_pt);
  // InputTree->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
  // InputTree->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta);
  // InputTree->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi);
  // InputTree->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt);


  int Nentries = InputTree->GetEntries();
  double SumGenWeights = 0.0;
  for (int i = 0; i < Nentries; i++) {
    InputTree->GetEntry(i);
    SumGenWeights += weight;
  }
  std::cout << "SumGenWeights (Total Entries): " << SumGenWeights << std::endl;
  
  // Events Loop
  for (int i = 0; i < Nentries; i++) {
    InputTree->GetEntry(i);
    if (i % 100000 == 0) {
      std::cout << "Processing event " << i << " / " << Nentries
                << "\r" << std::flush;
    }
    h_cutflow->Fill(0.5); // Total Events

    // HLT Selection
    if (lepton == "ele" && !HLT_Ele30_WPTight_Gsf) {
      continue;
    }
    if (lepton == "muon" && !HLT_IsoMu24) {
      continue;
    }

    h_cutflow->Fill(1.5); // After HLT trigger selection

    // FatJets correction (JEC and JER)
    JetCorrectionResult jec1 =
        applyJEC(ak8FatJetPt0, ak8FatJetEta0, ak8FatJetPhi0, ak8FatJetrawFactor0,
                 ak8FatJetCAmsoftdrop0, corrector);
    ak8FatJetPt0 = jec1.corrected_pt;
    ak8FatJetCAmsoftdrop0 = jec1.corrected_massSD;

    // JetCorrectionResult jer1 =
    //     applyJER(ak8FatJetPt0, ak8FatJetEta0, ak8FatJetPhi0, ak8FatJetCAmsoftdrop0,
    //              resolution_pt, resolution_pt_sf, rho, GenJetAK8_pt,
    //              GenJetAK8_eta, GenJetAK8_phi, nGenJetAK8);
    // ak8FatJetPt0 = jer1.corrected_pt;
    // ak8FatJetCAmsoftdrop0 = jer1.corrected_massSD;
    JetCorrectionResult jec2 =
        applyJEC(ak8FatJetPt1, ak8FatJetEta1, ak8FatJetPhi1, ak8FatJetrawFactor1,
                 ak8FatJetCAmsoftdrop1, corrector);
    ak8FatJetPt1 = jec2.corrected_pt;
    ak8FatJetCAmsoftdrop1 = jec2.corrected_massSD;

    // JetCorrectionResult jer2 =
    //     applyJER(ak8FatJetPt1, ak8FatJetEta1, ak8FatJetPhi1, ak8FatJetCAmsoftdrop1,
    //              resolution_pt, resolution_pt_sf, rho, GenJetAK8_pt,
    //              GenJetAK8_eta, GenJetAK8_phi, nGenJetAK8);
    // ak8FatJetPt1 = jer2.corrected_pt;
    // ak8FatJetCAmsoftdrop1 = jer2.corrected_massSD;


    // Selection
    if (year == "2022") {
      if (ak8FatJetPt0 < 270 || fabs(ak8FatJetEta0) > 2.4 || ak8FatJetMass0 < 50) {
          continue;
      }
      if (ak8FatJetPt1 > 250 && ak8FatJetMass1 > 50) {
        continue;
      }
      if (ElectronPt0 < 50 || ElectronPt1 > 30) {
        continue;
      }
      if (METPt < 50) {
        continue;
      }
    } else if (year == "2023") {
      if (ak8FatJetPt0 < 250 || fabs(ak8FatJetEta0) > 2.4 || ak8FatJetMass0 < 50) {
      continue;
      }
      h_cutflow->Fill(2.5); // After fatjet0 selection

      if (ak8FatJetPt1 > 200 && ak8FatJetMass1 > 50) {
        continue;
      }
      h_cutflow->Fill(3.5); // After fatjet1 selection
      if (lepton == "ele") {
        if (ElectronPt0 < 50 || ElectronPt1 > 30) {
          continue;
        }
      }
      else if (lepton == "muon") {
        if (MuonPt0 < 50 || MuonPt1 > 30) {
          continue;
        }
      } else {
        throw std::invalid_argument("Invalid lepton type: " + lepton);
      }

      h_cutflow->Fill(4.5); // After lepton kinematic selection

      if (METPt < 50) {
        continue;
      }
      h_cutflow->Fill(5.5); // After MET selection

    } else {
      throw std::invalid_argument("Invalid year: " + year);
    }
    

    double lepton_eta = (lepton == "ele") ? ElectronEta0 : MuonEta0;
    double lepton_phi = (lepton == "ele") ? ElectronPhi0 : MuonPhi0;

    double dR_LFJ = get_dR(lepton_eta, lepton_phi, ak8FatJetEta0, ak8FatJetPhi0);
    if (dR_LFJ < 1.5) {
      continue;
    }
    h_cutflow->Fill(6.5); // After lepton-fatjet deltaR selection

    double dR_JFJ_temp = -1;
    double dR_JL_temp = -1;
    double dR_JFJ = -1;
    double dR_JL = -1;
    int id_selected_jet = 0;
    for (int j = 0; j < 4; j++) {
      if (ak4JetPt[j] > 40){
        dR_JL_temp = get_dR(ak4JetEta[j], ak4JetPhi[j], lepton_eta, lepton_phi);
        dR_JFJ_temp = get_dR(ak4JetEta[j], ak4JetPhi[j], ak8FatJetEta0, ak8FatJetPhi0);
        if (dR_JL_temp > 0.4 && dR_JFJ_temp > 1.5) {
          dR_JL = dR_JL_temp;
          dR_JFJ = dR_JFJ_temp;
          id_selected_jet = j;
          break;
        }
      }
    }
    
    // weight
    weight = (weight / Double_t(SumGenWeights)) * xsec * param_dict.Lumi;
    double PU_weight = PU_Rew[(int)npu];
    if (PU_weight < 20.0) {
      weight = weight * PU_weight;
    }

    T_weight = weight;
    T_fatJet1_pt =  ak8FatJetPt0;
    T_fatJet1_eta = ak8FatJetEta0;
    T_fatJet1_phi = ak8FatJetPhi0;
    T_fatJet1_mass = ak8FatJetMass0;
    T_fatJet1_msd = ak8FatJetMsd0;
    T_fatJet1_CAmass = ak8FatJetCAmass0;
    T_fatJet1_CAmsoftdrop = ak8FatJetCAmsoftdrop0;
    T_fatJet1_CAmass_et = ak8FatJetCAmass_et0;
    T_fatJet1_CAmsoftdrop_et = ak8FatJetCAmsoftdrop_et0;
    T_fatJet1_CAmass_fatjet_et = ak8FatJetCAmass_fatjet_et0;

    T_fatJet1_GloParT_massVis = ak8FatJetCAglobalParT_massVisApplied0 * ak8FatJetCAmass0 * (1.0 - ak8FatJetrawFactor0);
    T_fatJet1_GloParT_massRes = ak8FatJetCAglobalParT_massResApplied0 * ak8FatJetCAmass0 * (1.0 - ak8FatJetrawFactor0);
    T_fatJet1_Tau3OverTau2 = ak8FatJetTau3OverTau20;
    T_fatJet1_ParticleNetLegacy_Xbb = ak8FatJetPNetXbbLegacy0;
    T_fatJet1_ParticleNetLegacy_XbbVsQCD = ak8FatJetPNetXbbvsQCDLegacy0;

    T_fatJet1_GloParT_Xbb = ak8FatJetParTXbb0;
    T_fatJet1_GloParT_XbbVsQCD = ak8FatJetParTXbbvsQCD0;
    T_fatJet1_GloParT_Xtauhtaue = ak8FatJetParTXtauhtaue0;

    T_fatJet2_pt =  ak8FatJetPt1;
    T_fatJet2_eta = ak8FatJetEta1;
    T_fatJet2_phi = ak8FatJetPhi1;
    T_fatJet2_mass = ak8FatJetMass1;
    T_fatJet2_msd = ak8FatJetMsd1;
    T_fatJet2_CAmass = ak8FatJetCAmass1;
    T_fatJet2_CAmsoftdrop = ak8FatJetCAmsoftdrop1;
    T_fatJet2_CAmass_et = ak8FatJetCAmass_et1;
    T_fatJet2_CAmsoftdrop_et = ak8FatJetCAmsoftdrop_et1;
    T_fatJet2_CAmass_fatjet_et = ak8FatJetCAmass_fatjet_et1;
    T_fatJet2_GloParT_massVis = ak8FatJetCAglobalParT_massVisApplied1 * ak8FatJetCAmass1 * (1.0 - ak8FatJetrawFactor1);
    T_fatJet2_GloParT_massRes = ak8FatJetCAglobalParT_massResApplied1 * ak8FatJetCAmass1 * (1.0 - ak8FatJetrawFactor1);

    T_MET = METPt;
    T_ele1_pt = ElectronPt0;
    T_ele1_eta = ElectronEta0;
    T_ele2_phi = ElectronPhi0;
    T_ele2_pt = ElectronPt1;
    T_ele2_eta = ElectronEta1;
    T_ele1_phi = ElectronPhi1;
    T_muon1_pt = MuonPt0;
    T_muon1_eta = MuonEta0;
    T_muon1_phi = MuonPhi0;
    T_muon2_pt = MuonPt1;
    T_muon2_eta = MuonEta1;
    T_muon2_phi = MuonPhi1;
    T_dR_LFJ = dR_LFJ;
    T_dR_JFJ = dR_JFJ;
    T_dR_JL = dR_JL;

    T_jet1_pt = ak4JetPt[id_selected_jet];
    T_jet1_eta = ak4JetEta[id_selected_jet];
    T_jet1_phi = ak4JetPhi[id_selected_jet];
    T_jet1_mass = ak4JetMass[id_selected_jet];
    T_jet1_btag = ak4JetbtagPNetB[id_selected_jet];
    // T_jet2_pt = ak4JetPt1;
    // T_jet2_eta = ak4JetEta1;
    // T_jet2_phi = ak4JetPhi1;
    // T_jet2_mass = ak4JetMass1;  
    // T_jet2_btag = ak4JetbtagPNetB1;

    outputTree->Fill();

  } // end event loop

  h_cutflow->Write();
  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
