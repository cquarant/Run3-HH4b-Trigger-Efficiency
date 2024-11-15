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
#include <bits/stdc++.h>
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

void histo_data(
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
  FactorizedJetCorrector *corrector_AK8 = new FactorizedJetCorrector(vPar_AK8);

  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  Float_t Lower_m[16] = {0,  5,   10,  20,  30,  40,  50,  60,
                         80, 100, 120, 150, 200, 250, 300, 350};
  Float_t Lower_pt[46] = {0,   10,  20,  30,  40,  50,  60,  70,  80,  90,
                          100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
                          200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
                          300, 320, 340, 360, 380, 400, 420, 440, 460, 480,
                          500, 550, 600, 700, 800, 1000};

  // probe FatJet 1 kinematics
  TH1D *_FatJet1_probe_pt =
      new TH1D("FatJet1_pt", "FatJet1_probe_pt", 200, 0, 1000);
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
  TH2D *_FatJet1_probe_Pt_Mass =
      new TH2D("FatJet1_probe_Pt_Mass", "FatJet1_probe_Pt_Mass", 45, Lower_pt,
               15, Lower_m);
  // probe FatJet 1 ParticleNet scores
  TH1D *_FatJet1_probe_PNet_QCD =
      new TH1D("FatJet1_probe_PNet_QCD", "FatJet1_probe_PNet_QCD", 100, 0, 1.0);
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
  TH2D *_FatJet1_tag_Pt_Mass = new TH2D(
      "FatJet1_tag_Pt_Mass", "FatJet1_tag_Pt_Mass", 45, Lower_pt, 15, Lower_m);
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

  // "both" FatJet Kinematics (missing ones)
  TH1D *_FatJet1_both_pt =
      new TH1D("FatJet1_both_pt", "FatJet1_both_pt", 200, 0, 1000);
  TH1D *_FatJet1_both_eta =
      new TH1D("FatJet1_both_eta", "FatJet1_both_eta", 100, -5, 5);
  TH1D *_FatJet1_both_phi =
      new TH1D("FatJet1_both_phi", "FatJet1_both_phi", 100, -5, 5);
  TH2D *_FatJet1_both_eta_phi = new TH2D(
      "FatJet1_both_eta_phi", "FatJet1_both_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet1_both_Mass =
      new TH1D("FatJet1_both_Mass", "FatJet1_both_Mass", 500, 0, 500);
  TH1D *_FatJet1_both_MassSD =
      new TH1D("FatJet1_both_MassSD", "FatJet1_both_MassSD", 500, 0, 500);
  TH2D *_FatJet1_both_Pt_Mass =
      new TH2D("FatJet1_both_Pt_Mass", "FatJet1_both_Pt_Mass", 45, Lower_pt, 15,
               Lower_m);
  // "both" FatJet 1 ParticleNet scores
  TH1D *_FatJet1_both_PNet_QCD =
      new TH1D("FatJet1_both_PNet_QCD", "FatJet1_both_PNet_QCD", 100, 0, 1.0);
  TH1D *_FatJet1_both_PNet_QCD0HF = new TH1D(
      "FatJet1_both_PNet_QCD0HF", "FatJet1_both_PNet_QCD0HF", 100, 0, 1.0);
  TH1D *_FatJet1_both_PNet_QCD1HF = new TH1D(
      "FatJet1_both_PNet_QCD1HF", "FatJet1_both_PNet_QCD1HF", 100, 0, 1.0);
  TH1D *_FatJet1_both_PNet_QCD2HF = new TH1D(
      "FatJet1_both_PNet_QCD2HF", "FatJet1_both_PNet_QCD2HF", 100, 0, 1.0);
  TH1D *_FatJet1_both_PNet_XbbVsQCD = new TH1D(
      "FatJet1_both_PNet_XbbVsQCD", "FatJet1_both_PNet_XbbVsQCD", 100, 0, 1.0);
  TH1D *_FatJet1_both_PNet_XccVsQCD = new TH1D(
      "FatJet1_both_PNet_XccVsQCD", "FatJet1_both_PNet_XccVsQCD", 100, 0, 1.0);
  TH1D *_FatJet1_both_PNet_XggVsQCD = new TH1D(
      "FatJet1_both_PNet_XggVsQCD", "FatJet1_both_PNet_XggVsQCD", 100, 0, 1.0);
  TH1D *_FatJet1_both_PNet_XqqVsQCD = new TH1D(
      "FatJet1_both_PNet_XqqVsQCD", "FatJet1_both_PNet_XqqVsQCD", 100, 0, 1.0);
  // tag FatJet 1 ParticleNetLegacy scores
  TH1D *_FatJet1_both_PNetLegacy_Xbb =
      new TH1D("FatJet1_both_PNetLegacy_Xbb", "FatJet1_both_PNetLegacy_Xbb",
               100, 0, 1.0);
  TH1D *_FatJet1_both_PNetLegacy_Xcc =
      new TH1D("FatJet1_both_PNetLegacy_Xcc", "FatJet1_both_PNetLegacy_Xcc",
               100, 0, 1.0);
  TH1D *_FatJet1_both_PNetLegacy_Xqq =
      new TH1D("FatJet1_both_PNetLegacy_Xqq", "FatJet1_both_PNetLegacy_Xqq",
               100, 0, 1.0);
  TH1D *_FatJet1_both_PNetLegacy_QCD =
      new TH1D("FatJet1_both_PNetLegacy_QCD", "FatJet1_both_PNetLegacy_QCD",
               100, 0, 1.0);
  TH1D *_FatJet1_both_PNetLegacy_QCDb =
      new TH1D("FatJet1_both_PNetLegacy_QCDb", "FatJet1_both_PNetLegacy_QCDb",
               100, 0, 1.0);
  TH1D *_FatJet1_both_PNetLegacy_QCDbb =
      new TH1D("FatJet1_both_PNetLegacy_QCDbb", "FatJet1_both_PNetLegacy_QCDbb",
               100, 0, 1.0);
  TH1D *_FatJet1_both_PNetLegacy_QCDothers =
      new TH1D("FatJet1_both_PNetLegacy_QCDothers",
               "FatJet1_both_PNetLegacy_QCDothers", 100, 0, 1.0);
  // "both" FatJet 1 GloParT scores
  TH1D *_FatJet1_both_GloParT_QCD0HF =
      new TH1D("FatJet1_both_GloParT_QCD0HF", "FatJet1_both_GloParT_QCD0HF",
               100, 0, 1.0);
  TH1D *_FatJet1_both_GloParT_QCD1HF =
      new TH1D("FatJet1_both_GloParT_QCD1HF", "FatJet1_both_GloParT_QCD1HF",
               100, 0, 1.0);
  TH1D *_FatJet1_both_GloParT_QCD2HF =
      new TH1D("FatJet1_both_GloParT_QCD2HF", "FatJet1_both_GloParT_QCD2HF",
               100, 0, 1.0);
  TH1D *_FatJet1_both_GloParT_Xbb = new TH1D(
      "FatJet1_both_GloParT_Xbb", "FatJet1_both_GloParT_Xbb", 100, 0, 1.0);
  TH1D *_FatJet1_both_GloParT_Xcc = new TH1D(
      "FatJet1_both_GloParT_Xcc", "FatJet1_both_GloParT_Xcc", 100, 0, 1.0);
  TH1D *_FatJet1_both_GloParT_Xqq = new TH1D(
      "FatJet1_both_GloParT_Xqq", "FatJet1_both_GloParT_Xqq", 100, 0, 1.0);
  TH1D *_FatJet1_both_GloParT_XbbVsQCD =
      new TH1D("FatJet1_both_GloParT_XbbVsQCD", "FatJet1_both_GloParT_XbbVsQCD",
               100, 0, 1.0);

  TH1D *_MET = new TH1D("MET", "MET", 100, 0, 500);
  TH1D *_Lep1_Pt = new TH1D("Lep1_Pt", "Lep1_Pt", 150, 0, 300);
  TH1D *_dR_LFJ = new TH1D("dR_LFJ", "dR_LFJ", 100, -1.0, 9.0);
  TH1D *_dR_J1FJ = new TH1D("dR_J1FJ", "dR_J1FJ", 100, -1.0, 9.0);
  TH1D *_dR_J2FJ = new TH1D("dR_J2FJ", "dR_J2FJ", 100, -1.0, 9.0);
  TH1D *_dR_JmaxL = new TH1D("dR_JmaxL", "dR_JmaxL", 100, -1.0, 9.0);

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

  // FatJet 1 kinematics
  Float_t FatJet1_pt;
  Float_t FatJet1_eta;
  Float_t FatJet1_phi;
  Float_t FatJet1_Mass;
  Float_t FatJet1_MassSD;
  Float_t FatJet1_rawFactor;
  Int_t FatJet1_hadronFlavour;
  Int_t FatJet1_nBHadrons;
  Int_t FatJet1_nCHadrons;
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

  // tag FatJet 1 kinematics
  Float_t FatJet2_pt;
  Float_t FatJet2_eta;
  Float_t FatJet2_phi;
  Float_t FatJet2_Mass;
  Float_t FatJet2_MassSD;
  Float_t FatJet2_rawFactor;
  Int_t FatJet2_nBHadrons;
  Int_t FatJet2_nCHadrons;

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

  // FatJet 1 kinematics
  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_mass", &FatJet1_Mass);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);
  InputTree->SetBranchAddress("fatJet1_hadronFlavour", &FatJet1_hadronFlavour);
  InputTree->SetBranchAddress("fatJet1_nBHadrons", &FatJet1_nBHadrons);
  InputTree->SetBranchAddress("fatJet1_nCHadrons", &FatJet1_nCHadrons);
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
  // GloParT scores
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
  InputTree->SetBranchAddress("fatJet2_nBHadrons", &FatJet2_nBHadrons);
  InputTree->SetBranchAddress("fatJet2_nCHadrons", &FatJet2_nCHadrons);

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
    if (channel_lower == "egamma" or channel_lower == "electron") {
      if (!(HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) == 11)) {
        continue;
      }
    } else if (channel_lower == "muon") {
      if (!(HLT_IsoMu27 && fabs(lep1_Id) == 13)) {
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
    if (!Certified)
      continue;

    // FatJets correction and selection
    if (FatJet1_pt > 0) {
      double Raw_FatJet1_probe_pt = FatJet1_pt * (1.0 - FatJet1_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_probe_pt);
      corrector_AK8->setJetEta(FatJet1_eta);
      corrector_AK8->setJetPhi(FatJet1_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet1_pt = Raw_FatJet1_probe_pt * corr;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * corr;
    }
    if (FatJet2_pt > 0) {
      double Raw_FatJet1_tag_pt = FatJet2_pt * (1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_tag_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      corrector_AK8->setJetPhi(FatJet2_phi);
      double corr = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet1_tag_pt * corr;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * corr;
    }

    if (FatJet1_pt < 250 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50)
      continue;

    /*
    // VBFTag veto
    if (isVBFtag) {
        continue;
    }
    */

    // Lepton selection or veto
    if (lep1_Pt < 50) {
      continue;
    }
    if (lep2_Pt > 30) {
      continue;
    }
    if (FatJet2_pt > 200 && FatJet2_MassSD > 50) {
      continue;
    }

    double dR_LFJ = get_dR(lep1_Eta, lep1_Phi, FatJet1_eta, FatJet1_phi);
    if (dR_LFJ < 1.5) {
      continue;
    }
    if (MET < 50) {
      continue;
    }
    if (FatJet1_MassSD < 50) {
      continue;
    }
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

    if (dR_JFJ_Max < 1.5) {
      continue;
    }
    if (dR_J1L <= 0.4 || dR_J2L <= 0.4) {
      continue;
    }
    if (dR_JmaxL > 3.5) {
      continue;
    }

    // Trigger Objects and Matchings
    // Matching 1st
    bool matched_TRG_1 = false;

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

    if (matched_to_AK8PFJet230_SoftDropMass40)
      matched_TRG_1 = true;

    // Matching 2nd
    bool matched_TRG_2 = false;
    // if (HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06) {
    //   matched_TRG_2=true;
    // }
    if (channel_lower == "egamma") {
      if ((HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06 &&
           abs(lep1_Id) == 11)) {
        matched_TRG_2 = true;
      }
    } else if (channel_lower == "muon") {
      if ((HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06 &&
           abs(lep1_Id) == 13)) {
        matched_TRG_2 = true;
      }
    } else {
      throw std::invalid_argument("Invalid channel: " + channel);
    }

    // Fill histograms
    // kinematics
    _FatJet1_probe_pt->Fill(FatJet1_pt);
    _FatJet1_probe_eta->Fill(FatJet1_eta);
    _FatJet1_probe_phi->Fill(FatJet1_phi);
    _FatJet1_probe_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
    _FatJet1_probe_Mass->Fill(FatJet1_Mass);
    _FatJet1_probe_MassSD->Fill(FatJet1_MassSD);
    _FatJet1_probe_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);
    // ParticleNet
    _FatJet1_probe_PNet_QCD->Fill(FatJet1PNet_QCD);
    _FatJet1_probe_PNet_QCD0HF->Fill(FatJet1PNet_QCD0HF);
    _FatJet1_probe_PNet_QCD1HF->Fill(FatJet1PNet_QCD1HF);
    _FatJet1_probe_PNet_QCD2HF->Fill(FatJet1PNet_QCD2HF);
    _FatJet1_probe_PNet_XbbVsQCD->Fill(FatJet1PNet_XbbVsQCD);
    _FatJet1_probe_PNet_XccVsQCD->Fill(FatJet1PNet_XccVsQCD);
    _FatJet1_probe_PNet_XggVsQCD->Fill(FatJet1PNet_XggVsQCD);
    _FatJet1_probe_PNet_XqqVsQCD->Fill(FatJet1PNet_XqqVsQCD);
    // ParticleNetLegacy
    _FatJet1_probe_PNetLegacy_Xbb->Fill(FatJet1PNetLegacy_Xbb);
    _FatJet1_probe_PNetLegacy_Xcc->Fill(FatJet1PNetLegacy_Xcc);
    _FatJet1_probe_PNetLegacy_Xqq->Fill(FatJet1PNetLegacy_Xqq);
    _FatJet1_probe_PNetLegacy_QCD->Fill(FatJet1PNetLegacy_QCD);
    _FatJet1_probe_PNetLegacy_QCDb->Fill(FatJet1PNetLegacy_QCDb);
    _FatJet1_probe_PNetLegacy_QCDbb->Fill(FatJet1PNetLegacy_QCDbb);
    _FatJet1_probe_PNetLegacy_QCDothers->Fill(FatJet1PNetLegacy_QCDothers);
    // GloParT
    _FatJet1_probe_GloParT_QCD0HF->Fill(FatJet1GloParT_QCD0HF);
    _FatJet1_probe_GloParT_QCD1HF->Fill(FatJet1GloParT_QCD1HF);
    _FatJet1_probe_GloParT_QCD2HF->Fill(FatJet1GloParT_QCD2HF);
    _FatJet1_probe_GloParT_Xbb->Fill(FatJet1GloParT_Xbb);
    _FatJet1_probe_GloParT_Xcc->Fill(FatJet1GloParT_Xcc);
    _FatJet1_probe_GloParT_Xqq->Fill(FatJet1GloParT_Xqq);
    _FatJet1_probe_GloParT_XbbVsQCD->Fill(FatJet1GloParT_XbbVsQCD);

    if (matched_TRG_2) {
      // kinematics
      _FatJet1_tag_pt->Fill(FatJet1_pt);
      _FatJet1_tag_eta->Fill(FatJet1_eta);
      _FatJet1_tag_phi->Fill(FatJet1_phi);
      _FatJet1_tag_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
      _FatJet1_tag_Mass->Fill(FatJet1_Mass);
      _FatJet1_tag_MassSD->Fill(FatJet1_MassSD);
      _FatJet1_tag_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);
      // ParticleNet
      _FatJet1_tag_PNet_QCD->Fill(FatJet1PNet_QCD);
      _FatJet1_tag_PNet_QCD0HF->Fill(FatJet1PNet_QCD0HF);
      _FatJet1_tag_PNet_QCD1HF->Fill(FatJet1PNet_QCD1HF);
      _FatJet1_tag_PNet_QCD2HF->Fill(FatJet1PNet_QCD2HF);
      _FatJet1_tag_PNet_XbbVsQCD->Fill(FatJet1PNet_XbbVsQCD);
      _FatJet1_tag_PNet_XccVsQCD->Fill(FatJet1PNet_XccVsQCD);
      _FatJet1_tag_PNet_XggVsQCD->Fill(FatJet1PNet_XggVsQCD);
      _FatJet1_tag_PNet_XqqVsQCD->Fill(FatJet1PNet_XqqVsQCD);
      // ParticleNetLegacy
      _FatJet1_tag_PNetLegacy_Xbb->Fill(FatJet1PNetLegacy_Xbb);
      _FatJet1_tag_PNetLegacy_Xcc->Fill(FatJet1PNetLegacy_Xcc);
      _FatJet1_tag_PNetLegacy_Xqq->Fill(FatJet1PNetLegacy_Xqq);
      _FatJet1_tag_PNetLegacy_QCD->Fill(FatJet1PNetLegacy_QCD);
      _FatJet1_tag_PNetLegacy_QCDb->Fill(FatJet1PNetLegacy_QCDb);
      _FatJet1_tag_PNetLegacy_QCDbb->Fill(FatJet1PNetLegacy_QCDbb);
      _FatJet1_tag_PNetLegacy_QCDothers->Fill(FatJet1PNetLegacy_QCDothers);
      // GloParT
      _FatJet1_tag_GloParT_QCD0HF->Fill(FatJet1GloParT_QCD0HF);
      _FatJet1_tag_GloParT_QCD1HF->Fill(FatJet1GloParT_QCD1HF);
      _FatJet1_tag_GloParT_QCD2HF->Fill(FatJet1GloParT_QCD2HF);
      _FatJet1_tag_GloParT_Xbb->Fill(FatJet1GloParT_Xbb);
      _FatJet1_tag_GloParT_Xcc->Fill(FatJet1GloParT_Xcc);
      _FatJet1_tag_GloParT_Xqq->Fill(FatJet1GloParT_Xqq);
      _FatJet1_tag_GloParT_XbbVsQCD->Fill(FatJet1GloParT_XbbVsQCD);
    }

    if (matched_TRG_1 && matched_TRG_2) {
      // kinematics
      _FatJet1_both_pt->Fill(FatJet1_pt);
      _FatJet1_both_eta->Fill(FatJet1_eta);
      _FatJet1_both_phi->Fill(FatJet1_phi);
      _FatJet1_both_eta_phi->Fill(FatJet1_eta, FatJet1_phi);
      _FatJet1_both_Mass->Fill(FatJet1_Mass);
      _FatJet1_both_MassSD->Fill(FatJet1_MassSD);
      _FatJet1_both_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD);
      // ParticleNet
      _FatJet1_both_PNet_QCD->Fill(FatJet1PNet_QCD);
      _FatJet1_both_PNet_QCD0HF->Fill(FatJet1PNet_QCD0HF);
      _FatJet1_both_PNet_QCD1HF->Fill(FatJet1PNet_QCD1HF);
      _FatJet1_both_PNet_QCD2HF->Fill(FatJet1PNet_QCD2HF);
      _FatJet1_both_PNet_XbbVsQCD->Fill(FatJet1PNet_XbbVsQCD);
      _FatJet1_both_PNet_XccVsQCD->Fill(FatJet1PNet_XccVsQCD);
      _FatJet1_both_PNet_XggVsQCD->Fill(FatJet1PNet_XggVsQCD);
      _FatJet1_both_PNet_XqqVsQCD->Fill(FatJet1PNet_XqqVsQCD);
      // ParticleNetLegacy
      _FatJet1_both_PNetLegacy_Xbb->Fill(FatJet1PNetLegacy_Xbb);
      _FatJet1_both_PNetLegacy_Xcc->Fill(FatJet1PNetLegacy_Xcc);
      _FatJet1_both_PNetLegacy_Xqq->Fill(FatJet1PNetLegacy_Xqq);
      _FatJet1_both_PNetLegacy_QCD->Fill(FatJet1PNetLegacy_QCD);
      _FatJet1_both_PNetLegacy_QCDb->Fill(FatJet1PNetLegacy_QCDb);
      _FatJet1_both_PNetLegacy_QCDbb->Fill(FatJet1PNetLegacy_QCDbb);
      _FatJet1_both_PNetLegacy_QCDothers->Fill(FatJet1PNetLegacy_QCDothers);
      // GloParT
      _FatJet1_both_GloParT_QCD0HF->Fill(FatJet1GloParT_QCD0HF);
      _FatJet1_both_GloParT_QCD1HF->Fill(FatJet1GloParT_QCD1HF);
      _FatJet1_both_GloParT_QCD2HF->Fill(FatJet1GloParT_QCD2HF);
      _FatJet1_both_GloParT_Xbb->Fill(FatJet1GloParT_Xbb);
      _FatJet1_both_GloParT_Xcc->Fill(FatJet1GloParT_Xcc);
      _FatJet1_both_GloParT_Xqq->Fill(FatJet1GloParT_Xqq);
      _FatJet1_both_GloParT_XbbVsQCD->Fill(FatJet1GloParT_XbbVsQCD);
    }

    _MET->Fill(MET);
    _Lep1_Pt->Fill(lep1_Pt);
    _dR_LFJ->Fill(dR_LFJ);
    _dR_J1FJ->Fill(dR_JFJ_Max);
    _dR_J2FJ->Fill(dR_JFJ_Min);
    _dR_JmaxL->Fill(dR_JmaxL);

  } // end event loop

  f->Write();
  std::cout << "Done. Written to " << output_path << std::endl;
}
