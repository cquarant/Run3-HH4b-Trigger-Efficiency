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


void tree_data_he(const std::string &year,      // 2022, 2023
                  const std::string &run_tag,   // Example: 2022C, 2023D
                  const std::string &channel,              // e.g. Muon, EGamma, JetMET
                  const std::string &sample_path, // path to the root file
                  const std::string &output_path, // path to the output root file
                  const std::string &jec_path_L2Relative,  // path to the AK4 JEC txt file
                  const std::string &jec_path_L2L3Residual, // path to the AK8 JEC txt file
                  const int begin_file_index=0,
                  const int end_file_index=999999999
) {

  // Load appropriate Good Lumi List based on era
  // TODO: add more eras if needed
  map<int, vector<pair<int, int>>> Good_Lumis;
  if (run_tag == "2022C") {
    Good_Lumis = {
#include "/afs/cern.ch/work/c/cquarant/Hbt2/CMSSW_16_0_0_pre3/src/TTbarBkgEstimation/Run3-HH4b-Trigger-Efficiency/GoodLumiList/GoodLumiList_Map_2022C.txt"
    };
  } else if (run_tag == "2022D") {
    Good_Lumis = {
#include "/afs/cern.ch/work/c/cquarant/Hbt2/CMSSW_16_0_0_pre3/src/TTbarBkgEstimation/Run3-HH4b-Trigger-Efficiency/GoodLumiList/GoodLumiList_Map_2022D.txt"
    };
  } else if (run_tag == "2022E") {
    Good_Lumis = {
#include "/afs/cern.ch/work/c/cquarant/Hbt2/CMSSW_16_0_0_pre3/src/TTbarBkgEstimation/Run3-HH4b-Trigger-Efficiency/GoodLumiList/GoodLumiList_Map_2022E.txt"
    };
  } else if (run_tag == "2022F") {
    Good_Lumis = {
#include "/afs/cern.ch/work/c/cquarant/Hbt2/CMSSW_16_0_0_pre3/src/TTbarBkgEstimation/Run3-HH4b-Trigger-Efficiency/GoodLumiList/GoodLumiList_Map_2022F.txt"
    };
  } else if (run_tag == "2022G") {
    Good_Lumis = {
#include "/afs/cern.ch/work/c/cquarant/Hbt2/CMSSW_16_0_0_pre3/src/TTbarBkgEstimation/Run3-HH4b-Trigger-Efficiency/GoodLumiList/GoodLumiList_Map_2022G.txt"
    };
  } else if (run_tag == "2023C") {
    Good_Lumis = {
#include "/afs/cern.ch/work/c/cquarant/Hbt2/CMSSW_16_0_0_pre3/src/TTbarBkgEstimation/Run3-HH4b-Trigger-Efficiency/GoodLumiList/GoodLumiList_Map_2023C.txt"
    };
  } else if (run_tag == "2023D") {
    Good_Lumis = {
#include "/afs/cern.ch/work/c/cquarant/Hbt2/CMSSW_16_0_0_pre3/src/TTbarBkgEstimation/Run3-HH4b-Trigger-Efficiency/GoodLumiList/GoodLumiList_Map_2023D.txt"
    };
  } else {
    throw std::runtime_error("Unsupported era for lumi dict: " + run_tag);
  }

  gSystem->Load("libFWCoreFWLite.so");

  if (year != "2022" && year != "2023") {
    throw std::invalid_argument("Invalid year: " + year);
  }

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

  
  // Outputs
  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  TH1F *h_cutflow = new TH1F("cutflow", "cutflow", 10, 0, 10);
  h_cutflow->GetXaxis()->SetBinLabel(1, "Total Events");
  h_cutflow->GetXaxis()->SetBinLabel(2, "HLT trigger selection");
  h_cutflow->GetXaxis()->SetBinLabel(3, "fatjet0 kinematic selection");   
  h_cutflow->GetXaxis()->SetBinLabel(4, "fatjet1 kinematic selection");   
  h_cutflow->GetXaxis()->SetBinLabel(5, "lepton kinematic selection");
  h_cutflow->GetXaxis()->SetBinLabel(6, "MET selection");
  h_cutflow->GetXaxis()->SetBinLabel(7, "lepton-FatJet dR cut");
  h_cutflow->GetXaxis()->SetBinLabel(8, "dR JL cut");
  h_cutflow->GetXaxis()->SetBinLabel(9, "dR JFJ cut (disabled)");

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
  Float_t T_MET, T_lep1_pt, T_lep1_eta, T_lep1_phi, T_dR_LFJ, T_dR_JFJ_max,
      T_dR_JFJ_min, T_dR_JmaxL;
  Float_t T_jet1_pt, T_jet1_eta, T_jet1_phi, T_jet1_mass, T_jet1_btag;
  Float_t T_jet2_pt, T_jet2_eta, T_jet2_phi, T_jet2_mass, T_jet2_btag;

  outputTree->Branch("weight", &T_weight, "weight/F");
  outputTree->Branch("fatJet1_pt", &T_fatJet1_pt, "fatJet1_pt/F");
  outputTree->Branch("fatJet1_eta", &T_fatJet1_eta, "fatJet1_eta/F");
  outputTree->Branch("fatJet1_phi", &T_fatJet1_phi, "fatJet1_phi/F");
  outputTree->Branch("fatJet1_msoftdrop", &T_fatJet1_CAmsoftdrop,
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
  outputTree->Branch("fatJet1_GloParT_XbbVsQCD", &T_fatJet1_GloParT_XbbVsQCD,
                     "fatJet1_GloParT_XbbVsQCD/F");
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
  outputTree->Branch("lep1_pt", &T_lep1_pt, "lep1_pt/F");
  outputTree->Branch("lep1_eta", &T_lep1_eta, "lep1_eta/F");
  outputTree->Branch("lep1_phi", &T_lep1_phi, "lep1_phi/F");
  outputTree->Branch("dR_LFJ", &T_dR_LFJ, "dR_LFJ/F");
  outputTree->Branch("dR_J1FJ", &T_dR_JFJ_max, "dR_J1FJ/F");
  outputTree->Branch("dR_J2FJ", &T_dR_JFJ_min, "dR_J2FJ/F");
  outputTree->Branch("dR_JmaxL", &T_dR_JmaxL, "dR_JmaxL/F");

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

  
  // TH1F *Events = (TH1F *)f1->Get("Events");
  // std::cout << "SumGenWeights: " << SumGenWeights << std::endl;
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
  
  int SumGenWeights = InputTree->GetEntries();
  std::cout << "SumGenWeights (Total Entries): " << SumGenWeights << std::endl;
  
  UInt_t run;
  UInt_t lumi;
  // Int_t npu;
  // Int_t isVBFtag;
  // Float_t rho;

  Long64_t HLT_Ele30_WPTight_Gsf;

  Float_t METPt;

  Double_t ElectronPt0;
  Double_t ElectronEta0;
  Double_t ElectronPhi0;
  Double_t ElectronPt1;
  Double_t ElectronEta1;
  Double_t ElectronPhi1;
  
  Double_t ak4JetPt0;
  Double_t ak4JetEta0;
  Double_t ak4JetPhi0;
  Double_t ak4JetMass0;
  Double_t ak4JetbtagPNetB0;
  Double_t ak4JetPt1;
  Double_t ak4JetEta1;
  Double_t ak4JetPhi1;
  Double_t ak4JetMass1;
  Double_t ak4JetbtagPNetB1;

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

  InputTree->SetBranchAddress("run", &run);
  InputTree->SetBranchAddress("luminosityBlock", &lumi);
  // InputTree->SetBranchAddress("nPU", &npu);

  // Electron triggers
  InputTree->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf);

  //Electron kinematics
  InputTree->SetBranchAddress("ElectronEta0", &ElectronEta0);
  InputTree->SetBranchAddress("ElectronEta1", &ElectronEta1);
  InputTree->SetBranchAddress("ElectronPhi0", &ElectronPhi0);
  InputTree->SetBranchAddress("ElectronPhi1", &ElectronPhi1);
  InputTree->SetBranchAddress("ElectronPt0", &ElectronPt0);
  InputTree->SetBranchAddress("ElectronPt1", &ElectronPt1);

  InputTree->SetBranchAddress("METPt", &METPt);

  // First two AK4 jets in pT order (PROBLEM: these are likely inside or close to the fat jet)
  // InputTree->SetBranchAddress("ak4JetEta0", &ak4JetEta0);
  // InputTree->SetBranchAddress("ak4JetEta1", &ak4JetEta1);
  // InputTree->SetBranchAddress("ak4JetPhi0", &ak4JetPhi0);
  // InputTree->SetBranchAddress("ak4JetPhi1", &ak4JetPhi1);
  // InputTree->SetBranchAddress("ak4JetPt0", &ak4JetPt0);
  // InputTree->SetBranchAddress("ak4JetPt1", &ak4JetPt1);
  // InputTree->SetBranchAddress("ak4JetMass0", &ak4JetMass0);
  // InputTree->SetBranchAddress("ak4JetMass1", &ak4JetMass1);
  // InputTree->SetBranchAddress("ak4JetbtagPNetB0", &ak4JetbtagPNetB0);
  // InputTree->SetBranchAddress("ak4JetbtagPNetB1", &ak4JetbtagPNetB1);
  InputTree->SetBranchAddress("AK4JetAwayEta0", &ak4JetEta0);
  InputTree->SetBranchAddress("AK4JetAwayEta1", &ak4JetEta1);
  InputTree->SetBranchAddress("AK4JetAwayPhi0", &ak4JetPhi0);
  InputTree->SetBranchAddress("AK4JetAwayPhi1", &ak4JetPhi1);
  InputTree->SetBranchAddress("AK4JetAwayPt0", &ak4JetPt0);
  InputTree->SetBranchAddress("AK4JetAwayPt1", &ak4JetPt1);
  InputTree->SetBranchAddress("AK4JetAwayMass0", &ak4JetMass0);
  InputTree->SetBranchAddress("AK4JetAwayMass1", &ak4JetMass1);
  InputTree->SetBranchAddress("AK4JetAwaybtagPNetB0", &ak4JetbtagPNetB0);
  InputTree->SetBranchAddress("AK4JetAwaybtagPNetB1", &ak4JetbtagPNetB1);
  

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
  InputTree->SetBranchAddress("ak8FatJetMsd1", &ak8FatJetMsd1);
  InputTree->SetBranchAddress("ak8FatJetCAmass1", &ak8FatJetCAmass1);
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

  // Events Loop
  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);
    if (i % 100000 == 0) {
      std::cout << "Processing event " << i << " / " << InputTree->GetEntries()
                << "\r" << std::flush;
    }
    h_cutflow->Fill(0.5); // Total Events

    // HLT Selection
    bool HLT_ele = (HLT_Ele30_WPTight_Gsf);
    if (!HLT_ele) {
      continue;
    }

    h_cutflow->Fill(1.5); // After HLT trigger selection

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

    // FatJets correction (JEC and JER)
    JetCorrectionResult jec1 =
        applyJEC(ak8FatJetPt0, ak8FatJetEta0, ak8FatJetPhi0, ak8FatJetrawFactor0,
                 ak8FatJetCAmsoftdrop0, corrector);
    ak8FatJetPt0 = jec1.corrected_pt;
    ak8FatJetCAmsoftdrop0 = jec1.corrected_massSD;

    JetCorrectionResult jec2 =
        applyJEC(ak8FatJetPt1, ak8FatJetEta1, ak8FatJetPhi1, ak8FatJetrawFactor1,
                 ak8FatJetCAmsoftdrop1, corrector);
    ak8FatJetPt1 = jec2.corrected_pt;
    ak8FatJetCAmsoftdrop1 = jec2.corrected_massSD;


    // Selection
    if (year == "2022") {
      if (ak8FatJetPt0 < 270 || fabs(ak8FatJetEta0) > 2.4 || ak8FatJetMsd0 < 50) {
          continue;
      }
      if (ak8FatJetPt1 > 250 && ak8FatJetMsd1 > 50) {
        continue;
      }
      if (ElectronPt0 < 50 || ElectronPt1 > 30) {
        continue;
      }
      if (METPt < 50) {
        continue;
      }
    } else if (year == "2023") {
      if (ak8FatJetPt0 < 250 || fabs(ak8FatJetEta0) > 2.5 || ak8FatJetMass0 < 50) {
      continue;
      }
      h_cutflow->Fill(2.5); // After fatjet0 kinematic selection
      if (ak8FatJetPt1 > 200 && ak8FatJetMass1 > 50) {
        continue;
      }
      h_cutflow->Fill(3.5); // After fatjet1 kinematic selection
      if (ElectronPt0 < 50 || ElectronPt1 > 30) {
        continue;
      }
      h_cutflow->Fill(4.5); // After lepton kinematic selection
      if (METPt < 50) {
        continue;
      }
      h_cutflow->Fill(5.5); // After MET selection
    } else {
      throw std::invalid_argument("Invalid year: " + year);
    }
    
    double dR_LFJ = get_dR(ElectronEta0, ElectronPhi0, ak8FatJetEta0, ak8FatJetPhi0);
    if (dR_LFJ < 1.5) {
      continue;
    }

    h_cutflow->Fill(6.5); // After lepton-FatJet dR cut
    
    double dR_J1FJ = -1;
    double dR_J1L = 10;
    if (ak4JetPt0 > 40) {
      dR_J1FJ = get_dR(ak4JetEta0, ak4JetPhi0, ak8FatJetEta0, ak8FatJetPhi0);
      dR_J1L = get_dR(ak4JetEta0, ak4JetPhi0, ElectronEta0, ElectronPhi0);
    }
    double dR_J2FJ = -1;
    double dR_J2L = 10;
    if (ak4JetPt1 > 40) {
      dR_J2FJ = get_dR(ak4JetEta1, ak4JetPhi1, ak8FatJetEta0, ak8FatJetPhi0);
      dR_J2L = get_dR(ak4JetEta1, ak4JetPhi1, ElectronEta0, ElectronPhi0);
    }
    

    double dR_JFJ_max = dR_J1FJ;
    double dR_JFJ_min = dR_J2FJ;
    double dR_JmaxL = dR_J1L;
    if (dR_J2FJ > dR_J1FJ) {
      dR_JFJ_max = dR_J2FJ;
      dR_JFJ_min = dR_J1FJ;
      dR_JmaxL = max(dR_J1L, dR_J2L);
    }

    if (dR_J1L <= 0.4 || dR_J2L <= 0.4) {
      continue;
    }
    h_cutflow->Fill(7.5); // After dR_JL cuts

    // if ( (dR_J1FJ > 0. && dR_J1FJ <= 1.5) || (dR_J2FJ > 0. && dR_J2FJ <= 1.5)) {
    //   continue;
    // }
    h_cutflow->Fill(8.5); // After dR_JFJ cuts


    T_weight = 1.0;
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
    T_lep1_pt = ElectronPt0;
    T_lep1_eta = ElectronEta0;
    T_lep1_phi = ElectronPhi0;
    T_dR_LFJ = dR_LFJ;
    T_dR_JFJ_max = dR_JFJ_max;
    T_dR_JFJ_min = dR_JFJ_min;
    T_dR_JmaxL = dR_JmaxL;

    T_jet1_pt = ak4JetPt0;
    T_jet1_eta = ak4JetEta0;
    T_jet1_phi = ak4JetPhi0;
    T_jet1_mass = ak4JetMass0;
    T_jet1_btag = ak4JetbtagPNetB0;
    T_jet2_pt = ak4JetPt1;
    T_jet2_eta = ak4JetEta1;
    T_jet2_phi = ak4JetPhi1;
    T_jet2_mass = ak4JetMass1;  
    T_jet2_btag = ak4JetbtagPNetB1;

    outputTree->Fill();

  } // end event loop

  h_cutflow->Write();
  f->Write();

  std::cout << "Done. Written to " << output_path << std::endl;
}
