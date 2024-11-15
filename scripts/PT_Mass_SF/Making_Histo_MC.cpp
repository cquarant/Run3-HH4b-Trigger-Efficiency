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
#include <bits/stdc++.h>

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
      } else if (name == "XSec_TTtoLNu2Q") {
        param_dict->XSec_TTtoLNu2Q = value;
      } else if (name == "XSec_TTto4Q") {
        param_dict->XSec_TTto4Q = value;
      } else if (name == "XSec_TTto2L2Nu") {
        param_dict->XSec_TTto2L2Nu = value;
      }
    }
  }
}

double getXSec(ParamDict *param_dict, std::string ttbar_type) {
  if (ttbar_type == "TTtoLNu2Q") {
    return param_dict->XSec_TTtoLNu2Q;
  } else if (ttbar_type == "TTto4Q") {
    return param_dict->XSec_TTto4Q;
  } else if (ttbar_type == "TTto2L2Nu") {
    return param_dict->XSec_TTto2L2Nu;
  } else {
    throw std::invalid_argument("Invalid ttbar_type: " + ttbar_type);
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

void Making_Histo_MC(
    const std::string &ttbar_type,   // default: TTtoLNu2Q
    const std::string &channel,      // Muon, EGamma, leptonic
    const std::string &sample_path,  // path to the root file
    const std::string &output_path,  // path to the output root file
    const std::string &pu_path,      // path to the pileup reweighting file
    const std::string &sf_path,      // path to the PT-Mass-SFs root file (not needed for this method)
    const std::string &param_path,   // path to the parameters file
    const std::string &jec_path_ak4, // path to the AK4 JEC txt file
    const std::string &jec_path_ak8  // path to the AK8 JEC txt file
) {
  gSystem->Load("libFWCoreFWLite.so");

  std::string channel_lower = to_lower(channel);
  std::cout << "TTBar Type: " << ttbar_type << std::endl;
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
  double xsec = getXSec(&param_dict, ttbar_type);

  // pu weight
  std::vector<double> PU_Rew = loadPUReweighting(pu_path);

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
  Float_t Lower_pt_N[9] = {230, 240, 250, 270, 300, 350, 500, 700, 1000};

  // probe histograms
  TH1D *_FatJet1_probe_pt = new TH1D("FatJet1_probe_pt", "FatJet1_probe_pt", 200, 0, 1000);
  TH1D *_FatJet1_probe_eta = new TH1D("FatJet1_probe_eta", "FatJet1_probe_eta", 100, -5, 5);
  TH1D *_FatJet1_probe_phi = new TH1D("FatJet1_probe_phi", "FatJet1_probe_phi", 100, -5, 5);
  TH2D *_FatJet1_probe_eta_phi =
      new TH2D("FatJet1_probe_eta_phi", "FatJet1_probe_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet1_probe_Mass = new TH1D("FatJet1_probe_Mass", "FatJet1_probe_Mass", 500, 0, 500);
  TH1D *_FatJet1_probe_MassSD =
      new TH1D("FatJet1_probe_MassSD", "FatJet1_probe_MassSD", 500, 0, 500);
  TH2D *_FatJet1_probe_Pt_Mass =
      new TH2D("FatJet1_probe_Pt_Mass", "FatJet1_probe_Pt_Mass", 45, Lower_pt, 15, Lower_m);

  // tag histograms
  TH1D *_FatJet1_tag_pt = new TH1D("FatJet1_tag_pt", "FatJet1_tag_pt", 200, 0, 1000);
  TH1D *_FatJet1_tag_eta = new TH1D("FatJet1_tag_eta", "FatJet1_tag_eta", 100, -5, 5);
  TH1D *_FatJet1_tag_phi = new TH1D("FatJet1_tag_phi", "FatJet1_tag_phi", 100, -5, 5);
  TH2D *_FatJet1_tag_eta_phi =
      new TH2D("FatJet1_tag_eta_phi", "FatJet1_tag_eta_phi", 100, -5, 5, 100, -5, 5);
  TH1D *_FatJet1_tag_Mass = new TH1D("FatJet1_tag_Mass", "FatJet1_tag_Mass", 500, 0, 500);
  TH1D *_FatJet1_tag_MassSD =
      new TH1D("FatJet1_tag_MassSD", "FatJet1_tag_MassSD", 500, 0, 500);
  TH2D *_FatJet1_tag_Pt_Mass =
      new TH2D("FatJet1_tag_Pt_Mass", "FatJet1_tag_Pt_Mass", 45, Lower_pt, 15, Lower_m);

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

  Float_t FatJet1_pt;
  Float_t FatJet1_eta;
  Float_t FatJet1_phi;
  Float_t FatJet1_Mass;
  Float_t FatJet1_MassSD;
  Float_t FatJet1_rawFactor;

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

  InputTree->SetBranchAddress("fatJet1_pt", &FatJet1_pt);
  InputTree->SetBranchAddress("fatJet1_eta", &FatJet1_eta);
  InputTree->SetBranchAddress("fatJet1_phi", &FatJet1_phi);
  InputTree->SetBranchAddress("fatJet1_mass", &FatJet1_Mass);
  InputTree->SetBranchAddress("fatJet1_msoftdrop", &FatJet1_MassSD);
  InputTree->SetBranchAddress("fatJet1_rawFactor", &FatJet1_rawFactor);

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
    } else {
      throw std::invalid_argument("Invalid channel: " + channel);
    }

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

    /*
       // Jet Smearing
       double res_pt_1;
       double res_pt_sf_1;
       JME::JetParameters JerPARAM_1 = {{JME::Binning::JetPt, FatJet1_pt},
    {JME::Binning::JetEta, FatJet1_eta},{JME::Binning::Rho, rho}};
       JME::JetParameters JerSFPARAM_1;
       JerSFPARAM_1.set(JME::Binning::JetPt,  FatJet1_pt);
       JerSFPARAM_1.set(JME::Binning::JetEta, FatJet1_eta);
       JerSFPARAM_1.set(JME::Binning::Rho, rho);
       res_pt_1 = resolution_pt_AK8.getResolution(JerPARAM_1);
       res_pt_sf_1 = resolution_pt_sf_AK8.getScaleFactor(JerSFPARAM_1);

       double res_pt_2;
       double res_pt_sf_2;
       JME::JetParameters JerPARAM_2 = {{JME::Binning::JetPt, FatJet2_pt},
    {JME::Binning::JetEta, FatJet2_eta},{JME::Binning::Rho, rho}};
       JME::JetParameters JerSFPARAM_2;
       JerSFPARAM_2.set(JME::Binning::JetPt,  FatJet2_pt);
       JerSFPARAM_2.set(JME::Binning::JetEta, FatJet2_eta);
       JerSFPARAM_2.set(JME::Binning::Rho, rho);
       res_pt_2 = resolution_pt_AK8.getResolution(JerPARAM_2);
       res_pt_sf_2 = resolution_pt_sf_AK8.getScaleFactor(JerSFPARAM_2);

       double SmearFactor_1   = 1;
       bool   GenJetMatched_1 = false;
       double SmearFactor_2   = 1;
       bool   GenJetMatched_2 = false;

       for(int nGJAK8=0;nGJAK8<nGenJetAK8; nGJAK8++)
         {
           if(!GenJetMatched_1 &&  sqrt(pow(FatJet1_eta-GenJetAK8_eta[nGJAK8],2)
    + pow(phi_dist(FatJet1_phi,GenJetAK8_phi[nGJAK8]),2)) < 0.2 &&
    (fabs(FatJet1_pt - GenJetAK8_pt[nGJAK8])/FatJet1_pt < 3*res_pt_1) )
             {
               SmearFactor_1 = 1.0 + (res_pt_sf_1 - 1.0) * (FatJet1_pt -
    GenJetAK8_pt[nGJAK8]) / FatJet1_pt; GenJetMatched_1 = true;
             }
           if(!GenJetMatched_2 &&  sqrt(pow(FatJet2_eta-GenJetAK8_eta[nGJAK8],2)
    + pow(phi_dist(FatJet2_phi,GenJetAK8_phi[nGJAK8]),2)) < 0.2 &&
    (fabs(FatJet2_pt - GenJetAK8_pt[nGJAK8])/FatJet2_pt < 3*res_pt_2) )
             {
               SmearFactor_2 = 1.0 + (res_pt_sf_2 - 1.0) * (FatJet2_pt -
    GenJetAK8_pt[nGJAK8]) / FatJet2_pt; GenJetMatched_2 = true;
             }
         }

    //     if(!GenJetMatched && res_pt_sf[nJ] > 1.0)
    //       {
    //        double sigma = res_pt[nJ] * sqrt(res_pt_sf[nJ]*res_pt_sf[nJ] - 1);
    //        normal_distribution<> d(0, sigma);
    //        SmearFactor = 1.0 + d(m_random_generator);
    //       }

         // Smear
         FatJet1_pt     = FatJet1_pt * SmearFactor_1;
         FatJet1_MassSD = FatJet1_MassSD * SmearFactor_1;
         FatJet2_pt     = FatJet2_pt * SmearFactor_2;
         FatJet2_MassSD = FatJet2_MassSD * SmearFactor_2;
    */

    // FatJets selection
    if (lep1_Pt < 50 || lep2_Pt > 30) {
      continue;
    }
    if (FatJet2_pt > 180 || FatJet1_pt < 160) {
      continue;
    }
    if (phi_dist(FatJet1_phi, lep1_Phi) < 2.0) {
      continue;
    }
    if (MET < 50) {
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

    if (dR_JFJ_Max < 0) {
      continue;
    }
    if (dR_J1L <= 0.4 || dR_J2L <= 0.4) {
      continue;
    }
    if (dR_JmaxL > 3.5) {
      continue;
    }

    // VBFTag veto
    //   if(isVBFtag) continue;

    // Lepton Selection or Veto
    //   if (fabs(lep1_Id) !=11 ) continue;

    // Trigger Objects and Matchings
    // Probe Matched
    bool Probe_Matched = false;

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

    if (matched_to_AK8PFJet230_SoftDropMass40) {
      Probe_Matched = true;
    }

    // weight
    weight = (weight / SumGenWeights) * xsec * param_dict.Lumi;
    double PU_weight = PU_Rew[(int)npu];
    weight = weight * PU_weight;

    // Fill histograms
    _FatJet1_probe_pt->Fill(FatJet1_pt, weight);
    _FatJet1_probe_eta->Fill(FatJet1_eta, weight);
    _FatJet1_probe_phi->Fill(FatJet1_phi, weight);
    _FatJet1_probe_eta_phi->Fill(FatJet1_eta, FatJet1_phi, weight);
    _FatJet1_probe_Mass->Fill(FatJet1_Mass, weight);
    _FatJet1_probe_MassSD->Fill(FatJet1_MassSD, weight);
    _FatJet1_probe_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD, weight);

    if (Probe_Matched) {
      _FatJet1_tag_pt->Fill(FatJet1_pt, weight);
      _FatJet1_tag_eta->Fill(FatJet1_eta, weight);
      _FatJet1_tag_phi->Fill(FatJet1_phi, weight);
      _FatJet1_tag_eta_phi->Fill(FatJet1_eta, FatJet1_phi, weight);
      _FatJet1_tag_Mass->Fill(FatJet1_Mass, weight);
      _FatJet1_tag_MassSD->Fill(FatJet1_MassSD, weight);
      _FatJet1_tag_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD, weight);
    }

  } // end event loop

  f->Write();
  
  std::cout << "Done. Written to " << output_path << std::endl;
}
