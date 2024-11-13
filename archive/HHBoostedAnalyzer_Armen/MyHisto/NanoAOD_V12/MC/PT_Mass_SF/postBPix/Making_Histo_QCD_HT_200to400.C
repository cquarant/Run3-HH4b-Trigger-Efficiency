#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
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
#include <cstdlib> // for std::getenv
#include <iostream>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <vector>

std::string getCMSSWBase() {
  const char *cmssw_base = std::getenv("CMSSW_BASE");
  if (!cmssw_base) {
    throw std::runtime_error("CMSSW_BASE environment variable not set! Did you "
                             "forget to run 'cmsenv'?");
  }
  return std::string(cmssw_base);
}

// D-phi
double phi_dist(double a, double b) {
  if (fabs(a - b) > 3.14159265) {
    return 6.2831853 - fabs(a - b);
  }
  return fabs(a - b);
}

bool inRange(int low, int high, int x) { return (low <= x && x <= high); }

double PU_Rew[100] = {
    0.750156,  1.89753,   1.16218,    1.17185,    1.24144,   1.24093,
    1.21528,   1.24569,   1.22439,    1.13186,    1.268,     1.22833,
    0.898481,  1.72457,   2.50901,    3.71031,    3.94382,   2.84716,
    2.23894,   2.00194,   1.94304,    1.9546,     2.13721,   2.41457,
    2.55297,   2.55506,   2.44414,    2.26349,    2.049,     1.82273,
    1.69234,   1.5774,    1.54434,    1.48068,    1.44632,   1.38943,
    1.33936,   1.29638,   1.26625,    1.26429,    1.26815,   1.27409,
    1.31266,   1.35538,   1.42214,    1.49389,    1.55454,   1.61395,
    1.63681,   1.62573,   1.59832,    1.54196,    1.43462,   1.30671,
    1.14609,   0.97302,   0.820247,   0.689621,   0.572241,  0.473172,
    0.392613,  0.32013,   0.260165,   0.210414,   0.182732,  0.14334,
    0.116381,  0.0954352, 0.0793652,  0.0667483,  0.0549356, 0.0456514,
    0.0380612, 0.0318324, 0.0281388,  0.0241728,  0.0210379, 0.0192466,
    0.0180967, 0.0183505, 0.0211177,  0.0249506,  0.023001,  0.022844,
    0.0209544, 0.018289,  0.0155916,  0.0123924,  0.0101449, 0.010468,
    0.006118,  0.0044397, 0.00357126, 0.00463679, 0.0842234, 0.0488443,
    0.0560634, 1,         1,          1};

// ******************************************

// #include "/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/MC/parameters_PostBPix.txt"
#include "HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/MC/parameters_PostBPix.txt"

void Making_Histo_QCD_HT_200to400() {
  gSystem->Load("libFWCoreFWLite.so");

  // ********************************************************* JEC
  std::string jec_path = getCMSSWBase() + "/src/JECs/";
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(JetCorrectorParameters(
      // "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      jec_path + "Summer23BPixPrompt23_V1_MC/Summer23BPixPrompt23_V1_MC_L2Relative_AK4PFPuppi.txt"));
  FactorizedJetCorrector *corrector = new FactorizedJetCorrector(vPar);

  vector<JetCorrectorParameters> vParAK8;
  vParAK8.push_back(JetCorrectorParameters(
      // "/afs/cern.ch/user/t/tumasyan/public/2023/JECs/"
      jec_path + "Summer23BPixPrompt23_V1_MC/Summer23BPixPrompt23_V1_MC_L2Relative_AK8PFPuppi.txt"));
  FactorizedJetCorrector *corrector_AK8 = new FactorizedJetCorrector(vParAK8);
  /*
  // ********************************************************* JER   (should be
  added when available) std::string resptstr =
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

  // *********************************************************
  */
  TFile *f = new TFile("Histograms_QCD_HT_200to400.root", "RECREATE");

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
  TH2D *_FatJet2_Pt_Mass_M = new TH2D("FatJet2_Pt_Mass_M", "FatJet2_Pt_Mass_M",
                                      45, Lower_pt, 15, Lower_m);

  Float_t Lower_pt_N[9] = {230, 240, 250, 270, 300, 350, 500, 700, 1000};
  TH2D *_FatJet2_Pt_MassN = new TH2D("FatJet2_Pt_MassN", "FatJet2_Pt_MassN", 8,
                                     Lower_pt_N, 15, Lower_m);
  TH2D *_FatJet2_Pt_MassN_M = new TH2D(
      "FatJet2_Pt_MassN_M", "FatJet2_Pt_MassN_M", 8, Lower_pt_N, 15, Lower_m);
  // TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/PostBPix/QCD_HT_200to400.root");
  TFile *f1 = new TFile("/eos/cms/store/group/phys_higgs/nonresonant_HH/bbbb/sixie/Run3Analysis/HH/HHTo4BNtupler/ArmenVersion_ICHEP2024/Data_2023/PostBPix/QCD_HT_200to400.root");

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

  Float_t FatJet3_pt;
  Float_t FatJet3_rawFactor;

  Int_t nGenJet;
  Float_t GenJet_eta[20];
  Float_t GenJet_phi[20];
  Float_t GenJet_pt[20];
  Int_t nGenJetAK8;
  Float_t GenJetAK8_eta[20];
  Float_t GenJetAK8_phi[20];
  Float_t GenJetAK8_pt[20];

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
  Float_t Trigger_Object_pt[20];
  Float_t Trigger_Object_eta[20];
  Float_t Trigger_Object_phi[20];
  Int_t Trigger_Object_bit[20];
  InputTree_TrgObj->SetBranchAddress("NTrigger_Objects", &NTrigger_Objects);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_pt", Trigger_Object_pt);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_eta", Trigger_Object_eta);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_phi", Trigger_Object_phi);
  InputTree_TrgObj->SetBranchAddress("Trigger_Object_bit", Trigger_Object_bit);

  // Events Loop
  for (int i = 0; i < InputTree->GetEntries(); i++) {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);

    // ********************************************************** HLT Selection
    //   if(HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06==0) continue;
    //   if(HLT_Ele32_WPTight_Gsf==0) continue;
    if (HLT_AK8PFJet230_SoftDropMass40 == 0)
      continue;

    // ********************************************************** FatJets
    // correction and selection
    if (FatJet1_pt > 0) {
      double Raw_FatJet1_pt = FatJet1_pt * (1.0 - FatJet1_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_pt);
      corrector_AK8->setJetEta(FatJet1_eta);
      corrector_AK8->setJetPhi(FatJet1_phi);
      double This_correction = corrector_AK8->getCorrection();
      FatJet1_pt = Raw_FatJet1_pt * This_correction;
      FatJet1_MassSD =
          FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
    }
    if (FatJet2_pt > 0) {
      double Raw_FatJet2_pt = FatJet2_pt * (1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet2_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      corrector_AK8->setJetPhi(FatJet2_phi);
      double This_correction = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt * This_correction;
      FatJet2_MassSD =
          FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
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
    // ********************************************************** FatJets
    // Selection
    if (FatJet3_pt > 150)
      continue;
    if (FatJet1_pt < 300 || fabs(FatJet1_eta) > 2.5 || FatJet1_MassSD < 80)
      continue;
    if (FatJet2_pt < 160)
      continue;
    if (phi_dist(FatJet1_phi, FatJet2_phi) < 2.5)
      continue;

    // ********************************************************** VBFTag veto
    //   if(isVBFtag) continue;

    // ********************************************************** Lepton
    // Selection or Veto
    //   if (fabs(lep1_Id) !=11 ) continue;

    // ********************************************************** Trigger
    // Objects and Matchings
    bool Tag_Matched = false;

    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if ((Trigger_Object_bit[itrg] & 4) == 4) {
        if (sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet1_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 100) {
          Tag_Matched = true;
          break;
        } else
          break;
      }
    if (!Tag_Matched)
      continue;

    // Probe Matched
    bool Probe_Matched = false;

    bool matched_to_AK8PFJet230_SoftDropMass40 = false;
    for (int itrg = 0; itrg < NTrigger_Objects; itrg++)
      if ((Trigger_Object_bit[itrg] & 4) == 4)
        if (sqrt(pow((FatJet2_eta - Trigger_Object_eta[itrg]), 2) +
                 pow(phi_dist(FatJet2_phi, Trigger_Object_phi[itrg]), 2)) <
                0.4 &&
            Trigger_Object_pt[itrg] > 100) {
          matched_to_AK8PFJet230_SoftDropMass40 = true;
          break;
        }

    if (matched_to_AK8PFJet230_SoftDropMass40)
      Probe_Matched = true;

    // ********************************************************** weight

    weight = (weight / SumGenWeights) * XSec_QCD_HT_200to400 * Lumi;
    double PU_weight = PU_Rew[(int)npu];
    weight = weight * PU_weight;

    // ********************************************************** Fill
    // Histograms

    _FatJet1_pt->Fill(FatJet1_pt, weight);
    _FatJet1_eta->Fill(FatJet1_eta, weight);
    _FatJet1_phi->Fill(FatJet1_phi, weight);
    _FatJet1_eta_phi->Fill(FatJet1_eta, FatJet1_phi, weight);
    _FatJet1_Mass->Fill(FatJet1_Mass, weight);
    _FatJet1_MassSD->Fill(FatJet1_MassSD, weight);
    _FatJet1_Pt_Mass->Fill(FatJet1_pt, FatJet1_MassSD, weight);

    _FatJet2_pt->Fill(FatJet2_pt, weight);
    _FatJet2_eta->Fill(FatJet2_eta, weight);
    _FatJet2_phi->Fill(FatJet2_phi, weight);
    _FatJet2_eta_phi->Fill(FatJet2_eta, FatJet2_phi, weight);
    _FatJet2_Mass->Fill(FatJet2_Mass, weight);
    _FatJet2_MassSD->Fill(FatJet2_MassSD, weight);
    _FatJet2_Pt_Mass->Fill(FatJet2_pt, FatJet2_MassSD, weight);
    _FatJet2_Pt_MassN->Fill(FatJet2_pt, FatJet2_MassSD, weight);

    if (Probe_Matched) {
      _FatJet2_Pt_Mass_M->Fill(FatJet2_pt, FatJet2_MassSD, weight);
      _FatJet2_Pt_MassN_M->Fill(FatJet2_pt, FatJet2_MassSD, weight);
    }

  } // end event loop

  f->Write();
}
