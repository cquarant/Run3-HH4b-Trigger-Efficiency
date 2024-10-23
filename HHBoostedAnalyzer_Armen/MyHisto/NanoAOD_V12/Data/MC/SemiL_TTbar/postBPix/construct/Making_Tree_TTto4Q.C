#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include "TCanvas.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <vector>
#include <map>
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// D-phi
 double phi_dist(double a, double b){
   if(fabs(a - b) > 3.14159265)
   {
    return 6.2831853 - fabs(a - b);
   }
   return fabs(a - b);
 }

bool inRange(int low, int high, int x)
{
    return (low <= x && x <= high);
}

double PU_Rew[100] = {0.750156,1.89753,1.16218,1.17185,1.24144,1.24093,1.21528,1.24569,1.22439,1.13186,1.268,1.22833,0.898481,1.72457,2.50901,3.71031,3.94382,2.84716,2.23894,2.00194,1.94304,1.9546,2.13721,2.41457,2.55297,2.55506,2.44414,2.26349,2.049,1.82273,1.69234,1.5774,1.54434,1.48068,1.44632,1.38943,1.33936,1.29638,1.26625,1.26429,1.26815,1.27409,1.31266,1.35538,1.42214,1.49389,1.55454,1.61395,1.63681,1.62573,1.59832,1.54196,1.43462,1.30671,1.14609,0.97302,0.820247,0.689621,0.572241,0.473172,0.392613,0.32013,0.260165,0.210414,0.182732,0.14334,0.116381,0.0954352,0.0793652,0.0667483,0.0549356,0.0456514,0.0380612,0.0318324,0.0281388,0.0241728,0.0210379,0.0192466,0.0180967,0.0183505,0.0211177,0.0249506,0.023001,0.022844,0.0209544,0.018289,0.0155916,0.0123924,0.0101449,0.010468,0.006118,0.0044397,0.00357126,0.00463679,0.0842234,0.0488443,0.0560634,1,1,1};
// ********************** PT-Mass-SFs **************************
  TFile* f_PT_Mass_SF       = new TFile("/afs/cern.ch/user/t/tumasyan/public/2023/Trigger_SFs/PT_Mass_SF_TTbar/PT_Mass_2dSF_PostBPix.root");
  TH2D * _Eff_Data    = (TH2D*)f_PT_Mass_SF->Get("Eff_Data_ETA0");
  TH2D * _Eff_MC      = (TH2D*)f_PT_Mass_SF->Get("Eff_MC_ETA0");
// ********************** BTG-SFs **************************
  TFile* f_BTG_SF       = new TFile("/afs/cern.ch/user/t/tumasyan/public/BTG_SFs_TXbb/TTbar/2023/BTG_SF_TXbb_2023_TTbar.root");
  TH1D * _BTG_Eff_Data  = (TH1D*)f_BTG_SF->Get("BTG_Eff_Data");
  TH1D * _BTG_Eff_MC    = (TH1D*)f_BTG_SF->Get("BTG_Eff_MC");
// ********************************************************************************

#include "/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/MC/parameters_PostBPix.txt"

void Making_Tree_TTto4Q()
{
gSystem->Load("libFWCoreFWLite.so");

// ********************************************************* JEC

vector<JetCorrectorParameters> vPar;
vPar.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23BPixPrompt23_V1_MC/Summer23BPixPrompt23_V1_MC_L2Relative_AK4PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK4 = new FactorizedJetCorrector(vPar);

vector<JetCorrectorParameters> vParAK8;
vParAK8.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23BPixPrompt23_V1_MC/Summer23BPixPrompt23_V1_MC_L2Relative_AK8PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK8 = new FactorizedJetCorrector(vParAK8);
/*
// ********************************************************* JER   (should be added when available)
std::string resptstr =    "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_PtResolution_AK4PFPuppi.txt";
std::string resptstr_sf = "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_SF_AK4PFPuppi.txt";
JME::JetResolution resolution_pt = JME::JetResolution(resptstr.c_str());
JME::JetResolutionScaleFactor resolution_pt_sf = JME::JetResolutionScaleFactor(resptstr_sf.c_str());

std::string resptstr_AK8 =    "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_PtResolution_AK8PFPuppi.txt";
std::string resptstr_sf_AK8 = "/afs/cern.ch/user/t/tumasyan/public/OldJEC_2022/JR_Winter22Run3_V1_MC/JR_Winter22Run3_V1_MC_SF_AK8PFPuppi.txt";
JME::JetResolution resolution_pt_AK8 = JME::JetResolution(resptstr_AK8.c_str());
JME::JetResolutionScaleFactor resolution_pt_sf_AK8 = JME::JetResolutionScaleFactor(resptstr_sf_AK8.c_str());

// *********************************************************
*/

 TFile *f = new TFile("Histograms_TTto4Q.root","RECREATE");


TTree *outputTree = new TTree("tree", "");

Float_t T_weight;
Float_t T_fatJet1_pt, T_fatJet1_eta, T_fatJet1_phi, T_fatJet1_msoftdrop, T_fatJet1_particleNet_XbbVsQCD, T_fatJet1_Tau3OverTau2, T_fatJet1_particleNet_Xbb_Legacy, T_fatJet1_particleNet_Xbb_Legacy_AN;
Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_msoftdrop;
Float_t T_MET, T_Lep1_Pt, T_Lep1_Eta, T_Lep1_Phi, T_Dr_LFJ, T_Dr_J1FJ, T_Dr_J2FJ, T_Dr_JmaxL;

outputTree->Branch("T_weight",      &T_weight,      "T_weight/F");
outputTree->Branch("T_fatJet1_pt",  &T_fatJet1_pt,  "T_fatJet1_pt/F");
outputTree->Branch("T_fatJet1_eta", &T_fatJet1_eta, "T_fatJet1_eta/F");
outputTree->Branch("T_fatJet1_phi", &T_fatJet1_phi, "T_fatJet1_phi/F");
outputTree->Branch("T_fatJet1_msoftdrop", &T_fatJet1_msoftdrop, "T_fatJet1_msoftdrop/F");
outputTree->Branch("T_fatJet1_particleNet_XbbVsQCD", &T_fatJet1_particleNet_XbbVsQCD, "T_fatJet1_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet1_Tau3OverTau2", &T_fatJet1_Tau3OverTau2, "T_fatJet1_Tau3OverTau2/F");
outputTree->Branch("T_fatJet1_particleNet_Xbb_Legacy", &T_fatJet1_particleNet_Xbb_Legacy, "T_fatJet1_particleNet_Xbb_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_Xbb_Legacy_AN", &T_fatJet1_particleNet_Xbb_Legacy_AN, "T_fatJet1_particleNet_Xbb_Legacy_AN/F");
outputTree->Branch("T_fatJet2_pt",  &T_fatJet2_pt,  "T_fatJet2_pt/F");
outputTree->Branch("T_fatJet2_eta", &T_fatJet2_eta, "T_fatJet2_eta/F");
outputTree->Branch("T_fatJet2_msoftdrop", &T_fatJet2_msoftdrop, "T_fatJet2_msoftdrop/F");
outputTree->Branch("T_MET",  &T_MET,  "T_MET/F");
outputTree->Branch("T_Lep1_Pt",  &T_Lep1_Pt,  "T_Lep1_Pt/F");
outputTree->Branch("T_Lep1_Eta",  &T_Lep1_Eta,  "T_Lep1_Eta/F");
outputTree->Branch("T_Lep1_Phi",  &T_Lep1_Phi,  "T_Lep1_Phi/F");
outputTree->Branch("T_Dr_LFJ",  &T_Dr_LFJ,  "T_Dr_LFJ/F");
outputTree->Branch("T_Dr_J1FJ",  &T_Dr_J1FJ,  "T_Dr_J1FJ/F");
outputTree->Branch("T_Dr_J2FJ",  &T_Dr_J2FJ,  "T_Dr_J2FJ/F");
outputTree->Branch("T_Dr_JmaxL",  &T_Dr_JmaxL,  "T_Dr_JmaxL/F");

TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/Legacy/PostBPix/TTto4Q.root");

TH1F  *NEvents = (TH1F*)f1->Get("NEvents");
double SumGenWeights = NEvents->GetBinContent(1);
TTree *InputTree = (TTree*)f1->Get("tree");

Float_t    weight;
UInt_t     run;
UInt_t     lumi;
Float_t    npu;
Int_t      isVBFtag;
Float_t    rho;

Bool_t     HLT_Ele32_WPTight_Gsf;
Bool_t     HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

Bool_t     HLT_Mu50;
Bool_t     HLT_IsoMu27;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06;

Bool_t     HLT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

Float_t    MET;

Float_t    lep1_Pt;
Float_t    lep1_Eta;
Float_t    lep1_Phi;
Int_t      lep1_Id;
Float_t    lep2_Pt;
Float_t    lep2_Eta;
Float_t    lep2_Phi;
Int_t      lep2_Id;

Float_t    Jet1_Pt;
Float_t    Jet1_Eta;
Float_t    Jet1_Phi;

Float_t    Jet2_Pt;
Float_t    Jet2_Eta;
Float_t    Jet2_Phi;

Float_t    FatJet1_pt;
Float_t    FatJet1_eta;
Float_t    FatJet1_phi;
Float_t    FatJet1_Mass;
Float_t    FatJet1_MassSD;
Float_t    FatJet1PNetMD_Xbb;
Float_t    FatJet1_HbbvsQCD;
Float_t    FatJet1_Tau3OverTau2;
Float_t    FatJet1_rawFactor;
Float_t    FatJet1PNetMD_Xbb_Legacy;
Float_t    FatJet1PNetMD_QCD_Legacy;

Float_t    FatJet2_pt;
Float_t    FatJet2_eta;
Float_t    FatJet2_phi;
Float_t    FatJet2_Mass;
Float_t    FatJet2_MassSD;
Float_t    FatJet2PNetMD_Xbb;
Float_t    FatJet2_Tau3OverTau2;
Float_t    FatJet2_rawFactor;
Float_t    FatJet2PNetMD_Xbb_Legacy;
Float_t    FatJet2PNetMD_QCD_Legacy;

Float_t    FatJet3_pt;
Float_t    FatJet3_rawFactor;

Int_t nGenJet;
Float_t GenJet_eta[20];
Float_t GenJet_phi[20];
Float_t GenJet_pt[20];
Int_t nGenJetAK8;
Float_t GenJetAK8_eta[20];
Float_t GenJetAK8_phi[20];
Float_t GenJetAK8_pt[20];

InputTree->SetBranchAddress("weight",&weight);
InputTree->SetBranchAddress("run",&run);
InputTree->SetBranchAddress("lumi",&lumi);
InputTree->SetBranchAddress("npu",&npu);
InputTree->SetBranchAddress("rho",&rho);
InputTree->SetBranchAddress("isVBFtag",&isVBFtag);

InputTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
InputTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("lep1Pt",&lep1_Pt);
InputTree->SetBranchAddress("lep1Eta",&lep1_Eta);
InputTree->SetBranchAddress("lep1Phi",&lep1_Phi);
InputTree->SetBranchAddress("lep1Id",&lep1_Id);
InputTree->SetBranchAddress("lep2Pt",&lep2_Pt);
InputTree->SetBranchAddress("lep2Eta",&lep2_Eta);
InputTree->SetBranchAddress("lep2Phi",&lep2_Phi);
InputTree->SetBranchAddress("lep2Id",&lep2_Id);
InputTree->SetBranchAddress("MET",&MET);

InputTree->SetBranchAddress("jet1Pt",&Jet1_Pt);
InputTree->SetBranchAddress("jet1Eta",&Jet1_Eta);
InputTree->SetBranchAddress("jet1Phi",&Jet1_Phi);
InputTree->SetBranchAddress("jet2Pt",&Jet2_Pt);
InputTree->SetBranchAddress("jet2Eta",&Jet2_Eta);
InputTree->SetBranchAddress("jet2Phi",&Jet2_Phi);

InputTree->SetBranchAddress("fatJet1_pt",&FatJet1_pt);
InputTree->SetBranchAddress("fatJet1_eta",&FatJet1_eta);
InputTree->SetBranchAddress("fatJet1_phi",&FatJet1_phi);
InputTree->SetBranchAddress("fatJet1_mass",&FatJet1_Mass);
InputTree->SetBranchAddress("fatJet1_msoftdrop",&FatJet1_MassSD);
InputTree->SetBranchAddress("fatJet1_rawFactor",&FatJet1_rawFactor);
InputTree->SetBranchAddress("fatJet1_Tau3OverTau2",&FatJet1_Tau3OverTau2);
InputTree->SetBranchAddress("fatJet1_particleNet_XbbVsQCD",&FatJet1PNetMD_Xbb);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xbb",&FatJet1PNetMD_Xbb_Legacy);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCD",&FatJet1PNetMD_QCD_Legacy);

InputTree->SetBranchAddress("fatJet2_pt",&FatJet2_pt);
InputTree->SetBranchAddress("fatJet2_eta",&FatJet2_eta);
InputTree->SetBranchAddress("fatJet2_phi",&FatJet2_phi);
InputTree->SetBranchAddress("fatJet2_mass",&FatJet2_Mass);
InputTree->SetBranchAddress("fatJet2_msoftdrop",&FatJet2_MassSD);
InputTree->SetBranchAddress("fatJet2_rawFactor",&FatJet2_rawFactor);
InputTree->SetBranchAddress("fatJet2_Tau3OverTau2",&FatJet2_Tau3OverTau2);
InputTree->SetBranchAddress("fatJet2_particleNet_XbbVsQCD",&FatJet2PNetMD_Xbb);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xbb",&FatJet2PNetMD_Xbb_Legacy);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCD",&FatJet2PNetMD_QCD_Legacy);

InputTree->SetBranchAddress("fatJet3_pt",&FatJet3_pt);
InputTree->SetBranchAddress("fatJet3_rawFactor",&FatJet3_rawFactor);

InputTree->SetBranchAddress("nGenJet",&nGenJet);
InputTree->SetBranchAddress("GenJet_eta",GenJet_eta);
InputTree->SetBranchAddress("GenJet_phi",GenJet_phi);
InputTree->SetBranchAddress("GenJet_pt", GenJet_pt);
InputTree->SetBranchAddress("nGenJetAK8",&nGenJetAK8);
InputTree->SetBranchAddress("GenJetAK8_eta",GenJetAK8_eta);
InputTree->SetBranchAddress("GenJetAK8_phi",GenJetAK8_phi);
InputTree->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt);

  // Events Loop
  for(int i=0;i<InputTree->GetEntries();i++)
   {
    InputTree->GetEntry(i);
// ********************************************************** HLT Selection
   if ( ! ( (HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) ==11) || (HLT_IsoMu27 && fabs(lep1_Id) ==13) ) )continue;

// ********************************************************** FatJets correction and selection
   if(FatJet1_pt > 0)
     {
      double Raw_FatJet1_pt = FatJet1_pt*(1.0 - FatJet1_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet1_pt);
      corrector_AK8->setJetEta(FatJet1_eta);
      corrector_AK8->setJetPhi(FatJet1_phi);
      double This_correction = corrector_AK8->getCorrection();
      FatJet1_pt     = Raw_FatJet1_pt * This_correction;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
     }
   if(FatJet2_pt > 0)
     {
      double Raw_FatJet2_pt = FatJet2_pt*(1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet2_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      corrector_AK8->setJetPhi(FatJet2_phi);
      double This_correction = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt* This_correction;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
     }

// ********************************************************** Jets / FatJets Smearing
/*
   // Jet Smearing
   double res_pt_1;
   double res_pt_sf_1;
   JME::JetParameters JerPARAM_1 = {{JME::Binning::JetPt, FatJet1_pt}, {JME::Binning::JetEta, FatJet1_eta},{JME::Binning::Rho, rho}};
   JME::JetParameters JerSFPARAM_1;
   JerSFPARAM_1.set(JME::Binning::JetPt,  FatJet1_pt);
   JerSFPARAM_1.set(JME::Binning::JetEta, FatJet1_eta);
   JerSFPARAM_1.set(JME::Binning::Rho, rho);
   res_pt_1 = resolution_pt_AK8.getResolution(JerPARAM_1);
   res_pt_sf_1 = resolution_pt_sf_AK8.getScaleFactor(JerSFPARAM_1);

   double res_pt_2;
   double res_pt_sf_2;
   JME::JetParameters JerPARAM_2 = {{JME::Binning::JetPt, FatJet2_pt}, {JME::Binning::JetEta, FatJet2_eta},{JME::Binning::Rho, rho}};
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
       if(!GenJetMatched_1 &&  sqrt(pow(FatJet1_eta-GenJetAK8_eta[nGJAK8],2) + pow(phi_dist(FatJet1_phi,GenJetAK8_phi[nGJAK8]),2)) < 0.2 && (fabs(FatJet1_pt - GenJetAK8_pt[nGJAK8])/FatJet1_pt < 3*res_pt_1) )
         {
           SmearFactor_1 = 1.0 + (res_pt_sf_1 - 1.0) * (FatJet1_pt - GenJetAK8_pt[nGJAK8]) / FatJet1_pt;
           GenJetMatched_1 = true;
         }
       if(!GenJetMatched_2 &&  sqrt(pow(FatJet2_eta-GenJetAK8_eta[nGJAK8],2) + pow(phi_dist(FatJet2_phi,GenJetAK8_phi[nGJAK8]),2)) < 0.2 && (fabs(FatJet2_pt - GenJetAK8_pt[nGJAK8])/FatJet2_pt < 3*res_pt_2) )
         {
           SmearFactor_2 = 1.0 + (res_pt_sf_2 - 1.0) * (FatJet2_pt - GenJetAK8_pt[nGJAK8]) / FatJet2_pt;
           GenJetMatched_2 = true;
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

    if(FatJet1_pt < 250 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50) continue;

// ********************************************************** VBFTag veto
//   if(isVBFtag) continue;

// ********************************************************** Lepton Selection or Veto

   if (lep1_Pt < 50) continue;
   if (lep2_Pt > 30) continue;
   if (FatJet2_pt > 200 && FatJet2_MassSD > 50) continue;
   double Dr_LFJ = sqrt ( pow(lep1_Eta - FatJet1_eta,2) + pow(phi_dist(lep1_Phi,FatJet1_phi),2) );
   if (Dr_LFJ < 1.5) continue;
//   if (FatJet1Tau3OverTau2 > 0.5) continue;
   if (MET < 50) continue;

   double Dr_J1FJ =-1;
   double Dr_J1L =10;
   if(Jet1_Pt > 40)
     {
      Dr_J1FJ = sqrt ( pow(Jet1_Eta - FatJet1_eta,2) + pow(phi_dist(Jet1_Phi,FatJet1_phi),2) );
      Dr_J1L  = sqrt ( pow(Jet1_Eta - lep1_Eta,2) + pow(phi_dist(Jet1_Phi,lep1_Phi),2) );
     }
   double Dr_J2FJ =-1;
   double Dr_J2L =10;
   if(Jet2_Pt > 40)
    {
      Dr_J2FJ = sqrt ( pow(Jet2_Eta - FatJet1_eta,2) + pow(phi_dist(Jet2_Phi,FatJet1_phi),2) );
      Dr_J2L  = sqrt ( pow(Jet2_Eta - lep1_Eta,2) + pow(phi_dist(Jet2_Phi,lep1_Phi),2) );
    }

   double Dr_JFJ_Max = Dr_J1FJ;
   double Dr_JFJ_Min = Dr_J2FJ;
   double Dr_JmaxL   = Dr_J1L;
   if(Dr_J2FJ > Dr_J1FJ)
     {
      Dr_JFJ_Max = Dr_J2FJ;
      Dr_JFJ_Min = Dr_J1FJ;
      Dr_JmaxL   = Dr_J2L;
     }

//  if(Dr_JFJ_Max < 1.5) continue;
  if( Dr_J1L <= 0.4 || Dr_J2L <= 0.4 ) continue;
//  if(Dr_JmaxL > 3.5 ) continue;

// ********************************************************** weight


   weight = (weight/SumGenWeights)*XSec_TTto4Q*Lumi;
   double PU_weight=PU_Rew[(int)npu];
   if (PU_weight<20.0)
      weight = weight*PU_weight;

// ********************************************************** Fill Histograms

  T_weight = weight;
  T_fatJet1_pt = FatJet1_pt;
  T_fatJet1_eta = FatJet1_eta;
  T_fatJet1_phi = FatJet1_phi;
  T_fatJet1_msoftdrop = FatJet1_MassSD;
  T_fatJet1_particleNet_XbbVsQCD = FatJet1PNetMD_Xbb;
  T_fatJet1_Tau3OverTau2 = FatJet1_Tau3OverTau2;
  T_fatJet1_particleNet_Xbb_Legacy = FatJet1PNetMD_Xbb_Legacy;
  double FatJet1PNetMD_Xbb_Legacy_AN = FatJet1PNetMD_Xbb_Legacy/(FatJet1PNetMD_Xbb_Legacy + FatJet1PNetMD_QCD_Legacy);
  T_fatJet1_particleNet_Xbb_Legacy_AN = FatJet1PNetMD_Xbb_Legacy_AN;
  T_fatJet2_pt = FatJet2_pt;
  T_fatJet2_eta = FatJet2_eta;
  T_fatJet2_msoftdrop = FatJet2_MassSD;

  T_MET      = MET;
  T_Lep1_Pt  = lep1_Pt;
  T_Lep1_Eta = lep1_Eta;
  T_Lep1_Phi = lep1_Phi;
  T_Dr_LFJ   = Dr_LFJ;
  T_Dr_J1FJ  = Dr_JFJ_Max;
  T_Dr_J2FJ  = Dr_JFJ_Min;
  T_Dr_JmaxL = Dr_JmaxL;

  outputTree->Fill();

 } // end event loop

f->Write();

}

