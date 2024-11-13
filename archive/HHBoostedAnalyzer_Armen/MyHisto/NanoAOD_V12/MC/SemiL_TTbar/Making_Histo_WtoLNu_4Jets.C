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

double PU_Rew[100] = {11.7389,9.08407,26.3901,42.8357,33.3515,18.6839,10.0213,5.11508,2.75789,1.62045,0.941633,0.579944,0.359025,0.229839,0.179189,0.179606,0.204167,0.233829,0.256693,0.270762,0.280817,0.293772,0.320734,0.371096,0.439664,0.508281,0.562052,0.59908,0.625454,0.646593,0.664557,0.679711,0.692435,0.703867,0.716255,0.732696,0.756506,0.790735,0.838193,0.901237,0.98116,1.07746,1.18749,1.30652,1.42829,1.54571,1.65169,1.73976,1.80481,1.84397,1.85694,1.84569,1.81444,1.76923,1.71722,1.66611,1.62358,1.59715,1.59424,1.62261,1.69115,1.81108,1.99775,2.27346,2.67174,3.24415,4.07131,5.28127,7.08028,9.80561,14.0171,20.6588,31.3499,48.9171,78.3945,128.949,217.711,377.687,674.814,1246.58,2394.24,4815.29,10222.6,23080.8,55715.2,144017,397685,1.16693e+06,3.61488e+06,1.17473e+07,3.98386e+07,1.4044e+08,5.13248e+08,1.94111e+09,7.58931e+09,3.06568e+10,1.27909e+11,5.5117e+11,2.45291e+12,1.12751e+13};

// ********************** PT-Mass-SFs **************************
  TFile* f_PT_Mass_SF       = new TFile("/afs/cern.ch/work/t/tumasyan/HHTo4B/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/MC/Trigger_SFs_Roots/PT_Mass_2dSF_TTbar.root");
  TH2D * _Eff_Data    = (TH2D*)f_PT_Mass_SF->Get("Eff_Data_PostEE");
  TH2D * _Eff_MC      = (TH2D*)f_PT_Mass_SF->Get("Eff_MC_PostEE");
// ********************** BTG-SFs **************************
  TFile* f_BTG_SF       = new TFile("/afs/cern.ch/work/t/tumasyan/HHTo4B/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/MC/Trigger_SFs_Roots/BTG_SF_TXbb_PostEE_TTbar.root");
  TH1D * _BTG_Eff_Data  = (TH1D*)f_BTG_SF->Get("BTG_Eff_Data");
  TH1D * _BTG_Eff_MC    = (TH1D*)f_BTG_SF->Get("BTG_Eff_MC");
// ********************************************************************************

#include "/afs/cern.ch/work/t/tumasyan/HHTo4B/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/MC/parameters_PostEE.txt"


void Making_Histo_WtoLNu_4Jets()
{
gSystem->Load("libFWCoreFWLite.so");

// ********************************************************* JEC

vector<JetCorrectorParameters> vPar;
vPar.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2022/Summer22EE_22Sep2023_V2_MC/Summer22EE_22Sep2023_V2_MC_L2Relative_AK4PFPuppi.txt"));
FactorizedJetCorrector*  corrector = new FactorizedJetCorrector(vPar);

vector<JetCorrectorParameters> vParAK8;
vParAK8.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2022/Summer22EE_22Sep2023_V2_MC/Summer22EE_22Sep2023_V2_MC_L2Relative_AK8PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK8 = new FactorizedJetCorrector(vParAK8);
/*
// ********************************************************* JER
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
 TFile *f = new TFile("Histograms_WtoLNu_4Jets.root","RECREATE");


TTree *outputTree = new TTree("tree", "");

Float_t T_weight;
Float_t T_fatJet1_pt, T_fatJet1_eta, T_fatJet1_phi, T_fatJet1_msoftdrop, T_fatJet1_particleNet_XbbVsQCD, T_fatJet1_Tau3OverTau2, T_fatJet1_particleNet_Xbb_Legacy, T_fatJet1_particleNet_Xbb_Legacy_AN;
Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_msoftdrop;
Float_t T_MET, T_Lep1_Pt, T_Dr_LFJ, T_Dr_J1FJ, T_Dr_J2FJ, T_Dr_JmaxL;

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
outputTree->Branch("T_Dr_LFJ",  &T_Dr_LFJ,  "T_Dr_LFJ/F");
outputTree->Branch("T_Dr_J1FJ",  &T_Dr_J1FJ,  "T_Dr_J1FJ/F");
outputTree->Branch("T_Dr_J2FJ",  &T_Dr_J2FJ,  "T_Dr_J2FJ/F");
outputTree->Branch("T_Dr_JmaxL",  &T_Dr_JmaxL,  "T_Dr_JmaxL/F");
TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2022/PostEE/WtoLNu_4Jets.root");

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
Bool_t     HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;

Bool_t     HLT_Mu50;
Bool_t     HLT_IsoMu27;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;

Bool_t     HLT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_AK8PFJet425_SoftDropMass40;
Bool_t     HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35;

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
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35);

InputTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
InputTree->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35);

InputTree->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass40", &HLT_AK8PFJet425_SoftDropMass40);
InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35);

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

// Trigger Objects
TTree *InputTree_TrgObj = (TTree*)f1->Get("tree_TrgObj");
Int_t NTrigger_Objects;
Float_t Trigger_Object_pt[20];
Float_t Trigger_Object_eta[20];
Float_t Trigger_Object_phi[20];
Int_t   Trigger_Object_bit[20];
InputTree_TrgObj->SetBranchAddress("NTrigger_Objects",   &NTrigger_Objects);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_pt",   Trigger_Object_pt);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_eta",  Trigger_Object_eta);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_phi",  Trigger_Object_phi);
InputTree_TrgObj->SetBranchAddress("Trigger_Object_bit",  Trigger_Object_bit);

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
      double This_correction = corrector_AK8->getCorrection();
      FatJet1_pt     = Raw_FatJet1_pt * This_correction;
      FatJet1_MassSD = FatJet1_MassSD * (1.0 - FatJet1_rawFactor) * This_correction;
     }
   if(FatJet2_pt > 0)
     {
      double Raw_FatJet2_pt = FatJet2_pt*(1.0 - FatJet2_rawFactor);
      corrector_AK8->setJetPt(Raw_FatJet2_pt);
      corrector_AK8->setJetEta(FatJet2_eta);
      double This_correction = corrector_AK8->getCorrection();
      FatJet2_pt = Raw_FatJet2_pt* This_correction;
      FatJet2_MassSD = FatJet2_MassSD * (1.0 - FatJet2_rawFactor) * This_correction;
     }
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
   if (MET < 50) continue;
//   if (FatJet1Tau3OverTau2 > 0.5) continue;
   if (FatJet1_MassSD < 50) continue;
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

   weight = (weight/SumGenWeights)*XSec_WtoLNu_4Jets*Lumi;
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
  T_Dr_LFJ   = Dr_LFJ;
  T_Dr_J1FJ  = Dr_JFJ_Max;
  T_Dr_J2FJ  = Dr_JFJ_Min;
  T_Dr_JmaxL = Dr_JmaxL;

  outputTree->Fill();

 } // end event loop

f->Write();

}
