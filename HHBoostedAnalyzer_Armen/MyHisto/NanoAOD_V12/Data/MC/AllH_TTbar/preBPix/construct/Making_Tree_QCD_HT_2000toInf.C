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

double PU_Rew[100] = {0.552875,1.03705,1.17684,1.09279,1.17336,1.17501,1.14858,1.11463,1.14098,1.17844,1.16723,1.2445,1.32702,1.52427,1.69805,1.7403,1.75334,1.8447,2.2354,2.77674,3.16751,3.10149,2.93969,2.90199,2.81243,2.735,2.60348,2.37099,2.15616,1.9852,1.84818,1.73102,1.66986,1.55718,1.49352,1.43964,1.39255,1.36398,1.34727,1.34113,1.34788,1.37159,1.40511,1.44573,1.50152,1.55852,1.60198,1.61822,1.60954,1.56008,1.47303,1.3643,1.21708,1.06146,0.908219,0.754511,0.617157,0.500328,0.406442,0.333317,0.276404,0.235751,0.203188,0.174879,0.161244,0.138147,0.123821,0.110808,0.0997031,0.0897595,0.0827855,0.0776173,0.0740599,0.0711833,0.0647596,0.0583041,0.0521346,0.0454088,0.0416007,0.0383446,0.0356971,0.0353548,0.0370903,0.0443115,0.0549878,0.0680183,0.0753875,0.0828112,0.103032,0.141383,0.306853,0.777128,1,1,1,1,1,1,1,1};

// ********************** PT-Mass-SFs **************************
  TFile* f_PT_Mass_SF   = new TFile("/afs/cern.ch/user/t/tumasyan/public/2023/Trigger_SFs/425/PT_Mass_2dSF_425_PreBPix.root");
  TH2D * _Eff_Data      = (TH2D*)f_PT_Mass_SF->Get("Eff_Data");
  TH2D * _Eff_MC        = (TH2D*)f_PT_Mass_SF->Get("Eff_MC");
// ********************************************************************************

#include "/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/MC/parameters_PreBPix.txt"

void Making_Tree_QCD_HT_2000toInf()
{
gSystem->Load("libFWCoreFWLite.so");

// ********************************************************* JEC

vector<JetCorrectorParameters> vPar;
vPar.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK4 = new FactorizedJetCorrector(vPar);

vector<JetCorrectorParameters> vParAK8;
vParAK8.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK8PFPuppi.txt"));
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
 TFile *f = new TFile("Histograms_QCD_HT_2000toInf.root","RECREATE");

TTree *outputTree = new TTree("tree", "");

Float_t T_weight;
Float_t T_fatJet1_pt, T_fatJet1_eta, T_fatJet1_phi, T_fatJet1_msoftdrop, T_fatJet1_particleNet_XbbVsQCD, T_fatJet1_Tau3OverTau2, T_fatJet1_particleNet_Xbb_Legacy, T_fatJet1_particleNet_Xbb_Legacy_AN;
Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_phi, T_fatJet2_msoftdrop, T_fatJet2_particleNet_XbbVsQCD, T_fatJet2_Tau3OverTau2, T_fatJet2_particleNet_Xbb_Legacy, T_fatJet2_particleNet_Xbb_Legacy_AN;
Float_t T_PTjj;

outputTree->Branch("T_weight",                           &T_weight,      "T_weight/F");

outputTree->Branch("T_fatJet1_pt",                        &T_fatJet1_pt,                        "T_fatJet1_pt/F");
outputTree->Branch("T_fatJet1_eta",                       &T_fatJet1_eta,                       "T_fatJet1_eta/F");
outputTree->Branch("T_fatJet1_phi",                       &T_fatJet1_phi,                       "T_fatJet1_phi/F");
outputTree->Branch("T_fatJet1_msoftdrop",                 &T_fatJet1_msoftdrop,                 "T_fatJet1_msoftdrop/F");
outputTree->Branch("T_fatJet1_particleNet_XbbVsQCD",      &T_fatJet1_particleNet_XbbVsQCD,      "T_fatJet1_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet1_Tau3OverTau2",              &T_fatJet1_Tau3OverTau2,              "T_fatJet1_Tau3OverTau2/F");
outputTree->Branch("T_fatJet1_particleNet_Xbb_Legacy",    &T_fatJet1_particleNet_Xbb_Legacy,    "T_fatJet1_particleNet_Xbb_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_Xbb_Legacy_AN", &T_fatJet1_particleNet_Xbb_Legacy_AN, "T_fatJet1_particleNet_Xbb_Legacy_AN/F");

outputTree->Branch("T_fatJet2_pt",                        &T_fatJet2_pt,                        "T_fatJet2_pt/F");
outputTree->Branch("T_fatJet2_eta",                       &T_fatJet2_eta,                       "T_fatJet2_eta/F");
outputTree->Branch("T_fatJet2_phi",                       &T_fatJet2_phi,                       "T_fatJet2_phi/F");
outputTree->Branch("T_fatJet2_msoftdrop",                 &T_fatJet2_msoftdrop,                 "T_fatJet2_msoftdrop/F");
outputTree->Branch("T_fatJet2_particleNet_XbbVsQCD",      &T_fatJet2_particleNet_XbbVsQCD,      "T_fatJet2_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet2_Tau3OverTau2",              &T_fatJet2_Tau3OverTau2,              "T_fatJet2_Tau3OverTau2/F");
outputTree->Branch("T_fatJet2_particleNet_Xbb_Legacy",    &T_fatJet2_particleNet_Xbb_Legacy,    "T_fatJet2_particleNet_Xbb_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_Xbb_Legacy_AN", &T_fatJet2_particleNet_Xbb_Legacy_AN, "T_fatJet2_particleNet_Xbb_Legacy_AN/F");

outputTree->Branch("T_PTjj",                              &T_PTjj,         "T_PTjj/F");

TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/Legacy/PreBPix/QCD_HT_2000toInf.root");

TH1F  *NEvents = (TH1F*)f1->Get("NEvents");
double SumGenWeights = NEvents->GetBinContent(1);
TTree *InputTree = (TTree*)f1->Get("tree");

Float_t    weight;
UInt_t     run;
UInt_t     lumi;
Float_t    npu;
Int_t      isVBFtag;
Float_t    rho;

Bool_t     HLT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_AK8PFJet425_SoftDropMass40;

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

InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass40", &HLT_AK8PFJet425_SoftDropMass40);

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
  for(int i=0;i<InputTree_TrgObj->GetEntries();i++)
   {
    InputTree->GetEntry(i);
    InputTree_TrgObj->GetEntry(i);

// ********************************************************** HLT Selection
   if(HLT_AK8PFJet425_SoftDropMass40==0) continue;

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
      if(FatJet1_pt < 450 || fabs(FatJet1_eta) > 2.4 || FatJet1_MassSD < 50 ) continue;
      if(FatJet2_pt < 450 || fabs(FatJet2_eta) > 2.4 || FatJet2_MassSD < 50 ) continue;

// ********************************************************** VBFTag veto
//   if(isVBFtag) continue;

// ********************************************************** Trigger Objects and Matchings

  // Matching 1st
  bool matched_TRG_1=false;

  bool matched_to_AK8PFJet230_SoftDropMass40_1 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if((Trigger_Object_bit[itrg] & 4) == 4)
       if(sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet1_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 100)
         {matched_to_AK8PFJet230_SoftDropMass40_1 = true; break;}

  bool matched_to_AK8PFJet250_1 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if(Trigger_Object_bit[itrg] == 1)
       if(sqrt(pow((FatJet1_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet1_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 425)
         {matched_to_AK8PFJet250_1 = true; break;}

  if (matched_to_AK8PFJet230_SoftDropMass40_1 && matched_to_AK8PFJet250_1) matched_TRG_1=true;

  // Matching 2nd
  bool matched_TRG_2=false;

  bool matched_to_AK8PFJet230_SoftDropMass40_2 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if((Trigger_Object_bit[itrg] & 4) == 4)
       if(sqrt(pow((FatJet2_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet2_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 100)
         {matched_to_AK8PFJet230_SoftDropMass40_2 = true; break;}

  bool matched_to_AK8PFJet250_2 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if(Trigger_Object_bit[itrg] == 1)
       if(sqrt(pow((FatJet2_eta - Trigger_Object_eta[itrg]),2) + pow(phi_dist(FatJet2_phi,Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 425)
         {matched_to_AK8PFJet250_2 = true; break;}

  if (matched_to_AK8PFJet230_SoftDropMass40_2 && matched_to_AK8PFJet250_2) matched_TRG_2=true;

  if (!(matched_TRG_1 || matched_TRG_2)) continue;

// ********************************************************** weight

   weight = (weight/SumGenWeights)*XSec_QCD_HT_2000toInf*Lumi;
   double PU_weight=PU_Rew[(int)npu];
   if (PU_weight<20.0)
      weight = weight*PU_weight;

   // Add Trigger PT_Mass-Scale Factor
   double Eff_Data_1 = 0;
   double Eff_MC_1   = 0;
   if(matched_TRG_1)
     {
      Int_t bin_PT_1     = _Eff_Data->GetXaxis()->FindBin(FatJet1_pt);
      Int_t bin_Mass_1   = _Eff_Data->GetYaxis()->FindBin(FatJet1_MassSD);
      Eff_Data_1 = 1.0;
      if(_Eff_Data->GetBinContent(bin_PT_1,bin_Mass_1) > 0)
         Eff_Data_1      = _Eff_Data->GetBinContent(bin_PT_1,bin_Mass_1);
      Eff_MC_1 =1.0;
      if(_Eff_MC->GetBinContent(bin_PT_1,bin_Mass_1) >0)
         Eff_MC_1        = _Eff_MC->GetBinContent(bin_PT_1,bin_Mass_1);
     }

   double Eff_Data_2 = 0;
   double Eff_MC_2   = 0;
   if(matched_TRG_2)
     {
      Int_t bin_PT_2     = _Eff_Data->GetXaxis()->FindBin(FatJet2_pt);
      Int_t bin_Mass_2   = _Eff_Data->GetYaxis()->FindBin(FatJet2_MassSD);

      Eff_Data_2 = 1.0;
      if(_Eff_Data->GetBinContent(bin_PT_2,bin_Mass_2) > 0)
         Eff_Data_2      = _Eff_Data->GetBinContent(bin_PT_2,bin_Mass_2);
      Eff_MC_2 =1.0;
      if(_Eff_MC->GetBinContent(bin_PT_2,bin_Mass_2) >0)
         Eff_MC_2        = _Eff_MC->GetBinContent(bin_PT_2,bin_Mass_2);
     }

   double Tot_Data = 1 - (1 - Eff_Data_1)*(1 - Eff_Data_2);
   double Tot_MC   = 1 - (1 - Eff_MC_1)*(1 - Eff_MC_2);
   double PT_Mass_SF = Tot_Data/Tot_MC;

   if(PT_Mass_SF > 0)
       weight = weight*PT_Mass_SF;

// ********************************************** PTJJ

   double FatJet1_px = FatJet1_pt*cos(FatJet1_phi), FatJet1_py = FatJet1_pt*sin(FatJet1_phi);
   double FatJet2_px = FatJet2_pt*cos(FatJet2_phi), FatJet2_py = FatJet2_pt*sin(FatJet2_phi);
   double PTjj = sqrt(pow((FatJet1_px +FatJet2_px),2) + pow((FatJet1_py +FatJet2_py),2));

// ********************************************************** Fill Histograms
  T_weight   = weight;

  T_fatJet1_pt                         = FatJet1_pt;
  T_fatJet1_eta                        = FatJet1_eta;
  T_fatJet1_phi                        = FatJet1_phi;
  T_fatJet1_msoftdrop                  = FatJet1_MassSD;
  T_fatJet1_particleNet_XbbVsQCD       = FatJet1PNetMD_Xbb;
  T_fatJet1_Tau3OverTau2               = FatJet1_Tau3OverTau2;
  T_fatJet1_particleNet_Xbb_Legacy     = FatJet1PNetMD_Xbb_Legacy;
  double FatJet1PNetMD_Xbb_Legacy_AN   = FatJet1PNetMD_Xbb_Legacy/(FatJet1PNetMD_Xbb_Legacy + FatJet1PNetMD_QCD_Legacy);
  T_fatJet1_particleNet_Xbb_Legacy_AN  = FatJet1PNetMD_Xbb_Legacy_AN;

  T_fatJet2_pt                         = FatJet2_pt;
  T_fatJet2_eta                        = FatJet2_eta;
  T_fatJet2_phi                        = FatJet2_phi;
  T_fatJet2_msoftdrop                  = FatJet2_MassSD;
  T_fatJet2_particleNet_XbbVsQCD       = FatJet2PNetMD_Xbb;
  T_fatJet2_Tau3OverTau2               = FatJet2_Tau3OverTau2;
  T_fatJet2_particleNet_Xbb_Legacy     = FatJet2PNetMD_Xbb_Legacy;
  double FatJet2PNetMD_Xbb_Legacy_AN   = FatJet2PNetMD_Xbb_Legacy/(FatJet2PNetMD_Xbb_Legacy + FatJet2PNetMD_QCD_Legacy);
  T_fatJet2_particleNet_Xbb_Legacy_AN  = FatJet2PNetMD_Xbb_Legacy_AN;

  T_PTjj = PTjj;

  outputTree->Fill();

 } // end event loop

f->Write();

}
