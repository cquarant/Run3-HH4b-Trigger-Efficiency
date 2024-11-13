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


// ******************************************
TString Run = "D";

void EGamma_Making_Histo_D()
{
 map<int,vector<pair<int,int>>> Good_Lumis
  {
   #include "/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/Data/GoodLumiList_Map_eraD.txt"
  };

// JEC
gSystem->Load("libFWCoreFWLite.so");
vector<JetCorrectorParameters> vPar_AK8;
vPar_AK8.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23BPixPrompt23_RunD_V1_DATA/Summer23BPixPrompt23_RunD_V1_DATA_L2Relative_AK8PFPuppi.txt"));
vPar_AK8.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23BPixPrompt23_RunD_V1_DATA/Summer23BPixPrompt23_RunD_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK8 = new FactorizedJetCorrector(vPar_AK8);

TFile *f = new TFile("EGamma_Histograms_"+Run+".root","RECREATE");

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

TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/Legacy/PostBPix/Run2023"+Run+"_EGamma.root");

TTree *InputTree = (TTree*)f1->Get("tree");
TH1F  *NEvents = (TH1F*)f1->Get("nPU_True");
double Tot_Events = NEvents->GetEntries();

UInt_t     run;
UInt_t     lumi;
Int_t      isVBFtag;

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

InputTree->SetBranchAddress("run",&run);
InputTree->SetBranchAddress("lumi",&lumi);
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

  // Events Loop
  for(int i=0;i<InputTree->GetEntries();i++)
   {
    InputTree->GetEntry(i);
// ********************************************************** HLT Selection
    if ( ! (HLT_Ele32_WPTight_Gsf && fabs(lep1_Id) ==11) ) continue;

// ********************************************************** JSON Certification
    bool Certified=false;
    for (pair<int, vector<pair<int,int>>> Run_Lumi : Good_Lumis)
      {
       if(Run_Lumi.first>run) break;
       if(Run_Lumi.first==run)
       for(auto LumiBlok : Run_Lumi.second)
         if ( inRange(LumiBlok.first,LumiBlok.second,lumi) ) { Certified=true; break; }
      }
    if(!Certified) continue;

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

// ********************************************************** Fill Histograms

  T_weight = 1.0;
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
