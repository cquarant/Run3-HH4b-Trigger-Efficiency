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
TString Run = "C";

void Making_Tree_C()
{
 map<int,vector<pair<int,int>>> Good_Lumis
  {
   #include "/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyHisto/NanoAOD_V12/Data/GoodLumiList_Map_eraC.txt"
  };

// JEC
gSystem->Load("libFWCoreFWLite.so");
vector<JetCorrectorParameters> vPar_AK8_Cv123;
vPar_AK8_Cv123.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2Relative_AK8PFPuppi.txt"));
vPar_AK8_Cv123.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK8_Cv123 = new FactorizedJetCorrector(vPar_AK8_Cv123);

vector<JetCorrectorParameters> vPar_AK8_Cv4;
vPar_AK8_Cv4.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2Relative_AK8PFPuppi.txt"));
vPar_AK8_Cv4.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2L3Residual_AK8PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK8_Cv4 = new FactorizedJetCorrector(vPar_AK8_Cv4);

gSystem->Load("libFWCoreFWLite.so");
vector<JetCorrectorParameters> vPar_AK4_Cv123;
vPar_AK4_Cv123.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2Relative_AK4PFPuppi.txt"));
vPar_AK4_Cv123.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2L3Residual_AK4PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK4_Cv123 = new FactorizedJetCorrector(vPar_AK4_Cv123);

vector<JetCorrectorParameters> vPar_AK4_Cv4;
vPar_AK4_Cv4.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2Relative_AK4PFPuppi.txt"));
vPar_AK4_Cv4.push_back(JetCorrectorParameters("/afs/cern.ch/user/t/tumasyan/public/2023/JECs/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2L3Residual_AK4PFPuppi.txt"));
FactorizedJetCorrector*  corrector_AK4_Cv4 = new FactorizedJetCorrector(vPar_AK4_Cv4);

// **************************************************** Output Tree Begin

TFile *f = new TFile("Tree_"+Run+".root","RECREATE");
TTree *outputTree = new TTree("tree", "");

Float_t T_weight;
outputTree->Branch("T_weight", &T_weight, "T_weight/F");

// FatJet1
Float_t T_fatJet1_pt, T_fatJet1_eta, T_fatJet1_phi, T_fatJet1_mass, T_fatJet1_msoftdrop, T_fatJet1_btagDDBvLV2, T_fatJet1_btagDeepB, T_fatJet1_particleNetWithMass_H4qvsQCD, T_fatJet1_particleNetWithMass_HbbvsQCD;
Float_t T_fatJet1_particleNetWithMass_QCD, T_fatJet1_particleNetWithMass_TvsQCD, T_fatJet1_particleNetWithMass_WvsQCD, T_fatJet1_particleNetWithMass_ZvsQCD, T_fatJet1_particleNet_QCD;
Float_t T_fatJet1_particleNet_QCD0HF, T_fatJet1_particleNet_QCD1HF, T_fatJet1_particleNet_QCD2HF, T_fatJet1_particleNet_XbbVsQCD, T_fatJet1_particleNet_XccVsQCD, T_fatJet1_particleNet_XggVsQCD;
Float_t T_fatJet1_particleNet_XqqVsQCD, T_fatJet1_particleNet_massCorr, T_fatJet1_n2b1, T_fatJet1_n3b1, T_fatJet1_Tau3OverTau2;
Bool_t  T_fatJet1_HasMuon, T_fatJet1_HasElectron, T_fatJet1_HasBJetMedium, T_fatJet1_OppositeHemisphereHasBJet;
Float_t T_fatJet1_subJetIdx1_pt, T_fatJet1_subJetIdx1_eta, T_fatJet1_subJetIdx1_phi, T_fatJet1_subJetIdx1_mass;
Float_t T_fatJet1_subJetIdx2_pt, T_fatJet1_subJetIdx2_eta, T_fatJet1_subJetIdx2_phi, T_fatJet1_subJetIdx2_mass;
Float_t T_fatJet1_subJets_alpha;
Float_t T_fatJet1_particleNet_Xbb_Legacy;
Float_t T_fatJet1_particleNet_Xcc_Legacy;
Float_t T_fatJet1_particleNet_Xqq_Legacy;
Float_t T_fatJet1_particleNet_QCD_Legacy;
Float_t T_fatJet1_particleNet_QCDb_Legacy;
Float_t T_fatJet1_particleNet_QCDbb_Legacy;
Float_t T_fatJet1_particleNet_QCDothers_Legacy;
Float_t T_fatJet1_particleNet_mass_Legacy;

outputTree->Branch("T_fatJet1_pt", &T_fatJet1_pt, "T_fatJet1_pt/F");
outputTree->Branch("T_fatJet1_eta", &T_fatJet1_eta, "T_fatJet1_eta/F");
outputTree->Branch("T_fatJet1_phi", &T_fatJet1_phi, "T_fatJet1_phi/F");
outputTree->Branch("T_fatJet1_mass", &T_fatJet1_mass, "T_fatJet1_mass/F");
outputTree->Branch("T_fatJet1_msoftdrop", &T_fatJet1_msoftdrop, "T_fatJet1_msoftdrop/F");
outputTree->Branch("T_fatJet1_btagDeepB", &T_fatJet1_btagDeepB, "T_fatJet1_btagDeepB/F");
outputTree->Branch("T_fatJet1_particleNetWithMass_H4qvsQCD", &T_fatJet1_particleNetWithMass_H4qvsQCD, "T_fatJet1_particleNetWithMass_H4qvsQCD/F");
outputTree->Branch("T_fatJet1_particleNetWithMass_HbbvsQCD", &T_fatJet1_particleNetWithMass_HbbvsQCD, "T_fatJet1_particleNetWithMass_HbbvsQCD/F");
outputTree->Branch("T_fatJet1_particleNetWithMass_QCD", &T_fatJet1_particleNetWithMass_QCD, "T_fatJet1_particleNetWithMass_QCD/F");
outputTree->Branch("T_fatJet1_particleNetWithMass_TvsQCD", &T_fatJet1_particleNetWithMass_TvsQCD, "T_fatJet1_particleNetWithMass_TvsQCD/F");
outputTree->Branch("T_fatJet1_particleNetWithMass_WvsQCD", &T_fatJet1_particleNetWithMass_WvsQCD, "T_fatJet1_particleNetWithMass_WvsQCD/F");
outputTree->Branch("T_fatJet1_particleNetWithMass_ZvsQCD", &T_fatJet1_particleNetWithMass_ZvsQCD, "T_fatJet1_particleNetWithMass_ZvsQCD/F");
outputTree->Branch("T_fatJet1_particleNet_QCD", &T_fatJet1_particleNet_QCD, "T_fatJet1_particleNet_QCD/F");
outputTree->Branch("T_fatJet1_particleNet_QCD0HF", &T_fatJet1_particleNet_QCD0HF, "T_fatJet1_particleNet_QCD0HF/F");
outputTree->Branch("T_fatJet1_particleNet_QCD1HF", &T_fatJet1_particleNet_QCD1HF, "T_fatJet1_particleNet_QCD1HF/F");
outputTree->Branch("T_fatJet1_particleNet_QCD2HF", &T_fatJet1_particleNet_QCD2HF, "T_fatJet1_particleNet_QCD2HF/F");
outputTree->Branch("T_fatJet1_particleNet_XbbVsQCD", &T_fatJet1_particleNet_XbbVsQCD, "T_fatJet1_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet1_particleNet_XccVsQCD", &T_fatJet1_particleNet_XccVsQCD, "T_fatJet1_particleNet_XccVsQCD/F");
outputTree->Branch("T_fatJet1_particleNet_XggVsQCD", &T_fatJet1_particleNet_XggVsQCD, "T_fatJet1_particleNet_XggVsQCD/F");
outputTree->Branch("T_fatJet1_particleNet_XqqVsQCD", &T_fatJet1_particleNet_XqqVsQCD, "T_fatJet1_particleNet_XqqVsQCD/F");
outputTree->Branch("T_fatJet1_particleNet_massCorr", &T_fatJet1_particleNet_massCorr, "T_fatJet1_particleNet_massCorr/F");
outputTree->Branch("T_fatJet1_n2b1", &T_fatJet1_n2b1, "T_fatJet1_n2b1/F");
outputTree->Branch("T_fatJet1_n3b1", &T_fatJet1_n3b1, "T_fatJet1_n3b1/F");
outputTree->Branch("T_fatJet1_Tau3OverTau2", &T_fatJet1_Tau3OverTau2, "T_fatJet1_Tau3OverTau2/F");
outputTree->Branch("T_fatJet1_HasMuon", &T_fatJet1_HasMuon, "T_fatJet1_HasMuon/O");
outputTree->Branch("T_fatJet1_HasElectron", &T_fatJet1_HasElectron, "T_fatJet1_HasElectron/O");
outputTree->Branch("T_fatJet1_HasBJetMedium", &T_fatJet1_HasBJetMedium, "T_fatJet1_HasBJetMedium/O");
outputTree->Branch("T_fatJet1_OppositeHemisphereHasBJet", &T_fatJet1_OppositeHemisphereHasBJet, "T_fatJet1_OppositeHemisphereHasBJet/O");
outputTree->Branch("T_fatJet1_subJetIdx1_pt", &T_fatJet1_subJetIdx1_pt, "T_fatJet1_subJetIdx1_pt/F");
outputTree->Branch("T_fatJet1_subJetIdx1_eta", &T_fatJet1_subJetIdx1_eta, "T_fatJet1_subJetIdx1_eta/F");
outputTree->Branch("T_fatJet1_subJetIdx1_phi", &T_fatJet1_subJetIdx1_phi, "T_fatJet1_subJetIdx1_phi/F");
outputTree->Branch("T_fatJet1_subJetIdx1_mass", &T_fatJet1_subJetIdx1_mass, "T_fatJet1_subJetIdx1_mass/F");
outputTree->Branch("T_fatJet1_subJetIdx2_pt", &T_fatJet1_subJetIdx2_pt, "T_fatJet1_subJetIdx2_pt/F");
outputTree->Branch("T_fatJet1_subJetIdx2_eta", &T_fatJet1_subJetIdx2_eta, "T_fatJet1_subJetIdx2_eta/F");
outputTree->Branch("T_fatJet1_subJetIdx2_phi", &T_fatJet1_subJetIdx2_phi, "T_fatJet1_subJetIdx2_phi/F");
outputTree->Branch("T_fatJet1_subJetIdx2_mass", &T_fatJet1_subJetIdx2_mass, "T_fatJet1_subJetIdx2_mass/F");
outputTree->Branch("T_fatJet1_subJets_alpha", &T_fatJet1_subJets_alpha, "T_fatJet1_subJets_alpha/F");
outputTree->Branch("T_fatJet1_particleNet_Xbb_Legacy", &T_fatJet1_particleNet_Xbb_Legacy, "T_fatJet1_particleNet_Xbb_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_Xcc_Legacy", &T_fatJet1_particleNet_Xcc_Legacy, "T_fatJet1_particleNet_Xcc_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_Xqq_Legacy", &T_fatJet1_particleNet_Xqq_Legacy, "T_fatJet1_particleNet_Xqq_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_QCD_Legacy", &T_fatJet1_particleNet_QCD_Legacy, "T_fatJet1_particleNet_QCD_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_QCDb_Legacy", &T_fatJet1_particleNet_QCDb_Legacy, "T_fatJet1_particleNet_QCDb_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_QCDbb_Legacy", &T_fatJet1_particleNet_QCDbb_Legacy, "T_fatJet1_particleNet_QCDbb_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_QCDothers_Legacy", &T_fatJet1_particleNet_QCDothers_Legacy, "T_fatJet1_particleNet_QCDothers_Legacy/F");
outputTree->Branch("T_fatJet1_particleNet_mass_Legacy", &T_fatJet1_particleNet_mass_Legacy, "T_fatJet1_particleNet_mass_Legacy/F");

// FatJet2
Float_t T_fatJet2_pt, T_fatJet2_eta, T_fatJet2_phi, T_fatJet2_mass, T_fatJet2_msoftdrop, T_fatJet2_btagDDBvLV2, T_fatJet2_btagDeepB, T_fatJet2_particleNetWithMass_H4qvsQCD, T_fatJet2_particleNetWithMass_HbbvsQCD;
Float_t T_fatJet2_particleNetWithMass_QCD, T_fatJet2_particleNetWithMass_TvsQCD, T_fatJet2_particleNetWithMass_WvsQCD, T_fatJet2_particleNetWithMass_ZvsQCD, T_fatJet2_particleNet_QCD;
Float_t T_fatJet2_particleNet_QCD0HF, T_fatJet2_particleNet_QCD1HF, T_fatJet2_particleNet_QCD2HF, T_fatJet2_particleNet_XbbVsQCD, T_fatJet2_particleNet_XccVsQCD, T_fatJet2_particleNet_XggVsQCD;
Float_t T_fatJet2_particleNet_XqqVsQCD, T_fatJet2_particleNet_massCorr, T_fatJet2_n2b1, T_fatJet2_n3b1, T_fatJet2_Tau3OverTau2;
Bool_t  T_fatJet2_HasMuon, T_fatJet2_HasElectron, T_fatJet2_HasBJetMedium, T_fatJet2_OppositeHemisphereHasBJet;
Float_t T_fatJet2_subJetIdx1_pt, T_fatJet2_subJetIdx1_eta, T_fatJet2_subJetIdx1_phi, T_fatJet2_subJetIdx1_mass;
Float_t T_fatJet2_subJetIdx2_pt, T_fatJet2_subJetIdx2_eta, T_fatJet2_subJetIdx2_phi, T_fatJet2_subJetIdx2_mass;
Float_t T_fatJet2_subJets_alpha;
Float_t T_fatJet2_particleNet_Xbb_Legacy;
Float_t T_fatJet2_particleNet_Xcc_Legacy;
Float_t T_fatJet2_particleNet_Xqq_Legacy;
Float_t T_fatJet2_particleNet_QCD_Legacy;
Float_t T_fatJet2_particleNet_QCDb_Legacy;
Float_t T_fatJet2_particleNet_QCDbb_Legacy;
Float_t T_fatJet2_particleNet_QCDothers_Legacy;
Float_t T_fatJet2_particleNet_mass_Legacy;

outputTree->Branch("T_fatJet2_pt", &T_fatJet2_pt, "T_fatJet2_pt/F");
outputTree->Branch("T_fatJet2_eta", &T_fatJet2_eta, "T_fatJet2_eta/F");
outputTree->Branch("T_fatJet2_phi", &T_fatJet2_phi, "T_fatJet2_phi/F");
outputTree->Branch("T_fatJet2_mass", &T_fatJet2_mass, "T_fatJet2_mass/F");
outputTree->Branch("T_fatJet2_msoftdrop", &T_fatJet2_msoftdrop, "T_fatJet2_msoftdrop/F");
outputTree->Branch("T_fatJet2_btagDeepB", &T_fatJet2_btagDeepB, "T_fatJet2_btagDeepB/F");
outputTree->Branch("T_fatJet2_particleNetWithMass_H4qvsQCD", &T_fatJet2_particleNetWithMass_H4qvsQCD, "T_fatJet2_particleNetWithMass_H4qvsQCD/F");
outputTree->Branch("T_fatJet2_particleNetWithMass_HbbvsQCD", &T_fatJet2_particleNetWithMass_HbbvsQCD, "T_fatJet2_particleNetWithMass_HbbvsQCD/F");
outputTree->Branch("T_fatJet2_particleNetWithMass_QCD", &T_fatJet2_particleNetWithMass_QCD, "T_fatJet2_particleNetWithMass_QCD/F");
outputTree->Branch("T_fatJet2_particleNetWithMass_TvsQCD", &T_fatJet2_particleNetWithMass_TvsQCD, "T_fatJet2_particleNetWithMass_TvsQCD/F");
outputTree->Branch("T_fatJet2_particleNetWithMass_WvsQCD", &T_fatJet2_particleNetWithMass_WvsQCD, "T_fatJet2_particleNetWithMass_WvsQCD/F");
outputTree->Branch("T_fatJet2_particleNetWithMass_ZvsQCD", &T_fatJet2_particleNetWithMass_ZvsQCD, "T_fatJet2_particleNetWithMass_ZvsQCD/F");
outputTree->Branch("T_fatJet2_particleNet_QCD", &T_fatJet2_particleNet_QCD, "T_fatJet2_particleNet_QCD/F");
outputTree->Branch("T_fatJet2_particleNet_QCD0HF", &T_fatJet2_particleNet_QCD0HF, "T_fatJet2_particleNet_QCD0HF/F");
outputTree->Branch("T_fatJet2_particleNet_QCD1HF", &T_fatJet2_particleNet_QCD1HF, "T_fatJet2_particleNet_QCD1HF/F");
outputTree->Branch("T_fatJet2_particleNet_QCD2HF", &T_fatJet2_particleNet_QCD2HF, "T_fatJet2_particleNet_QCD2HF/F");
outputTree->Branch("T_fatJet2_particleNet_XbbVsQCD", &T_fatJet2_particleNet_XbbVsQCD, "T_fatJet2_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet2_particleNet_XccVsQCD", &T_fatJet2_particleNet_XccVsQCD, "T_fatJet2_particleNet_XccVsQCD/F");
outputTree->Branch("T_fatJet2_particleNet_XggVsQCD", &T_fatJet2_particleNet_XggVsQCD, "T_fatJet2_particleNet_XggVsQCD/F");
outputTree->Branch("T_fatJet2_particleNet_XqqVsQCD", &T_fatJet2_particleNet_XqqVsQCD, "T_fatJet2_particleNet_XqqVsQCD/F");
outputTree->Branch("T_fatJet2_particleNet_massCorr", &T_fatJet2_particleNet_massCorr, "T_fatJet2_particleNet_massCorr/F");
outputTree->Branch("T_fatJet2_n2b1", &T_fatJet2_n2b1, "T_fatJet2_n2b1/F");
outputTree->Branch("T_fatJet2_n3b1", &T_fatJet2_n3b1, "T_fatJet2_n3b1/F");
outputTree->Branch("T_fatJet2_Tau3OverTau2", &T_fatJet2_Tau3OverTau2, "T_fatJet2_Tau3OverTau2/F");
outputTree->Branch("T_fatJet2_HasMuon", &T_fatJet2_HasMuon, "T_fatJet2_HasMuon/O");
outputTree->Branch("T_fatJet2_HasElectron", &T_fatJet2_HasElectron, "T_fatJet2_HasElectron/O");
outputTree->Branch("T_fatJet2_HasBJetMedium", &T_fatJet2_HasBJetMedium, "T_fatJet2_HasBJetMedium/O");
outputTree->Branch("T_fatJet2_OppositeHemisphereHasBJet", &T_fatJet2_OppositeHemisphereHasBJet, "T_fatJet2_OppositeHemisphereHasBJet/O");
outputTree->Branch("T_fatJet2_subJetIdx1_pt", &T_fatJet2_subJetIdx1_pt, "T_fatJet2_subJetIdx1_pt/F");
outputTree->Branch("T_fatJet2_subJetIdx1_eta", &T_fatJet2_subJetIdx1_eta, "T_fatJet2_subJetIdx1_eta/F");
outputTree->Branch("T_fatJet2_subJetIdx1_phi", &T_fatJet2_subJetIdx1_phi, "T_fatJet2_subJetIdx1_phi/F");
outputTree->Branch("T_fatJet2_subJetIdx1_mass", &T_fatJet2_subJetIdx1_mass, "T_fatJet2_subJetIdx1_mass/F");
outputTree->Branch("T_fatJet2_subJetIdx2_pt", &T_fatJet2_subJetIdx2_pt, "T_fatJet2_subJetIdx2_pt/F");
outputTree->Branch("T_fatJet2_subJetIdx2_eta", &T_fatJet2_subJetIdx2_eta, "T_fatJet2_subJetIdx2_eta/F");
outputTree->Branch("T_fatJet2_subJetIdx2_phi", &T_fatJet2_subJetIdx2_phi, "T_fatJet2_subJetIdx2_phi/F");
outputTree->Branch("T_fatJet2_subJetIdx2_mass", &T_fatJet2_subJetIdx2_mass, "T_fatJet2_subJetIdx2_mass/F");
outputTree->Branch("T_fatJet2_subJets_alpha", &T_fatJet2_subJets_alpha, "T_fatJet2_subJets_alpha/F");
outputTree->Branch("T_fatJet2_particleNet_Xbb_Legacy", &T_fatJet2_particleNet_Xbb_Legacy, "T_fatJet2_particleNet_Xbb_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_Xcc_Legacy", &T_fatJet2_particleNet_Xcc_Legacy, "T_fatJet2_particleNet_Xcc_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_Xqq_Legacy", &T_fatJet2_particleNet_Xqq_Legacy, "T_fatJet2_particleNet_Xqq_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_QCD_Legacy", &T_fatJet2_particleNet_QCD_Legacy, "T_fatJet2_particleNet_QCD_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_QCDb_Legacy", &T_fatJet2_particleNet_QCDb_Legacy, "T_fatJet2_particleNet_QCDb_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_QCDbb_Legacy", &T_fatJet2_particleNet_QCDbb_Legacy, "T_fatJet2_particleNet_QCDbb_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_QCDothers_Legacy", &T_fatJet2_particleNet_QCDothers_Legacy, "T_fatJet2_particleNet_QCDothers_Legacy/F");
outputTree->Branch("T_fatJet2_particleNet_mass_Legacy", &T_fatJet2_particleNet_mass_Legacy, "T_fatJet2_particleNet_mass_Legacy/F");

// FatJet3
Float_t T_fatJet3_pt, T_fatJet3_eta, T_fatJet3_phi, T_fatJet3_mass, T_fatJet3_msoftdrop, T_fatJet3_particleNet_XbbVsQCD, T_fatJet3_particleNet_Xbb_Legacy;

outputTree->Branch("T_fatJet3_pt", &T_fatJet3_pt, "T_fatJet3_pt/F");
outputTree->Branch("T_fatJet3_eta", &T_fatJet3_eta, "T_fatJet3_eta/F");
outputTree->Branch("T_fatJet3_phi", &T_fatJet3_phi, "T_fatJet3_phi/F");
outputTree->Branch("T_fatJet3_mass", &T_fatJet3_mass, "T_fatJet3_mass/F");
outputTree->Branch("T_fatJet3_msoftdrop", &T_fatJet3_msoftdrop, "T_fatJet3_msoftdrop/F");
outputTree->Branch("T_fatJet3_particleNet_XbbVsQCD", &T_fatJet3_particleNet_XbbVsQCD, "T_fatJet3_particleNet_XbbVsQCD/F");
outputTree->Branch("T_fatJet3_particleNet_Xbb_Legacy", &T_fatJet3_particleNet_Xbb_Legacy, "T_fatJet3_particleNet_Xbb_Legacy/F");

// Other variables
Int_t T_NJets,T_nBTaggedJets;
Float_t T_MET, T_alpha_FJ12, T_ptFJ12, T_etaFJ12, T_phiFJ12, T_massFJ12, T_detaFJ12, T_dphiFJ12, T_dRFJ12, T_FJ1_PtOverMFJ12, T_FJ2_PtOverMFJ12, T_FJ1_PtOverMSDFJ1, T_FJ2_PtOverMSDFJ2;
Float_t T_ptFJ1OverptFJ2, T_massFJ1OvermassFJ2, T_ptak4J1, T_ptak4J2, T_dRak4J1FJ1, T_dRak4J1FJ2, T_dRak4J2FJ1, T_dRak4J2FJ2;

outputTree->Branch("T_NJets", &T_NJets, "T_NJets/I");
outputTree->Branch("T_nBTaggedJets", &T_nBTaggedJets, "T_nBTaggedJets/I");
outputTree->Branch("T_MET", &T_MET, "T_MET/F");
outputTree->Branch("T_alpha_FJ12", &T_alpha_FJ12, "T_alpha_FJ12/F");
outputTree->Branch("T_ptFJ12", &T_ptFJ12, "T_ptFJ12/F");
outputTree->Branch("T_etaFJ12", &T_etaFJ12, "T_etaFJ12/F");
outputTree->Branch("T_phiFJ12", &T_phiFJ12, "T_phiFJ12/F");
outputTree->Branch("T_massFJ12", &T_massFJ12, "T_massFJ12/F");
outputTree->Branch("T_detaFJ12", &T_detaFJ12, "T_detaFJ12/F");
outputTree->Branch("T_dphiFJ12", &T_dphiFJ12, "T_dphiFJ12/F");
outputTree->Branch("T_dRFJ12", &T_dRFJ12, "T_dRFJ12/F");
outputTree->Branch("T_FJ1_PtOverMFJ12", &T_FJ1_PtOverMFJ12, "T_FJ1_PtOverMFJ12/F");
outputTree->Branch("T_FJ2_PtOverMFJ12", &T_FJ2_PtOverMFJ12, "T_FJ2_PtOverMFJ12/F");
outputTree->Branch("T_FJ1_PtOverMSDFJ1", &T_FJ1_PtOverMSDFJ1, "T_FJ1_PtOverMSDFJ1/F");
outputTree->Branch("T_FJ2_PtOverMSDFJ2", &T_FJ2_PtOverMSDFJ2, "T_FJ2_PtOverMSDFJ2/F");
outputTree->Branch("T_ptFJ1OverptFJ2", &T_ptFJ1OverptFJ2, "T_ptFJ1OverptFJ2/F");
outputTree->Branch("T_massFJ1OvermassFJ2", &T_massFJ1OvermassFJ2, "T_massFJ1OvermassFJ2/F");
outputTree->Branch("T_ptak4J1", &T_ptak4J1, "T_ptak4J1/F");
outputTree->Branch("T_ptak4J2", &T_ptak4J2, "T_ptak4J2/F");
outputTree->Branch("T_dRak4J1FJ1", &T_dRak4J1FJ1, "T_dRak4J1FJ1/F");
outputTree->Branch("T_dRak4J1FJ2", &T_dRak4J1FJ2, "T_dRak4J1FJ2/F");
outputTree->Branch("T_dRak4J2FJ1", &T_dRak4J2FJ1, "T_dRak4J2FJ1/F");
outputTree->Branch("T_dRak4J2FJ2", &T_dRak4J2FJ2, "T_dRak4J2FJ2/F");

// **************************************************** Output Tree END

TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/Legacy/PreBPix/Run2023"+Run+"_JetMET.root");

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
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06;
Bool_t     HLT_AK8PFJet230_SoftDropMass40;
Bool_t     HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;

Int_t      NJets;
Int_t      nBTaggedJets;
Float_t    ak4J_pt[2];
Float_t    ak4J_eta[2];
Float_t    ak4J_phi[2];
Float_t    ak4J_rawFactor[2];
Float_t    MET;
Float_t    lep1_Pt;
Float_t    lep1_Eta;
Float_t    lep1_Phi;
Int_t      lep1_Id;

Float_t fatJet_pt[3];
Float_t fatJet_eta[3];
Float_t fatJet_phi[3];
Float_t fatJet_mass[3];
Float_t fatJet_msoftdrop[3];
Float_t fatJet_btagDDBvLV2[3];
Float_t fatJet_btagDDCvBV2[3];
Float_t fatJet_btagDDCvLV2[3];
Float_t fatJet_btagDeepB[3];
Float_t fatJet_particleNetWithMass_H4qvsQCD[3];
Float_t fatJet_particleNetWithMass_HbbvsQCD[3];
Float_t fatJet_particleNetWithMass_HccvsQCD[3];
Float_t fatJet_particleNetWithMass_QCD[3];
Float_t fatJet_particleNetWithMass_TvsQCD[3];
Float_t fatJet_particleNetWithMass_WvsQCD[3];
Float_t fatJet_particleNetWithMass_ZvsQCD[3];
Float_t fatJet_particleNet_QCD[3];
Float_t fatJet_particleNet_QCD0HF[3];
Float_t fatJet_particleNet_QCD1HF[3];
Float_t fatJet_particleNet_QCD2HF[3];
Float_t fatJet_particleNet_XbbVsQCD[3];
Float_t fatJet_particleNet_XccVsQCD[3];
Float_t fatJet_particleNet_XggVsQCD[3];
Float_t fatJet_particleNet_XqqVsQCD[3];
Float_t fatJet_particleNet_massCorr[3];
Float_t fatJet_n2b1[3];
Float_t fatJet_n3b1[3];
Float_t fatJet_Tau3OverTau2[3];
Float_t fatJet_rawFactor[3];
Bool_t  fatJet_HasMuon[3];
Bool_t  fatJet_HasElectron[3];
Bool_t  fatJet_HasBJetMedium[3];
Bool_t  fatJet_OppositeHemisphereHasBJet[3];
Float_t fatJet_subJetIdx1_pt[3];
Float_t fatJet_subJetIdx1_eta[3];
Float_t fatJet_subJetIdx1_phi[3];
Float_t fatJet_subJetIdx1_mass[3];
Float_t fatJet_subJetIdx2_pt[3];
Float_t fatJet_subJetIdx2_eta[3];
Float_t fatJet_subJetIdx2_phi[3];
Float_t fatJet_subJetIdx2_mass[3];
Float_t fatJet_particleNet_Xbb_Legacy[3];
Float_t fatJet_particleNet_Xcc_Legacy[3];
Float_t fatJet_particleNet_Xqq_Legacy[3];
Float_t fatJet_particleNet_QCD_Legacy[3];
Float_t fatJet_particleNet_QCDb_Legacy[3];
Float_t fatJet_particleNet_QCDbb_Legacy[3];
Float_t fatJet_particleNet_QCDothers_Legacy[3];
Float_t fatJet_particleNet_mass_Legacy[3];

InputTree->SetBranchAddress("run",&run);
InputTree->SetBranchAddress("lumi",&lumi);
InputTree->SetBranchAddress("isVBFtag",&isVBFtag);

InputTree->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("HLT_Mu50", &HLT_Mu50);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40);
InputTree->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06);

InputTree->SetBranchAddress("lep1Pt",&lep1_Pt);
InputTree->SetBranchAddress("lep1Eta",&lep1_Eta);
InputTree->SetBranchAddress("lep1Phi",&lep1_Phi);
InputTree->SetBranchAddress("lep1Id",&lep1_Id);
InputTree->SetBranchAddress("MET",&MET);
InputTree->SetBranchAddress("NJets",&NJets);
InputTree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
InputTree->SetBranchAddress("jet1Pt",&ak4J_pt[0]);
InputTree->SetBranchAddress("jet1Eta",&ak4J_eta[0]);
InputTree->SetBranchAddress("jet1Phi",&ak4J_phi[0]);
InputTree->SetBranchAddress("jet1_rawFactor",&ak4J_rawFactor[0]);
InputTree->SetBranchAddress("jet2Pt",&ak4J_pt[1]);
InputTree->SetBranchAddress("jet2Eta",&ak4J_eta[1]);
InputTree->SetBranchAddress("jet2Phi",&ak4J_phi[1]);
InputTree->SetBranchAddress("jet2_rawFactor",&ak4J_rawFactor[1]);

InputTree->SetBranchAddress("fatJet1_pt",&fatJet_pt[0]);
InputTree->SetBranchAddress("fatJet1_eta",&fatJet_eta[0]);
InputTree->SetBranchAddress("fatJet1_phi",&fatJet_phi[0]);
InputTree->SetBranchAddress("fatJet1_mass",&fatJet_mass[0]);
InputTree->SetBranchAddress("fatJet1_msoftdrop",&fatJet_msoftdrop[0]);
InputTree->SetBranchAddress("fatJet1_btagDDBvLV2",&fatJet_btagDDBvLV2[0]);
InputTree->SetBranchAddress("fatJet1_btagDDCvBV2",&fatJet_btagDDCvBV2[0]);
InputTree->SetBranchAddress("fatJet1_btagDDCvLV2",&fatJet_btagDDCvLV2[0]);
InputTree->SetBranchAddress("fatJet1_btagDeepB",&fatJet_btagDeepB[0]);
InputTree->SetBranchAddress("fatJet1_particleNetWithMass_H4qvsQCD",&fatJet_particleNetWithMass_H4qvsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNetWithMass_HbbvsQCD",&fatJet_particleNetWithMass_HbbvsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNetWithMass_HccvsQCD",&fatJet_particleNetWithMass_HccvsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNetWithMass_QCD",&fatJet_particleNetWithMass_QCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNetWithMass_TvsQCD",&fatJet_particleNetWithMass_TvsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNetWithMass_WvsQCD",&fatJet_particleNetWithMass_WvsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNetWithMass_ZvsQCD",&fatJet_particleNetWithMass_ZvsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_QCD",&fatJet_particleNet_QCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_QCD0HF",&fatJet_particleNet_QCD0HF[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_QCD1HF",&fatJet_particleNet_QCD1HF[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_QCD2HF",&fatJet_particleNet_QCD2HF[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_XbbVsQCD",&fatJet_particleNet_XbbVsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_XccVsQCD",&fatJet_particleNet_XccVsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_XggVsQCD",&fatJet_particleNet_XggVsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_XqqVsQCD",&fatJet_particleNet_XqqVsQCD[0]);
InputTree->SetBranchAddress("fatJet1_particleNet_massCorr",&fatJet_particleNet_massCorr[0]);
InputTree->SetBranchAddress("fatJet1_n2b1",&fatJet_n2b1[0]);
InputTree->SetBranchAddress("fatJet1_n3b1",&fatJet_n3b1[0]);
InputTree->SetBranchAddress("fatJet1_Tau3OverTau2",&fatJet_Tau3OverTau2[0]);
InputTree->SetBranchAddress("fatJet1_rawFactor",&fatJet_rawFactor[0]);
InputTree->SetBranchAddress("fatJet1_HasMuon",&fatJet_HasMuon[0]);
InputTree->SetBranchAddress("fatJet1_HasElectron",&fatJet_HasElectron[0]);
InputTree->SetBranchAddress("fatJet1_HasBJetMedium",&fatJet_HasBJetMedium[0]);
InputTree->SetBranchAddress("fatJet1_OppositeHemisphereHasBJet",&fatJet_OppositeHemisphereHasBJet[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx1_pt",&fatJet_subJetIdx1_pt[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx1_eta",&fatJet_subJetIdx1_eta[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx1_phi",&fatJet_subJetIdx1_phi[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx1_mass",&fatJet_subJetIdx1_mass[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx2_pt",&fatJet_subJetIdx2_pt[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx2_eta",&fatJet_subJetIdx2_eta[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx2_phi",&fatJet_subJetIdx2_phi[0]);
InputTree->SetBranchAddress("fatJet1_subJetIdx2_mass",&fatJet_subJetIdx2_mass[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xbb",&fatJet_particleNet_Xbb_Legacy[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xcc",&fatJet_particleNet_Xcc_Legacy[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_Xqq",&fatJet_particleNet_Xqq_Legacy[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCD",&fatJet_particleNet_QCD_Legacy[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDb",&fatJet_particleNet_QCDb_Legacy[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDbb",&fatJet_particleNet_QCDbb_Legacy[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_QCDothers",&fatJet_particleNet_QCDothers_Legacy[0]);
InputTree->SetBranchAddress("fatJet1_particleNetLegacy_mass",&fatJet_particleNet_mass_Legacy[0]);

InputTree->SetBranchAddress("fatJet2_pt",&fatJet_pt[1]);
InputTree->SetBranchAddress("fatJet2_eta",&fatJet_eta[1]);
InputTree->SetBranchAddress("fatJet2_phi",&fatJet_phi[1]);
InputTree->SetBranchAddress("fatJet2_mass",&fatJet_mass[1]);
InputTree->SetBranchAddress("fatJet2_msoftdrop",&fatJet_msoftdrop[1]);
InputTree->SetBranchAddress("fatJet2_btagDDBvLV2",&fatJet_btagDDBvLV2[1]);
InputTree->SetBranchAddress("fatJet2_btagDDCvBV2",&fatJet_btagDDCvBV2[1]);
InputTree->SetBranchAddress("fatJet2_btagDDCvLV2",&fatJet_btagDDCvLV2[1]);
InputTree->SetBranchAddress("fatJet2_btagDeepB",&fatJet_btagDeepB[1]);
InputTree->SetBranchAddress("fatJet2_particleNetWithMass_H4qvsQCD",&fatJet_particleNetWithMass_H4qvsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNetWithMass_HbbvsQCD",&fatJet_particleNetWithMass_HbbvsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNetWithMass_HccvsQCD",&fatJet_particleNetWithMass_HccvsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNetWithMass_QCD",&fatJet_particleNetWithMass_QCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNetWithMass_TvsQCD",&fatJet_particleNetWithMass_TvsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNetWithMass_WvsQCD",&fatJet_particleNetWithMass_WvsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNetWithMass_ZvsQCD",&fatJet_particleNetWithMass_ZvsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_QCD",&fatJet_particleNet_QCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_QCD0HF",&fatJet_particleNet_QCD0HF[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_QCD1HF",&fatJet_particleNet_QCD1HF[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_QCD2HF",&fatJet_particleNet_QCD2HF[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_XbbVsQCD",&fatJet_particleNet_XbbVsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_XccVsQCD",&fatJet_particleNet_XccVsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_XggVsQCD",&fatJet_particleNet_XggVsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_XqqVsQCD",&fatJet_particleNet_XqqVsQCD[1]);
InputTree->SetBranchAddress("fatJet2_particleNet_massCorr",&fatJet_particleNet_massCorr[1]);
InputTree->SetBranchAddress("fatJet2_n2b1",&fatJet_n2b1[1]);
InputTree->SetBranchAddress("fatJet2_n3b1",&fatJet_n3b1[1]);
InputTree->SetBranchAddress("fatJet2_Tau3OverTau2",&fatJet_Tau3OverTau2[1]);
InputTree->SetBranchAddress("fatJet2_rawFactor",&fatJet_rawFactor[1]);
InputTree->SetBranchAddress("fatJet2_HasMuon",&fatJet_HasMuon[1]);
InputTree->SetBranchAddress("fatJet2_HasElectron",&fatJet_HasElectron[1]);
InputTree->SetBranchAddress("fatJet2_HasBJetMedium",&fatJet_HasBJetMedium[1]);
InputTree->SetBranchAddress("fatJet2_OppositeHemisphereHasBJet",&fatJet_OppositeHemisphereHasBJet[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx1_pt",&fatJet_subJetIdx1_pt[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx1_eta",&fatJet_subJetIdx1_eta[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx1_phi",&fatJet_subJetIdx1_phi[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx1_mass",&fatJet_subJetIdx1_mass[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx2_pt",&fatJet_subJetIdx2_pt[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx2_eta",&fatJet_subJetIdx2_eta[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx2_phi",&fatJet_subJetIdx2_phi[1]);
InputTree->SetBranchAddress("fatJet2_subJetIdx2_mass",&fatJet_subJetIdx2_mass[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xbb",&fatJet_particleNet_Xbb_Legacy[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xcc",&fatJet_particleNet_Xcc_Legacy[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_Xqq",&fatJet_particleNet_Xqq_Legacy[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCD",&fatJet_particleNet_QCD_Legacy[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCDb",&fatJet_particleNet_QCDb_Legacy[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCDbb",&fatJet_particleNet_QCDbb_Legacy[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_QCDothers",&fatJet_particleNet_QCDothers_Legacy[1]);
InputTree->SetBranchAddress("fatJet2_particleNetLegacy_mass",&fatJet_particleNet_mass_Legacy[1]);

InputTree->SetBranchAddress("fatJet3_pt",&fatJet_pt[2]);
InputTree->SetBranchAddress("fatJet3_eta",&fatJet_eta[2]);
InputTree->SetBranchAddress("fatJet3_phi",&fatJet_phi[2]);
InputTree->SetBranchAddress("fatJet3_mass",&fatJet_mass[2]);
InputTree->SetBranchAddress("fatJet3_msoftdrop",&fatJet_msoftdrop[2]);
InputTree->SetBranchAddress("fatJet3_btagDDBvLV2",&fatJet_btagDDBvLV2[2]);
InputTree->SetBranchAddress("fatJet3_btagDDCvBV2",&fatJet_btagDDCvBV2[2]);
InputTree->SetBranchAddress("fatJet3_btagDDCvLV2",&fatJet_btagDDCvLV2[2]);
InputTree->SetBranchAddress("fatJet3_btagDeepB",&fatJet_btagDeepB[2]);
InputTree->SetBranchAddress("fatJet3_particleNetWithMass_H4qvsQCD",&fatJet_particleNetWithMass_H4qvsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNetWithMass_HbbvsQCD",&fatJet_particleNetWithMass_HbbvsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNetWithMass_HccvsQCD",&fatJet_particleNetWithMass_HccvsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNetWithMass_QCD",&fatJet_particleNetWithMass_QCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNetWithMass_TvsQCD",&fatJet_particleNetWithMass_TvsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNetWithMass_WvsQCD",&fatJet_particleNetWithMass_WvsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNetWithMass_ZvsQCD",&fatJet_particleNetWithMass_ZvsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_QCD",&fatJet_particleNet_QCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_QCD0HF",&fatJet_particleNet_QCD0HF[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_QCD1HF",&fatJet_particleNet_QCD1HF[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_QCD2HF",&fatJet_particleNet_QCD2HF[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_XbbVsQCD",&fatJet_particleNet_XbbVsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_XccVsQCD",&fatJet_particleNet_XccVsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_XggVsQCD",&fatJet_particleNet_XggVsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_XqqVsQCD",&fatJet_particleNet_XqqVsQCD[2]);
InputTree->SetBranchAddress("fatJet3_particleNet_massCorr",&fatJet_particleNet_massCorr[2]);
InputTree->SetBranchAddress("fatJet3_n2b1",&fatJet_n2b1[2]);
InputTree->SetBranchAddress("fatJet3_n3b1",&fatJet_n3b1[2]);
InputTree->SetBranchAddress("fatJet3_Tau3OverTau2",&fatJet_Tau3OverTau2[2]);
InputTree->SetBranchAddress("fatJet3_rawFactor",&fatJet_rawFactor[2]);
InputTree->SetBranchAddress("fatJet3_HasMuon",&fatJet_HasMuon[2]);
InputTree->SetBranchAddress("fatJet3_HasElectron",&fatJet_HasElectron[2]);
InputTree->SetBranchAddress("fatJet3_HasBJetMedium",&fatJet_HasBJetMedium[2]);
InputTree->SetBranchAddress("fatJet3_OppositeHemisphereHasBJet",&fatJet_OppositeHemisphereHasBJet[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx1_pt",&fatJet_subJetIdx1_pt[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx1_eta",&fatJet_subJetIdx1_eta[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx1_phi",&fatJet_subJetIdx1_phi[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx1_mass",&fatJet_subJetIdx1_mass[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx2_pt",&fatJet_subJetIdx2_pt[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx2_eta",&fatJet_subJetIdx2_eta[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx2_phi",&fatJet_subJetIdx2_phi[2]);
InputTree->SetBranchAddress("fatJet3_subJetIdx2_mass",&fatJet_subJetIdx2_mass[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_Xbb",&fatJet_particleNet_Xbb_Legacy[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_Xcc",&fatJet_particleNet_Xcc_Legacy[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_Xqq",&fatJet_particleNet_Xqq_Legacy[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_QCD",&fatJet_particleNet_QCD_Legacy[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_QCDb",&fatJet_particleNet_QCDb_Legacy[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_QCDbb",&fatJet_particleNet_QCDbb_Legacy[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_QCDothers",&fatJet_particleNet_QCDothers_Legacy[2]);
InputTree->SetBranchAddress("fatJet3_particleNetLegacy_mass",&fatJet_particleNet_mass_Legacy[2]);

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
   if(HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06==0) continue;

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

    if (run<367661) continue;

// ********************************************************** AK4 b-Jets correction
   if(run <367765)
    {
     for(int i=0; i<2; i++)
      if(ak4J_pt[i] > 40)
       {
        double Raw_Jet_pt = ak4J_pt[i]*(1.0 - ak4J_rawFactor[i]);
        corrector_AK4_Cv123->setJetPt(Raw_Jet_pt);
        corrector_AK4_Cv123->setJetEta(ak4J_eta[i]);
        corrector_AK4_Cv123->setJetPhi(ak4J_phi[i]);
        double This_correction = corrector_AK4_Cv123->getCorrection();
        ak4J_pt[i]        = Raw_Jet_pt * This_correction;
       }
     }

   if(run >=367765)
    {
     for(int i=0; i<2; i++)
      if(ak4J_pt[i] > 40)
       {
        double Raw_Jet_pt = ak4J_pt[i]*(1.0 - ak4J_rawFactor[i]);
        corrector_AK4_Cv4->setJetPt(Raw_Jet_pt);
        corrector_AK4_Cv4->setJetEta(ak4J_eta[i]);
        corrector_AK4_Cv4->setJetPhi(ak4J_phi[i]);
        double This_correction = corrector_AK4_Cv4->getCorrection();
        ak4J_pt[i]        = Raw_Jet_pt * This_correction;
       }
     }
// ********************************************************** FatJets correction
   if(run <367765)
    {
     for(int i=0; i<3; i++)
      if(fatJet_pt[i] > 50)
       {
        double Raw_fatJet_pt = fatJet_pt[i]*(1.0 - fatJet_rawFactor[i]);
        corrector_AK8_Cv123->setJetPt(Raw_fatJet_pt);
        corrector_AK8_Cv123->setJetEta(fatJet_eta[i]);
        corrector_AK8_Cv123->setJetPhi(fatJet_phi[i]);
        double This_correction = corrector_AK8_Cv123->getCorrection();
        fatJet_pt[i]        = Raw_fatJet_pt * This_correction;
        fatJet_mass[i]      = fatJet_mass[i] * (1.0 - fatJet_rawFactor[i]) * This_correction;
        fatJet_msoftdrop[i] = fatJet_msoftdrop[i] * (1.0 - fatJet_rawFactor[i]) * This_correction;
       }
     }

   if(run >=367765)
    {
     for(int i=0; i<3; i++)
      if(fatJet_pt[i] > 50)
       {
        double Raw_fatJet_pt = fatJet_pt[i]*(1.0 - fatJet_rawFactor[i]);
        corrector_AK8_Cv4->setJetPt(Raw_fatJet_pt);
        corrector_AK8_Cv4->setJetEta(fatJet_eta[i]);
        corrector_AK8_Cv4->setJetPhi(fatJet_phi[i]);
        double This_correction = corrector_AK8_Cv4->getCorrection();
        fatJet_pt[i]        = Raw_fatJet_pt * This_correction;
        fatJet_mass[i]      = fatJet_mass[i] * (1.0 - fatJet_rawFactor[i]) * This_correction;
        fatJet_msoftdrop[i] = fatJet_msoftdrop[i] * (1.0 - fatJet_rawFactor[i]) * This_correction;
       }
     }
// ********************************************************** Ordering FatJets by Xbb-score
    for(int ii=0;ii<3; ii++)
      if(fatJet_pt[ii] > 230)
         fatJet_particleNet_Xbb_Legacy[ii] = fatJet_particleNet_Xbb_Legacy[ii]/(fatJet_particleNet_Xbb_Legacy[ii] + fatJet_particleNet_QCD_Legacy[ii]);

    int Idx_1=-10;
    int Idx_2=-10;
    int Idx_3=-10;
    double maxXbb_1=-1000;
    double maxXbb_2=-1000;

    for(int ii=0;ii<3; ii++)
      if(fatJet_pt[ii] > 230 && fatJet_particleNet_Xbb_Legacy[ii] > maxXbb_1)
       {
         maxXbb_1 = fatJet_particleNet_Xbb_Legacy[ii];
         Idx_1 =ii; // 1st Xbb tagged
       }

    for(int ii=0;ii<3; ii++)
      if(fatJet_pt[ii] > 230 && ii!=Idx_1 && fatJet_particleNet_Xbb_Legacy[ii] > maxXbb_2)
       {
         maxXbb_2 = fatJet_particleNet_Xbb_Legacy[ii];
         Idx_2 =ii;  // 2nd Xbb tagged
       }

    for(int ii=0;ii<3; ii++)
      if(fatJet_pt[ii] > 230 && ii!=Idx_1 && ii!=Idx_2)
         Idx_3 =ii; // 3nrd Xbb tagged

    if (Idx_1 < 0 || Idx_1 >2 || Idx_2 < 0 || Idx_2 > 2) continue;
    if (Idx_3 < 0) Idx_3=2;

    if(fatJet_pt[Idx_1] < 250 || fabs(fatJet_eta[Idx_1]) > 2.4 || fatJet_msoftdrop[Idx_1] < 50) continue;
    if(fatJet_pt[Idx_2] < 250 || fabs(fatJet_eta[Idx_2]) > 2.4 || fatJet_msoftdrop[Idx_2] < 50) continue;

    if(fatJet_particleNet_Xbb_Legacy[Idx_1] < 0.8) continue;

// ********************************************************** VBFTag veto
//   if(isVBFtag) continue;

// ********************************************************** Lepton Selection or Veto
//   if (fabs(lep1_Id) !=11 ) continue;

// ********************************************************** Trigger Objects and Matchings
  // Matching 1st
  bool matched_TRG_1=false;
  bool matched_to_AK8PFJet230_SoftDropMass40_1 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if((Trigger_Object_bit[itrg] & 4) == 4)
       if(sqrt(pow((fatJet_eta[Idx_1] - Trigger_Object_eta[itrg]),2) + pow(phi_dist(fatJet_phi[Idx_1],Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 100)
         {matched_to_AK8PFJet230_SoftDropMass40_1 = true; break;}
  if (matched_to_AK8PFJet230_SoftDropMass40_1) matched_TRG_1=true;

  // Matching 2nd
  bool matched_TRG_2=false;
  bool matched_to_AK8PFJet230_SoftDropMass40_2 = false;
  for(int itrg=0;itrg<NTrigger_Objects;itrg++)
    if((Trigger_Object_bit[itrg] & 4) == 4)
       if(sqrt(pow((fatJet_eta[Idx_2] - Trigger_Object_eta[itrg]),2) + pow(phi_dist(fatJet_phi[Idx_2],Trigger_Object_phi[itrg]),2)) < 0.4 && Trigger_Object_pt[itrg] > 100)
         {matched_to_AK8PFJet230_SoftDropMass40_2 = true; break;}
  if (matched_to_AK8PFJet230_SoftDropMass40_2) matched_TRG_2=true;

  if (!(matched_TRG_1 || matched_TRG_2)) continue;

// ********************************************** Final Tree Variables
   T_weight = 1.0;

   T_fatJet1_pt = fatJet_pt[Idx_1];
   T_fatJet1_eta = fatJet_eta[Idx_1];
   T_fatJet1_phi = fatJet_phi[Idx_1];
   T_fatJet1_mass = fatJet_mass[Idx_1];
   T_fatJet1_msoftdrop = fatJet_msoftdrop[Idx_1];
   T_fatJet1_btagDDBvLV2 = fatJet_btagDDBvLV2[Idx_1];
   T_fatJet1_btagDeepB = fatJet_btagDeepB[Idx_1];
   T_fatJet1_particleNetWithMass_H4qvsQCD = fatJet_particleNetWithMass_H4qvsQCD[Idx_1];
   T_fatJet1_particleNetWithMass_HbbvsQCD = fatJet_particleNetWithMass_HbbvsQCD[Idx_1];
   T_fatJet1_particleNetWithMass_QCD = fatJet_particleNetWithMass_QCD[Idx_1];
   T_fatJet1_particleNetWithMass_TvsQCD = fatJet_particleNetWithMass_TvsQCD[Idx_1];
   T_fatJet1_particleNetWithMass_WvsQCD = fatJet_particleNetWithMass_WvsQCD[Idx_1];
   T_fatJet1_particleNetWithMass_ZvsQCD = fatJet_particleNetWithMass_ZvsQCD[Idx_1];
   T_fatJet1_particleNet_QCD = fatJet_particleNet_QCD[Idx_1];
   T_fatJet1_particleNet_QCD0HF = fatJet_particleNet_QCD0HF[Idx_1];
   T_fatJet1_particleNet_QCD1HF = fatJet_particleNet_QCD1HF[Idx_1];
   T_fatJet1_particleNet_QCD2HF = fatJet_particleNet_QCD2HF[Idx_1];
   T_fatJet1_particleNet_XbbVsQCD = fatJet_particleNet_XbbVsQCD[Idx_1];
   T_fatJet1_particleNet_XccVsQCD = fatJet_particleNet_XccVsQCD[Idx_1];
   T_fatJet1_particleNet_XggVsQCD = fatJet_particleNet_XggVsQCD[Idx_1];
   T_fatJet1_particleNet_XqqVsQCD = fatJet_particleNet_XqqVsQCD[Idx_1];
   T_fatJet1_particleNet_massCorr = fatJet_mass[Idx_1]*fatJet_particleNet_massCorr[Idx_1];
   T_fatJet1_n2b1 = fatJet_n2b1[Idx_1];
   T_fatJet1_n3b1 = fatJet_n3b1[Idx_1];
   T_fatJet1_Tau3OverTau2 = fatJet_Tau3OverTau2[Idx_1];
   T_fatJet1_HasMuon = fatJet_HasMuon[Idx_1];
   T_fatJet1_HasElectron = fatJet_HasElectron[Idx_1];
   T_fatJet1_HasBJetMedium = fatJet_HasBJetMedium[Idx_1];
   T_fatJet1_OppositeHemisphereHasBJet = fatJet_OppositeHemisphereHasBJet[Idx_1];
   T_fatJet1_subJetIdx1_pt   = fatJet_subJetIdx1_pt[Idx_1];
   T_fatJet1_subJetIdx1_eta  = fatJet_subJetIdx1_eta[Idx_1];
   T_fatJet1_subJetIdx1_phi  = fatJet_subJetIdx1_phi[Idx_1];
   T_fatJet1_subJetIdx1_mass = fatJet_subJetIdx1_mass[Idx_1];
   T_fatJet1_subJetIdx2_pt   = fatJet_subJetIdx2_pt[Idx_1];
   T_fatJet1_subJetIdx2_eta  = fatJet_subJetIdx2_eta[Idx_1];
   T_fatJet1_subJetIdx2_phi  = fatJet_subJetIdx2_phi[Idx_1];
   T_fatJet1_subJetIdx2_mass = fatJet_subJetIdx2_mass[Idx_1];
   T_fatJet1_particleNet_Xbb_Legacy = fatJet_particleNet_Xbb_Legacy[Idx_1];
   T_fatJet1_particleNet_Xcc_Legacy = fatJet_particleNet_Xcc_Legacy[Idx_1];
   T_fatJet1_particleNet_Xqq_Legacy = fatJet_particleNet_Xqq_Legacy[Idx_1];
   T_fatJet1_particleNet_QCD_Legacy = fatJet_particleNet_QCD_Legacy[Idx_1];
   T_fatJet1_particleNet_QCDb_Legacy = fatJet_particleNet_QCDb_Legacy[Idx_1];
   T_fatJet1_particleNet_QCDbb_Legacy = fatJet_particleNet_QCDbb_Legacy[Idx_1];
   T_fatJet1_particleNet_QCDothers_Legacy = fatJet_particleNet_QCDothers_Legacy[Idx_1];
   T_fatJet1_particleNet_mass_Legacy = fatJet_particleNet_mass_Legacy[Idx_1];

   T_fatJet2_pt = fatJet_pt[Idx_2];
   T_fatJet2_eta = fatJet_eta[Idx_2];
   T_fatJet2_phi = fatJet_phi[Idx_2];
   T_fatJet2_mass = fatJet_mass[Idx_2];
   T_fatJet2_msoftdrop = fatJet_msoftdrop[Idx_2];
   T_fatJet2_btagDDBvLV2 = fatJet_btagDDBvLV2[Idx_2];
   T_fatJet2_btagDeepB = fatJet_btagDeepB[Idx_2];
   T_fatJet2_particleNetWithMass_H4qvsQCD = fatJet_particleNetWithMass_H4qvsQCD[Idx_2];
   T_fatJet2_particleNetWithMass_HbbvsQCD = fatJet_particleNetWithMass_HbbvsQCD[Idx_2];
   T_fatJet2_particleNetWithMass_QCD = fatJet_particleNetWithMass_QCD[Idx_2];
   T_fatJet2_particleNetWithMass_TvsQCD = fatJet_particleNetWithMass_TvsQCD[Idx_2];
   T_fatJet2_particleNetWithMass_WvsQCD = fatJet_particleNetWithMass_WvsQCD[Idx_2];
   T_fatJet2_particleNetWithMass_ZvsQCD = fatJet_particleNetWithMass_ZvsQCD[Idx_2];
   T_fatJet2_particleNet_QCD = fatJet_particleNet_QCD[Idx_2];
   T_fatJet2_particleNet_QCD0HF = fatJet_particleNet_QCD0HF[Idx_2];
   T_fatJet2_particleNet_QCD1HF = fatJet_particleNet_QCD1HF[Idx_2];
   T_fatJet2_particleNet_QCD2HF = fatJet_particleNet_QCD2HF[Idx_2];
   T_fatJet2_particleNet_XbbVsQCD = fatJet_particleNet_XbbVsQCD[Idx_2];
   T_fatJet2_particleNet_XccVsQCD = fatJet_particleNet_XccVsQCD[Idx_2];
   T_fatJet2_particleNet_XggVsQCD = fatJet_particleNet_XggVsQCD[Idx_2];
   T_fatJet2_particleNet_XqqVsQCD = fatJet_particleNet_XqqVsQCD[Idx_2];
   T_fatJet2_particleNet_massCorr = fatJet_mass[Idx_2]*fatJet_particleNet_massCorr[Idx_2];
   T_fatJet2_n2b1 = fatJet_n2b1[Idx_2];
   T_fatJet2_n3b1 = fatJet_n3b1[Idx_2];
   T_fatJet2_Tau3OverTau2 = fatJet_Tau3OverTau2[Idx_2];
   T_fatJet2_HasMuon = fatJet_HasMuon[Idx_2];
   T_fatJet2_HasElectron = fatJet_HasElectron[Idx_2];
   T_fatJet2_HasBJetMedium = fatJet_HasBJetMedium[Idx_2];
   T_fatJet2_OppositeHemisphereHasBJet = fatJet_OppositeHemisphereHasBJet[Idx_2];
   T_fatJet2_subJetIdx1_pt   = fatJet_subJetIdx1_pt[Idx_2];
   T_fatJet2_subJetIdx1_eta  = fatJet_subJetIdx1_eta[Idx_2];
   T_fatJet2_subJetIdx1_phi  = fatJet_subJetIdx1_phi[Idx_2];
   T_fatJet2_subJetIdx1_mass = fatJet_subJetIdx1_mass[Idx_2];
   T_fatJet2_subJetIdx2_pt   = fatJet_subJetIdx2_pt[Idx_2];
   T_fatJet2_subJetIdx2_eta  = fatJet_subJetIdx2_eta[Idx_2];
   T_fatJet2_subJetIdx2_phi  = fatJet_subJetIdx2_phi[Idx_2];
   T_fatJet2_subJetIdx2_mass = fatJet_subJetIdx2_mass[Idx_2];
   T_fatJet2_particleNet_Xbb_Legacy = fatJet_particleNet_Xbb_Legacy[Idx_2];
   T_fatJet2_particleNet_Xcc_Legacy = fatJet_particleNet_Xcc_Legacy[Idx_2];
   T_fatJet2_particleNet_Xqq_Legacy = fatJet_particleNet_Xqq_Legacy[Idx_2];
   T_fatJet2_particleNet_QCD_Legacy = fatJet_particleNet_QCD_Legacy[Idx_2];
   T_fatJet2_particleNet_QCDb_Legacy = fatJet_particleNet_QCDb_Legacy[Idx_2];
   T_fatJet2_particleNet_QCDbb_Legacy = fatJet_particleNet_QCDbb_Legacy[Idx_2];
   T_fatJet2_particleNet_QCDothers_Legacy = fatJet_particleNet_QCDothers_Legacy[Idx_2];
   T_fatJet2_particleNet_mass_Legacy = fatJet_particleNet_mass_Legacy[Idx_2];

   T_fatJet3_pt = fatJet_pt[Idx_3];
   T_fatJet3_eta = fatJet_eta[Idx_3];
   T_fatJet3_phi = fatJet_phi[Idx_3];
   T_fatJet3_mass = fatJet_mass[Idx_3];
   T_fatJet3_msoftdrop = fatJet_msoftdrop[Idx_3];
   T_fatJet3_particleNet_XbbVsQCD = fatJet_particleNet_XbbVsQCD[Idx_3];

   // Other variables
   T_NJets = NJets;
   T_nBTaggedJets = nBTaggedJets;
   T_MET = MET;
   // jj-system
   TLorentzVector p4_fatJet1,p4_fatJet2;
   p4_fatJet1.SetPtEtaPhiM(fatJet_pt[Idx_1],fatJet_eta[Idx_1],fatJet_phi[Idx_1],fatJet_msoftdrop[Idx_1]);
   p4_fatJet2.SetPtEtaPhiM(fatJet_pt[Idx_2],fatJet_eta[Idx_2],fatJet_phi[Idx_2],fatJet_msoftdrop[Idx_2]);
   TLorentzVector p4_System = p4_fatJet1 + p4_fatJet2;
   T_ptFJ12   = p4_System.Pt();
   T_etaFJ12  = p4_System.Eta();
   T_phiFJ12  = p4_System.Phi();
   T_massFJ12 = p4_System.M();
   T_detaFJ12 = fabs(fatJet_eta[Idx_1] - fatJet_eta[Idx_2]);
   T_dphiFJ12 = phi_dist(fatJet_phi[Idx_1],fatJet_phi[Idx_2]);
   T_dRFJ12   = sqrt(T_detaFJ12*T_detaFJ12 + T_dphiFJ12*T_dphiFJ12);
   T_FJ1_PtOverMFJ12    = T_fatJet1_pt/T_massFJ12;
   T_FJ2_PtOverMFJ12    = T_fatJet2_pt/T_massFJ12;
   T_FJ1_PtOverMSDFJ1   = T_fatJet1_pt/T_fatJet1_msoftdrop;
   T_FJ2_PtOverMSDFJ2   = T_fatJet2_pt/T_fatJet2_msoftdrop;
   T_ptFJ1OverptFJ2     = T_fatJet1_pt/T_fatJet2_pt;
   T_massFJ1OvermassFJ2 = T_fatJet1_msoftdrop/T_fatJet2_msoftdrop;
   T_ptak4J1     = ak4J_pt[0];
   T_ptak4J2     = ak4J_pt[1];
   T_dRak4J1FJ1 = -99.0;
   T_dRak4J1FJ2 = -99.0;
   if(ak4J_pt[0] >40)
    {
     T_dRak4J1FJ1 = sqrt( pow((T_fatJet1_eta - ak4J_eta[0]),2) + pow(phi_dist(T_fatJet1_phi,ak4J_phi[0]),2)  );
     T_dRak4J1FJ2 = sqrt( pow((T_fatJet2_eta - ak4J_eta[0]),2) + pow(phi_dist(T_fatJet2_phi,ak4J_phi[0]),2)  );
    }
   T_dRak4J2FJ1 = -99.0;
   T_dRak4J2FJ2 = -99.0;
   if(ak4J_pt[1] >40)
    {
     T_dRak4J2FJ1 = sqrt( pow((T_fatJet1_eta - ak4J_eta[1]),2) + pow(phi_dist(T_fatJet1_phi,ak4J_phi[1]),2)  );
     T_dRak4J2FJ2 = sqrt( pow((T_fatJet2_eta - ak4J_eta[1]),2) + pow(phi_dist(T_fatJet2_phi,ak4J_phi[1]),2)  );
    }

   // Mixing Angle hh
   p4_fatJet1.Boost(-p4_System.BoostVector());
   p4_fatJet2.Boost(-p4_System.BoostVector());
   TVector3 P3_fatJet1(p4_fatJet1.Px(),p4_fatJet1.Py(),p4_fatJet1.Pz());
   TVector3 P3_fatJet2(p4_fatJet2.Px(),p4_fatJet2.Py(),p4_fatJet2.Pz());
   TVector3 System_Direction(p4_System.Px(),p4_System.Py(),p4_System.Pz());
   double alpha_FatJets = min(P3_fatJet1.Angle(System_Direction),P3_fatJet2.Angle(System_Direction));
   T_alpha_FJ12 = alpha_FatJets;

   // Mixing Angle h1
   TLorentzVector p4_fatJet1_subJet1,p4_fatJet1_subJet2;
   p4_fatJet1_subJet1.SetPtEtaPhiM(T_fatJet1_subJetIdx1_pt, T_fatJet1_subJetIdx1_eta, T_fatJet1_subJetIdx1_phi, T_fatJet1_subJetIdx1_mass);
   p4_fatJet1_subJet2.SetPtEtaPhiM(T_fatJet1_subJetIdx2_pt, T_fatJet1_subJetIdx2_eta, T_fatJet1_subJetIdx2_phi, T_fatJet1_subJetIdx2_mass);
   TLorentzVector p4_fatJet1_subJets_system = p4_fatJet1_subJet1 + p4_fatJet1_subJet2;
   p4_fatJet1_subJet1.Boost(-p4_fatJet1_subJets_system.BoostVector());
   p4_fatJet1_subJet2.Boost(-p4_fatJet1_subJets_system.BoostVector());
   TVector3 P3_fatJet1_SubJ1(p4_fatJet1_subJet1.Px(), p4_fatJet1_subJet1.Py(), p4_fatJet1_subJet1.Pz());
   TVector3 P3_fatJet1_SubJ2(p4_fatJet1_subJet2.Px(), p4_fatJet1_subJet2.Py(), p4_fatJet1_subJet2.Pz());
   TVector3 fatJet1_SubJets_System_Direction(p4_fatJet1_subJets_system.Px(), p4_fatJet1_subJets_system.Py(), p4_fatJet1_subJets_system.Pz());
   double fatJet1_subJets_alpha = min(P3_fatJet1_SubJ1.Angle(fatJet1_SubJets_System_Direction), P3_fatJet1_SubJ2.Angle(fatJet1_SubJets_System_Direction));
   T_fatJet1_subJets_alpha = fatJet1_subJets_alpha;

   // Mixing Angle h2
   TLorentzVector p4_fatJet2_subJet1,p4_fatJet2_subJet2;
   p4_fatJet2_subJet1.SetPtEtaPhiM(T_fatJet2_subJetIdx1_pt, T_fatJet2_subJetIdx1_eta, T_fatJet2_subJetIdx1_phi, T_fatJet2_subJetIdx1_mass);
   p4_fatJet2_subJet2.SetPtEtaPhiM(T_fatJet2_subJetIdx2_pt, T_fatJet2_subJetIdx2_eta, T_fatJet2_subJetIdx2_phi, T_fatJet2_subJetIdx2_mass);
   TLorentzVector p4_fatJet2_subJets_system = p4_fatJet2_subJet1 + p4_fatJet2_subJet2;
   p4_fatJet2_subJet1.Boost(-p4_fatJet2_subJets_system.BoostVector());
   p4_fatJet2_subJet2.Boost(-p4_fatJet2_subJets_system.BoostVector());
   TVector3 P3_fatJet2_SubJ1(p4_fatJet2_subJet1.Px(), p4_fatJet2_subJet1.Py(), p4_fatJet2_subJet1.Pz());
   TVector3 P3_fatJet2_SubJ2(p4_fatJet2_subJet2.Px(), p4_fatJet2_subJet2.Py(), p4_fatJet2_subJet2.Pz());
   TVector3 fatJet2_SubJets_System_Direction(p4_fatJet2_subJets_system.Px(), p4_fatJet2_subJets_system.Py(), p4_fatJet2_subJets_system.Pz());
   double fatJet2_subJets_alpha = min(P3_fatJet2_SubJ1.Angle(fatJet2_SubJets_System_Direction), P3_fatJet2_SubJ2.Angle(fatJet2_SubJets_System_Direction));
   T_fatJet2_subJets_alpha = fatJet2_subJets_alpha;

// ********************************************************** Fill the Tree
   outputTree->Fill();

 } // end event loop

f->Write();

}


