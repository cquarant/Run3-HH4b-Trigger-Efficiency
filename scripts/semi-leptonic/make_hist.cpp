#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

void make_hist(const std::string &input_path, // path to the input root file
               const std::string &output_path // path to the output root file
) {
  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  TH1D *_FatJet1_pt = new TH1D("FatJet1_pt", "FatJet1_pt", 200, 0, 1200);
  TH1D *_FatJet1_eta = new TH1D("FatJet1_eta", "FatJet1_eta", 100, -5, 5);
  TH1D *_FatJet1_MassSD =
      new TH1D("FatJet1_MassSD", "FatJet1_MassSD", 500, 0, 500);
  TH1D *_FatJet1_ParticleNetLegacy_XbbVsQCD =
      new TH1D("FatJet1_ParticleNetLegacy_XbbVsQCD",
               "FatJet1_ParticleNetLegacy_XbbVsQCD", 100, 0.0, 1.0);
  TH1D *_FatJet1_GloParT_XbbVsQCD = new TH1D(
      "FatJet1_GloParT_XbbVsQCD", "FatJet1_GloParT_XbbVsQCD", 100, 0.0, 1.0);
  TH1D *_FatJet1_Tau3OverTau2 =
      new TH1D("FatJet1_Tau3OverTau2", "FatJet1_Tau3OverTau2", 100, 0.0, 1.0);
  TH1D *_FatJet2_pt = new TH1D("FatJet2_pt", "FatJet2_pt", 200, 0, 1200);
  TH1D *_FatJet2_eta = new TH1D("FatJet2_eta", "FatJet2_eta", 100, -5, 5);
  TH1D *_FatJet2_MassSD =
      new TH1D("FatJet2_MassSD", "FatJet2_MassSD", 500, 0, 500);

  TH1D *_MET = new TH1D("MET", "MET", 100, 0, 500);
  TH1D *_lep1_Pt = new TH1D("lep1_Pt", "lep1_Pt", 300, 0, 300);
  TH1D *_dR_LFJ = new TH1D("dR_LFJ", "dR_LFJ", 100, 0, 10);
  TH1D *_dR_J1FJ = new TH1D("dR_J1FJ", "dR_J1FJ", 100, 0, 10);
  TH1D *_dR_J2FJ = new TH1D("dR_J2FJ", "dR_J2FJ", 100, 0, 10);
  TH1D *_dR_JmaxL = new TH1D("dR_JmaxL", "dR_JmaxL", 100, 0, 10);

  TFile *f0 = TFile::Open(input_path.c_str());
  TTree *ntuples = (TTree *)f0->Get("tree");

  Float_t weight;
  Float_t fatJet1_pt, fatJet1_eta, fatJet1_phi, fatJet1_msoftdrop;
  Float_t fatJet1_ParticleNetLegacy_XbbVsQCD;
  Float_t fatJet1_GloParT_XbbVsQCD;
  Float_t fatJet1_Tau3OverTau2;
  Float_t fatJet2_pt, fatJet2_eta, fatJet2_msoftdrop;
  Float_t MET, lep1_Pt, dR_LFJ, dR_J1FJ, dR_J2FJ, dR_JmaxL;

  ntuples->SetBranchAddress("weight", &weight);
  ntuples->SetBranchAddress("fatJet1_pt", &fatJet1_pt);
  ntuples->SetBranchAddress("fatJet1_eta", &fatJet1_eta);
  ntuples->SetBranchAddress("fatJet1_msoftdrop", &fatJet1_msoftdrop);
  ntuples->SetBranchAddress("fatJet1_ParticleNetLegacy_XbbVsQCD",
                            &fatJet1_ParticleNetLegacy_XbbVsQCD);
  ntuples->SetBranchAddress("fatJet1_Tau3OverTau2", &fatJet1_Tau3OverTau2);
  ntuples->SetBranchAddress("fatJet1_GloParT_XbbVsQCD", &fatJet1_GloParT_XbbVsQCD);
  ntuples->SetBranchAddress("fatJet2_pt", &fatJet2_pt);
  ntuples->SetBranchAddress("fatJet2_eta", &fatJet2_eta);
  ntuples->SetBranchAddress("fatJet2_msoftdrop", &fatJet2_msoftdrop);
  ntuples->SetBranchAddress("MET", &MET);
  ntuples->SetBranchAddress("lep1_Pt", &lep1_Pt);
  ntuples->SetBranchAddress("dR_LFJ", &dR_LFJ);
  ntuples->SetBranchAddress("dR_J1FJ", &dR_J1FJ);
  ntuples->SetBranchAddress("dR_J2FJ", &dR_J2FJ);
  ntuples->SetBranchAddress("dR_JmaxL", &dR_JmaxL);

  double pt_MIN = 250.0;

  for (int i = 0; i < ntuples->GetEntries(); i++) {
    ntuples->GetEntry(i);

    if (dR_J1FJ < 0.0 && dR_J2FJ < 0.0) {
      continue;
    }
    if (dR_JmaxL > 3.5) {
      continue;
    }

    _FatJet1_pt->Fill(fatJet1_pt, weight);
    _FatJet1_eta->Fill(fatJet1_eta, weight);
    _FatJet1_MassSD->Fill(fatJet1_msoftdrop, weight);
    _FatJet1_ParticleNetLegacy_XbbVsQCD->Fill(
        fatJet1_ParticleNetLegacy_XbbVsQCD, weight);
    _FatJet1_GloParT_XbbVsQCD->Fill(fatJet1_GloParT_XbbVsQCD, weight);
    _FatJet1_Tau3OverTau2->Fill(fatJet1_Tau3OverTau2, weight);
    _FatJet2_pt->Fill(fatJet2_pt, weight);
    _FatJet2_eta->Fill(fatJet2_eta, weight);
    _FatJet2_MassSD->Fill(fatJet2_msoftdrop, weight);
    _MET->Fill(MET, weight);
    _lep1_Pt->Fill(lep1_Pt, weight);
    _dR_LFJ->Fill(dR_LFJ, weight);
    _dR_J1FJ->Fill(dR_J1FJ, weight);
    _dR_J2FJ->Fill(dR_J2FJ, weight);
    _dR_JmaxL->Fill(dR_JmaxL, weight);
  }

  f->Write();
}
