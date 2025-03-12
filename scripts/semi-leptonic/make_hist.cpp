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
#include "TXbb.h"

#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

void make_hist(
    const std::string &input_path,         // path to the input root file
    const std::string &output_path,        // path to the output root file
    const std::string &sf_path_tau32 = "", // path to the Tau32 SF root file
    const std::string &sf_path_txbb = ""   // path to the Xbb SF root file
) {
  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  std::cout << "Input Path: " << input_path << std::endl;
  std::cout << "Output Path: " << output_path << std::endl;

  bool apply_sf_tau32 = false;
  bool apply_sf_txbb = false;
  TH1D *_SF_Tau32 = nullptr;
  TH1D *_SF_Xbb = nullptr;
  TFile *f_Tau3toTau2 = nullptr;
  TFile *f_Xbb = nullptr;
  if (sf_path_tau32 != "") {
      f_Tau3toTau2 = new TFile(sf_path_tau32.c_str());
      _SF_Tau32 = (TH1D *)f_Tau3toTau2->Get("SF");
      if (!_SF_Tau32) {
        throw std::runtime_error("Could not load Tau32 SF histogram");
      }
      std::cout << "Tau32 SF Path: " << sf_path_tau32 << std::endl;
      apply_sf_tau32 = true;
    }

  if (sf_path_txbb != "") {
    f_Xbb = new TFile(sf_path_txbb.c_str());
    _SF_Xbb = (TH1D *)f_Xbb->Get("SF");
    if (!_SF_Xbb) {
      throw std::runtime_error("Could not load Xbb SF histogram");
    }
    std::cout << "Xbb SF Path: " << sf_path_txbb << std::endl;
    apply_sf_txbb = true;
  }

  TH1D *_FatJet1_pt = new TH1D("FatJet1_pt", "FatJet1_pt", 200, 0, 1200);
  TH1D *_FatJet1_eta = new TH1D("FatJet1_eta", "FatJet1_eta", 100, -5, 5);
  TH1D *_FatJet1_MassSD =
      new TH1D("FatJet1_MassSD", "FatJet1_MassSD", 500, 0, 500);
  TH1D *_FatJet1_ParticleNetLegacy_XbbVsQCD =
      new TH1D("FatJet1_ParticleNetLegacy_XbbVsQCD",
               "FatJet1_ParticleNetLegacy_XbbVsQCD", N_TXbb, 0.0, 1.0);
  TH1D *_FatJet1_GloParT_XbbVsQCD = new TH1D(
      "FatJet1_GloParT_XbbVsQCD", "FatJet1_GloParT_XbbVsQCD", N_TXbb, 0.0, 1.0);
  TH1D *_FatJet1_Tau3OverTau2 =
      new TH1D("FatJet1_Tau3OverTau2", "FatJet1_Tau3OverTau2", N_TXbb, 0.0, 1.0);
  TH1D *_FatJet2_pt = new TH1D("FatJet2_pt", "FatJet2_pt", 200, 0, 1200);
  TH1D *_FatJet2_eta = new TH1D("FatJet2_eta", "FatJet2_eta", 100, -5, 5);
  TH1D *_FatJet2_MassSD =
      new TH1D("FatJet2_MassSD", "FatJet2_MassSD", 500, 0, 500);

  TH1D *_MET = new TH1D("MET", "MET", 100, 0, 500);
  TH1D *_lep1_pt = new TH1D("lep1_pt", "lep1_pt", 300, 0, 300);
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
  Float_t MET, lep1_pt, dR_LFJ, dR_J1FJ, dR_J2FJ, dR_JmaxL;

  ntuples->SetBranchAddress("weight", &weight);
  ntuples->SetBranchAddress("fatJet1_pt", &fatJet1_pt);
  ntuples->SetBranchAddress("fatJet1_eta", &fatJet1_eta);
  ntuples->SetBranchAddress("fatJet1_msoftdrop", &fatJet1_msoftdrop);
  ntuples->SetBranchAddress("fatJet1_ParticleNetLegacy_XbbVsQCD",
                            &fatJet1_ParticleNetLegacy_XbbVsQCD);
  ntuples->SetBranchAddress("fatJet1_Tau3OverTau2", &fatJet1_Tau3OverTau2);
  ntuples->SetBranchAddress("fatJet1_GloParT_XbbVsQCD",
                            &fatJet1_GloParT_XbbVsQCD);
  ntuples->SetBranchAddress("fatJet2_pt", &fatJet2_pt);
  ntuples->SetBranchAddress("fatJet2_eta", &fatJet2_eta);
  ntuples->SetBranchAddress("fatJet2_msoftdrop", &fatJet2_msoftdrop);
  ntuples->SetBranchAddress("MET", &MET);
  ntuples->SetBranchAddress("lep1_pt", &lep1_pt);
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

    // adding Xbb SF
    Float_t sf_TXbb = 1.0;
    if (apply_sf_txbb) {
      Int_t bin_TXbb = _SF_Xbb->GetXaxis()->FindBin(fatJet1_GloParT_XbbVsQCD);
      sf_TXbb = _SF_Xbb->GetBinContent(bin_TXbb);
      if (sf_TXbb <= 0) {
        sf_TXbb = 1.0;
      }
    }
    weight = weight * sf_TXbb;

    Float_t sf_tau32 = 1.0;
    if (apply_sf_tau32) {
      Int_t bin_tau32 = _SF_Tau32->GetXaxis()->FindBin(fatJet1_Tau3OverTau2);
      sf_tau32 = _SF_Tau32->GetBinContent(bin_tau32);
      if (sf_tau32 <= 0) {
        sf_tau32 = 1.0;
      }
    }
    weight = weight * sf_tau32;

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
    _lep1_pt->Fill(lep1_pt, weight);
    _dR_LFJ->Fill(dR_LFJ, weight);
    _dR_J1FJ->Fill(dR_J1FJ, weight);
    _dR_J2FJ->Fill(dR_J2FJ, weight);
    _dR_JmaxL->Fill(dR_JmaxL, weight);
  }

  f0->Close();  // Close input file
  f->cd();      // Make sure we're writing to the output file
  
  // Write each histogram explicitly
  _FatJet1_pt->Write();
  _FatJet1_eta->Write();
  _FatJet1_MassSD->Write();
  _FatJet1_ParticleNetLegacy_XbbVsQCD->Write();
  _FatJet1_GloParT_XbbVsQCD->Write();
  _FatJet1_Tau3OverTau2->Write();
  _FatJet2_pt->Write();
  _FatJet2_eta->Write();
  _FatJet2_MassSD->Write();
  _MET->Write();
  _lep1_pt->Write();
  _dR_LFJ->Write();
  _dR_J1FJ->Write();
  _dR_J2FJ->Write();
  _dR_JmaxL->Write();
  
  f->Write();
  
  delete f;
}
