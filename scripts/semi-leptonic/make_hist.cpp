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

#define GLOPART_MASS_MIN 0
#define GLOPART_MASS_MAX 350
#define GLOPART_MASS_STEP 1


void updateHistErr(TH1D* hist, Float_t value, Float_t err, Float_t weight) {
    // Find the bin corresponding to the value
    int bin = hist->FindBin(value);
    
    // Update the bin error using quadrature sum
    double curr_err = hist->GetBinError(bin);
    double new_err = sqrt(pow(curr_err, 2) + pow(weight * err, 2));
    hist->SetBinError(bin, new_err);
}

void make_hist(
    const std::string &input_path,         // path to the input root file
    const std::string &output_path,        // path to the output root file
    const std::string &sf_path_tau32 = "", // path to the Tau32 SF root file
    const std::string &sf_path_txbb = ""  // path to the Xbb SF root file
) {
  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  std::cout << "Input Path: " << input_path << std::endl;
  std::cout << "Output Path: " << output_path << std::endl;

  bool apply_sf_tau32 = false;
  bool apply_sf_TXbb = false;
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
    apply_sf_TXbb = true;
  }

  int N_GloParT_Mass = (GLOPART_MASS_MAX - GLOPART_MASS_MIN) / GLOPART_MASS_STEP;
  TH1D *_FatJet1_pt = new TH1D("FatJet1_pt", "FatJet1_pt", 200, 0, 1200);
  TH1D *_FatJet1_eta = new TH1D("FatJet1_eta", "FatJet1_eta", 100, -5, 5);
  TH1D *_FatJet1_MassSD =
      new TH1D("FatJet1_MassSD", "FatJet1_MassSD", 500, 0, 500);
  TH1D *_FatJet1_GloParT_MassVis =
      new TH1D("FatJet1_GloParT_MassVis", "FatJet1_GloParT_MassVis", 
               N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
  TH1D *_FatJet1_GloParT_MassRes =
      new TH1D("FatJet1_GloParT_MassRes", "FatJet1_GloParT_MassRes", 
               N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
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
  TH1D *_FatJet2_GloParT_MassVis =
      new TH1D("FatJet2_GloParT_MassVis", "FatJet2_GloParT_MassVis", 
               N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
  TH1D *_FatJet2_GloParT_MassRes =
      new TH1D("FatJet2_GloParT_MassRes", "FatJet2_GloParT_MassRes", 
               N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);

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
  Float_t fatJet1_GloParT_massVis, fatJet1_GloParT_massRes;
  Float_t fatJet1_ParticleNetLegacy_XbbVsQCD;
  Float_t fatJet1_GloParT_XbbVsQCD;
  Float_t fatJet1_Tau3OverTau2;
  Float_t fatJet2_pt, fatJet2_eta, fatJet2_msoftdrop;
  Float_t fatJet2_GloParT_massVis, fatJet2_GloParT_massRes;
  Float_t MET, lep1_pt, dR_LFJ, dR_J1FJ, dR_J2FJ, dR_JmaxL;

  ntuples->SetBranchAddress("weight", &weight);
  ntuples->SetBranchAddress("fatJet1_pt", &fatJet1_pt);
  ntuples->SetBranchAddress("fatJet1_eta", &fatJet1_eta);
  ntuples->SetBranchAddress("fatJet1_msoftdrop", &fatJet1_msoftdrop);
  ntuples->SetBranchAddress("fatJet1_GloParT_massVis", &fatJet1_GloParT_massVis);
  ntuples->SetBranchAddress("fatJet1_GloParT_massRes", &fatJet1_GloParT_massRes);
  ntuples->SetBranchAddress("fatJet1_ParticleNetLegacy_XbbVsQCD",
                            &fatJet1_ParticleNetLegacy_XbbVsQCD);
  ntuples->SetBranchAddress("fatJet1_Tau3OverTau2", &fatJet1_Tau3OverTau2);
  ntuples->SetBranchAddress("fatJet1_GloParT_XbbVsQCD",
                            &fatJet1_GloParT_XbbVsQCD);
  ntuples->SetBranchAddress("fatJet2_pt", &fatJet2_pt);
  ntuples->SetBranchAddress("fatJet2_eta", &fatJet2_eta);
  ntuples->SetBranchAddress("fatJet2_msoftdrop", &fatJet2_msoftdrop);
  ntuples->SetBranchAddress("fatJet2_GloParT_massVis", &fatJet2_GloParT_massVis);
  ntuples->SetBranchAddress("fatJet2_GloParT_massRes", &fatJet2_GloParT_massRes);
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
    Float_t sf_TXbb_err = 0.0;
    if (apply_sf_TXbb) {
      Int_t bin_TXbb = _SF_Xbb->GetXaxis()->FindBin(fatJet1_GloParT_XbbVsQCD);
      sf_TXbb = _SF_Xbb->GetBinContent(bin_TXbb);
      sf_TXbb_err = _SF_Xbb->GetBinError(bin_TXbb);
      if (sf_TXbb <= 0 || sf_TXbb_err < 0) {
        sf_TXbb = 1.0;
        sf_TXbb_err = 0.0;
      }
    }

    Float_t sf_tau32 = 1.0;
    Float_t sf_tau32_err = 0.0;
    if (apply_sf_tau32) {
      Int_t bin_tau32 = _SF_Tau32->GetXaxis()->FindBin(fatJet1_Tau3OverTau2);
      sf_tau32 = _SF_Tau32->GetBinContent(bin_tau32);
      sf_tau32_err = _SF_Tau32->GetBinError(bin_tau32);
      if (sf_tau32 <= 0 || sf_tau32_err < 0) {
        sf_tau32 = 1.0;
        sf_tau32_err = 0.0;
      }
    }

    weight = weight * sf_TXbb * sf_tau32;
    Float_t sf_err_sq = 0.0;
    if (apply_sf_TXbb) {
        // partial derivative of weight w.r.t. sf_TXbb
        sf_err_sq += pow(sf_TXbb_err * sf_tau32, 2);
    }
    if (apply_sf_tau32) {
      // partial derivative of weight w.r.t. sf_tau32
        sf_err_sq += pow(sf_TXbb * sf_tau32_err, 2);
    }
    Float_t sf_err = sqrt(sf_err_sq);

    _FatJet1_pt->Fill(fatJet1_pt, weight);
    _FatJet1_eta->Fill(fatJet1_eta, weight);
    _FatJet1_MassSD->Fill(fatJet1_msoftdrop, weight);
    _FatJet1_GloParT_MassVis->Fill(fatJet1_GloParT_massVis, weight);
    _FatJet1_GloParT_MassRes->Fill(fatJet1_GloParT_massRes, weight);
    _FatJet1_ParticleNetLegacy_XbbVsQCD->Fill(
        fatJet1_ParticleNetLegacy_XbbVsQCD, weight);
    _FatJet1_GloParT_XbbVsQCD->Fill(fatJet1_GloParT_XbbVsQCD, weight);
    _FatJet1_Tau3OverTau2->Fill(fatJet1_Tau3OverTau2, weight);
    _FatJet2_pt->Fill(fatJet2_pt, weight);
    _FatJet2_eta->Fill(fatJet2_eta, weight);
    _FatJet2_MassSD->Fill(fatJet2_msoftdrop, weight);
    _FatJet2_GloParT_MassVis->Fill(fatJet2_GloParT_massVis, weight);
    _FatJet2_GloParT_MassRes->Fill(fatJet2_GloParT_massRes, weight);
    _MET->Fill(MET, weight);
    _lep1_pt->Fill(lep1_pt, weight);
    _dR_LFJ->Fill(dR_LFJ, weight);
    _dR_J1FJ->Fill(dR_J1FJ, weight);
    _dR_J2FJ->Fill(dR_J2FJ, weight);
    _dR_JmaxL->Fill(dR_JmaxL, weight);

    if (apply_sf_TXbb || apply_sf_tau32) {
      updateHistErr(_FatJet1_pt, fatJet1_pt, sf_err, weight);
      updateHistErr(_FatJet1_eta, fatJet1_eta, sf_err, weight);
      updateHistErr(_FatJet1_MassSD, fatJet1_msoftdrop, sf_err, weight);
      updateHistErr(_FatJet1_GloParT_MassVis, fatJet1_GloParT_massVis, sf_err, weight);
      updateHistErr(_FatJet1_GloParT_MassRes, fatJet1_GloParT_massRes, sf_err, weight);
      updateHistErr(_FatJet1_ParticleNetLegacy_XbbVsQCD, fatJet1_ParticleNetLegacy_XbbVsQCD, sf_err, weight);
      updateHistErr(_FatJet1_GloParT_XbbVsQCD, fatJet1_GloParT_XbbVsQCD, sf_err, weight);
      updateHistErr(_FatJet1_Tau3OverTau2, fatJet1_Tau3OverTau2, sf_err, weight);
      updateHistErr(_FatJet2_pt, fatJet2_pt, sf_err, weight);
      updateHistErr(_FatJet2_eta, fatJet2_eta, sf_err, weight);
      updateHistErr(_FatJet2_MassSD, fatJet2_msoftdrop, sf_err, weight);
      updateHistErr(_FatJet2_GloParT_MassVis, fatJet2_GloParT_massVis, sf_err, weight);
      updateHistErr(_FatJet2_GloParT_MassRes, fatJet2_GloParT_massRes, sf_err, weight);
      updateHistErr(_MET, MET, sf_err, weight);
      updateHistErr(_lep1_pt, lep1_pt, sf_err, weight);
      updateHistErr(_dR_LFJ, dR_LFJ, sf_err, weight);
      updateHistErr(_dR_J1FJ, dR_J1FJ, sf_err, weight);
      updateHistErr(_dR_J2FJ, dR_J2FJ, sf_err, weight);
      updateHistErr(_dR_JmaxL, dR_JmaxL, sf_err, weight);
    }
  }

  f0->Close();  // Close input file
  
  if (apply_sf_TXbb || apply_sf_tau32) {
    f->cd();      // Make sure we're writing to the output file
    _FatJet1_pt->Write();
    _FatJet1_eta->Write();
    _FatJet1_MassSD->Write();
    _FatJet1_GloParT_MassVis->Write();
    _FatJet1_GloParT_MassRes->Write();
    _FatJet1_ParticleNetLegacy_XbbVsQCD->Write();
    _FatJet1_GloParT_XbbVsQCD->Write();
    _FatJet1_Tau3OverTau2->Write();
    _FatJet2_pt->Write();
    _FatJet2_eta->Write();
    _FatJet2_MassSD->Write();
    _FatJet2_GloParT_MassVis->Write();
    _FatJet2_GloParT_MassRes->Write();
    _MET->Write();
    _lep1_pt->Write();
    _dR_LFJ->Write();
    _dR_J1FJ->Write();
    _dR_J2FJ->Write();
    _dR_JmaxL->Write();
  }
  f->Write();

  // Close scale factor files
  if (f_Tau3toTau2) {
    f_Tau3toTau2->Close();
    delete f_Tau3toTau2;
  }
  if (f_Xbb) {
    f_Xbb->Close();
    delete f_Xbb;
  }
  
  delete f;
}
