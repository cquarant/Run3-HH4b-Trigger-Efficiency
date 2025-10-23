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

#define GLOPART_MASS_MIN 90
#define GLOPART_MASS_MAX 160
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
    const std::string &sf_path_txbb = "",  // path to the TXbb SF root file
    const std::string &sf_path_ptjj = ""   // path to the ptjj SF root file
) {
  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  std::cout << "Input Path: " << input_path << std::endl;
  std::cout << "Output Path: " << output_path << std::endl;

  bool apply_sf_tau32 = !sf_path_tau32.empty();
  bool apply_sf_TXbb = !sf_path_txbb.empty();
  bool apply_sf_ptjj = !sf_path_ptjj.empty();
  
  TH1D *_SF_Tau32 = nullptr;
  TH1D *_SF_Xbb = nullptr;
  TH1D *_SF_ptjj = nullptr;
  TFile *f_Tau3toTau2 = nullptr;
  TFile *f_Xbb = nullptr;
  TFile *f_ptjj = nullptr;
  
  if (apply_sf_tau32) {
    f_Tau3toTau2 = new TFile(sf_path_tau32.c_str());
    _SF_Tau32 = (TH1D *)f_Tau3toTau2->Get("SF");
    if (!_SF_Tau32) {
      throw std::runtime_error("Could not load Tau32 SF histogram");
    }
    std::cout << "Tau32 SF Path: " << sf_path_tau32 << std::endl;
  }

  if (apply_sf_TXbb) {
    f_Xbb = new TFile(sf_path_txbb.c_str());
    _SF_Xbb = (TH1D *)f_Xbb->Get("SF");
    if (!_SF_Xbb) {
      throw std::runtime_error("Could not load Xbb SF histogram");
    }
    std::cout << "Xbb SF Path: " << sf_path_txbb << std::endl;
  }

  if (apply_sf_ptjj) {
    f_ptjj = new TFile(sf_path_ptjj.c_str());
    _SF_ptjj = (TH1D *)f_ptjj->Get("SF");
    if (!_SF_ptjj) {
      throw std::runtime_error("Could not load ptjj SF histogram");
    }
    std::cout << "ptjj SF Path: " << sf_path_ptjj << std::endl;
  }

  int N_GloParT_Mass = (GLOPART_MASS_MAX - GLOPART_MASS_MIN) / GLOPART_MASS_STEP;
  
  // FatJet1 histograms
  TH1D *_FatJet1_pt = new TH1D("FatJet1_pt", "FatJet1_pt", 200, 0, 1200);
  TH1D *_FatJet1_eta = new TH1D("FatJet1_eta", "FatJet1_eta", 100, -5, 5);
  TH1D *_FatJet1_MassSD = new TH1D("FatJet1_MassSD", "FatJet1_MassSD", 500, 0, 500);
  TH1D *_FatJet1_GloParT_MassVis = new TH1D("FatJet1_GloParT_MassVis", "FatJet1_GloParT_MassVis", 
                                             N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
  TH1D *_FatJet1_GloParT_MassRes = new TH1D("FatJet1_GloParT_MassRes", "FatJet1_GloParT_MassRes", 
                                             N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
  TH1D *_FatJet1_GloParT_XbbVsQCD = new TH1D("FatJet1_GloParT_XbbVsQCD", "FatJet1_GloParT_XbbVsQCD", 
                                              N_TXbb, 0.0, 1.0);
  TH1D *_FatJet1_Tau3OverTau2 = new TH1D("FatJet1_Tau3OverTau2", "FatJet1_Tau3OverTau2", 
                                          100, 0.0, 1.0);
  
  // FatJet2 histograms
  TH1D *_FatJet2_pt = new TH1D("FatJet2_pt", "FatJet2_pt", 200, 0, 1200);
  TH1D *_FatJet2_eta = new TH1D("FatJet2_eta", "FatJet2_eta", 100, -5, 5);
  TH1D *_FatJet2_MassSD = new TH1D("FatJet2_MassSD", "FatJet2_MassSD", 500, 0, 500);
  TH1D *_FatJet2_GloParT_MassVis = new TH1D("FatJet2_GloParT_MassVis", "FatJet2_GloParT_MassVis", 
                                             N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
  TH1D *_FatJet2_GloParT_MassRes = new TH1D("FatJet2_GloParT_MassRes", "FatJet2_GloParT_MassRes", 
                                             N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
  TH1D *_FatJet2_GloParT_XbbVsQCD = new TH1D("FatJet2_GloParT_XbbVsQCD", "FatJet2_GloParT_XbbVsQCD", 
                                              N_TXbb, 0.0, 1.0);
  TH1D *_FatJet2_Tau3OverTau2 = new TH1D("FatJet2_Tau3OverTau2", "FatJet2_Tau3OverTau2", 
                                          100, 0.0, 1.0);

  // pTjj histogram
  TH1D *_pTjj = new TH1D("pTjj", "pTjj", 200, 0, 1200);
  
  // Open and set up input file
  TFile *f0 = TFile::Open(input_path.c_str());
  TTree *ntuples = (TTree *)f0->Get("tree");

  // Declare variables to be read from tree
  Float_t weight;
  Float_t fatJet1_pt, fatJet1_eta, fatJet1_phi, fatJet1_msoftdrop;
  Float_t fatJet1_GloParT_massVis, fatJet1_GloParT_massRes;
  Float_t fatJet1_GloParT_XbbVsQCD, fatJet1_Tau3OverTau2;
  Float_t fatJet2_pt, fatJet2_eta, fatJet2_msoftdrop;
  Float_t fatJet2_GloParT_massVis, fatJet2_GloParT_massRes;
  Float_t fatJet2_GloParT_XbbVsQCD, fatJet2_Tau3OverTau2;
  Float_t pTjj;
  
  // Set branch addresses
  ntuples->SetBranchAddress("weight", &weight);
  ntuples->SetBranchAddress("fatJet1_pt", &fatJet1_pt);
  ntuples->SetBranchAddress("fatJet1_eta", &fatJet1_eta);
  ntuples->SetBranchAddress("fatJet1_msoftdrop", &fatJet1_msoftdrop);
  ntuples->SetBranchAddress("fatJet1_GloParT_massVis", &fatJet1_GloParT_massVis);
  ntuples->SetBranchAddress("fatJet1_GloParT_massRes", &fatJet1_GloParT_massRes);
  ntuples->SetBranchAddress("fatJet1_Tau3OverTau2", &fatJet1_Tau3OverTau2);
  ntuples->SetBranchAddress("fatJet1_GloParT_XbbVsQCD", &fatJet1_GloParT_XbbVsQCD);
  
  ntuples->SetBranchAddress("fatJet2_pt", &fatJet2_pt);
  ntuples->SetBranchAddress("fatJet2_eta", &fatJet2_eta);
  ntuples->SetBranchAddress("fatJet2_msoftdrop", &fatJet2_msoftdrop);
  ntuples->SetBranchAddress("fatJet2_GloParT_massVis", &fatJet2_GloParT_massVis);
  ntuples->SetBranchAddress("fatJet2_GloParT_massRes", &fatJet2_GloParT_massRes);
  ntuples->SetBranchAddress("fatJet2_Tau3OverTau2", &fatJet2_Tau3OverTau2);
  ntuples->SetBranchAddress("fatJet2_GloParT_XbbVsQCD", &fatJet2_GloParT_XbbVsQCD);
  ntuples->SetBranchAddress("pTjj", &pTjj);

  // Loop over all events
  for (int i = 0; i < ntuples->GetEntries(); i++) {
    ntuples->GetEntry(i);

    // Apply selection cuts
    if (fatJet1_pt <= 450 || fabs(fatJet1_eta) >= 2.5 || fatJet1_msoftdrop <= 50 || 
        fatJet1_GloParT_XbbVsQCD <= 0.1 || fatJet1_Tau3OverTau2 >= 0.46) continue;
    if (fatJet2_pt <= 450 || fabs(fatJet2_eta) >= 2.5 || fatJet2_msoftdrop <= 50 || 
        fatJet2_GloParT_XbbVsQCD <= 0.1 || fatJet2_Tau3OverTau2 >= 0.46) continue;

    // Mass cuts
    // if (fatJet1_GloParT_massVis < GLOPART_MASS_MIN || fatJet1_GloParT_massVis > GLOPART_MASS_MAX) continue;
    if (fatJet2_GloParT_massVis < GLOPART_MASS_MIN || fatJet2_GloParT_massVis > GLOPART_MASS_MAX) continue;

    // Get scale factors and their uncertainties
    Float_t sf_TXbb = 1.0, sf_TXbb_err = 0.0;
    Float_t sf_tau32 = 1.0, sf_tau32_err = 0.0;
    Float_t sf_ptjj = 1.0, sf_ptjj_err = 0.0;
    
    if (apply_sf_TXbb) {
      Int_t bin_TXbb = _SF_Xbb->GetXaxis()->FindBin(fatJet1_GloParT_XbbVsQCD);
      sf_TXbb = _SF_Xbb->GetBinContent(bin_TXbb);
      sf_TXbb_err = _SF_Xbb->GetBinError(bin_TXbb);
      if (sf_TXbb <= 0 || sf_TXbb_err < 0) {
        sf_TXbb = 1.0;
        sf_TXbb_err = 0.0;
      }
    }

    if (apply_sf_tau32) {
      Int_t bin_tau32 = _SF_Tau32->GetXaxis()->FindBin(fatJet1_Tau3OverTau2);
      sf_tau32 = _SF_Tau32->GetBinContent(bin_tau32);
      sf_tau32_err = _SF_Tau32->GetBinError(bin_tau32);
      if (sf_tau32 <= 0 || sf_tau32_err < 0) {
        sf_tau32 = 1.0;
        sf_tau32_err = 0.0;
      }
    }

    if (apply_sf_ptjj) {
      Int_t bin_ptjj = _SF_ptjj->GetXaxis()->FindBin(pTjj);
      sf_ptjj = _SF_ptjj->GetBinContent(bin_ptjj);
      sf_ptjj_err = _SF_ptjj->GetBinError(bin_ptjj);
      // std::cout << "pTjj: " << pTjj << ", bin: " << bin_ptjj << ", SF: " << sf_ptjj << " +/- " << sf_ptjj_err << std::endl;
    }

    // Apply scale factors to weight
    weight = weight * sf_TXbb * sf_tau32 * sf_ptjj;
    
    // Calculate combined systematic uncertainty using quadrature sum
    Float_t sf_err_sq = 0.0;
    if (apply_sf_TXbb) {
      sf_err_sq += pow(sf_TXbb_err * sf_tau32 * sf_ptjj, 2);
    }
    if (apply_sf_tau32) {
      sf_err_sq += pow(sf_TXbb * sf_tau32_err * sf_ptjj, 2);
    }
    if (apply_sf_ptjj) {
      sf_err_sq += pow(sf_TXbb * sf_tau32 * sf_ptjj_err, 2);
    }
    Float_t sf_err = sqrt(sf_err_sq);

    // Fill histograms with nominal weights
    _FatJet1_pt->Fill(fatJet1_pt, weight);
    _FatJet1_eta->Fill(fatJet1_eta, weight);
    _FatJet1_MassSD->Fill(fatJet1_msoftdrop, weight);
    _FatJet1_GloParT_MassVis->Fill(fatJet1_GloParT_massVis, weight);
    _FatJet1_GloParT_MassRes->Fill(fatJet1_GloParT_massRes, weight);
    _FatJet1_GloParT_XbbVsQCD->Fill(fatJet1_GloParT_XbbVsQCD, weight);
    _FatJet1_Tau3OverTau2->Fill(fatJet1_Tau3OverTau2, weight);
    
    _FatJet2_pt->Fill(fatJet2_pt, weight);
    _FatJet2_eta->Fill(fatJet2_eta, weight);
    _FatJet2_MassSD->Fill(fatJet2_msoftdrop, weight);
    _FatJet2_GloParT_MassVis->Fill(fatJet2_GloParT_massVis, weight);
    _FatJet2_GloParT_MassRes->Fill(fatJet2_GloParT_massRes, weight);
    _FatJet2_GloParT_XbbVsQCD->Fill(fatJet2_GloParT_XbbVsQCD, weight);
    _FatJet2_Tau3OverTau2->Fill(fatJet2_Tau3OverTau2, weight);

    _pTjj->Fill(pTjj, weight);

    // Update histogram errors with systematic uncertainties
    if (apply_sf_TXbb || apply_sf_tau32 || apply_sf_ptjj) {
      updateHistErr(_FatJet1_pt, fatJet1_pt, sf_err, weight);
      updateHistErr(_FatJet1_eta, fatJet1_eta, sf_err, weight);
      updateHistErr(_FatJet1_MassSD, fatJet1_msoftdrop, sf_err, weight);
      updateHistErr(_FatJet1_GloParT_MassVis, fatJet1_GloParT_massVis, sf_err, weight);
      updateHistErr(_FatJet1_GloParT_MassRes, fatJet1_GloParT_massRes, sf_err, weight);
      updateHistErr(_FatJet1_GloParT_XbbVsQCD, fatJet1_GloParT_XbbVsQCD, sf_err, weight);
      updateHistErr(_FatJet1_Tau3OverTau2, fatJet1_Tau3OverTau2, sf_err, weight);
      
      updateHistErr(_FatJet2_pt, fatJet2_pt, sf_err, weight);
      updateHistErr(_FatJet2_eta, fatJet2_eta, sf_err, weight);
      updateHistErr(_FatJet2_MassSD, fatJet2_msoftdrop, sf_err, weight);
      updateHistErr(_FatJet2_GloParT_MassVis, fatJet2_GloParT_massVis, sf_err, weight);
      updateHistErr(_FatJet2_GloParT_MassRes, fatJet2_GloParT_massRes, sf_err, weight);
      updateHistErr(_FatJet2_GloParT_XbbVsQCD, fatJet2_GloParT_XbbVsQCD, sf_err, weight);
      updateHistErr(_FatJet2_Tau3OverTau2, fatJet2_Tau3OverTau2, sf_err, weight);

      updateHistErr(_pTjj, pTjj, sf_err, weight);
    }
  }

  // Close input file
  f0->Close();

  // Write histograms to output file
  if (apply_sf_TXbb || apply_sf_tau32 || apply_sf_ptjj) {
    f->cd();  // Make sure we're writing to the output file
    _FatJet1_pt->Write();
    _FatJet1_eta->Write();
    _FatJet1_MassSD->Write();
    _FatJet1_GloParT_MassVis->Write();
    _FatJet1_GloParT_MassRes->Write();
    _FatJet1_GloParT_XbbVsQCD->Write();
    _FatJet1_Tau3OverTau2->Write();
    
    _FatJet2_pt->Write();
    _FatJet2_eta->Write();
    _FatJet2_MassSD->Write();
    _FatJet2_GloParT_MassVis->Write();
    _FatJet2_GloParT_MassRes->Write();
    _FatJet2_GloParT_XbbVsQCD->Write();
    _FatJet2_Tau3OverTau2->Write();

    _pTjj->Write();
  }
  
  f->Write();
  std::cout << "Wrote histograms to output file." << std::endl;
  
  // Clean up scale factor files
  if (f_Tau3toTau2) {
    f_Tau3toTau2->Close();
    delete f_Tau3toTau2;
  }
  if (f_Xbb) {
    f_Xbb->Close();
    delete f_Xbb;
  }
  if (f_ptjj) {
    f_ptjj->Close();
    delete f_ptjj;
  }
  
  // Close and save output file
  f->Close();
  delete f;
}