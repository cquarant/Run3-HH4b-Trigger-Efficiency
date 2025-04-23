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
  
  // Define the variations based on whether any scale factors are applied
  std::vector<std::string> variations = {"nominal"};
  
  // Only add up/down variations if any scale factors are applied
  bool apply_any_sf = apply_sf_tau32 || apply_sf_TXbb || apply_sf_ptjj;
  if (apply_any_sf) {
    variations.push_back("up");
    variations.push_back("down");
  }
  
  // Create map to store all histograms for each variation
  std::map<std::string, std::map<std::string, TH1D*>> hists;
  
  // Initialize histograms for all variations
  for (const auto& var : variations) {
    // FatJet1 histograms
    hists[var]["FatJet1_pt"] = new TH1D(("FatJet1_pt_" + var).c_str(), 
                                      "FatJet1_pt", 200, 0, 1200);
    hists[var]["FatJet1_eta"] = new TH1D(("FatJet1_eta_" + var).c_str(), 
                                       "FatJet1_eta", 100, -5, 5);
    hists[var]["FatJet1_MassSD"] = new TH1D(("FatJet1_MassSD_" + var).c_str(), 
                                          "FatJet1_MassSD", 500, 0, 500);
    hists[var]["FatJet1_GloParT_MassVis"] = new TH1D(("FatJet1_GloParT_MassVis_" + var).c_str(), 
                                                   "FatJet1_GloParT_MassVis", 
                                                   N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
    hists[var]["FatJet1_GloParT_MassRes"] = new TH1D(("FatJet1_GloParT_MassRes_" + var).c_str(), 
                                                   "FatJet1_GloParT_MassRes", 
                                                   N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
    hists[var]["FatJet1_GloParT_XbbVsQCD"] = new TH1D(("FatJet1_GloParT_XbbVsQCD_" + var).c_str(), 
                                                    "FatJet1_GloParT_XbbVsQCD", N_TXbb, 0.0, 1.0);
    hists[var]["FatJet1_Tau3OverTau2"] = new TH1D(("FatJet1_Tau3OverTau2_" + var).c_str(), 
                                                "FatJet1_Tau3OverTau2", 100, 0.0, 1.0);
    
    // FatJet 2 Histograms
    hists[var]["FatJet2_pt"] = new TH1D(("FatJet2_pt_" + var).c_str(), 
                                      "FatJet2_pt", 200, 0, 1200);
    hists[var]["FatJet2_eta"] = new TH1D(("FatJet2_eta_" + var).c_str(), 
                                       "FatJet2_eta", 100, -5, 5);
    hists[var]["FatJet2_MassSD"] = new TH1D(("FatJet2_MassSD_" + var).c_str(), 
                                          "FatJet2_MassSD", 500, 0, 500);
    hists[var]["FatJet2_GloParT_MassVis"] = new TH1D(("FatJet2_GloParT_MassVis_" + var).c_str(), 
                                                   "FatJet2_GloParT_MassVis", 
                                                   N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
    hists[var]["FatJet2_GloParT_MassRes"] = new TH1D(("FatJet2_GloParT_MassRes_" + var).c_str(), 
                                                   "FatJet2_GloParT_MassRes", 
                                                   N_GloParT_Mass, GLOPART_MASS_MIN, GLOPART_MASS_MAX);
    hists[var]["FatJet2_GloParT_XbbVsQCD"] = new TH1D(("FatJet2_GloParT_XbbVsQCD_" + var).c_str(), 
                                                    "FatJet2_GloParT_XbbVsQCD", N_TXbb, 0.0, 1.0);
    hists[var]["FatJet2_Tau3OverTau2"] = new TH1D(("FatJet2_Tau3OverTau2_" + var).c_str(), 
                                                "FatJet2_Tau3OverTau2", 100, 0.0, 1.0);
                                                
    // Enable proper error calculation for all histograms
    for (auto& hist_pair : hists[var]) {
      hist_pair.second->Sumw2();
    }
  }
  
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

  // Loop over all events
  for (int i = 0; i < ntuples->GetEntries(); i++) {
    ntuples->GetEntry(i);

    // Apply selection cuts
    if (fatJet1_pt <= 450 || fabs(fatJet1_eta) >= 2.5 || fatJet1_msoftdrop <= 50 || 
        fatJet1_GloParT_XbbVsQCD <= 0.1 || fatJet1_Tau3OverTau2 >= 0.46) continue;
    if (fatJet2_pt <= 450 || fabs(fatJet2_eta) >= 2.5 || fatJet2_msoftdrop <= 50 || 
        fatJet2_GloParT_XbbVsQCD <= 0.1 || fatJet2_Tau3OverTau2 >= 0.46) continue;

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
      Int_t bin_ptjj = _SF_ptjj->GetXaxis()->FindBin(fatJet1_pt);
      sf_ptjj = _SF_ptjj->GetBinContent(bin_ptjj);
      sf_ptjj_err = _SF_ptjj->GetBinError(bin_ptjj);
      if (sf_ptjj <= 0 || sf_ptjj_err < 0) {
        sf_ptjj = 1.0;
        sf_ptjj_err = 0.0;
      }
    }

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
    
    // Calculate weights
    Float_t weight_nominal = weight * sf_TXbb * sf_tau32 * sf_ptjj;
    
    // Fill nominal histograms
    hists["nominal"]["FatJet1_pt"]->Fill(fatJet1_pt, weight_nominal);
    hists["nominal"]["FatJet1_eta"]->Fill(fatJet1_eta, weight_nominal);
    hists["nominal"]["FatJet1_MassSD"]->Fill(fatJet1_msoftdrop, weight_nominal);
    hists["nominal"]["FatJet1_GloParT_MassVis"]->Fill(fatJet1_GloParT_massVis, weight_nominal);
    hists["nominal"]["FatJet1_GloParT_MassRes"]->Fill(fatJet1_GloParT_massRes, weight_nominal);
    hists["nominal"]["FatJet1_GloParT_XbbVsQCD"]->Fill(fatJet1_GloParT_XbbVsQCD, weight_nominal);
    hists["nominal"]["FatJet1_Tau3OverTau2"]->Fill(fatJet1_Tau3OverTau2, weight_nominal);
    
    hists["nominal"]["FatJet2_pt"]->Fill(fatJet2_pt, weight_nominal);
    hists["nominal"]["FatJet2_eta"]->Fill(fatJet2_eta, weight_nominal);
    hists["nominal"]["FatJet2_MassSD"]->Fill(fatJet2_msoftdrop, weight_nominal);
    hists["nominal"]["FatJet2_GloParT_MassVis"]->Fill(fatJet2_GloParT_massVis, weight_nominal);
    hists["nominal"]["FatJet2_GloParT_MassRes"]->Fill(fatJet2_GloParT_massRes, weight_nominal);
    hists["nominal"]["FatJet2_GloParT_XbbVsQCD"]->Fill(fatJet2_GloParT_XbbVsQCD, weight_nominal);
    hists["nominal"]["FatJet2_Tau3OverTau2"]->Fill(fatJet2_Tau3OverTau2, weight_nominal);
    
    // Fill up/down variations if any scale factors are applied
    if (apply_any_sf) {
      // weight_up = weight * (sf + sf_err) = weight_nominal + weight * sf_err
      // weight_dn = weight * (sf - sf_err) = weight_nominal - weight * sf_err
      Float_t weight_up = weight_nominal + weight * sf_err;
      Float_t weight_down = weight_nominal - weight * sf_err;
      
      // Prevent negative weights in down variation
      if (weight_down < 0) weight_down = 0;
      
      // Up variation
      hists["up"]["FatJet1_pt"]->Fill(fatJet1_pt, weight_up);
      hists["up"]["FatJet1_eta"]->Fill(fatJet1_eta, weight_up);
      hists["up"]["FatJet1_MassSD"]->Fill(fatJet1_msoftdrop, weight_up);
      hists["up"]["FatJet1_GloParT_MassVis"]->Fill(fatJet1_GloParT_massVis, weight_up);
      hists["up"]["FatJet1_GloParT_MassRes"]->Fill(fatJet1_GloParT_massRes, weight_up);
      hists["up"]["FatJet1_GloParT_XbbVsQCD"]->Fill(fatJet1_GloParT_XbbVsQCD, weight_up);
      hists["up"]["FatJet1_Tau3OverTau2"]->Fill(fatJet1_Tau3OverTau2, weight_up);
      
      hists["up"]["FatJet2_pt"]->Fill(fatJet2_pt, weight_up);
      hists["up"]["FatJet2_eta"]->Fill(fatJet2_eta, weight_up);
      hists["up"]["FatJet2_MassSD"]->Fill(fatJet2_msoftdrop, weight_up);
      hists["up"]["FatJet2_GloParT_MassVis"]->Fill(fatJet2_GloParT_massVis, weight_up);
      hists["up"]["FatJet2_GloParT_MassRes"]->Fill(fatJet2_GloParT_massRes, weight_up);
      hists["up"]["FatJet2_GloParT_XbbVsQCD"]->Fill(fatJet2_GloParT_XbbVsQCD, weight_up);
      hists["up"]["FatJet2_Tau3OverTau2"]->Fill(fatJet2_Tau3OverTau2, weight_up);
      
      // Down variation
      hists["down"]["FatJet1_pt"]->Fill(fatJet1_pt, weight_down);
      hists["down"]["FatJet1_eta"]->Fill(fatJet1_eta, weight_down);
      hists["down"]["FatJet1_MassSD"]->Fill(fatJet1_msoftdrop, weight_down);
      hists["down"]["FatJet1_GloParT_MassVis"]->Fill(fatJet1_GloParT_massVis, weight_down);
      hists["down"]["FatJet1_GloParT_MassRes"]->Fill(fatJet1_GloParT_massRes, weight_down);
      hists["down"]["FatJet1_GloParT_XbbVsQCD"]->Fill(fatJet1_GloParT_XbbVsQCD, weight_down);
      hists["down"]["FatJet1_Tau3OverTau2"]->Fill(fatJet1_Tau3OverTau2, weight_down);
      
      hists["down"]["FatJet2_pt"]->Fill(fatJet2_pt, weight_down);
      hists["down"]["FatJet2_eta"]->Fill(fatJet2_eta, weight_down);
      hists["down"]["FatJet2_MassSD"]->Fill(fatJet2_msoftdrop, weight_down);
      hists["down"]["FatJet2_GloParT_MassVis"]->Fill(fatJet2_GloParT_massVis, weight_down);
      hists["down"]["FatJet2_GloParT_MassRes"]->Fill(fatJet2_GloParT_massRes, weight_down);
      hists["down"]["FatJet2_GloParT_XbbVsQCD"]->Fill(fatJet2_GloParT_XbbVsQCD, weight_down);
      hists["down"]["FatJet2_Tau3OverTau2"]->Fill(fatJet2_Tau3OverTau2, weight_down);
    }
  }

  // Close input file
  f0->Close();

  // Force write each histogram to deal with potential bugs
  if (apply_any_sf) {
    f->cd();  // Make sure we're writing to the output file
    
    // Write nominal histograms explicitly
    hists["nominal"]["FatJet1_pt"]->Write();
    hists["nominal"]["FatJet1_eta"]->Write();
    hists["nominal"]["FatJet1_MassSD"]->Write();
    hists["nominal"]["FatJet1_GloParT_MassVis"]->Write();
    hists["nominal"]["FatJet1_GloParT_MassRes"]->Write();
    hists["nominal"]["FatJet1_GloParT_XbbVsQCD"]->Write();
    hists["nominal"]["FatJet1_Tau3OverTau2"]->Write();
    
    hists["nominal"]["FatJet2_pt"]->Write();
    hists["nominal"]["FatJet2_eta"]->Write();
    hists["nominal"]["FatJet2_MassSD"]->Write();
    hists["nominal"]["FatJet2_GloParT_MassVis"]->Write();
    hists["nominal"]["FatJet2_GloParT_MassRes"]->Write();
    hists["nominal"]["FatJet2_GloParT_XbbVsQCD"]->Write();
    hists["nominal"]["FatJet2_Tau3OverTau2"]->Write();
    
    // Write up variation histograms explicitly
    hists["up"]["FatJet1_pt"]->Write();
    hists["up"]["FatJet1_eta"]->Write();
    hists["up"]["FatJet1_MassSD"]->Write();
    hists["up"]["FatJet1_GloParT_MassVis"]->Write();
    hists["up"]["FatJet1_GloParT_MassRes"]->Write();
    hists["up"]["FatJet1_GloParT_XbbVsQCD"]->Write();
    hists["up"]["FatJet1_Tau3OverTau2"]->Write();
    
    hists["up"]["FatJet2_pt"]->Write();
    hists["up"]["FatJet2_eta"]->Write();
    hists["up"]["FatJet2_MassSD"]->Write();
    hists["up"]["FatJet2_GloParT_MassVis"]->Write();
    hists["up"]["FatJet2_GloParT_MassRes"]->Write();
    hists["up"]["FatJet2_GloParT_XbbVsQCD"]->Write();
    hists["up"]["FatJet2_Tau3OverTau2"]->Write();
    
    // Write down variation histograms explicitly
    hists["down"]["FatJet1_pt"]->Write();
    hists["down"]["FatJet1_eta"]->Write();
    hists["down"]["FatJet1_MassSD"]->Write();
    hists["down"]["FatJet1_GloParT_MassVis"]->Write();
    hists["down"]["FatJet1_GloParT_MassRes"]->Write();
    hists["down"]["FatJet1_GloParT_XbbVsQCD"]->Write();
    hists["down"]["FatJet1_Tau3OverTau2"]->Write();
    
    hists["down"]["FatJet2_pt"]->Write();
    hists["down"]["FatJet2_eta"]->Write();
    hists["down"]["FatJet2_MassSD"]->Write();
    hists["down"]["FatJet2_GloParT_MassVis"]->Write();
    hists["down"]["FatJet2_GloParT_MassRes"]->Write();
    hists["down"]["FatJet2_GloParT_XbbVsQCD"]->Write();
    hists["down"]["FatJet2_Tau3OverTau2"]->Write();
  }
  
  // Call Write() on the file to ensure everything is written
  f->Write();
  std::cout << "Wrote histograms to output file." << std::endl;
  
  // Clean up
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