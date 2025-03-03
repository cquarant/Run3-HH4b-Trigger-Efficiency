#include "TGaxis.h"
#include "TRandom.h"
#include "TStyle.h"

#include "TCanvas.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>
#include "kfact.h"
#include "TXbb.h"


void SF_TXbb(const std::string &year,      // 2022, 2023
             const std::string &path_data,    // path to the data root file
             const std::string &path_QCD,     // path to the QCD root file
             const std::string &path_VV,      // path to the VV root file
             const std::string &path_VJ,      // path to the VJ root file
             const std::string &path_TTbar,   // path to the TTbar root file
             const std::string &path_ttHto2B, // path to the ttHto2B root file
             const std::string &output_path   // path to the output root file
) {
  // Setup output file
  TFile *f = new TFile(output_path.c_str(), "RECREATE");
  
  // Define binning
  Float_t Bound[14] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
                       0.94, 0.97, 0.99, 1.00};
  Float_t BoundZoomed[6] = {0.8, 0.9, 0.94, 0.97, 0.99, 1.00};
  
  // Create histograms
  TH1D *_SF = new TH1D("SF", "SF", 13, Bound);
  TH1D *hint_68 = new TH1D("hint_68", "", 370, 0.0, 1.0);
  TH1D *_SF_zoomed = new TH1D("SF_zoomed", "SF_zoomed", 5, BoundZoomed);
  TH1D *hint_68_zoomed = new TH1D("hint_68_zoomed", "", 100, 0.8, 1.0);
  
  // Style settings
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  
  // Load input histograms
  TString JET = "FatJet1";
  TString VAR = "_GloParT_XbbVsQCD";
  TString variable = JET + VAR;
  
  // Open input files
  TFile *f_Data = new TFile(path_data.c_str());
  TFile *f_QCD = new TFile(path_QCD.c_str());
  TFile *f_VV = new TFile(path_VV.c_str());
  TFile *f_VJ = new TFile(path_VJ.c_str());
  TFile *f_TTbar = new TFile(path_TTbar.c_str());
  TFile *f_ttHto2B = new TFile(path_ttHto2B.c_str());
  
  // Get histograms
  TH1D *_Data_var = (TH1D *)f_Data->Get(variable);
  TH1D *_QCD_var = (TH1D *)f_QCD->Get(variable);
  TH1D *_TTbar_var = (TH1D *)f_TTbar->Get(variable);
  TH1D *_VV_var = (TH1D *)f_VV->Get(variable);
  TH1D *_VJ_var = (TH1D *)f_VJ->Get(variable);
  TH1D *_ttHto2B_var = (TH1D *)f_ttHto2B->Get(variable);
  
  // Apply k-factor
  Float_t kfact = (year == "2023") ? KFACT_2023 : KFACT_2022;
  _QCD_var->Scale(kfact);
  _TTbar_var->Scale(kfact);
  _VJ_var->Scale(kfact);
  _VV_var->Scale(kfact);
  _ttHto2B_var->Scale(kfact);
  
  // Subtract backgrounds from data
  for (int iB = 1; iB <= _Data_var->GetSize(); ++iB) {
    _Data_var->SetBinContent(
        iB, _Data_var->GetBinContent(iB) - _VJ_var->GetBinContent(iB) -
                _VV_var->GetBinContent(iB) - _QCD_var->GetBinContent(iB)
                - _ttHto2B_var->GetBinContent(iB));
  }
  
  // Create rebinned histograms
  TH1D *_My_MC = new TH1D("My_MC", "My_MC", 13, Bound);
  TH1D *_My_Data = new TH1D("My_Data", "My_Data", 13, Bound);
  TH1D *_My_MC_zoomed = new TH1D("My_MC_zoomed", "My_MC_zoomed", 5, BoundZoomed);
  TH1D *_My_Data_zoomed = new TH1D("My_Data_zoomed", "My_Data_zoomed", 5, BoundZoomed);
  
  // Initialize with zeros
  for (int i = 1; i <= _My_MC->GetSize(); i++) _My_MC->SetBinContent(i, 0.0);
  for (int i = 1; i <= _My_Data->GetSize(); i++) _My_Data->SetBinContent(i, 0.0);
  for (int i = 1; i <= _My_MC_zoomed->GetSize(); i++) _My_MC_zoomed->SetBinContent(i, 0.0);
  for (int i = 1; i <= _My_Data_zoomed->GetSize(); i++) _My_Data_zoomed->SetBinContent(i, 0.0);
  
  // Fill rebinned histograms
  for (int i = 1; i <= _TTbar_var->GetSize() - 2; i++) {
    Float_t bin_value = i / 100.0 - 0.005;
    
    // Full range
    Int_t bin = _My_MC->GetXaxis()->FindBin(bin_value);
    _My_MC->SetBinContent(bin, _My_MC->GetBinContent(bin) + _TTbar_var->GetBinContent(i));
    _My_Data->SetBinContent(bin, _My_Data->GetBinContent(bin) + _Data_var->GetBinContent(i));
    
    // Zoomed range (TXbb > 0.8)
    if (bin_value >= 0.8) {
      Int_t bin_zoomed = _My_MC_zoomed->GetXaxis()->FindBin(bin_value);
      _My_MC_zoomed->SetBinContent(bin_zoomed, _My_MC_zoomed->GetBinContent(bin_zoomed) + _TTbar_var->GetBinContent(i));
      _My_Data_zoomed->SetBinContent(bin_zoomed, _My_Data_zoomed->GetBinContent(bin_zoomed) + _Data_var->GetBinContent(i));
    }
  }
  
  // Apply rebinning to original histograms
  int reb = 5 * (N_TXbb / 100);
  _Data_var->Rebin(reb);
  _VJ_var->Rebin(reb);
  _VV_var->Rebin(reb);
  _QCD_var->Rebin(reb);
  _TTbar_var->Rebin(reb);
  _ttHto2B_var->Rebin(reb);
  
  // ---- Create full range plot ----
  TCanvas *c1 = new TCanvas("SF_TXbb", "SF_TXbb", 100, 52, 800, 686);
  c1->SetGridy();
  
  // Setup upper pad
  TPad *c1_1 = new TPad("c1_1", "c1_1", 0.1, 0.35, 0.9, 0.95);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetLogy();
  c1_1->SetBottomMargin(0.15);
  c1_1->SetLeftMargin(0.15);
  c1_1->SetGridy();
  
  // Draw MC and Data histograms
  _My_MC->SetTitle("");
  _My_MC->GetXaxis()->SetTitle("FJet - gParT_XbbVsQCD");
  _My_MC->GetYaxis()->SetTitle("Events");
  _My_MC->SetFillColor(kOrange - 2);
  _My_MC->SetLineWidth(1);
  
  double ymin = 10e-2;
  double ymax = 15 * (_My_MC->GetMaximum());
  _My_MC->SetMaximum(ymax);
  _My_MC->SetMinimum(ymin);
  _My_MC->Draw("HIST");
  
  _My_Data->SetMarkerColor(1);
  _My_Data->SetMarkerStyle(20);
  _My_Data->SetMarkerSize(1.1);
  _My_Data->SetLineColor(1);
  _My_Data->SetLineWidth(3);
  _My_Data->Draw("E same");
  
  // Add legend
  TLegend *leg = new TLegend(0.58, 0.8, 0.78, 0.88);
  leg->SetTextSize(0.03);
  leg->AddEntry(_My_Data, "Data (bkg. subtracted)", "p");
  leg->AddEntry(_My_MC, "TTbar", "f");
  leg->Draw("same");
  
  c1->cd();
  
  // Setup lower pad
  TPad *c1_2 = new TPad("lower", "pad", 0.1, 0.05, 0.9, 0.35);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetBottomMargin(0.2);
  c1_2->SetLeftMargin(0.15);
  c1_2->SetGridy();
  
  // Create and draw SF ratio
  TH1D *_SF_pt = (TH1D *)_My_Data->Clone();
  TH1D *_MC_pt = (TH1D *)_My_MC->Clone();
  _MC_pt->Sumw2();
  _SF_pt->Sumw2();
  _SF_pt->Divide(_MC_pt);
  
  _SF_pt->SetTitle("");
  _SF_pt->SetMarkerColor(1);
  _SF_pt->SetMarkerStyle(20);
  _SF_pt->SetMarkerSize(1.2);
  _SF_pt->SetLineColor(1);
  _SF_pt->GetYaxis()->SetTitle("Scale factor");
  _SF_pt->SetMaximum(1.5);
  _SF_pt->SetMinimum(0.5);
  _SF_pt->Draw("e1");
  
  // Fit and draw confidence interval
  TF1 *Mypol = new TF1("Mypol", "pol3", 0.01, 1.0);
  _SF_pt->Fit("Mypol", "R");
  
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68, 0.68);
  hint_68->SetFillColor(kYellow);
  hint_68->Draw("e3same");
  _SF_pt->Draw("e1same");
  
  // Add legend
  TLegend *leg1 = new TLegend(0.35, 0.6, 0.55, 0.85);
  leg1->SetTextSize(0.07);
  leg1->AddEntry(Mypol, "pol3", "l");
  leg1->AddEntry(hint_68, "1 #sigma", "f");
  leg1->Draw("same");
  
  // Fill the output SF histogram
  for (int i = 1; i <= _SF_pt->GetNbinsX(); i++) {
    _SF->SetBinContent(i, _SF_pt->GetBinContent(i));
    _SF->SetBinError(i, _SF_pt->GetBinError(i));
  }
  
  // ---- Create zoomed plot (TXbb > 0.8) ----
  TCanvas *c2 = new TCanvas("SF_TXbb_zoomed", "TXbb > 0.8", 100, 52, 800, 686);
  c2->SetGridy();
  
  // Setup upper pad
  TPad *c2_1 = new TPad("c2_1", "c2_1", 0.1, 0.35, 0.9, 0.95);
  c2_1->Draw();
  c2_1->cd();
  c2_1->SetLogy();
  c2_1->SetBottomMargin(0.15);
  c2_1->SetLeftMargin(0.15);
  c2_1->SetGridy();
  
  // Draw MC and Data histograms
  _My_MC_zoomed->SetTitle("");
  _My_MC_zoomed->GetXaxis()->SetTitle("FJet - gParT_XbbVsQCD (Zoomed >0.8)");
  _My_MC_zoomed->GetYaxis()->SetTitle("Events");
  _My_MC_zoomed->SetFillColor(kOrange - 2);
  _My_MC_zoomed->SetLineWidth(1);
  
  double ymin_zoomed = 10e-2;
  double ymax_zoomed = 15 * (_My_MC_zoomed->GetMaximum());
  _My_MC_zoomed->SetMaximum(ymax_zoomed);
  _My_MC_zoomed->SetMinimum(ymin_zoomed);
  _My_MC_zoomed->Draw("HIST");
  
  _My_Data_zoomed->SetMarkerColor(1);
  _My_Data_zoomed->SetMarkerStyle(20);
  _My_Data_zoomed->SetMarkerSize(1.1);
  _My_Data_zoomed->SetLineColor(1);
  _My_Data_zoomed->SetLineWidth(3);
  _My_Data_zoomed->Draw("E same");
  
  // Add legend
  TLegend *leg_zoomed = new TLegend(0.58, 0.8, 0.78, 0.88);
  leg_zoomed->SetTextSize(0.03);
  leg_zoomed->AddEntry(_My_Data_zoomed, "Data (bkg. subtracted)", "p");
  leg_zoomed->AddEntry(_My_MC_zoomed, "TTbar", "f");
  leg_zoomed->Draw("same");
  
  c2->cd();
  
  // Setup lower pad
  TPad *c2_2 = new TPad("lower_zoomed", "pad_zoomed", 0.1, 0.05, 0.9, 0.35);
  c2_2->Draw();
  c2_2->cd();
  c2_2->SetBottomMargin(0.2);
  c2_2->SetLeftMargin(0.15);
  c2_2->SetGridy();
  
  // Create and draw SF ratio
  TH1D *_SF_pt_zoomed = (TH1D *)_My_Data_zoomed->Clone();
  TH1D *_MC_pt_zoomed = (TH1D *)_My_MC_zoomed->Clone();
  _MC_pt_zoomed->Sumw2();
  _SF_pt_zoomed->Sumw2();
  _SF_pt_zoomed->Divide(_MC_pt_zoomed);
  
  _SF_pt_zoomed->SetTitle("");
  _SF_pt_zoomed->SetMarkerColor(1);
  _SF_pt_zoomed->SetMarkerStyle(20);
  _SF_pt_zoomed->SetMarkerSize(1.2);
  _SF_pt_zoomed->SetLineColor(1);
  _SF_pt_zoomed->GetYaxis()->SetTitle("Scale factor");
  _SF_pt_zoomed->SetMaximum(1.5);
  _SF_pt_zoomed->SetMinimum(0.5);
  _SF_pt_zoomed->Draw("e1");
  
  // Fit and draw confidence interval
  TF1 *Mypol_zoomed = new TF1("Mypol_zoomed", "pol3", 0.8, 1.0);
  _SF_pt_zoomed->Fit("Mypol_zoomed", "R");
  
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68_zoomed, 0.68);
  hint_68_zoomed->SetFillColor(kYellow);
  hint_68_zoomed->Draw("e3same");
  _SF_pt_zoomed->Draw("e1same");
  
  // Add legend
  TLegend *leg1_zoomed = new TLegend(0.35, 0.6, 0.55, 0.85);
  leg1_zoomed->SetTextSize(0.07);
  leg1_zoomed->AddEntry(Mypol_zoomed, "pol3", "l");
  leg1_zoomed->AddEntry(hint_68_zoomed, "1 #sigma", "f");
  leg1_zoomed->Draw("same");
  
  // Fill the output SF histogram
  for (int i = 1; i <= _SF_pt_zoomed->GetNbinsX(); i++) {
    _SF_zoomed->SetBinContent(i, _SF_pt_zoomed->GetBinContent(i));
    _SF_zoomed->SetBinError(i, _SF_pt_zoomed->GetBinError(i));
  }
  
  // Save everything to output file
  f->cd();
  c1->Write("SF_TXbb_canvas");
  c2->Write("SF_TXbb_zoomed_canvas");
  _SF->Write();
  _SF_zoomed->Write();
  hint_68->Write();
  hint_68_zoomed->Write();
  f->Close();
  
  // Close input files
  f_Data->Close();
  f_QCD->Close();
  f_VV->Close();
  f_VJ->Close();
  f_TTbar->Close();
  f_ttHto2B->Close();
}