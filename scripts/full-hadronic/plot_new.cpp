#include "TCanvas.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TStyle.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "kfact.h"
#include "TXbb.h"

#define PLOT_CHI_2 false

void plot(const std::string &year,         // 2022, 2023
          const std::string &path_data,    // path to the data root file
          const std::string &path_QCD,     // path to the QCD root file
          const std::string &path_VV,      // path to the VV root file
          const std::string &path_VJ,      // path to the VJ root file
          const std::string &path_TTbar,   // path to the TTbar root file
          const std::string &path_ttHto2B, // path to the ttHto2B root file
          const std::string &output_path,  // path to the output root file
          const std::string &var,          // variable to plot
          const std::string &var_label,    // jet to plot
          const std::string &TXbb_bin=""   // custom binning (comma separated)
) {

  TString variable = var;
  TString XTitle = var_label;

  gStyle->SetOptFile(0);
  gStyle->SetOptStat("m");

  TCanvas *c1 = new TCanvas("c1", "c1", 100, 50, 850, 750);
  gStyle->SetOptStat(0);
  c1->Range(0, 0, 1, 1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.08);
  c1->SetBottomMargin(0.13);
  c1->SetFrameFillStyle(0);
  c1->SetFrameLineStyle(0);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderSize(10);

  TFile *f_Data = new TFile(path_data.c_str());
  TFile *f_QCD = new TFile(path_QCD.c_str());
  TFile *f_VV = new TFile(path_VV.c_str());
  TFile *f_VJ = new TFile(path_VJ.c_str());
  TFile *f_TTbar = new TFile(path_TTbar.c_str());
  TFile *f_ttHto2B = new TFile(path_ttHto2B.c_str());

  // Get nominal histograms - first try with "_nominal" suffix for TTbar
  TH1D *_Data_var = (TH1D *)f_Data->Get(variable+ "_nominal");
  TH1D *_QCD_var = (TH1D *)f_QCD->Get(variable+ "_nominal");
  TH1D *_TTbar_var = (TH1D *)f_TTbar->Get(variable + "_nominal");
  TH1D *_VV_var = (TH1D *)f_VV->Get(variable+ "_nominal");
  TH1D *_VJ_var = (TH1D *)f_VJ->Get(variable+ "_nominal");
  TH1D *_ttHto2B_var = (TH1D *)f_ttHto2B->Get(variable+ "_nominal");

  // If nominal TTbar doesn't exist with suffix, try without
  if (!_TTbar_var) _TTbar_var = (TH1D *)f_TTbar->Get(variable);

  // Check for TTbar up/down variations
  TH1D *_TTbar_var_up = (TH1D *)f_TTbar->Get(variable + "_up");
  TH1D *_TTbar_var_down = (TH1D *)f_TTbar->Get(variable + "_down");

  // Check if we have systematic variations for TTbar
  bool has_TTbar_syst = _TTbar_var_up && _TTbar_var_down;

  // Enable proper error calculation for all histograms
  _Data_var->Sumw2();
  _QCD_var->Sumw2();
  _TTbar_var->Sumw2();
  _VV_var->Sumw2();
  _VJ_var->Sumw2();
  _ttHto2B_var->Sumw2();

  // Apply k-factor scaling
  Float_t kfact;
  if (year == "2023") {
      kfact = KFACT_2023;
  } else if (year == "2022") {
      kfact = KFACT_2022;
  } else {
    throw std::invalid_argument("Invalid year: " + year);
  }
  
  _QCD_var->Scale(kfact);
  _TTbar_var->Scale(kfact);
  _VJ_var->Scale(kfact);
  _VV_var->Scale(kfact);
  _ttHto2B_var->Scale(kfact);

  // Scale TTbar up/down variations if they exist
  if (has_TTbar_syst) {
    _TTbar_var_up->Scale(kfact);
    _TTbar_var_down->Scale(kfact);
  }

  // Handle custom binning for XbbVsQCD
  if (variable.Contains("Xbb") && !TXbb_bin.empty()) {
    // Parse custom binning string
    std::vector<Double_t> bins;
    std::stringstream ss(TXbb_bin);
    std::string value;
    
    while (std::getline(ss, value, ',')) {
        bins.push_back(std::stod(value));
    }
    
    // Create array from vector
    const int nBins = bins.size() - 1;
    Double_t* bins_array = bins.data();
    
    // Rebin histograms
    TH1D* rebinned_Data = (TH1D*)_Data_var->Rebin(nBins, "rebinned_Data", bins_array);
    TH1D* rebinned_QCD = (TH1D*)_QCD_var->Rebin(nBins, "rebinned_QCD", bins_array);
    TH1D* rebinned_VV = (TH1D*)_VV_var->Rebin(nBins, "rebinned_VV", bins_array);
    TH1D* rebinned_VJ = (TH1D*)_VJ_var->Rebin(nBins, "rebinned_VJ", bins_array);
    TH1D* rebinned_TTbar = (TH1D*)_TTbar_var->Rebin(nBins, "rebinned_TTbar", bins_array);
    TH1D* rebinned_ttHto2B = (TH1D*)_ttHto2B_var->Rebin(nBins, "rebinned_ttHto2B", bins_array);
    
    // Also rebin TTbar up/down variations if they exist
    TH1D* rebinned_TTbar_up = nullptr;
    TH1D* rebinned_TTbar_down = nullptr;
    
    if (has_TTbar_syst) {
      rebinned_TTbar_up = (TH1D*)_TTbar_var_up->Rebin(nBins, "rebinned_TTbar_up", bins_array);
      rebinned_TTbar_down = (TH1D*)_TTbar_var_down->Rebin(nBins, "rebinned_TTbar_down", bins_array);
      
      _TTbar_var_up = rebinned_TTbar_up;
      _TTbar_var_down = rebinned_TTbar_down;
    }
    
    // Replace original histograms with rebinned ones
    _Data_var = rebinned_Data;
    _QCD_var = rebinned_QCD;
    _VV_var = rebinned_VV;
    _VJ_var = rebinned_VJ;
    _TTbar_var = rebinned_TTbar;
    _ttHto2B_var = rebinned_ttHto2B;
  } else {
    int reb = 5 * (N_TXbb / 100);

    if (variable.Contains("GloParT") && variable.Contains("Mass")) {
      reb = 2 * (N_TXbb / 100);
    }
    if (variable.Contains("pTjj")) {
      reb = 2 * (N_TXbb / 100);
    }

    _Data_var->Rebin(reb);
    _VJ_var->Rebin(reb);
    _VV_var->Rebin(reb);
    _QCD_var->Rebin(reb);
    _TTbar_var->Rebin(reb);
    _ttHto2B_var->Rebin(reb);
    
    // Rebin TTbar up/down variations if they exist
    if (has_TTbar_syst) {
      _TTbar_var_up->Rebin(reb);
      _TTbar_var_down->Rebin(reb);
    }
  }
  
  // Create stacked histograms for nominal values
  TH1D *_QCD_stack = (TH1D *)_QCD_var->Clone("QCD_stack");
  TH1D *_ttHto2B_stack = (TH1D *)_ttHto2B_var->Clone("ttHto2B_stack");
  TH1D *_VV_stack = (TH1D *)_VV_var->Clone("VV_stack");
  TH1D *_VJ_stack = (TH1D *)_VJ_var->Clone("VJ_stack");
  TH1D *_TTbar_stack = (TH1D *)_TTbar_var->Clone("TTbar_stack");
  
  _ttHto2B_stack->Add(_QCD_stack);  // ttHto2B on top of QCD
  _VV_stack->Add(_ttHto2B_stack);   // VV on top of ttHto2B
  _VJ_stack->Add(_VV_stack);        // VJ on top of VV
  _TTbar_stack->Add(_VJ_stack);     // TTbar on top of VJ

  // Create up/down stacks for systematic variations if TTbar up/down exist
  TH1D *_TTbar_stack_up = nullptr;
  TH1D *_TTbar_stack_down = nullptr;
  
  if (has_TTbar_syst) {
    // Calculate the difference between nominal TTbar and up/down variations
    TH1D *_TTbar_diff_up = (TH1D *)_TTbar_var_up->Clone("TTbar_diff_up");
    TH1D *_TTbar_diff_down = (TH1D *)_TTbar_var_down->Clone("TTbar_diff_down");
    
    _TTbar_diff_up->Add(_TTbar_var, -1.0);    // up - nominal
    _TTbar_diff_down->Add(_TTbar_var, -1.0);  // down - nominal
    
    // Create the total stack with variations
    _TTbar_stack_up = (TH1D *)_TTbar_stack->Clone("TTbar_stack_up");
    _TTbar_stack_down = (TH1D *)_TTbar_stack->Clone("TTbar_stack_down");
    
    _TTbar_stack_up->Add(_TTbar_diff_up);    // Total + (up - nominal)
    _TTbar_stack_down->Add(_TTbar_diff_down); // Total + (down - nominal)
  }

  // Set up the plotting canvas
  TPad *c1_1 = new TPad("c1_1", "c1_1", 0.01, 0.04, 0.75, 0.9);
  c1_1->Draw();
  c1_1->cd();
  c1_1->Range(-93.75, -1.992728, 531.25, 5.643644);
  c1_1->SetFillColor(0);
  c1_1->SetBorderMode(0);
  c1_1->SetBorderSize(10);
  c1_1->SetLogy();
  c1_1->SetTickx(1);
  c1_1->SetTicky(1);
  c1_1->SetLeftMargin(0.15);
  c1_1->SetRightMargin(0.05);
  c1_1->SetTopMargin(0.01);
  c1_1->SetBottomMargin(0.3);
  c1_1->SetFrameFillStyle(0);
  c1_1->SetFrameLineStyle(0);
  c1_1->SetFrameLineWidth(2);
  c1_1->SetFrameBorderMode(0);
  c1_1->SetFrameBorderSize(10);
  c1_1->SetFrameFillStyle(0);
  c1_1->SetFrameLineStyle(0);
  c1_1->SetFrameLineWidth(2);
  c1_1->SetFrameBorderMode(0);
  c1_1->SetFrameBorderSize(10);

  TString YTitle_10 = "Events";

  _TTbar_stack->GetXaxis()->SetLabelFont(42);
  _TTbar_stack->GetXaxis()->SetLabelOffset(0.15);
  _TTbar_stack->GetXaxis()->SetTitle(XTitle);
  _TTbar_stack->GetXaxis()->SetTitleSize(0.05);
  _TTbar_stack->GetXaxis()->SetTitleOffset(3.3);
  _TTbar_stack->GetXaxis()->SetTitleFont(42);
  _TTbar_stack->GetYaxis()->SetNdivisions(505);
  _TTbar_stack->GetYaxis()->SetLabelFont(42);
  _TTbar_stack->GetYaxis()->SetLabelSize(0.04);
  _TTbar_stack->GetYaxis()->SetLabelOffset(0.01);
  _TTbar_stack->GetYaxis()->SetTitle(YTitle_10);
  _TTbar_stack->GetYaxis()->SetTitleOffset(1.5);
  _TTbar_stack->GetYaxis()->SetTitleSize(0.05);
  _TTbar_stack->GetYaxis()->SetTitleFont(42);
  _TTbar_stack->SetTitle("");

  double ymin = 10e-2;
  double ymax = 5 * (_TTbar_stack->GetMaximum());
  
  // Drawing order needs to match the stacking order (drawing in reverse):
  _TTbar_stack->Draw("HIST");
  _TTbar_stack->SetFillColor(kOrange - 2);
  _TTbar_stack->SetLineWidth(1);
  _TTbar_stack->SetLineStyle(1);
  _TTbar_stack->SetMaximum(ymax);
  _TTbar_stack->SetMinimum(ymin);
  
  _VJ_stack->SetFillColor(kGreen - 3);
  _VJ_stack->SetLineWidth(1);
  _VJ_stack->SetLineStyle(1);
  _VJ_stack->SetMaximum(ymax);
  _VJ_stack->SetMinimum(ymin);
  _VJ_stack->Draw("HIST same");
  
  _VV_stack->SetFillColor(kGreen);
  _VV_stack->SetLineWidth(1);
  _VV_stack->SetLineStyle(1);
  _VV_stack->SetMaximum(ymax);
  _VV_stack->SetMinimum(ymin);
  _VV_stack->Draw("HIST same");
  
  _ttHto2B_stack->SetFillColor(kBlue);
  _ttHto2B_stack->SetLineWidth(1);
  _ttHto2B_stack->SetLineStyle(1);
  _ttHto2B_stack->SetMaximum(ymax);
  _ttHto2B_stack->SetMinimum(ymin);
  _ttHto2B_stack->Draw("HIST same");
  
  _QCD_stack->Draw("HIST same");
  _QCD_stack->SetFillColor(43);
  _QCD_stack->SetLineWidth(1);
  _QCD_stack->SetLineStyle(1);
  _QCD_stack->SetMaximum(ymax);
  _QCD_stack->SetMinimum(ymin);

  // Create the error band histogram
  TH1D *bkgdErr = (TH1D *)_TTbar_stack->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3004); // More visible hatching pattern
  bkgdErr->SetFillColor(kBlack); // Black is more visible
  bkgdErr->SetLineColor(kBlack);
  bkgdErr->SetLineWidth(1);
  bkgdErr->SetMarkerStyle(0);
  bkgdErr->SetMarkerSize(0);

  // Calculate total errors for the MC histogram including systematics
  for (int iB = 1; iB <= _TTbar_stack->GetNbinsX(); ++iB) {
    float eStat = _TTbar_stack->GetBinError(iB);
    float eSyst = 0.0;
    
    // Add systematic uncertainties if available
    if (has_TTbar_syst) {
      float bin_nominal = _TTbar_stack->GetBinContent(iB);
      float bin_up = _TTbar_stack_up->GetBinContent(iB);
      float bin_down = _TTbar_stack_down->GetBinContent(iB);
      
      // Take the largest deviation as the systematic uncertainty
      float diff_up = fabs(bin_up - bin_nominal);
      float diff_down = fabs(bin_down - bin_nominal);
      eSyst = (diff_up > diff_down) ? diff_up : diff_down;
    }
    
    // Combine statistical and systematic uncertainties in quadrature
    float Err = TMath::Sqrt(eStat * eStat + eSyst * eSyst);
    bkgdErr->SetBinError(iB, Err);
    bkgdErr->SetBinContent(iB, _TTbar_stack->GetBinContent(iB));
  }

  // Draw the error band on the main plot
  bkgdErr->Draw("E2 same"); // E2 draws error bands

  // Setup data points with proper error bars
  _Data_var->SetMarkerColor(kBlack);
  _Data_var->SetMarkerStyle(20); // Filled circle
  _Data_var->SetMarkerSize(1.1);
  _Data_var->SetLineColor(kBlack);
  _Data_var->SetLineWidth(2);
  
  // Draw data with error bars
  _Data_var->Draw("E same"); // E draws error bars

  c1_1->Modified();
  c1->cd();

  // Calculate chi-square
  float chi2 = 0;
  for (int iB = 1; iB <= _Data_var->GetSize(); ++iB) {
    float xData = _Data_var->GetBinContent(iB);
    float xMC = _TTbar_stack->GetBinContent(iB);
    if (xMC > 1e-8) {
      float diff2 = (xData - xMC) * (xData - xMC);
      chi2 += diff2 / xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  // Create and set up ratio plot
  TH1D *ratioH = (TH1D *)_Data_var->Clone("ratioH");
  TH1D *ratioErrH = (TH1D *)bkgdErr->Clone("ratioErrH");
  
  // Setup ratio histogram style
  ratioH->SetMarkerColor(kBlack);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(kBlack);
  ratioH->SetLineWidth(2);
  
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.02);
  ratioH->GetXaxis()->SetLabelSize(0.15);
  ratioH->GetXaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetTickLength(0.07);

  ratioH->GetYaxis()->SetRangeUser(0.5, 1.5); // Adjust range to better show errors
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetYaxis()->SetTitle("Data / MC");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.01);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleOffset(0.38);
  ratioH->GetYaxis()->SetTickLength(0.02);
  ratioH->GetYaxis()->SetTitleSize(0.19);

  // Calculate ratio and errors properly
  for (int iB = 1; iB <= ratioH->GetNbinsX(); ++iB) {
    float dataVal = _Data_var->GetBinContent(iB);
    float dataErr = _Data_var->GetBinError(iB);
    float mcVal = _TTbar_stack->GetBinContent(iB);
    float mcErr = bkgdErr->GetBinError(iB);
    
    // Set the error band to be centered at 1.0
    ratioErrH->SetBinContent(iB, 1.0);
    
    // Calculate relative error for the MC prediction
    if (mcVal > 1e-6) {
      ratioErrH->SetBinError(iB, mcErr / mcVal);
    } else {
      ratioErrH->SetBinError(iB, 0);
    }
    
    // Calculate the data/MC ratio with proper error propagation
    if (mcVal > 1e-6) {
      float ratio = dataVal / mcVal;
      
      // Propagate errors - use data error for simplicity
      float ratioErr = 0;
      if (dataVal > 0) {
        ratioErr = dataErr / mcVal;
      }
      
      ratioH->SetBinContent(iB, ratio);
      ratioH->SetBinError(iB, ratioErr);
    } else {
      ratioH->SetBinContent(iB, 0);
      ratioH->SetBinError(iB, 0);
    }
  }

  // print the ratio
  std::cout << "Data / MC: ";
  for (int i = 1; i <= ratioH->GetNbinsX(); i++) {
    std::cout << ratioH->GetBinContent(i) << ", ";
  }
  std::cout << std::endl;
  std::cout << "Data / MC error: ";
  for (int i = 1; i <= ratioH->GetNbinsX(); i++) {
    std::cout << ratioH->GetBinError(i) << ", ";
  }
  std::cout << std::endl;

  // Set up ratio pad
  TPad *c1_2 = new TPad("lower", "pad", 0.01, 0.11, 0.75, 0.28);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetFillColor(0);
  c1_2->SetBorderMode(0);
  c1_2->SetBorderSize(10);
  c1_2->SetGridy();
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);
  c1_2->SetLeftMargin(0.15);
  c1_2->SetRightMargin(0.05);
  c1_2->SetTopMargin(0.1);
  c1_2->SetBottomMargin(0.15);
  c1_2->SetFrameFillStyle(0);
  c1_2->SetFrameLineStyle(0);
  c1_2->SetFrameLineWidth(2);
  c1_2->SetFrameBorderMode(0);
  c1_2->SetFrameBorderSize(10);
  c1_2->SetFrameFillStyle(0);
  c1_2->SetFrameLineStyle(0);
  c1_2->SetFrameLineWidth(2);
  c1_2->SetFrameBorderMode(0);
  c1_2->SetFrameBorderSize(10);

  // Make the error band in the ratio plot more visible
  ratioErrH->SetFillStyle(3004);
  ratioErrH->SetFillColor(kBlack);
  ratioErrH->SetLineColor(kBlack);
  ratioErrH->SetLineWidth(1);
  ratioErrH->SetMarkerStyle(0);
  ratioErrH->SetMarkerSize(0);

  // Draw ratio plot with error bars
  ratioH->Draw("E1"); // E1 draws error bars
  ratioH->SetTitle("");
  
  // Draw error band
  ratioErrH->Draw("E2 same");
  
  // Draw reference line at y=1 
  TLine *line = new TLine(ratioH->GetXaxis()->GetXmin(), 1.0, 
                         ratioH->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(kRed);
  line->SetLineStyle(2); // Dashed
  line->SetLineWidth(2);
  line->Draw();

  c1_2->Modified();
  c1_2->RedrawAxis();
  c1->cd();

  // Add legend
  TLegend *leg = new TLegend(0.75, 0.7, 0.95, 0.90);
  leg->SetTextSize(0.03);
  leg->AddEntry(_Data_var, "Data", "p");
  leg->AddEntry(_TTbar_stack, "TTbar", "f");
  leg->AddEntry(_VJ_stack, "V+Jets", "f");
  leg->AddEntry(_VV_stack, "VV", "f");
  leg->AddEntry(_ttHto2B_stack, "ttHto2B", "f");
  leg->AddEntry(_QCD_stack, "QCD", "f");
  leg->AddEntry(bkgdErr, "MC stat.+syst. unc.", "f");
  leg->Draw("same");

  // Add chi-square text
  TPaveText *t1 = new TPaveText(0.16, 0.92, 0.3, 0.97);
  t1->SetFillColor(0);
  t1->SetTextColor(1);
  t1->SetLineWidth(2);
  t1->SetTextFont(22);
  t1->SetTextSize(0.032);
  if (PLOT_CHI_2) {
    t1->AddText(Form("#chi^{2} = %.3f", chi2));
  }
  t1->Draw("same");

  // Add SF status text
  TPaveText *t11 = new TPaveText(0.2, 0.88, 0.3, 9);
  t11->SetFillColor(0);
  t11->SetTextColor(1);
  t11->SetLineWidth(2);
  t11->SetTextFont(22);
  t11->SetTextSize(0.035);
  if (has_TTbar_syst) {
    t11->AddText("p_{T} SFs ON");
  } else {
    t11->AddText("p_{T} SFs OFF");
  }
  t11->Draw("same");

  // Save canvas
  c1->SaveAs(output_path.c_str());
  
  // Clean up memory
  if (has_TTbar_syst) {
    delete _TTbar_stack_up;
    delete _TTbar_stack_down;
  }
}