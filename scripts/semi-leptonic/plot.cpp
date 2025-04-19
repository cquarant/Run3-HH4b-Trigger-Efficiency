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

#define PLOT_CHI_2 true

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

  TH1D *_Data_var = (TH1D *)f_Data->Get(variable);
  TH1D *_QCD_var = (TH1D *)f_QCD->Get(variable);
  TH1D *_TTbar_var = (TH1D *)f_TTbar->Get(variable);
  TH1D *_VV_var = (TH1D *)f_VV->Get(variable);
  TH1D *_VJ_var = (TH1D *)f_VJ->Get(variable);
  TH1D *_ttHto2B_var = (TH1D *)f_ttHto2B->Get(variable);

  // Enable proper error calculation for all histograms
  _Data_var->Sumw2();
  _QCD_var->Sumw2();
  _TTbar_var->Sumw2();
  _VV_var->Sumw2();
  _VJ_var->Sumw2();
  _ttHto2B_var->Sumw2();

  // Float_t kfact = 0.9547;
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

  // chech if var contains "Xbb"
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
  }
  

  _ttHto2B_var->Add(_QCD_var);  // ttHto2B on top of QCD
  _VV_var->Add(_ttHto2B_var);   // VV on top of ttHto2B
  _VJ_var->Add(_VV_var);        // VJ on top of VV
  _TTbar_var->Add(_VJ_var);     // TTbar on top of VJ

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

  _TTbar_var->GetXaxis()->SetLabelFont(42);
  _TTbar_var->GetXaxis()->SetLabelOffset(0.15);
  _TTbar_var->GetXaxis()->SetTitle(XTitle);
  _TTbar_var->GetXaxis()->SetTitleSize(0.05);
  _TTbar_var->GetXaxis()->SetTitleOffset(3.3);
  _TTbar_var->GetXaxis()->SetTitleFont(42);
  _TTbar_var->GetYaxis()->SetNdivisions(505);
  _TTbar_var->GetYaxis()->SetLabelFont(42);
  _TTbar_var->GetYaxis()->SetLabelSize(0.04);
  _TTbar_var->GetYaxis()->SetLabelOffset(0.01);
  _TTbar_var->GetYaxis()->SetTitle(YTitle_10);
  _TTbar_var->GetYaxis()->SetTitleOffset(1.5);
  _TTbar_var->GetYaxis()->SetTitleSize(0.05);
  _TTbar_var->GetYaxis()->SetTitleFont(42);
  _TTbar_var->SetTitle("");

  double ymin = 10e-2;
  double ymax = 5 * (_TTbar_var->GetMaximum());
  
  // Drawing order needs to match the stacking order (drawing in reverse):
  _TTbar_var->Draw("HIST");
  _TTbar_var->SetFillColor(kOrange - 2);
  _TTbar_var->SetLineWidth(1);
  _TTbar_var->SetLineStyle(1);
  _TTbar_var->SetMaximum(ymax);
  _TTbar_var->SetMinimum(ymin);
  
  _VJ_var->SetFillColor(kGreen - 3);
  _VJ_var->SetLineWidth(1);
  _VJ_var->SetLineStyle(1);
  _VJ_var->SetMaximum(ymax);
  _VJ_var->SetMinimum(ymin);
  _VJ_var->Draw("HIST same");
  
  _VV_var->SetFillColor(kGreen);
  _VV_var->SetLineWidth(1);
  _VV_var->SetLineStyle(1);
  _VV_var->SetMaximum(ymax);
  _VV_var->SetMinimum(ymin);
  _VV_var->Draw("HIST same");
  
  _ttHto2B_var->SetFillColor(kBlue);
  _ttHto2B_var->SetLineWidth(1);
  _ttHto2B_var->SetLineStyle(1);
  _ttHto2B_var->SetMaximum(ymax);
  _ttHto2B_var->SetMinimum(ymin);
  _ttHto2B_var->Draw("HIST same");
  
  _QCD_var->Draw("HIST same");
  _QCD_var->SetFillColor(43);
  _QCD_var->SetLineWidth(1);
  _QCD_var->SetLineStyle(1);
  _QCD_var->SetMaximum(ymax);
  _QCD_var->SetMinimum(ymin);

  TH1D *bkgdErr = (TH1D *)_TTbar_var->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3004); // More visible hatching pattern
  bkgdErr->SetFillColor(kBlack); // Black is more visible
  bkgdErr->SetLineColor(kBlack);
  bkgdErr->SetLineWidth(1);
  bkgdErr->SetMarkerStyle(0);
  bkgdErr->SetMarkerSize(0);

  // Calculate total errors for the MC histogram
  for (int iB = 1; iB <= _TTbar_var->GetNbinsX(); ++iB) {
    float eStat = _TTbar_var->GetBinError(iB);
    // You might want to add systematic uncertainties here if available
    float Err = TMath::Sqrt(eStat * eStat);
    bkgdErr->SetBinError(iB, Err);
    bkgdErr->SetBinContent(iB, _TTbar_var->GetBinContent(iB));
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

  float chi2 = 0;
  for (int iB = 1; iB <= _Data_var->GetSize(); ++iB) {
    float xData = _Data_var->GetBinContent(iB);
    float xMC = _TTbar_var->GetBinContent(iB);
    if (xMC > 1e-8) {
      float diff2 = (xData - xMC) * (xData - xMC);
      chi2 += diff2 / xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

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
    float mcVal = _TTbar_var->GetBinContent(iB);
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

  TLegend *leg = new TLegend(0.86, 0.4, 0.96, 0.75);
  leg = new TLegend(0.75, 0.7, 0.95, 0.90);
  leg->SetTextSize(0.03);
  leg->AddEntry(_Data_var, "Data", "p");
  leg->AddEntry(_TTbar_var, "TTbar", "f");
  leg->AddEntry(_VJ_var, "V+Jets", "f");
  leg->AddEntry(_VV_var, "VV", "f");
  leg->AddEntry(_ttHto2B_var, "ttHto2B", "f");
  leg->AddEntry(_QCD_var, "QCD", "f");
  leg->AddEntry(bkgdErr, "MC stat. unc.", "f");
  leg->Draw("same");

  TPaveText *t1 = new TPaveText(0.16, 0.92, 0.3, 0.97);
  t1->SetFillColor(0);
  t1->SetTextColor(1);
  t1->SetLineWidth(2);
  t1->SetTextFont(22);
  t1->SetTextSize(0.032);
  // t1->AddText("                                                           CMS "
  //             "(2023 C-D EGamma + Muon) ~ 21.7 fb^{-1} #it{preliminary}");
  if (PLOT_CHI_2) {
    t1->AddText(Form("#chi^{2} = %.3f", chi2));
  }
  t1->Draw("same");

  TPaveText *t11 = new TPaveText(0.2, 0.88, 0.3, 9);
  t11->SetFillColor(0);
  t11->SetTextColor(1);
  t11->SetLineWidth(2);
  t11->SetTextFont(22);
  t11->SetTextSize(0.035);
  t11->AddText("p_{T} SFs OFF");

  // TImage *img = TImage::Create();
  // img->FromPad(c1);
  // img->WriteImage(output_path.c_str());
  c1->SaveAs(output_path.c_str());
}