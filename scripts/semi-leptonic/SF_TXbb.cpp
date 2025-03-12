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
#define SF_min 0
#define SF_max 2.0


void SF_TXbb(const std::string &year,      // 2022, 2023
             const std::string &path_data,    // path to the data root file
             const std::string &path_QCD,     // path to the QCD root file
             const std::string &path_VV,      // path to the VV root file
             const std::string &path_VJ,      // path to the VJ root file
             const std::string &path_TTbar,   // path to the TTbar root file
             const std::string &path_ttHto2B, // path to the ttHto2B root file
             const std::string &output_path   // path to the output root file
) {
  int reb = 5 * (N_TXbb / 100);

  TFile *f = new TFile(output_path.c_str(), "RECREATE");
  TH1D *_SF = new TH1D("SF", "SF", N_TXbb, 0.0, 1.0);
  _SF->Rebin(reb);
  TH1D *hint_68 = new TH1D("hint_68", "", 370, 0.18, 1.0);
  
  // Define custom bins for zoomed SF
  const int nBins_zoomed = 5;
  Double_t bins_zoomed[nBins_zoomed + 1] = {0.8, 0.9, 0.94, 0.97, 0.99, 1.00};

  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);

  TString PlotTitle = " ";

  TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle, 100, 52, 800, 686);
  c1->SetGridy();
  c1->Range(0, 0, 1, 1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetFrameFillStyle(0);
  c1->SetFrameLineStyle(0);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderSize(10);

  TPad *c1_1 = new TPad("c1_1", "c1_1", 0.1, 0.05, 0.9, 0.95);
  c1_1->Draw();
  c1_1->cd();
  c1_1->Range(-93.75, -1.992728, 531.25, 5.643644);
  c1_1->SetFillColor(0);
  c1_1->SetBorderMode(0);
  c1_1->SetBorderSize(10);
  c1_1->SetTickx(1);
  c1_1->SetLogy();
  c1_1->SetTicky(1);
  c1_1->SetBottomMargin(0.35);
  c1_1->SetLeftMargin(0.2);
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
  c1_1->SetGridy();

  TString JET = "FatJet1";
  TString VAR = "_GloParT_XbbVsQCD";

  TString variable = JET + VAR;
  TString XTitle = "FatJet TXbb";

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

  //  Float_t kfact=0.92;
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

  for (int iB = 1; iB <= _Data_var->GetSize(); ++iB) {
    _Data_var->SetBinContent(
        iB, _Data_var->GetBinContent(iB) - _VJ_var->GetBinContent(iB) -
                _VV_var->GetBinContent(iB) - _QCD_var->GetBinContent(iB) - 
                _ttHto2B_var->GetBinContent(iB));
  }

  _Data_var->Rebin(reb);
  _VJ_var->Rebin(reb);
  _VV_var->Rebin(reb);
  _QCD_var->Rebin(reb);
  _TTbar_var->Rebin(reb);
  _ttHto2B_var->Rebin(reb);

  _TTbar_var->GetXaxis()->SetLabelFont(42);
  _TTbar_var->GetXaxis()->SetLabelOffset(0.15);
  _TTbar_var->GetXaxis()->SetTitle(XTitle);
  _TTbar_var->GetXaxis()->SetTitleSize(0.05);
  _TTbar_var->GetXaxis()->SetTitleOffset(3.3);
  _TTbar_var->GetXaxis()->SetTitleFont(42);
  _TTbar_var->GetXaxis()->SetNdivisions(505);
  _TTbar_var->GetYaxis()->SetNdivisions(505);
  _TTbar_var->GetYaxis()->SetLabelFont(42);
  _TTbar_var->GetYaxis()->SetLabelSize(0.04);
  _TTbar_var->GetYaxis()->SetLabelOffset(0.01);
  _TTbar_var->GetYaxis()->SetTitle("Events");
  _TTbar_var->GetYaxis()->SetTitleOffset(1.5);
  _TTbar_var->GetYaxis()->SetTitleSize(0.05);
  _TTbar_var->GetYaxis()->SetTitleFont(42);
  _TTbar_var->SetTitle("");

  double ymin = 10e-2;
  double ymax = 5 * (_TTbar_var->GetMaximum());

  _TTbar_var->Draw("HIST");
  _TTbar_var->SetFillColor(kOrange - 2);
  _TTbar_var->SetLineWidth(1);
  _TTbar_var->SetLineStyle(1);
  _TTbar_var->SetMaximum(ymax);
  _TTbar_var->SetMinimum(ymin);

  _Data_var->Draw("E same");
  _Data_var->SetMarkerColor(1);
  _Data_var->SetMarkerStyle(20);
  _Data_var->SetMarkerSize(1.1);
  _Data_var->SetLineColor(1);
  _Data_var->SetLineWidth(3);
  _Data_var->SetMaximum(ymax);
  _Data_var->SetMinimum(ymin);

  TLegend *leg = new TLegend(0.2, 0.4, 0.3, 0.75);
  leg = new TLegend(0.22, 0.78, 0.5, 0.88);
  leg->SetTextSize(0.03);
  leg->AddEntry(_Data_var, "Data (bkg. subtracted)", "p");
  leg->AddEntry(_TTbar_var, "TTbar", "f");
  leg->Draw("same");

  c1_1->Modified();
  c1->cd();

  TH1D *_SF_pt = (TH1D *)_Data_var->Clone();
  TH1D *_MC_pt = (TH1D *)_TTbar_var->Clone();

  _MC_pt->Sumw2();
  _SF_pt->Sumw2();
  _SF_pt->Divide(_MC_pt);
  _SF_pt->SetMarkerStyle(20);

  _SF_pt->SetMarkerColor(1);
  _SF_pt->SetMarkerStyle(20);
  _SF_pt->SetMarkerSize(1.2);
  _SF_pt->SetLineColor(1);
  _SF_pt->GetXaxis()->SetLabelFont(42);
  _SF_pt->GetXaxis()->SetLabelOffset(0.02);
  _SF_pt->GetXaxis()->SetLabelSize(0.2);
  _SF_pt->GetXaxis()->SetNdivisions(505);
  _SF_pt->GetXaxis()->SetTickLength(0.07);
  _SF_pt->GetYaxis()->SetRangeUser(SF_min, SF_max);
  _SF_pt->GetYaxis()->SetNdivisions(505);
  _SF_pt->GetYaxis()->SetTitle("Scale factor");
  _SF_pt->GetYaxis()->SetLabelFont(42);
  _SF_pt->GetYaxis()->SetLabelOffset(0.1 / 5);
  _SF_pt->GetYaxis()->SetLabelSize(0.13);
  _SF_pt->GetYaxis()->SetTitleOffset(0.38);
  _SF_pt->GetYaxis()->SetTickLength(0.02);
  _SF_pt->GetYaxis()->SetTitleSize(0.16);
  _SF_pt->GetYaxis()->SetLabelOffset(0.01);
  _SF_pt->SetMarkerSize(1.2);

  TPad *c1_2 = new TPad("lower", "pad", 0.1, 0.18, 0.9, 0.35);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetFillColor(0);
  c1_2->SetBorderMode(0);
  c1_2->SetBorderSize(10);
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);
  c1_2->SetBottomMargin(0.2);
  c1_2->SetLeftMargin(0.2);
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
  c1_2->SetGridy();

  _SF_pt->Draw("e1");
  _SF_pt->SetTitle("");

  cout << "SF_pt: ";
  for (int i = 1; i <= _SF_pt->GetNbinsX(); i++)
    cout << _SF_pt->GetBinContent(i) << ", ";
  cout << endl;

  TF1 *Mypol = new TF1("Mypol", "pol8", 0.1, 1.0);
  _SF_pt->Fit("Mypol", "R");

  TH1D *hint = new TH1D("hint", "", 185, 0.1, 1.0);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68, 0.68);

  hint->SetFillColor(kCyan);
  hint_68->SetFillColor(kYellow);

  hint->Draw("e3same");
  hint_68->Draw("e3same");
  _SF_pt->Draw("e1same");

  TLegend *leg1 = new TLegend(0.25, 0.3, 0.45, 0.75);
  leg1->SetTextSize(0.15);
  leg1->AddEntry(Mypol, "pol8", "l");
  leg1->AddEntry(hint_68, "1#sigma", "f");
  leg1->AddEntry(hint, "2#sigma", "f");
  leg1->Draw("same");

  c1_2->Modified();
  c1_2->RedrawAxis();
  c1->cd();

  // Eff Histrograms to root
  for (int i = 1; i <= _SF_pt->GetNbinsX(); i++) {
    _SF->SetBinContent(i, _SF_pt->GetBinContent(i));
    _SF->SetBinError(i, _SF_pt->GetBinError(i));
  }
  
  // Create a new canvas for the zoomed SF
  TCanvas *c2 = new TCanvas("c2", "Zoomed SF", 800, 600);
  c2->SetGridy();
  
  // Create zoomed SF histogram with custom bins
  TH1D *_SF_zoomed = new TH1D("SF_zoomed", "SF_zoomed", nBins_zoomed, bins_zoomed);
  TH1D *hint_68_zoomed = new TH1D("hint_68_zoomed", "", 100, 0.8, 1.0);
  
  // Fill the zoomed SF histogram from the fit function
  for (int i = 1; i <= nBins_zoomed; i++) {
    double binCenter = _SF_zoomed->GetBinCenter(i);
    double value = Mypol->Eval(binCenter);
    _SF_zoomed->SetBinContent(i, value);
    
    // Estimate error from the fit (can be improved)
    double error = 0.0;
    if (binCenter >= 0.1 && binCenter <= 1.0) {
      // Find the closest bin in the original SF histogram
      int closestBin = _SF_pt->FindBin(binCenter);
      error = _SF_pt->GetBinError(closestBin);
    }
    _SF_zoomed->SetBinError(i, error);
  }
  
  _SF_zoomed->SetMarkerColor(kBlue);
  _SF_zoomed->SetMarkerStyle(20);
  _SF_zoomed->SetMarkerSize(1.2);
  _SF_zoomed->SetLineColor(kBlue);
  _SF_zoomed->GetXaxis()->SetTitle("FatJet TXbb");
  _SF_zoomed->GetYaxis()->SetTitle("Scale factor");
  _SF_zoomed->GetYaxis()->SetRangeUser(SF_min, SF_max);
  _SF_zoomed->SetTitle("Zoomed Scale Factor");
  _SF_zoomed->Draw("e1");
  
  // Fit the zoomed region with a polynomial
  TF1 *Mypol_zoomed = new TF1("Mypol_zoomed", "pol3", 0.8, 1.0);
  _SF_zoomed->Fit("Mypol_zoomed", "R");
  
  // Get confidence intervals for the zoomed fit
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68_zoomed, 0.68);
  hint_68_zoomed->SetFillColor(kYellow);
  hint_68_zoomed->Draw("e3same");
  _SF_zoomed->Draw("e1same");
  
  TLegend *leg2 = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(_SF_zoomed, "SF (zoomed)", "p");
  leg2->AddEntry(Mypol_zoomed, "pol3 fit", "l");
  leg2->AddEntry(hint_68_zoomed, "1#sigma confidence", "f");
  leg2->Draw();
  
  c2->Modified();
  c2->Update();

  f->cd();
  _SF_zoomed->Write();
  hint_68_zoomed->Write();
  Mypol->Write();
  Mypol_zoomed->Write();
  f->Write();
}