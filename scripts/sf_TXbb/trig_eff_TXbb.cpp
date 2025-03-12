#include "TCanvas.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TXbb.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

#define N_BINS_FINE 500
#define POLY_FIT "pol4"

void set_pad_style(TPad *pad) {
  pad->SetFillColor(0);
  pad->SetBorderMode(0);
  pad->SetBorderSize(2);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.04);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.12);
  pad->SetFrameFillStyle(0);
  pad->SetFrameBorderMode(0);
  pad->SetFrameBorderSize(2);
}

void set_histo_style_1D(TH1D *hist, const std::string &x_axis_title,
                        const std::string &y_axis_title, float min, float max) {
  hist->SetTitle("");

  // X-axis style
  hist->GetXaxis()->SetTitle(x_axis_title.c_str());
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetTitleFont(42);

  // Y-axis style
  hist->GetYaxis()->SetTitle(y_axis_title.c_str());
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleFont(42);

  // Set range
  hist->SetMaximum(max);
  hist->SetMinimum(min);

  // Add minor ticks
  hist->GetXaxis()->SetNdivisions(510);
  hist->GetYaxis()->SetNdivisions(510);
}

void trig_eff_TXbb(const std::string &mc_path, const std::string &data_path,
                   const std::string &tagger_name,
                   const std::string &output_root_path,
                   const std::string &figure_path,
                   const std::string &TXbb_bin = "") {

  // const int rebin = 4 * (N_TXbb / 100);
  const int rebin = 4 * (N_TXbb / 100);

  // Open files
  TFile *f = new TFile(output_root_path.c_str(), "RECREATE");
  TFile *fData = new TFile(data_path.c_str());
  TFile *fMC = new TFile(mc_path.c_str());

  // Define variables
  TString tag_var = "ProbeJet_all_" + tagger_name + "_XbbVsQCD";
  TString probe_var = "ProbeJet_pass_" + tagger_name + "_XbbVsQCD";

  // Get histograms
  TH1D *_mc_all = (TH1D *)fMC->Get(tag_var);
  TH1D *_mc_pass = (TH1D *)fMC->Get(probe_var);
  TH1D *_data_all = (TH1D *)fData->Get(tag_var);
  TH1D *_data_pass = (TH1D *)fData->Get(probe_var);

  // Check if histograms exist
  if (!_mc_all || !_mc_pass || !_data_all || !_data_pass) {
    std::cerr << "Error: could not find histograms in input files" << std::endl;
    return;
  }

  // int nBins = 10;
  // Double_t bins_array[11] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
  // 0.9, 1.0}; TH1D* mc_all_reb = (TH1D*)_mc_all->Rebin(nBins, "mc_all",
  // bins_array); TH1D* mc_pass_reb = (TH1D*)_mc_pass->Rebin(nBins, "mc_pass",
  // bins_array); TH1D* data_all_reb = (TH1D*)_data_all->Rebin(nBins,
  // "data_all", bins_array); TH1D* data_pass_reb =
  // (TH1D*)_data_pass->Rebin(nBins, "data_pass", bins_array);

  // _mc_all = mc_all_reb;
  // _mc_pass = mc_pass_reb;
  // _data_all = data_all_reb;
  // _data_pass = data_pass_reb;

  // Rebin histograms
  // Check if custom binning is provided
  if (!TXbb_bin.empty()) {
    // Parse custom binning string
    std::vector<Double_t> bins;
    std::stringstream ss(TXbb_bin);
    std::string value;

    while (std::getline(ss, value, ',')) {
      bins.push_back(std::stod(value));
    }

    // Create array from vector
    const int nBins = bins.size() - 1;
    Double_t *bins_array = bins.data();

    // Rebin histograms with custom binning
    TH1D *mc_all_reb = (TH1D *)_mc_all->Rebin(nBins, "mc_all_reb", bins_array);
    TH1D *mc_pass_reb =
        (TH1D *)_mc_pass->Rebin(nBins, "mc_pass_reb", bins_array);
    TH1D *data_all_reb =
        (TH1D *)_data_all->Rebin(nBins, "data_all_reb", bins_array);
    TH1D *data_pass_reb =
        (TH1D *)_data_pass->Rebin(nBins, "data_pass_reb", bins_array);

    // Replace original histograms with rebinned ones
    _mc_all = mc_all_reb;
    _mc_pass = mc_pass_reb;
    _data_all = data_all_reb;
    _data_pass = data_pass_reb;
  } else {
    // Use default rebinning if no custom bins provided
    _mc_all->Rebin(rebin);
    _mc_pass->Rebin(rebin);
    _data_all->Rebin(rebin);
    _data_pass->Rebin(rebin);
  }

  // Calculate efficiencies
  TH1D *_eff_mc = (TH1D *)_mc_pass->Clone("eff_MC");
  _eff_mc->Sumw2();
  _eff_mc->Divide(_mc_all);

  TH1D *_eff_data = (TH1D *)_data_pass->Clone("eff_data");
  _eff_data->Sumw2();
  _eff_data->Divide(_data_all);

  // Calculate scale factors
  TH1D *_SF_TXbb = (TH1D *)_eff_data->Clone("SF_TXbb");
  _SF_TXbb->Sumw2();
  _SF_TXbb->Divide(_eff_mc);

  // Calculate additional histograms as requested
  // 1. pass(data) / pass(MC)
  TH1D *_numerator_data_over_MC =
      (TH1D *)_data_pass->Clone("numerator_data_over_MC");
  _numerator_data_over_MC->Sumw2();
  _numerator_data_over_MC->Divide(_mc_pass);

  // 2. all(data) / all(MC)
  TH1D *_denominator_data_over_MC =
      (TH1D *)_data_all->Clone("denominator_data_over_MC");
  _denominator_data_over_MC->Sumw2();
  _denominator_data_over_MC->Divide(_mc_all);

  // Set global style
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.2f");

  // Create canvas and pads for the main efficiency plot
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);

  set_pad_style(pad1);
  set_pad_style(pad2);
  pad1->SetBottomMargin(0.05);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);

  pad1->Draw();
  pad2->Draw();

  // Draw upper pad (efficiency)
  pad1->cd();
  set_histo_style_1D(_eff_data, "", "Efficiency", 0.0, 1.0);

  _eff_data->SetMarkerStyle(20); // Filled circle
  _eff_data->SetMarkerColor(kRed);
  _eff_data->SetLineColor(kRed);
  _eff_data->SetMarkerSize(1.0);
  _eff_data->Draw("EP");

  _eff_mc->SetMarkerStyle(23); // Filled triangle down
  _eff_mc->SetMarkerColor(kBlue);
  _eff_mc->SetLineColor(kBlue);
  _eff_mc->SetMarkerSize(1.0);
  _eff_mc->Draw("EP SAME");

  // Add legend
  TLegend *legend = new TLegend(0.2, 0.7, 0.4, 0.82);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(_eff_data, "Data", "ep");
  legend->AddEntry(_eff_mc, "MC", "ep");
  legend->Draw();

  // Draw lower pad (ratio)
  pad2->cd();
  set_histo_style_1D(_SF_TXbb, "T_{Xbb}", "Data / MC", 0.3, 1.7);
  _SF_TXbb->SetMarkerStyle(20);
  _SF_TXbb->GetXaxis()->SetLabelSize(0.08);
  _SF_TXbb->GetXaxis()->SetTitleSize(0.08);
  _SF_TXbb->GetXaxis()->SetTitleOffset(0.9);
  _SF_TXbb->GetYaxis()->SetLabelSize(0.08);
  _SF_TXbb->GetYaxis()->SetTitleSize(0.08);
  _SF_TXbb->GetYaxis()->SetTitleOffset(0.5);
  _SF_TXbb->Draw("EP");

  // Add a poly fit to the scale factor
  TF1 *sf_fit = new TF1("SF_TXbb_fit", POLY_FIT, _SF_TXbb->GetXaxis()->GetXmin(), _SF_TXbb->GetXaxis()->GetXmax());
  _SF_TXbb->Fit(sf_fit, "R0");
  // Evaluate the fit function on a fine grid
  TH1D *_SF_TXbb_fit_eval = new TH1D("SF_TXbb_fit_eval", "SF_TXbb_fit_eval", N_BINS_FINE, _SF_TXbb->GetXaxis()->GetXmin(), _SF_TXbb->GetXaxis()->GetXmax());
  for (int i = 1; i <= N_BINS_FINE; i++) {
      double x = _SF_TXbb_fit_eval->GetBinCenter(i);
      double y = sf_fit->Eval(x);
      _SF_TXbb_fit_eval->SetBinContent(i, y);
      // Set small errors to make the histogram plottable
      _SF_TXbb_fit_eval->SetBinError(i, 0.001);
    }


  // Add grid lines
  pad1->SetGrid(0, 1);
  pad2->SetGrid(0, 1);

  // Save main canvas
  c->SaveAs(figure_path.c_str());

  // Create second canvas for numerator and denominator ratios with SF below
  TCanvas *c2 = new TCanvas("c2", "Ratios", 800, 800);
  c2->cd();

  TPad *ratiopad1 = new TPad("ratiopad1", "ratiopad1", 0, 0.35, 1, 1.0);
  TPad *ratiopad2 = new TPad("ratiopad2", "ratiopad2", 0, 0.0, 1, 0.35);

  set_pad_style(ratiopad1);
  set_pad_style(ratiopad2);
  ratiopad1->SetBottomMargin(0.05);
  ratiopad2->SetTopMargin(0.05);
  ratiopad2->SetBottomMargin(0.3);

  ratiopad1->Draw();
  ratiopad2->Draw();

  // Set range for numerator/denominator ratios
  float ratio_min = 0.0;
  float ratio_max = 5.0;

  // Check if mc_path contains "QCD" and adjust range if needed
  if (mc_path.find("QCD") != std::string::npos) {
    ratio_min = 0.0;
    ratio_max = 0.06;
  }

  // Draw upper pad (numerator and denominator ratios)
  ratiopad1->cd();
  set_histo_style_1D(_numerator_data_over_MC, "", "Ratio", ratio_min,
                     ratio_max);

  _numerator_data_over_MC->SetMarkerStyle(22); // Triangle
  _numerator_data_over_MC->SetMarkerColor(kGreen + 2);
  _numerator_data_over_MC->SetLineColor(kGreen + 2);
  _numerator_data_over_MC->SetMarkerSize(1.2);
  _numerator_data_over_MC->Draw("EP");

  _denominator_data_over_MC->SetMarkerStyle(23); // Diamond
  _denominator_data_over_MC->SetMarkerColor(kMagenta + 1);
  _denominator_data_over_MC->SetLineColor(kMagenta + 1);
  _denominator_data_over_MC->SetMarkerSize(1.2);
  _denominator_data_over_MC->Draw("EP SAME");

  // Add legend for numerator and denominator
  TLegend *legend2 = new TLegend(0.2, 0.7, 0.5, 0.82);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(_numerator_data_over_MC, "Numerator", "ep");
  legend2->AddEntry(_denominator_data_over_MC, "Denominator", "ep");
  legend2->Draw();

  // Draw lower pad (SF)
  ratiopad2->cd();

  // Clone SF for the second plot to set different axis properties
  TH1D *_SF_TXbb2 = (TH1D *)_SF_TXbb->Clone("SF_TXbb2");
  set_histo_style_1D(_SF_TXbb2, "T_{Xbb}", "Data / MC", 0.3, 1.7);
  _SF_TXbb2->SetMarkerStyle(20);
  _SF_TXbb2->GetXaxis()->SetLabelSize(0.08);
  _SF_TXbb2->GetXaxis()->SetTitleSize(0.08);
  _SF_TXbb2->GetXaxis()->SetTitleOffset(0.9);
  _SF_TXbb2->GetYaxis()->SetLabelSize(0.08);
  _SF_TXbb2->GetYaxis()->SetTitleSize(0.08);
  _SF_TXbb2->GetYaxis()->SetTitleOffset(0.5);
  _SF_TXbb2->Draw("EP");

  // Add grid lines
  ratiopad1->SetGrid(0, 1);
  ratiopad2->SetGrid(0, 1);

  // Save second canvas
  c2->SaveAs(
      (figure_path.substr(0, figure_path.find_last_of('.')) + "_ratios.pdf")
          .c_str());

  // Save histograms and fit to ROOT file
  f->cd();
  _eff_data->Write();
  _eff_mc->Write();
  _SF_TXbb->Write();
  sf_fit->Write("SF_TXbb_fit");
  _SF_TXbb_fit_eval->Write();
  _numerator_data_over_MC->Write();
  _denominator_data_over_MC->Write();
  f->Write();
  f->Close();

  // Clean up
  delete legend;
  delete legend2;
  delete sf_fit;
  delete _SF_TXbb_fit_eval;
  delete pad1;
  delete pad2;
  delete ratiopad1;
  delete ratiopad2;
  delete _SF_TXbb2;
  delete c;
  delete c2;
  delete f;
  delete fData;
  delete fMC;
}