// #include "scripts/HttStyles.cc"
#include "TGaxis.h"
#include "TRandom.h"
#include "TStyle.h"

#include "TCanvas.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>

// #define x_min 100.0  // Minimum for mSD
// #define x_max 350.0  // Adjust based on your preferred maximum
// #define y_min 300.0  // Minimum for pT
// #define y_max 1000.0 // Adjust based on your preferred maximum
// #define x_min 50.0   // Minimum for mSD
// #define x_max 350.0  // Adjust based on your preferred maximum
// #define y_min 280.0  // Minimum for pT
// #define y_max 1000.0 // Adjust based on your preferred maximum
#define x_min 50.0     // Minimum for mSD
#define x_max 350.0    // Minimum for mSD
#define y_min 250.0    // Minimum for pT
#define y_max 1000.0   // Maximum for pT
#define SET_LOG true

double get_dynamic_marker_size(TH2D *hist, double min_size = 0.05,
                               double max_size = 1.0) {
  // Bin dimensions
  int nbinsx = hist->GetNbinsX();
  int nbinsy = hist->GetNbinsY();

  // Pad dimensions
  Double_t width = gPad->GetWw() * gPad->GetAbsWNDC();
  Double_t height = gPad->GetWh() * gPad->GetAbsHNDC();

  // Smallest bin size in pixels
  double smallest_bin_width = width / nbinsx;
  double smallest_bin_height = height / nbinsy;
  double smallest_dimension =
      TMath::Min(smallest_bin_width, smallest_bin_height);

  // Scale marker size based on bin size
  double scale_factor = smallest_dimension / 100.0; // adjustable
  double marker_size = scale_factor * max_size;

  // Clamp between min and max size
  return TMath::Max(min_size, TMath::Min(marker_size, max_size));
}

void set_histo_style(TH2D *hist, std::string plot_title,
                     const std::string &x_axis_title,
                     const std::string &y_axis_title,
                     float hist_min, float hist_max) {
  hist->SetMarkerStyle(23);
  // hist->SetMarkerSize(0.05);
  // hist->SetMarkerSize(1.0);
  hist->SetMarkerSize(0.9);
  // double marker_size = get_dynamic_marker_size(hist, 0.05, 1.0);
  // hist->SetMarkerSize(marker_size);
  hist->SetTitle(plot_title.c_str());
  // x axis settings
  hist->GetXaxis()->SetTitle(x_axis_title.c_str());
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelOffset(0.02);
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.2);
  // hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetMoreLogLabels(kTRUE); // Show more labels in log scale
  hist->GetXaxis()->SetNoExponent(kTRUE);    // Don't use exponential notation
  hist->GetXaxis()->SetNdivisions(510);      // primary:5, secondary:10

  // y axis settings
  hist->GetYaxis()->SetTitle(y_axis_title.c_str());
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetNoExponent(kTRUE); // Don't use exponential notation
  hist->GetYaxis()->SetNdivisions(510);   // primary:5, secondary:10

  // z axis
  hist->GetZaxis()->SetTitle("Efficiency");

  if (SET_LOG) {
    hist->GetYaxis()->SetMoreLogLabels(kTRUE); // Show more labels in log scale
    gPad->SetLogx(1);
    gPad->SetLogy(1);
  }

  hist->SetMaximum(hist_max);
  hist->SetMinimum(hist_min);
}

void trig_eff_mass_pt(const std::string &hist_mc_path,
                      const std::string &hist_data_path,
                      const std::string &output_path,
                      const std::string &figure_mc_path,
                      const std::string &figure_data_path,
                      const std::string &figure_sf_path) { // Added new argument
  int rebPT = 1;
  int rebM = 1;

  TFile *f = new TFile(output_path.c_str(), "RECREATE");

  // Float_t bins_pt[9] = {300, 350, 400, 450, 500, 600, 700, 850, 1000};
  // int num_pt_bins = 8;

  // Float_t bins_m[8] = {60, 90, 120, 150, 180, 210, 240, 300};
  // int num_m_bins = 7;
  // Float_t bins_pt[46] = {0.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,
  //                        80.0,  90.0,  100.0, 110.0, 120.0, 130.0, 140.0, 150.0,
  //                        160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0,
  //                        240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 320.0,
  //                        340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0,
  //                        500.0, 550.0, 600.0, 700.0, 800.0, 1000.0};
  // int num_pt_bins = 45;

  // Float_t bins_m[16] = {0.0,  5.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,
  //                       80.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 350.0};
  // int num_m_bins = 15;
  Float_t bins_pt[9] = {250, 275, 300, 350, 400, 450, 500, 600, 1000};
  int num_pt_bins = 8;

  Float_t bins_m[10] = {50, 60, 80, 100, 120, 150, 200, 250, 300, 350};
  int num_m_bins = 9;

  TH2D *_eff_data = new TH2D("Eff_Data", "Eff_Data", num_m_bins, bins_m,
                             num_pt_bins, bins_pt);
  TH2D *_eff_mc =
      new TH2D("Eff_MC", "Eff_MC", num_m_bins, bins_m, num_pt_bins, bins_pt);
  TH2D *_sf =
      new TH2D("SF", "Scale Factors", num_m_bins, bins_m, num_pt_bins, bins_pt);

  _eff_data->Rebin2D(rebPT, rebM);
  _eff_mc->Rebin2D(rebPT, rebM);
  _sf->Rebin2D(rebPT, rebM);

  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.3f");
  std::string plot_title = "";

  TCanvas *c1 = new TCanvas("Data", "Data", 100, 52, 1200, 800);
  TCanvas *c2 = new TCanvas("MC", "MC", 100, 52, 1200, 800);
  TCanvas *c3 =
      new TCanvas("SF", "Scale Factor", 100, 52, 1200, 800); // Added new canvas

  for (auto c : {c1, c2, c3}) { // Added c3 to the loop
    c->Range(0, 0, 1, 1);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(10);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFrameFillStyle(0);
    c->SetFrameLineStyle(0);
    c->SetFrameLineWidth(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderSize(10);
    c->SetBottomMargin(0.2);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.15);

    if (SET_LOG) {
      c->SetLogx(1);
      c->SetLogy(1);
    }
  }

  TFile *fData = new TFile(hist_data_path.c_str());
  TFile *fMC = new TFile(hist_mc_path.c_str());

  TH2D *_mc_all = (TH2D *)fMC->Get("ProbeJet_all_Mass_Pt");
  TH2D *_mc_pass = (TH2D *)fMC->Get("ProbeJet_pass_Mass_Pt");
  TH2D *_data_all = (TH2D *)fData->Get("ProbeJet_all_Mass_Pt");
  TH2D *_data_pass = (TH2D *)fData->Get("ProbeJet_pass_Mass_Pt");

  if (!_mc_all || !_mc_pass || !_data_all || !_data_pass) {
      std::cerr << "Error: could not find histograms in input files" << std::endl;
      return;
  }

  // Enable automatic error calculation for all histograms
  _mc_all->Sumw2();
  _mc_pass->Sumw2();
  _data_all->Sumw2();
  _data_pass->Sumw2();
  _eff_data->Sumw2();
  _eff_mc->Sumw2();
  _sf->Sumw2();

  // Apply rebinning
  _mc_all->Rebin2D(rebPT, rebM);
  _mc_pass->Rebin2D(rebPT, rebM);
  _data_all->Rebin2D(rebPT, rebM);
  _data_pass->Rebin2D(rebPT, rebM);

  // Calculate efficiencies with automatic error propagation
  _eff_data->Divide(_data_pass, _data_all, 1.0, 1.0, "B");  // "B" for binomial errors
  _eff_mc->Divide(_mc_pass, _mc_all, 1.0, 1.0, "B");

  // Calculate scale factors with automatic error propagation
  _sf->Divide(_eff_data, _eff_mc);

  // Draw Data efficiency
  std::string x_axis_title = "FatJet m_{SD} [GeV]";
  std::string y_axis_title = "FatJet p_{T} [GeV]";

  set_histo_style(_eff_data, plot_title, x_axis_title, y_axis_title, 0, 1.2);
  c1->cd();
  _eff_data->GetXaxis()->SetRangeUser(x_min, x_max);
  _eff_data->GetYaxis()->SetRangeUser(y_min, y_max);
  _eff_data->Draw("colz text e");
  c1->SaveAs(figure_data_path.c_str());

  // Draw MC efficiency
  set_histo_style(_eff_mc, plot_title, x_axis_title, y_axis_title, 0, 1.2);
  c2->cd();
  _eff_mc->GetXaxis()->SetRangeUser(x_min, x_max);
  _eff_mc->GetYaxis()->SetRangeUser(y_min, y_max);
  _eff_mc->Draw("colz text e");
  c2->SaveAs(figure_mc_path.c_str());

  // Draw Scale Factors
  set_histo_style(_sf, plot_title, x_axis_title, y_axis_title, 0.5, 1.5);
  c3->cd();
  _sf->GetXaxis()->SetRangeUser(x_min, x_max);
  _sf->GetYaxis()->SetRangeUser(y_min, y_max);
  _sf->Draw("colz text e");
  c3->SaveAs(figure_sf_path.c_str());

  // Write to output file
  f->cd();
  _eff_data->Write();
  _eff_mc->Write();
  _sf->Write();
  f->Close();

  // Clean up
  delete f;
  delete fData;
  delete fMC;
  delete c1;
  delete c2;
  delete c3;
}