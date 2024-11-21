// #include "scripts/HttStyles.cc"
#include "TGaxis.h"
#include "TRandom.h"
#include "TStyle.h"

#include "TCanvas.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>

void set_histo_style(TH2D* hist, std::string plot_title, const std::string& x_axis_title, const std::string& y_axis_title) {
  hist->SetMarkerStyle(23);
  hist->SetTitle(plot_title.c_str());
  // x axis
  hist->GetXaxis()->SetTitle(x_axis_title.c_str());
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelOffset(0.02);
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetTitleFont(42);
  // hist->GetXaxis()->SetNdivisions(505);
  // y axis
  hist->GetYaxis()->SetTitle(y_axis_title.c_str());
  // hist->GetYaxis()->SetNdivisions(510);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetTitleOffset(1);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetTitleFont(42);
  // z axis
  hist->GetZaxis()->SetTitle("Efficiency");

  hist->SetMaximum(1.0);
  hist->SetMinimum(0.0);
}

void trig_eff_mass_pt(const std::string &hist_mc_path,
                      const std::string &hist_data_path,
                      const std::string &output_path,
                      const std::string &figure_mc_path,
                      const std::string &figure_data_path) {
  int rebPT = 1;
  int rebM = 1;

  TFile *f = new TFile(output_path.c_str(), "RECREATE");
  
  Float_t bins_m[13] = {0,   20,  40,  60,  80,  100, 120,
                        140, 160, 180, 200, 220, 240};
  int num_m_bins = 12;
  Float_t bins_pt[10] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
  int num_pt_bins = 9;
  TH2D *_eff_data = new TH2D("Eff_Data", "Eff_Data", num_m_bins, bins_m, num_pt_bins, bins_pt);
  TH2D *_eff_mc = new TH2D("Eff_MC", "Eff_MC", num_m_bins, bins_m, num_pt_bins, bins_pt);

  _eff_data->Rebin2D(rebPT, rebM);
  _eff_mc->Rebin2D(rebPT, rebM);

  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.2f");

  // std::string plot_title = "HLT_AK8PFJet230_SoftDropMass40 (|#eta| < 2.5)";
  std::string plot_title = "";

  TCanvas *c1 = new TCanvas("Data", "Data", 100, 52, 1200, 800);
  TCanvas *c2 = new TCanvas("MC", "MC", 100, 52, 1200, 800);
  
  for (auto c : {c1, c2}) {
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
}

  TFile *fData = new TFile(hist_data_path.c_str());
  TFile *fMC = new TFile(hist_mc_path.c_str());

  TH2D *_mc_tag = (TH2D *)fMC->Get("FatJet1_tag_Mass_Pt");
  TH2D *_mc_probe = (TH2D *)fMC->Get("FatJet1_probe_Mass_Pt");
  TH2D *_data_tag = (TH2D *)fData->Get("FatJet1_tag_Mass_Pt");
  TH2D *_data_probe = (TH2D *)fData->Get("FatJet1_probe_Mass_Pt");

  // check not null
  if (!_mc_tag || !_mc_probe || !_data_tag || !_data_probe) {
    std::cerr << "Error: could not find histograms in input files" << std::endl;
    return;
  }

  _mc_tag->Rebin2D(rebPT, rebM);
  _mc_probe->Rebin2D(rebPT, rebM);
  _mc_probe->Sumw2();
  _mc_probe->Divide(_mc_tag);

  _data_tag->Rebin2D(rebPT, rebM);
  _data_probe->Rebin2D(rebPT, rebM);
  _data_probe->Sumw2();
  _data_probe->Divide(_data_tag);


  // Draw Data
  std::string x_axis_title = "FatJet m_{SD} [GeV]";
  std::string y_axis_title = "FatJet p_{T} [GeV]";
  set_histo_style(_data_probe, plot_title, x_axis_title, y_axis_title);
  _data_probe->SetMaximum(1.2);
  _data_probe->SetMinimum(0.0);
  c1->cd();
  _data_probe->Draw("colz text");
  c1->SaveAs(figure_data_path.c_str());

  // Draw MC
  set_histo_style(_mc_probe, plot_title, x_axis_title, y_axis_title);
  _mc_probe->SetMaximum(1.2);
  _mc_probe->SetMinimum(0.0);
  c2->cd();
  _mc_probe->Draw("colz text");
  c2->SaveAs(figure_mc_path.c_str());


  // Fill SFs Histo
  for (int i = 0; i <= _data_probe->GetNbinsX(); i++) {
    for (int j = 0; j <= _data_probe->GetNbinsY(); j++) {
      _eff_data->SetBinContent(i, j, _data_probe->GetBinContent(i, j));
      _eff_data->SetBinError(i, j, _data_probe->GetBinError(i, j));
      _eff_mc->SetBinContent(i, j, _mc_probe->GetBinContent(i, j));
      _eff_mc->SetBinError(i, j, _mc_probe->GetBinError(i, j));
    }
  }

  f->Write();
}