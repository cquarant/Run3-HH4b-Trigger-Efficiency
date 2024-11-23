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
  hist->SetMarkerSize(0.05);
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
  // hist->GetYaxis()->SetMoreLogLabels(kTRUE); // Show more labels in log scale
  hist->GetYaxis()->SetNoExponent(kTRUE);    // Don't use exponential notation
  hist->GetYaxis()->SetNdivisions(510);      // primary:5, secondary:10
  
  // z axis
  hist->GetZaxis()->SetTitle("Efficiency");

  // gPad->SetLogx(1);
  // gPad->SetLogy(1);

  hist->SetMaximum(1.2);
  hist->SetMinimum(0.0);
}

void trig_eff_mass_pt(const std::string &hist_mc_path,
                      const std::string &hist_data_path,
                      const std::string &output_path,
                      const std::string &figure_mc_path,
                      const std::string &figure_data_path,
                      const std::string &figure_sf_path) {  // Added new argument
    int rebPT = 1;
    int rebM = 1;

    TFile *f = new TFile(output_path.c_str(), "RECREATE");
    
    // Float_t bins_pt[9] = {300, 350, 400, 450, 500, 600, 700, 850, 1000};
    // int num_pt_bins = 8;

    // Float_t bins_m[8] = {60, 90, 120, 150, 180, 210, 240, 300};
    // int num_m_bins = 7;
    Float_t bins_pt[46] = {0.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,
                         80.0,  90.0,  100.0, 110.0, 120.0, 130.0, 140.0, 150.0,
                         160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0,
                         240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 320.0,
                         340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0,
                         500.0, 550.0, 600.0, 700.0, 800.0, 1000.0};
    int num_pt_bins = 45;

    Float_t bins_m[16] = {0.0,  5.0,   10.0,  20.0,  30.0,  40.0,  50.0,  60.0,
                          80.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 350.0};
    int num_m_bins = 15;

    TH2D *_eff_data = new TH2D("Eff_Data", "Eff_Data", num_m_bins, bins_m, num_pt_bins, bins_pt);
    TH2D *_eff_mc = new TH2D("Eff_MC", "Eff_MC", num_m_bins, bins_m, num_pt_bins, bins_pt);
    TH2D *_sf = new TH2D("SF", "Scale Factors", num_m_bins, bins_m, num_pt_bins, bins_pt);

    _eff_data->Rebin2D(rebPT, rebM);
    _eff_mc->Rebin2D(rebPT, rebM);
    _sf->Rebin2D(rebPT, rebM);

    gStyle->SetOptFile(0);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.3f");
    std::string plot_title = "";

    TCanvas *c1 = new TCanvas("Data", "Data", 100, 52, 1200, 800);
    TCanvas *c2 = new TCanvas("MC", "MC", 100, 52, 1200, 800);
    TCanvas *c3 = new TCanvas("SF", "Scale Factor", 100, 52, 1200, 800);  // Added new canvas
    
    for (auto c : {c1, c2, c3}) {  // Added c3 to the loop
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
        // c->SetLogx(1);
        // c->SetLogy(1);
    }

    TFile *fData = new TFile(hist_data_path.c_str());
    TFile *fMC = new TFile(hist_mc_path.c_str());

    TH2D *_mc_tag = (TH2D *)fMC->Get("FatJet1_tag_Mass_Pt");
    TH2D *_mc_probe = (TH2D *)fMC->Get("FatJet1_probe_Mass_Pt");
    TH2D *_data_tag = (TH2D *)fData->Get("FatJet1_tag_Mass_Pt");
    TH2D *_data_probe = (TH2D *)fData->Get("FatJet1_probe_Mass_Pt");

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
    c1->cd();
    _data_probe->Draw("colz text e");
    c1->SaveAs(figure_data_path.c_str());

    // Draw MC
    set_histo_style(_mc_probe, plot_title, x_axis_title, y_axis_title);
    c2->cd();
    _mc_probe->Draw("colz text e");
    c2->SaveAs(figure_mc_path.c_str());

    // Fill histograms and compute scale factors
    for (int i = 0; i <= _data_probe->GetNbinsX(); i++) {
        for (int j = 0; j <= _data_probe->GetNbinsY(); j++) {
            _eff_data->SetBinContent(i, j, _data_probe->GetBinContent(i, j));
            _eff_data->SetBinError(i, j, _data_probe->GetBinError(i, j));
            _eff_mc->SetBinContent(i, j, _mc_probe->GetBinContent(i, j));
            _eff_mc->SetBinError(i, j, _mc_probe->GetBinError(i, j));
            
            // Calculate scale factors
            double data_eff = _data_probe->GetBinContent(i, j);
            double mc_eff = _mc_probe->GetBinContent(i, j);
            double data_err = _data_probe->GetBinError(i, j);
            double mc_err = _mc_probe->GetBinError(i, j);
            
            if (mc_eff > 0) {
                double sf = data_eff / mc_eff;
                // Error propagation for ratio
                double sf_err = sf * sqrt(pow(data_err/data_eff, 2) + pow(mc_err/mc_eff, 2));
                _sf->SetBinContent(i, j, sf);
                _sf->SetBinError(i, j, sf_err);
            }
        }
    }

    // Draw Scale Factors
    set_histo_style(_sf, plot_title, x_axis_title, y_axis_title);
    c3->cd();
    _sf->Draw("colz text e");
    c3->SaveAs(figure_sf_path.c_str());

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