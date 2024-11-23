#include "TGaxis.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>

void set_histo_style_1D(TH1D* hist, std::string plot_title, const std::string& x_axis_title, const std::string& y_axis_title, float min, float max) {
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.2);
    hist->SetLineColor(1);
    hist->SetMarkerColor(1);
    hist->SetTitle(plot_title.c_str());

    hist->GetXaxis()->SetTitle(x_axis_title.c_str());
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelOffset(0.02);
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleFont(42);

    hist->GetYaxis()->SetTitle(y_axis_title.c_str());
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleOffset(1);
    hist->GetYaxis()->SetTitleFont(42);
    hist->SetMaximum(max);
    hist->SetMinimum(min);
}


void trig_eff_TXbb(const std::string& data_path, const std::string& mc_path,
                   const std::string& output_root_path, const std::string& figure_mc_path,
                   const std::string& figure_data_path, const std::string& figure_sf_path) {
    
    // Define the new binning
    const int nBins = 20;
    Double_t xbb_bins[21] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0
    };

    TFile* f = new TFile(output_root_path.c_str(), "RECREATE");
    TFile* fData = new TFile(data_path.c_str());
    TFile* fMC = new TFile(mc_path.c_str());

    TString tag_var = "FatJet1_tag_GloParT_XbbVsQCD";
    TString probe_var = "FatJet1_probe_GloParT_XbbVsQCD";

    // Get original histograms
    TH1D* _mc_tag_orig = (TH1D*)fMC->Get(tag_var);
    TH1D* _mc_probe_orig = (TH1D*)fMC->Get(probe_var);
    TH1D* _data_tag_orig = (TH1D*)fData->Get(tag_var);
    TH1D* _data_probe_orig = (TH1D*)fData->Get(probe_var);

    // Check not null
    if (!_mc_tag_orig || !_mc_probe_orig || !_data_tag_orig || !_data_probe_orig) {
        std::cerr << "Error: could not find histograms in input files" << std::endl;
        return;
    }

    // Create new histograms with desired binning
    TH1D* _mc_tag = new TH1D("mc_tag", "MC Tag", nBins, xbb_bins);
    TH1D* _mc_probe = new TH1D("mc_probe", "MC Probe", nBins, xbb_bins);
    TH1D* _data_tag = new TH1D("data_tag", "Data Tag", nBins, xbb_bins);
    TH1D* _data_probe = new TH1D("data_probe", "Data Probe", nBins, xbb_bins);

    // Fill new histograms
    for (int i = 1; i <= _mc_tag_orig->GetNbinsX(); i++) {
        double x = _mc_tag_orig->GetBinCenter(i);
        double content = _mc_tag_orig->GetBinContent(i);
        double error = _mc_tag_orig->GetBinError(i);
        _mc_tag->Fill(x, content);
        int newBin = _mc_tag->FindBin(x);
        _mc_tag->SetBinError(newBin, error);
    }

    for (int i = 1; i <= _mc_probe_orig->GetNbinsX(); i++) {
        double x = _mc_probe_orig->GetBinCenter(i);
        double content = _mc_probe_orig->GetBinContent(i);
        double error = _mc_probe_orig->GetBinError(i);
        _mc_probe->Fill(x, content);
        int newBin = _mc_probe->FindBin(x);
        _mc_probe->SetBinError(newBin, error);
    }

    for (int i = 1; i <= _data_tag_orig->GetNbinsX(); i++) {
        double x = _data_tag_orig->GetBinCenter(i);
        double content = _data_tag_orig->GetBinContent(i);
        double error = _data_tag_orig->GetBinError(i);
        _data_tag->Fill(x, content);
        int newBin = _data_tag->FindBin(x);
        _data_tag->SetBinError(newBin, error);
    }

    for (int i = 1; i <= _data_probe_orig->GetNbinsX(); i++) {
        double x = _data_probe_orig->GetBinCenter(i);
        double content = _data_probe_orig->GetBinContent(i);
        double error = _data_probe_orig->GetBinError(i);
        _data_probe->Fill(x, content);
        int newBin = _data_probe->FindBin(x);
        _data_probe->SetBinError(newBin, error);
    }

    // Calculate efficiencies
    TH1D* _eff_mc = (TH1D*)_mc_probe->Clone("Eff_MC");
    _eff_mc->Sumw2();
    _eff_mc->Divide(_mc_tag);

    TH1D* _eff_data = (TH1D*)_data_probe->Clone("Eff_Data");
    _eff_data->Sumw2();
    _eff_data->Divide(_data_tag);

    // Set styles and draw histograms
    gStyle->SetOptFile(0);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.2f");

    TCanvas* c1 = new TCanvas("Data", "Data", 100, 52, 1200, 800);
    TCanvas* c2 = new TCanvas("MC", "MC", 100, 52, 1200, 800);
    TCanvas* c3 = new TCanvas("ScaleFactor", "ScaleFactor", 100, 52, 1200, 800);

    for (auto c : {c1, c2, c3}) {
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
        c->SetBottomMargin(0.15);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.15);
    }

    set_histo_style_1D(_eff_data, "", "T_{Xbb}", "Efficiency", 0.5, 1.0);
    c1->cd();
    _eff_data->Draw("E1");
    c1->SaveAs(figure_data_path.c_str());

    set_histo_style_1D(_eff_mc, "", "T_{Xbb}", "Efficiency", 0.5, 1.0);
    c2->cd();
    _eff_mc->Draw("E1");
    c2->SaveAs(figure_mc_path.c_str());

    // Compute scale factors
    TH1D* _SF_TXbb = (TH1D*)_eff_data->Clone("SF_TXbb");
    _SF_TXbb->Sumw2();
    _SF_TXbb->Divide(_eff_mc);

    set_histo_style_1D(_SF_TXbb, "", "T_{Xbb}", "Data / MC", 0.5, 1.5);
    c3->cd();
    _SF_TXbb->Draw("E1");
    c3->SaveAs(figure_sf_path.c_str());

    // Save histograms to ROOT file
    f->cd();
    _eff_data->Write();
    _eff_mc->Write();
    _SF_TXbb->Write();
    f->Write();
    f->Close();

    // Clean up
    delete _mc_tag;
    delete _mc_probe;
    delete _data_tag;
    delete _data_probe;
    delete _eff_mc;
    delete _eff_data;
    delete _SF_TXbb;
    delete c1;
    delete c2;
    delete c3;
    delete f;
    delete fData;
    delete fMC;
}