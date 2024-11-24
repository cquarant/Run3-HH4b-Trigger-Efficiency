#include "TGaxis.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>

void set_pad_style(TPad* pad) {
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

void set_histo_style_1D(TH1D* hist, const std::string& x_axis_title, const std::string& y_axis_title, float min, float max) {
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

void trig_eff_TXbb(const std::string& data_path, const std::string& mc_path,
                   const std::string& output_root_path, const std::string& figure_path) {
    
    // Define the binning
    const int nBins = 20;
    Double_t xbb_bins[21] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0
    };

    // Open files
    TFile* f = new TFile(output_root_path.c_str(), "RECREATE");
    TFile* fData = new TFile(data_path.c_str());
    TFile* fMC = new TFile(mc_path.c_str());

    // Define variables
    TString tag_var = "FatJet1_tag_GloParT_XbbVsQCD";
    TString probe_var = "FatJet1_probe_GloParT_XbbVsQCD";

    // Get original histograms
    TH1D* _mc_tag_orig = (TH1D*)fMC->Get(tag_var);
    TH1D* _mc_probe_orig = (TH1D*)fMC->Get(probe_var);
    TH1D* _data_tag_orig = (TH1D*)fData->Get(tag_var);
    TH1D* _data_probe_orig = (TH1D*)fData->Get(probe_var);

    // Check if histograms exist
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

    // Calculate scale factors
    TH1D* _SF_TXbb = (TH1D*)_eff_data->Clone("SF_TXbb");
    _SF_TXbb->Sumw2();
    _SF_TXbb->Divide(_eff_mc);

    // Set global style
    gStyle->SetOptFile(0);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.2f");

    // Create canvas and pads
    TCanvas* c = new TCanvas("c", "c", 800, 800);
    c->cd();
    
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
    
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
    
    _eff_data->SetMarkerStyle(20);  // Filled circle
    _eff_data->SetMarkerColor(kRed);
    _eff_data->SetLineColor(kRed);
    _eff_data->SetMarkerSize(1.0);
    _eff_data->Draw("EP");
    
    _eff_mc->SetMarkerStyle(23);  // Filled triangle down
    _eff_mc->SetMarkerColor(kBlue);
    _eff_mc->SetLineColor(kBlue);
    _eff_mc->SetMarkerSize(1.0);
    _eff_mc->Draw("EP SAME");
    
    // Add legend
    TLegend* legend = new TLegend(0.7, 0.2, 0.90, 0.35);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(_eff_data, "Data", "ep");
    legend->AddEntry(_eff_mc, "MC TTbar", "ep");
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
    
    // Add horizontal line at 1
    TLine* line = new TLine(0, 1, 1, 1);
    line->SetLineStyle(2);
    line->Draw();
    
    // Add grid lines
    pad1->SetGrid(0, 1);
    pad2->SetGrid(0, 1);
    
    // Save canvas
    c->SaveAs(figure_path.c_str());
    
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
    delete legend;
    delete line;
    delete pad1;
    delete pad2;
    delete c;
    delete f;
    delete fData;
    delete fMC;
}