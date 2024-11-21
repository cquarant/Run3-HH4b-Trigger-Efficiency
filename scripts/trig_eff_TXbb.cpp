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
    int reb = 4;

    TFile* f = new TFile(output_root_path.c_str(), "RECREATE");
    TFile* fData = new TFile(data_path.c_str());
    TFile* fMC = new TFile(mc_path.c_str());

    TString tag_var = "FatJet1_tag_GloParT_XbbVsQCD";
    TString probe_var = "FatJet1_probe_GloParT_XbbVsQCD";

    TH1D* _mc_tag = (TH1D*)fMC->Get(tag_var);
    TH1D* _mc_probe = (TH1D*)fMC->Get(probe_var);
    TH1D* _data_tag = (TH1D*)fData->Get(tag_var);
    TH1D* _data_probe = (TH1D*)fData->Get(probe_var);

    // Check not null
    if (!_mc_tag || !_mc_probe || !_data_tag || !_data_probe) {
        std::cerr << "Error: could not find histograms in input files" << std::endl;
        return;
    }

    // Rebin histograms
    _mc_tag->Rebin(reb);
    _mc_probe->Rebin(reb);
    _mc_probe->Sumw2();
    _mc_probe->Divide(_mc_tag);

    _data_tag->Rebin(reb);
    _data_probe->Rebin(reb);
    _data_probe->Sumw2();
    _data_probe->Divide(_data_tag);

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


    set_histo_style_1D(_data_probe, "", "T_{Xbb}", "Efficiency", 0.5, 1.0);
    c1->cd();
    _data_probe->Draw("colz");
    c1->SaveAs(figure_data_path.c_str());

    set_histo_style_1D(_mc_probe, "", "T_{Xbb}", "Efficiency", 0.5, 1.0);
    c2->cd();
    _mc_probe->Draw("colz");
    c2->SaveAs(figure_mc_path.c_str());

    // Compute scale factors
    TH1D* _SF_TXbb = (TH1D*)_data_probe->Clone("SF_TXbb");
    _SF_TXbb->Sumw2();
    _SF_TXbb->Divide(_mc_probe);

    set_histo_style_1D(_SF_TXbb, "", "T_{Xbb}", "Data / MC", 0.5, 1.5);
    c3->cd();
    _SF_TXbb->Draw("colz");
    c3->SaveAs(figure_sf_path.c_str());

    // Save histograms to ROOT file
    f->cd();
    _SF_TXbb->Write();
    f->Write();
    f->Close();

    // Clean up
    delete fData;
    delete fMC;
    delete c1;
    delete c2;
    delete c3;
}