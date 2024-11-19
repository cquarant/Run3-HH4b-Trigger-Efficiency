#include "TGaxis.h"
#include "TStyle.h"
#include "TRandom.h"

#include "TCanvas.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>

void set_histo_style_1D(TH1D* hist, const std::string& x_axis_title, const std::string& y_axis_title) {
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.2);
    hist->SetLineColor(1);
    hist->SetMarkerColor(1);

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
    hist->GetYaxis()->SetTitleOffset(1.5);
    hist->GetYaxis()->SetTitleFont(42);
    hist->SetMaximum(1.4);
    hist->SetMinimum(0.6);
}

void trig_eff_TXbb(const std::string& data_path, const std::string& mc_path,
                   const std::string& output_path, const std::string& figure_data_path,
                   const std::string& figure_mc_path) {
    int reb = 4;

    TFile* f = new TFile(output_path.c_str(), "RECREATE");
    TFile* fData = new TFile(data_path.c_str());
    TFile* fMC = new TFile(mc_path.c_str());

    TString probe_var = "FatJet1_probe_GloParT_XbbVsQCD";
    TString tag_var = "FatJet1_tag_GloParT_XbbVsQCD";

    TH1D* _mc_probe = (TH1D*)fMC->Get(probe_var);
    TH1D* _mc_tag = (TH1D*)fMC->Get(tag_var);
    TH1D* _data_probe = (TH1D*)fData->Get(probe_var);
    TH1D* _data_tag = (TH1D*)fData->Get(tag_var);

    // Check not null
    if (!_mc_probe || !_mc_tag || !_data_probe || !_data_tag) {
        std::cerr << "Error: could not find histograms in input files" << std::endl;
        return;
    }

    // Rebin histograms
    _mc_probe->Rebin(reb);
    _mc_tag->Rebin(reb);
    _mc_tag->Sumw2();
    _mc_tag->Divide(_mc_probe);

    _data_probe->Rebin(reb);
    _data_tag->Rebin(reb);
    _data_tag->Sumw2();
    _data_tag->Divide(_data_probe);

    // Set styles and draw histograms
    TCanvas* c1 = new TCanvas("Data", "Data", 100, 52, 1200, 800);
    set_histo_style_1D(_data_tag, "Probe - X_{bb}", "Efficiency");
    c1->cd();
    _data_tag->Draw("E1");
    c1->SaveAs(figure_data_path.c_str());

    TCanvas* c2 = new TCanvas("MC", "MC", 100, 52, 1200, 800);
    set_histo_style_1D(_mc_tag, "Probe - X_{bb}", "Efficiency");
    c2->cd();
    _mc_tag->Draw("E1");
    c2->SaveAs(figure_mc_path.c_str());

    // Compute scale factors
    TH1D* _SF_TXbb = (TH1D*)_data_tag->Clone("SF_TXbb");
    _SF_TXbb->Sumw2();
    _SF_TXbb->Divide(_mc_tag);

    TCanvas* c3 = new TCanvas("ScaleFactor", "ScaleFactor", 100, 52, 1200, 800);
    set_histo_style_1D(_SF_TXbb, "Probe - T_{Xbb}", "Data / MC");
    c3->cd();
    _SF_TXbb->Draw("E1");
    c3->SaveAs(output_path.c_str());

    // Save histograms
    f->cd();
    _SF_TXbb->Write();
    f->Write();
    f->Close();
}
