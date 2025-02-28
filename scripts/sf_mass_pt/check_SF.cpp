#include "TCanvas.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TStyle.h"
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <math.h>

double GetMassAveragedSF(TH2D* _eff_data, TH2D* _eff_mc, double ProbeJet_pt) {
    int nMassBins = _eff_data->GetNbinsX();
    Int_t bin_pt = _eff_data->GetYaxis()->FindBin(ProbeJet_pt);
    
    double sum_eff_data = 0.0;
    double sum_eff_mc = 0.0;
    double total_weight = 0.0;
    
    for (int mass_bin = 1; mass_bin <= nMassBins; mass_bin++) {
        double eff_data_bin = _eff_data->GetBinContent(mass_bin, bin_pt);
        double eff_mc_bin = _eff_mc->GetBinContent(mass_bin, bin_pt);
        
        if (eff_data_bin > 0 && eff_mc_bin > 0) {
            double weight = eff_mc_bin;
            sum_eff_data += eff_data_bin * weight;
            sum_eff_mc += eff_mc_bin * weight;
            total_weight += weight;
        }
    }
    
    double SF_mass_pt = 1.0;
    if (total_weight > 0) {
        double avg_eff_data = sum_eff_data / total_weight;
        double avg_eff_mc = sum_eff_mc / total_weight;
        if (avg_eff_mc > 0) {
            SF_mass_pt = avg_eff_data / avg_eff_mc;
        }
    }
    
    return SF_mass_pt;
}

void check_SF() {
    int reb = 4;
    TString variable = "ProbeJet_all_pt";
    TString XTitle = "Probe Jet p_{T} [GeV]";

    gStyle->SetOptFile(0);
    gStyle->SetOptStat("m");

    // Canvas setup
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

    // Open files
    TFile *f_Data = new TFile("output/Histograms_2022EE_data_QCD.root");
    TFile *f_MC = new TFile("output/Histograms_2022EE_MC_QCD.root");
    TFile *f_SF = new TFile("output/efficiency_mass_pt_2022EE_QCD.root");

    // check if files are open
    if (!f_Data->IsOpen() || !f_MC->IsOpen() || !f_SF->IsOpen()) {
        std::cout << "Error: files not found" << std::endl;
        return;
    }

    // Get efficiency histograms
    TH2D* _eff_data = (TH2D*)f_SF->Get("Eff_Data");  // Adjust name if needed
    TH2D* _eff_mc = (TH2D*)f_SF->Get("Eff_MC");      // Adjust name if needed
    // Check if histograms exist
    if (!_eff_data || !_eff_mc) {
        std::cerr << "Error: could not find histograms in input files" << std::endl;
        return;
    }

    // Get input histograms
    TH1D *_Data_var = (TH1D *)f_Data->Get(variable);
    TH1D *_QCD_var = (TH1D *)f_MC->Get(variable);
    // Check if histograms exist
    if (!_Data_var || !_QCD_var) {
        std::cerr << "Error: could not find histograms in input files" << std::endl;
        return;
    }

    // Create histogram for MC after SF application
    TH1D *_QCD_var_withSF = (TH1D*)_QCD_var->Clone("_QCD_var_withSF");
    
    // Apply mass-averaged SF to MC
    for (int bin = 1; bin <= _QCD_var_withSF->GetNbinsX(); bin++) {
        double pt = _QCD_var_withSF->GetBinCenter(bin);
        double SF = GetMassAveragedSF(_eff_data, _eff_mc, pt);
        double content = _QCD_var_withSF->GetBinContent(bin);
        double error = _QCD_var_withSF->GetBinError(bin);
        
        _QCD_var_withSF->SetBinContent(bin, content * SF);
        _QCD_var_withSF->SetBinError(bin, error * SF);
    }

    _QCD_var_withSF->Scale(0.01);  // Scale factor for visualization

    // Upper pad for distributions
    TPad *c1_1 = new TPad("c1_1", "c1_1", 0.01, 0.04, 0.75, 0.9);
    c1_1->Draw();
    c1_1->cd();
    c1_1->SetLogy();
    c1_1->SetLeftMargin(0.15);
    c1_1->SetRightMargin(0.05);
    c1_1->SetTopMargin(0.0);
    c1_1->SetBottomMargin(0.3);

    // Rebin and set styles
    _Data_var->Rebin(reb);
    _QCD_var_withSF->Rebin(reb);
    
    _QCD_var_withSF->GetXaxis()->SetTitle(XTitle);
    _QCD_var_withSF->GetXaxis()->SetTitleSize(0.05);
    _QCD_var_withSF->GetXaxis()->SetTitleOffset(3.3);
    _QCD_var_withSF->GetYaxis()->SetTitle("Events");
    _QCD_var_withSF->GetYaxis()->SetTitleOffset(1.5);
    _QCD_var_withSF->SetTitle("");

    // Set y-axis range
    double ymin = 10e-2;
    double ymax = 5 * (_QCD_var_withSF->GetMaximum());
    _QCD_var_withSF->SetMaximum(ymax);
    _QCD_var_withSF->SetMinimum(ymin);

    // Draw distributions
    _QCD_var_withSF->Draw("HIST");
    _QCD_var_withSF->SetFillColor(43);
    _QCD_var_withSF->SetLineWidth(1);

    _Data_var->Draw("E same");
    _Data_var->SetMarkerColor(1);
    _Data_var->SetMarkerStyle(20);
    _Data_var->SetMarkerSize(1.1);
    _Data_var->SetLineColor(1);
    _Data_var->SetLineWidth(3);

    // Calculate chi2
    float chi2 = 0;
    for (int iB = 1; iB <= _Data_var->GetSize(); ++iB) {
        float xData = _Data_var->GetBinContent(iB);
        float xMC = _QCD_var_withSF->GetBinContent(iB);
        if (xMC > 1e-1) {
            float diff2 = (xData - xMC) * (xData - xMC);
            chi2 += diff2 / xMC;
        }
    }
    std::cout << "Chi2 = " << chi2 << std::endl;

    // Create ratio plot
    TH1D *ratioH = (TH1D *)_Data_var->Clone("ratioH");
    ratioH->Divide(_QCD_var_withSF);
    
    // Lower pad for ratio
    TPad *c1_2 = new TPad("c1_2", "pad", 0.01, 0.11, 0.75, 0.28);
    c1_2->Draw();
    c1_2->cd();
    c1_2->SetGridy();
    c1_2->SetLeftMargin(0.15);
    c1_2->SetRightMargin(0.05);
    c1_2->SetTopMargin(0.1);
    c1_2->SetBottomMargin(0.15);

    ratioH->SetTitle("");
    ratioH->GetYaxis()->SetRangeUser(0.2, 2.0);
    ratioH->GetYaxis()->SetTitle("Data / MC");
    ratioH->GetYaxis()->SetTitleSize(0.19);
    ratioH->GetYaxis()->SetTitleOffset(0.38);
    ratioH->SetMarkerStyle(20);
    ratioH->Draw("EP");

    // Add legend
    c1->cd();
    TLegend *leg = new TLegend(0.75, 0.7, 0.95, 0.90);
    leg->SetTextSize(0.03);
    leg->AddEntry(_Data_var, "Data", "p");
    leg->AddEntry(_QCD_var_withSF, "QCD (w/ SF)", "f");
    leg->Draw("same");

    // Add CMS text
    TPaveText *t1 = new TPaveText(0.16, 0.92, 0.3, 0.97);
    t1->SetFillColor(0);
    t1->SetTextFont(22);
    t1->SetTextSize(0.032);
    t1->AddText("#it{CMS 2023 PreBPix + PostBPix preliminary ~ 217 pb^{-1}}");
    t1->Draw("same");

    // Save plot
    c1->SaveAs("sf_check/" + variable + "_withSF.png");
}