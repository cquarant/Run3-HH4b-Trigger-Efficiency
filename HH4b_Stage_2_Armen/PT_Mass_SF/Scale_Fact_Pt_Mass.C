#include "HttStyles.cc"
#include "TGaxis.h"
#include "TRandom.h"
#include "TStyle.h"

#include "TCanvas.h"
#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>

void Scale_Fact_Pt_Mass() {
  int reb = 2;
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  SetStyle();

  TString PlotTitle = "           HLT_AK8PFJet250    ( |#eta| < 1.4 )";

  TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle, 100, 52, 800, 686);
  //  TCanvas *c1 = new TCanvas("", "",100,52,800,686);
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

  // ------------>Primitives in pad: c1_1
  TPad *c1_1 = new TPad("c1_1", "c1_1", 0.1, 0.05, 0.9, 0.95);
  c1_1->Draw();
  c1_1->cd();
  c1_1->Range(-93.75, -1.992728, 531.25, 5.643644);
  c1_1->SetFillColor(0);
  c1_1->SetBorderMode(0);
  c1_1->SetBorderSize(10);
  c1_1->SetTickx(1);
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
  c1_1->SetTitle("");
  c1_1->SetGridy();

  TString ETA_REG = "Eta1";

  if (ETA_REG == "Eta2") {
    PlotTitle = "       HLT_AK8PFJet250    ( 1.4 < |#eta| < 2.5 )";
  }

  TFile *fData = new TFile("Data/TTbar/Data_Trig_Eff_Histo_JetMET_" + ETA_REG + ".root");
  TFile *fMC = new TFile("MC/QCD_Trig_Eff_Histo_" + ETA_REG + ".root");

  TH2D *_MC_Probe_pt = (TH2D *)fMC->Get("Pt_Mass");
  TH2D *_MC_ProbeM_pt = (TH2D *)fMC->Get("Pt_Mass_M");
  TH2D *_Data_Probe_pt = (TH2D *)fData->Get("Pt_Mass");
  TH2D *_Data_ProbeM_pt = (TH2D *)fData->Get("Pt_Mass_M");

  _Data_Probe_pt->Rebin2D(reb);
  _Data_ProbeM_pt->Rebin2D(reb);
  _Data_ProbeM_pt->Sumw2();
  _Data_ProbeM_pt->Divide(_Data_Probe_pt);
  _Data_ProbeM_pt->SetTitle(PlotTitle);
  _Data_ProbeM_pt->GetXaxis()->SetLabelFont(42);
  _Data_ProbeM_pt->GetXaxis()->SetLabelOffset(0.01);
  _Data_ProbeM_pt->GetXaxis()->SetTitle("P_{T}-probe");
  _Data_ProbeM_pt->GetXaxis()->SetTitleSize(0.045);
  _Data_ProbeM_pt->GetXaxis()->SetTitleOffset(1.2);
  _Data_ProbeM_pt->GetXaxis()->SetTitleFont(42);
  _Data_ProbeM_pt->GetYaxis()->SetNdivisions(510);
  _Data_ProbeM_pt->GetXaxis()->SetNdivisions(505);
  _Data_ProbeM_pt->GetYaxis()->SetLabelFont(42);
  _Data_ProbeM_pt->GetYaxis()->SetLabelSize(0.04);
  _Data_ProbeM_pt->GetYaxis()->SetLabelOffset(0.01);
  _Data_ProbeM_pt->GetYaxis()->SetTitle("Mass-Probe");
  _Data_ProbeM_pt->GetYaxis()->SetTitleOffset(1.5);
  _Data_ProbeM_pt->GetYaxis()->SetTitleSize(0.045);
  _Data_ProbeM_pt->GetYaxis()->SetTitleFont(42);
  _Data_ProbeM_pt->GetXaxis()->SetTitle("P_{T}-Probe");

  _MC_Probe_pt->Rebin2D(reb);
  _MC_ProbeM_pt->Rebin2D(reb);
  _MC_ProbeM_pt->Sumw2();
  _MC_ProbeM_pt->Divide(_MC_Probe_pt);

  _Data_ProbeM_pt->Divide(_MC_ProbeM_pt);
  _Data_ProbeM_pt->GetZaxis()->SetTitle("Efficiency");
  _Data_ProbeM_pt->SetMaximum(2.1);
  _Data_ProbeM_pt->SetMinimum(0.0);
  // _Data_ProbeM_pt->Draw("colz text");
  _Data_ProbeM_pt->Draw("col z");

  c1->Modified();
  c1->cd();
  c1->SaveAs("Mass_PT/Eff_PT_Eta_" + ETA_REG + ".png");
}