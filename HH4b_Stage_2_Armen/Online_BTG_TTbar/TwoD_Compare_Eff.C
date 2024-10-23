#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void TwoD_Compare_Eff()
{

// *******************************************************************

  int rebPT=2;
  int rebM=1;
// *******************************************************************

  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  SetStyle();
  gStyle->SetPaintTextFormat("1.2f");

  TString PlotTitle="HLT_AK8PFJet250_SoftDropMass40     (|#eta| < 1.4)";

  TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle,100,52,900,686);
//  c1->SetGridy();
  c1->Range(0,0,1,1);
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
  c1->SetBottomMargin(0.2);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);

  TString ETA_REG="";

  if(ETA_REG=="Eta2") PlotTitle="HLT_AK8PFJet250_SoftDropMass40     ( 1.4 < |#eta| < 2.5 )";

  TFile* fData = new TFile("Data/Histograms_Data.root");
  TFile* fMC = new TFile("MC/Histograms_TTbar.root");

  TH2D * _MC_Probe  = (TH2D*)fMC->Get("FatJet1_Pt_Mass");
  TH2D * _MC_ProbeM = (TH2D*)fMC->Get("FatJet2_Pt_Mass");
  TH2D * _Data_Probe  = (TH2D*)fData->Get("FatJet1_Pt_Mass");
  TH2D * _Data_ProbeM = (TH2D*)fData->Get("FatJet2_Pt_Mass");

// ******************************************************************* MC
  _MC_Probe->Rebin2D(rebPT,rebM);
  _MC_ProbeM->Rebin2D(rebPT,rebM);
  _MC_ProbeM->Sumw2();
  _MC_ProbeM->Divide(_MC_Probe);

  _Data_Probe ->Rebin2D(rebPT,rebM);
  _Data_ProbeM->Rebin2D(rebPT,rebM);
  _Data_ProbeM->Sumw2();
  _Data_ProbeM->Divide(_Data_Probe);

//  _Data_ProbeM=_MC_ProbeM;

  _Data_ProbeM->SetMarkerStyle(23);
//  _Data_ProbeM->SetLineColor(2);
//  _Data_ProbeM->SetMarkerColor(2);
  _Data_ProbeM->GetXaxis()->SetLabelFont(42);
  _Data_ProbeM->GetXaxis()->SetLabelOffset(0.02);
  _Data_ProbeM->GetXaxis()->SetTitle("Probe-P_{T}");
  _Data_ProbeM->GetXaxis()->SetTitleSize(0.045);
  _Data_ProbeM->GetXaxis()->SetTitleOffset(1.2);
  _Data_ProbeM->GetXaxis()->SetTitleFont(42);
  _Data_ProbeM->GetYaxis()->SetNdivisions(510);
  _Data_ProbeM->GetXaxis()->SetNdivisions(505);
  _Data_ProbeM->GetYaxis()->SetLabelFont(42);
  _Data_ProbeM->GetYaxis()->SetLabelSize(0.04);
  _Data_ProbeM->GetYaxis()->SetLabelOffset(0.01);
  _Data_ProbeM->GetYaxis()->SetTitle("Probe-MassSD");
  _Data_ProbeM->GetYaxis()->SetTitleOffset(1.5);
  _Data_ProbeM->GetYaxis()->SetTitleSize(0.045);
  _Data_ProbeM->GetYaxis()->SetTitleFont(42);
  _Data_ProbeM->SetTitle(PlotTitle);

  _Data_ProbeM->GetZaxis()->SetTitle("Efficiency");
//  _Data_ProbeM->SetTitle("2.3 < #eta < 4.7");

  _Data_ProbeM->Divide(_MC_ProbeM);


  _Data_ProbeM->SetMaximum(1.8);
//  _Data_ProbeM->SetMaximum(1.0);
  _Data_ProbeM->SetMinimum(0.0);
  _Data_ProbeM->Draw("colz");
//  _Data_ProbeM->Draw("colz text");

   c1->Modified();
   c1->cd();
   c1->SaveAs("Mass_PT/Eff_PT_Eta_"+ETA_REG+".png");

}