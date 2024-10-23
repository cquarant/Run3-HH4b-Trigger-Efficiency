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

  int rebPT=1;
  int rebM=1;

// *******************************************************************
  TFile *f = new TFile("PT_Mass_2dSF.root","RECREATE");
  Float_t Lower_m[13]   = {20,30,40,50,60,80,100,120,150,200,250,300,350};
  Float_t Lower_pt[13]  = {300,320,350,370,400,450,500,550,600,700,800,1000,1200};
  TH2D *_Eff_Data        = new TH2D("Eff_Data","Eff_Data",12,Lower_pt,12,Lower_m);
  TH2D *_Eff_MC          = new TH2D("Eff_MC","Eff_MC",12,Lower_pt,12,Lower_m);

  _Eff_Data->Rebin2D(rebPT,rebM);
  _Eff_MC->Rebin2D(rebPT,rebM);

  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  SetStyle();
  gStyle->SetPaintTextFormat("1.2f");

  TString PlotTitle="HLT_AK8PFJet425_SoftDropMass40     (|#eta| < 2.5)";

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

//  TFile* fData = new TFile("Data/PreBPix/Data_Trig_Eff_Histo_JetMET.root");
//  TFile* fMC = new TFile("MC/PreBPix/QCD_Trig_Eff_Histo.root");
  TFile* fData = new TFile("Data/PostBPix/Data_Trig_Eff_Histo_JetMET.root");
  TFile* fMC = new TFile("MC/PostBPix/QCD_Trig_Eff_Histo.root");

  TH2D * _MC_Tag  = (TH2D*)fMC->Get("FatJet1_Pt_Mass");
  TH2D * _MC_TagM = (TH2D*)fMC->Get("FatJet1_Pt_Mass_M");
  TH2D * _MC_Probe  = (TH2D*)fMC->Get("FatJet2_Pt_Mass");
  TH2D * _MC_ProbeM = (TH2D*)fMC->Get("FatJet2_Pt_Mass_M");
  TH2D * _Data_Tag  = (TH2D*)fData->Get("FatJet1_Pt_Mass");
  TH2D * _Data_TagM = (TH2D*)fData->Get("FatJet1_Pt_Mass_M");
  TH2D * _Data_Probe  = (TH2D*)fData->Get("FatJet2_Pt_Mass");
  TH2D * _Data_ProbeM = (TH2D*)fData->Get("FatJet2_Pt_Mass_M");

// ******************************************************************* MC
  _MC_Tag->Rebin2D(rebPT,rebM);
  _MC_TagM->Rebin2D(rebPT,rebM);
  _MC_Probe->Rebin2D(rebPT,rebM);
  _MC_ProbeM->Rebin2D(rebPT,rebM);

  _MC_Probe->Add(_MC_Tag);
  _MC_ProbeM->Add(_MC_TagM);

  _MC_ProbeM->Sumw2();
  _MC_ProbeM->Divide(_MC_Probe);

  _Data_Probe ->Rebin2D(rebPT,rebM);
  _Data_ProbeM->Rebin2D(rebPT,rebM);
  _Data_Tag ->Rebin2D(rebPT,rebM);
  _Data_TagM->Rebin2D(rebPT,rebM);

  _Data_Probe->Add(_Data_Tag);
  _Data_ProbeM->Add(_Data_TagM);

  _Data_TagM->Sumw2();
  _Data_ProbeM->Divide(_Data_Probe);

//  _Data_ProbeM=_MC_ProbeM;     // For efficiency
//  _MC_ProbeM=_Data_ProbeM;     // For efficiency

  _Data_ProbeM->SetMarkerStyle(23);
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
  _Data_ProbeM->GetZaxis()->SetTitle("Data / MC");


//  _Data_ProbeM=_MC_ProbeM;
//  _Data_ProbeM->Divide(_MC_ProbeM);


  _Data_ProbeM->SetMaximum(1.8);
//  _Data_ProbeM->SetMaximum(1.0);
  _Data_ProbeM->SetMinimum(0.0);
//  _Data_ProbeM->Draw("colz");
  _Data_ProbeM->Draw("colz text");

   c1->Modified();
   c1->cd();
   c1->SaveAs("Mass_PT/Eff_PT.png");


// ************************************* Fill SFs Histo 


   for(int i=0;i<=_Data_ProbeM->GetNbinsX();i++)
    for(int j=0;j<=_Data_ProbeM->GetNbinsY();j++)
    {
     _Eff_Data->SetBinContent(i,j,_Data_ProbeM->GetBinContent(i,j));
     _Eff_MC->SetBinContent(i,j,_MC_ProbeM->GetBinContent(i,j));
    }


  f->Write();

}