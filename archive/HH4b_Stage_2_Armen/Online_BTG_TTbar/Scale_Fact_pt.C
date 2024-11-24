#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void Scale_Fact_pt()
{
  int reb=4;
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  SetStyle();

// **********************************************************************
  TFile *f = new TFile("BTG_SF.root","RECREATE");
  TH1D *_BTG_Eff_Data   = new TH1D("BTG_Eff_Data","BTG_Eff_Data",100,0.0,1.0);
  TH1D *_BTG_Eff_MC     = new TH1D("BTG_Eff_MC","BTG_Eff_MC",100,0.0,1.0);
  _BTG_Eff_MC->Rebin(reb);
  _BTG_Eff_Data->Rebin(reb);
// **********************************************************************

  TString PlotTitle="    AK8PFJet230_SoftDropMass40_PNetBB0p06 ";
  TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle,100,52,800,686);
//  TCanvas *c1 = new TCanvas("", "",100,52,800,686);
  c1->SetGridy();
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

  // ------------>Primitives in pad: c1_1
  TPad *c1_1 = new TPad("c1_1", "c1_1",0.1,0.05,0.9,0.95);
  c1_1->Draw();
  c1_1->cd();
  c1_1->Range(-93.75,-1.992728,531.25,5.643644);
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


//  TString erA= "PreBPix";
  TString erA= "PostBPix";
//  TString erA= "All";

//  TString probe_var  = "FatJet2_HbbvsQCD";
//  TString probeM_var = "FatJet3_HbbvsQCD";
//  TString probe_var  = "FatJet2PNetMD_Xbb";
//  TString probeM_var = "FatJet3PNetMD_Xbb";
  TString probe_var  = "FatJet2PNetMD_Xbb_Legacy_AN";
  TString probeM_var = "FatJet3PNetMD_Xbb_Legacy_AN";

//  TFile* fData = new TFile("Data/Old_With_SK8_HLT/Data_"+erA+".root");
//  TFile* fMC = new TFile("MC/Old_With_SK8_HLT/Histograms_TTtoLNu2Q_"+erA+".root");
  TFile* fData = new TFile("Data/Data_"+erA+".root");
  TFile* fMC = new TFile("MC/Histograms_TTtoLNu2Q_"+erA+".root");

  TH1D * _MC_Probe_pt  = (TH1D*)fMC->Get(probe_var);
  TH1D * _MC_ProbeM_pt = (TH1D*)fMC->Get(probeM_var);
  TH1D * _Data_Probe_pt  = (TH1D*)fData->Get(probe_var);
  TH1D * _Data_ProbeM_pt = (TH1D*)fData->Get(probeM_var);

  _MC_Probe_pt->Rebin(reb);
  _MC_ProbeM_pt->Rebin(reb);
  _MC_ProbeM_pt->Sumw2();
  _MC_ProbeM_pt->Divide(_MC_Probe_pt);
  _MC_ProbeM_pt->SetTitle(PlotTitle);
  _MC_ProbeM_pt->SetMarkerStyle(23);
  _MC_ProbeM_pt->SetMarkerSize(1.3);
  _MC_ProbeM_pt->SetLineColor(kBlue);
  _MC_ProbeM_pt->SetMarkerColor(kBlue);
  _MC_ProbeM_pt->SetMaximum(1.0);
  _MC_ProbeM_pt->GetXaxis()->SetLabelFont(42);
  _MC_ProbeM_pt->GetXaxis()->SetLabelOffset(0.1);
  _MC_ProbeM_pt->GetXaxis()->SetTitle("P_{T}-probe");
  _MC_ProbeM_pt->GetXaxis()->SetTitleSize(0.045);
  _MC_ProbeM_pt->GetXaxis()->SetTitleOffset(3.5);
  _MC_ProbeM_pt->GetXaxis()->SetTitleFont(42);
  _MC_ProbeM_pt->GetYaxis()->SetNdivisions(510);
  _MC_ProbeM_pt->GetXaxis()->SetNdivisions(505);
  _MC_ProbeM_pt->GetYaxis()->SetLabelFont(42);
  _MC_ProbeM_pt->GetYaxis()->SetLabelSize(0.04);
  _MC_ProbeM_pt->GetYaxis()->SetLabelOffset(0.01);
  _MC_ProbeM_pt->GetYaxis()->SetTitle("Efficiency");
  _MC_ProbeM_pt->GetYaxis()->SetTitleOffset(1.5);
  _MC_ProbeM_pt->GetYaxis()->SetTitleSize(0.045);
  _MC_ProbeM_pt->GetYaxis()->SetTitleFont(42);
  _MC_ProbeM_pt->GetXaxis()->SetTitle("Probe - X_{bb}");
  _MC_ProbeM_pt->Draw("CP");

  _Data_Probe_pt->Rebin(reb);
  _Data_ProbeM_pt->Rebin(reb);
  _Data_ProbeM_pt->Sumw2();
  _Data_ProbeM_pt->Divide(_Data_Probe_pt);
  _Data_ProbeM_pt->SetMarkerStyle(22);
  _Data_ProbeM_pt->SetMarkerSize(1.3);
  _Data_ProbeM_pt->SetLineColor(kRed);
  _Data_ProbeM_pt->SetMarkerColor(kRed);
  _Data_ProbeM_pt->Draw("CP same");

// ******************************* Legend
  TLegend * leg = new TLegend(0.6,0.5,0.8,0.6);
  SetLegendStyle(leg);
  leg->SetTextSize(0.03);
//  leg->AddEntry(_Data_ProbeM_pt,"Data 2023 "+erA,"p");
//  leg->AddEntry(_MC_ProbeM_pt,"MC TTbar "+erA,"p");
//  leg->AddEntry(_Data_ProbeM_pt,"Data 2023 ","p");
  leg->AddEntry(_Data_ProbeM_pt,"Data ","p");
  leg->AddEntry(_MC_ProbeM_pt,"MC TTbar ","p");
  leg->Draw("same");

  c1_1->Modified();
  c1->cd();



// ********************************************************************* SF
   TH1D * _SF_pt = (TH1D*) _Data_ProbeM_pt->Clone();
   TH1D * _MC_pt = (TH1D*) _MC_ProbeM_pt->Clone();

   _MC_pt->Sumw2();
   _SF_pt->Sumw2();
   _SF_pt->Divide(_MC_pt);
   _SF_pt->SetMarkerStyle(20);

   _SF_pt->SetMarkerColor(1);
   _SF_pt->SetMarkerStyle(20);
   _SF_pt->SetMarkerSize(1.2);
   _SF_pt->SetLineColor(1);
   _SF_pt->GetXaxis()->SetLabelFont(42);
   _SF_pt->GetXaxis()->SetLabelOffset(0.02);
   _SF_pt->GetXaxis()->SetLabelSize(0.2);
   _SF_pt->GetXaxis()->SetNdivisions(505);
   _SF_pt->GetXaxis()->SetTickLength(0.07);
   _SF_pt->GetYaxis()->SetRangeUser(0.6,1.4); // ********************************************
   _SF_pt->GetYaxis()->SetNdivisions(505);
   _SF_pt->GetYaxis()->SetTitle("Data / MC");
   _SF_pt->GetYaxis()->SetLabelFont(42);
   _SF_pt->GetYaxis()->SetLabelOffset(0.1/5);
   _SF_pt->GetYaxis()->SetLabelSize(0.13);
   _SF_pt->GetYaxis()->SetTitleOffset(0.38);
   _SF_pt->GetYaxis()->SetTickLength(0.02);
   _SF_pt->GetYaxis()->SetTitleSize(0.16);
   _SF_pt->GetYaxis()->SetLabelOffset(0.01);
   _SF_pt->SetMarkerSize(1.2);




    TPad * c1_2 = new TPad("lower", "pad",0.1,0.18,0.9,0.35);
    c1_2->Draw();
    c1_2->cd();
    c1_2->SetFillColor(0);
    c1_2->SetBorderMode(0);
    c1_2->SetBorderSize(10);
    c1_2->SetTickx(1);
    c1_2->SetTicky(1);
    c1_2->SetBottomMargin(0.2);
    c1_2->SetLeftMargin(0.2);
    c1_2->SetFrameFillStyle(0);
    c1_2->SetFrameLineStyle(0);
    c1_2->SetFrameLineWidth(2);
    c1_2->SetFrameBorderMode(0);
    c1_2->SetFrameBorderSize(10);
    c1_2->SetFrameFillStyle(0);
    c1_2->SetFrameLineStyle(0);
    c1_2->SetFrameLineWidth(2);
    c1_2->SetFrameBorderMode(0);
    c1_2->SetFrameBorderSize(10);
    c1_2->SetGridy();

    _SF_pt->Draw("e1");
    _SF_pt->SetTitle("");


     for(int i=1;i<=30;i++)
       cout<<_SF_pt->GetBinContent(i)<<endl;

    TF1 *Mypol   = new TF1("Mypol","pol1",0.8,1.0);
    Mypol->SetLineColor(2);
    _SF_pt->Fit("Mypol","R");

   TH1D * hint   = new TH1D("hint","",100,0.8,1.0);
   TH1D * hint_68 = new TH1D("hint_68","",500,0.8,1.0);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68,0.68);

   hint   ->SetFillColor(kCyan);
   hint_68->SetFillColor(kYellow);
   hint_68->SetMarkerSize(0.0);

//   hint    ->Draw("e3same");
   hint_68 ->Draw("e3same");
   _SF_pt  ->Draw("e1same");

    TLegend * leg2 = new TLegend(0.8,0.7,0.99,0.85);
    SetLegendStyle(leg2);
    leg2->SetTextSize(0.1);
    leg2->AddEntry(hint_68,"#sigma","f");
//    leg2->AddEntry(hint,"2#sigma","f");
    leg2->Draw("same");


    c1_2->Modified();
    c1_2->RedrawAxis();
    c1->cd();


// ****************************************************************************
// Eff Histrograms to root

   for(int i=1;i<=_Data_ProbeM_pt->GetNbinsX();i++)
    {
     _BTG_Eff_Data->SetBinContent(i,_Data_ProbeM_pt->GetBinContent(i));
     _BTG_Eff_Data->SetBinError(i,_Data_ProbeM_pt->GetBinError(i));
     _BTG_Eff_MC->SetBinContent(i,_MC_ProbeM_pt->GetBinContent(i));
     _BTG_Eff_MC->SetBinError(i,_MC_ProbeM_pt->GetBinError(i));
    }

  f->Write();

// ****************************************************************************


}