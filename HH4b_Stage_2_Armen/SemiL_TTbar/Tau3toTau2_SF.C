#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void Tau3toTau2_SF()
{

  int reb =5;

// **********************************************************************
  TFile *f = new TFile("Tau3OverTau2_SF.root","RECREATE");
  TH1D *_SF   = new TH1D("SF","SF",100,0.0,1.0);
  _SF->Rebin(reb);
// **********************************************************************
   TH1D * hint_68 = new TH1D("hint_68","",370,0.18,0.5);


  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);

  TString PlotTitle=" ";

  TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle,100,52,800,686);
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
  c1_1->SetLogy();
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
  c1_1->SetGridy();

  TString Era = "";
//  Era="PostEE";

  TString JET = "FatJet1";
//  TString JET = "HBTG";
//  TString VAR = "_pt";
//  TString VAR = "_eta";
//  TString VAR = "_MassSD";
  TString VAR = "_Tau3OverTau2";
//  TString VAR = "PNet_QCD";
//  TString VAR = "PNetHbbvsQCD";
//  TString VAR = "PNetMD_Xbb";
//  TString VAR = "PNetMD_Xbb_Modified_as_in_AN_Xbb";

  TString variable = JET+VAR;
//  variable="MET";
//  variable="Dr_LFJ";
//  variable="Lep1_Pt";
  TString XTitle = "FJet - #tau_{3 }/#tau_{2}";

/*
  TFile* f_Data   = new TFile("Data/"+Era+"Histograms_Data.root");
  TFile* f_QCD    = new TFile("MC/"+Era+"Histograms_QCD.root");
  TFile* f_TTbar  = new TFile("MC/"+Era+"Histograms_TTbar.root");
  TFile* f_VV     = new TFile("MC/"+Era+"Histograms_VV.root");
  TFile* f_VJ     = new TFile("MC/"+Era+"Histograms_VJ.root");
*/

  TFile* f_Data   = new TFile("NTuples/Full_2023/Histograms_Data.root");
  TFile* f_QCD    = new TFile("NTuples/Full_2023/Histograms_QCD.root");
  TFile* f_VV     = new TFile("NTuples/Full_2023/Histograms_VV.root");
  TFile* f_VJ     = new TFile("NTuples/Full_2023/Histograms_VJ.root");
  TFile* f_WtoLNu_4Jets     = new TFile("NTuples/Full_2023/Histograms_WtoLNu_4Jets.root");

  TFile* f_TTbar  = new TFile("NTuples/Full_2023/Histograms_TTbar.root");


  TH1D * _Data_var   = (TH1D*)f_Data->Get(variable);

  TH1D * _QCD_var    = (TH1D*)f_QCD->Get(variable);
  TH1D * _TTbar_var  = (TH1D*)f_TTbar->Get(variable);
  TH1D * _VV_var     = (TH1D*)f_VV->Get(variable);
  TH1D * _VJ_var     = (TH1D*)f_VJ->Get(variable);

  TH1D * _WtoLNu_4Jets     = (TH1D*)f_WtoLNu_4Jets->Get(variable);
  _VJ_var->Add(_WtoLNu_4Jets);


  Float_t kfact=0.92;
  _QCD_var->Scale(kfact);
  _TTbar_var->Scale(kfact);
  _VJ_var->Scale(kfact);
  _VV_var->Scale(kfact);

    for (int iB=1; iB<=_Data_var->GetSize(); ++iB) {
      _Data_var->SetBinContent(iB,_Data_var->GetBinContent(iB) - _VJ_var->GetBinContent(iB) - _VV_var->GetBinContent(iB) - _QCD_var->GetBinContent(iB));
    }



  _Data_var->Rebin(reb);
  _VJ_var->Rebin(reb);
  _VV_var->Rebin(reb);
  _QCD_var->Rebin(reb);
  _TTbar_var->Rebin(reb);

// ********************************************************************* Data _vs MC
  _TTbar_var->GetXaxis()->SetLabelFont(42);
  _TTbar_var->GetXaxis()->SetLabelOffset(0.15);
  _TTbar_var->GetXaxis()->SetTitle(XTitle);
  _TTbar_var->GetXaxis()->SetTitleSize(0.05);
  _TTbar_var->GetXaxis()->SetTitleOffset(3.3);
  _TTbar_var->GetXaxis()->SetTitleFont(42);
  _TTbar_var->GetXaxis()->SetNdivisions(505);
  _TTbar_var->GetYaxis()->SetNdivisions(505);
  _TTbar_var->GetYaxis()->SetLabelFont(42);
  _TTbar_var->GetYaxis()->SetLabelSize(0.04);
  _TTbar_var->GetYaxis()->SetLabelOffset(0.01);
  _TTbar_var->GetYaxis()->SetTitle("Events");
  _TTbar_var->GetYaxis()->SetTitleOffset(1.5);
  _TTbar_var->GetYaxis()->SetTitleSize(0.05);
  _TTbar_var->GetYaxis()->SetTitleFont(42);
  _TTbar_var->SetTitle("");
//  _TTbar_var->Scale(Offline_BTG_SF);

  double ymin=10e-2;
  double ymax=5*(_TTbar_var->GetMaximum());

  _TTbar_var->Draw("HIST");
  _TTbar_var->SetFillColor(kOrange-2);
  _TTbar_var->SetLineWidth(1);
  _TTbar_var->SetLineStyle(1);
  _TTbar_var->SetMaximum(ymax);
  _TTbar_var->SetMinimum(ymin);

  _Data_var->Draw("E same");
  _Data_var->SetMarkerColor(1);
  _Data_var->SetMarkerStyle(20);
  _Data_var->SetMarkerSize(1.1);
  _Data_var->SetLineColor(1);
  _Data_var->SetLineWidth(3);
  _Data_var->SetMaximum(ymax);
  _Data_var->SetMinimum(ymin);

    TLegend * leg = new TLegend(0.2,0.4,0.3,0.75);
    leg = new TLegend(0.22,0.78,0.5,0.88);
    SetLegendStyle(leg);
    leg->SetTextSize(0.03);
    leg->AddEntry(_Data_var,"Data (bkg. subtracted)","p");
    leg->AddEntry(_TTbar_var,"TTbar","f");
    leg->Draw("same");

   c1_1->Modified();
   c1->cd();

// ********************************************************************* SF

   TH1D * _SF_pt = (TH1D*) _Data_var->Clone();
   TH1D * _MC_pt = (TH1D*) _TTbar_var->Clone();

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
   _SF_pt->GetYaxis()->SetRangeUser(0.3,1.7); // ********************************************
   _SF_pt->GetYaxis()->SetNdivisions(505);
   _SF_pt->GetYaxis()->SetTitle("Scale factor");
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


    TF1 *Mypol   = new TF1("Mypol","pol5",0.1,1.0);
    _SF_pt->Fit("Mypol","R");

   TH1D * hint   = new TH1D("hint","",185,0.1,1.0);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68,0.68);

    hint   ->SetFillColor(kCyan);
    hint_68->SetFillColor(kYellow);

  hint    ->Draw("e3same");
  hint_68 ->Draw("e3same");
  _SF_pt  ->Draw("e1same");


    TLegend * leg1 = new TLegend(0.25,0.3,0.45,0.75);
    SetLegendStyle(leg1);
    leg1->SetTextSize(0.15);
    leg1->AddEntry(Mypol,"pol5","l");
    leg1->AddEntry(hint_68,"1#sigma","f");
    leg1->AddEntry(hint,"2#sigma","f");
    leg1->Draw("same");


    c1_2->Modified();
    c1_2->RedrawAxis();
    c1->cd();

//*********************************************************************
// Eff Histrograms to root

   for(int i=1;i<=_SF_pt->GetNbinsX();i++)
    {
     _SF->SetBinContent(i,_SF_pt->GetBinContent(i));
     _SF->SetBinError(i,_SF_pt->GetBinError(i));

    }

  f->Write();
//*********************************************************************

}