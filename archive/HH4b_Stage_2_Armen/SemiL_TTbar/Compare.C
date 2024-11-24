#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"
#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void Compare()
{

  int reb =2;
  TString JET = "FatJet2";
//  TString VAR = "_pt";
//  TString VAR = "_eta";
//  TString VAR = "_MassSD";
//  TString VAR = "PNet_QCD";
//  TString VAR = "PNetHbbvsQCD";
  TString VAR = "PNetMD_Xbb";
//  TString VAR = "_Tau3OverTau2";


  TString XTitle = JET+"_"+VAR;
  TString variable = JET+VAR;
//  variable = "MET";

  gStyle->SetOptFile(0);
  gStyle->SetOptStat("m");
  SetStyle();

  TCanvas *c1 = new TCanvas("c1", "c1",100,52,1200,550);
  gStyle->SetOptStat(0);
  c1->Range(0,0,1,1);
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

  TFile* f_QCD      = new TFile("MC/WithoutSF/Histograms_QCD.root");
  TFile* f_QCD_SF   = new TFile("MC/Histograms_QCD.root");
  TFile* f_TTbar    = new TFile("MC/WithoutSF/Histograms_TTbar.root");
  TFile* f_TTbar_SF = new TFile("MC/Histograms_TTbar.root");
  TFile* f_VV       = new TFile("MC/WithoutSF/Histograms_VV.root");
  TFile* f_VV_SF    = new TFile("MC/Histograms_VV.root");
  TFile* f_Data     = new TFile("Data/Histograms_Data.root");

  TH1D * _QCD = (TH1D*)f_QCD->Get(variable);
  TH1D * _QCD_SF = (TH1D*)f_QCD_SF->Get(variable);
  TH1D * _TTbar = (TH1D*)f_TTbar->Get(variable);
  TH1D * _TTbar_SF = (TH1D*)f_TTbar_SF->Get(variable);
  TH1D * _VV = (TH1D*)f_VV->Get(variable);
  TH1D * _VV_SF = (TH1D*)f_VV_SF->Get(variable);
  TH1D * _Data_var_10 = (TH1D*)f_Data->Get(variable);
  TH1D * _Data_var_11 = (TH1D*)f_Data->Get(variable);


// **************************************************************** Prescale
//   _QCD->Scale(0.01);
//   _QCD_SF->Scale(0.01);
//   _TTbar->Scale(0.01);
//   _TTbar_SF->Scale(0.01);
//   _VV->Scale(0.01);
//   _VV_SF->Scale(0.01);

//  _Data_var_10->Scale(1.00);
//  _Data_var_11->Scale(1.00);

// ****************************************************************
  _TTbar->Add(_VV);
  _QCD->Add(_TTbar);
  _TTbar_SF->Add(_VV_SF);
  _QCD_SF->Add(_TTbar_SF);

// ****************************************************************** Uncorrceted
// ------------>Primitives in pad: c1_1
   TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.04,0.46,0.9);
   c1_1->Draw();
   c1_1->cd();
   c1_1->Range(-93.75,-1.992728,531.25,5.643644);
   c1_1->SetFillColor(0);
   c1_1->SetBorderMode(0);
   c1_1->SetBorderSize(10);
   c1_1->SetLogy();
   c1_1->SetTickx(1);
   c1_1->SetTicky(1);
   c1_1->SetLeftMargin(0.15);
   c1_1->SetRightMargin(0.1);
   c1_1->SetTopMargin(0.0);
   c1_1->SetBottomMargin(0.3);
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


  _QCD->Rebin(reb);
  _TTbar->Rebin(reb);
  _VV->Rebin(reb);
  _Data_var_10->Rebin(reb);
//  TString YTitle_10 = "Events / xxx GeV";
  TString YTitle_10 = "Events";

  _QCD->GetXaxis()->SetLabelFont(42);
  _QCD->GetXaxis()->SetLabelOffset(0.15);
  _QCD->GetXaxis()->SetTitle(XTitle);
  _QCD->GetXaxis()->SetTitleSize(0.05);
  _QCD->GetXaxis()->SetTitleOffset(3.3);
  _QCD->GetXaxis()->SetTitleFont(42);
  _QCD->GetYaxis()->SetNdivisions(505);
  _QCD->GetYaxis()->SetLabelFont(42);
  _QCD->GetYaxis()->SetLabelSize(0.04);
  _QCD->GetYaxis()->SetLabelOffset(0.01);
  _QCD->GetYaxis()->SetTitle(YTitle_10);
  _QCD->GetYaxis()->SetTitleOffset(1.5);
  _QCD->GetYaxis()->SetTitleSize(0.05);
  _QCD->GetYaxis()->SetTitleFont(42);
  _QCD->SetTitle("");

  double ymin=10e-1;
  double ymax=5*(_QCD->GetMaximum());

  _QCD->SetFillColor(43);
  _QCD->SetLineWidth(1);
  _QCD->SetLineStyle(1);
  _QCD->SetMaximum(ymax);
  _QCD->SetMinimum(ymin);
  _QCD->Draw("HIST");

  _TTbar->SetFillColor(kOrange);
  _TTbar->SetLineWidth(1);
  _TTbar->SetLineStyle(1);
  _TTbar->SetMaximum(ymax);
  _TTbar->SetMinimum(ymin);
  _TTbar->Draw("HIST same");

  _VV->SetFillColor(kGreen);
  _VV->SetLineWidth(1);
  _VV->SetLineStyle(1);
  _VV->SetMaximum(ymax);
  _VV->SetMinimum(ymin);
  _VV->Draw("HIST same");

  _Data_var_10->Draw("E same");
  _Data_var_10->SetMarkerColor(1);
  _Data_var_10->SetMarkerStyle(20);
  _Data_var_10->SetMarkerSize(1.1);
  _Data_var_10->SetLineColor(1);
  _Data_var_10->SetLineWidth(3);
  _Data_var_10->SetMaximum(ymax);
  _Data_var_10->SetMinimum(ymin);

   c1_1->Modified();
   c1->cd();


// ************************* Errors  ratioHis
    TH1D * bkgdErr_10 = (TH1D*)_QCD->Clone("bkgdErr");
    bkgdErr_10->SetFillStyle(3013);
    bkgdErr_10->SetFillColor(1);
    bkgdErr_10->SetMarkerStyle(21);
    bkgdErr_10->SetMarkerSize(0);

    for (int iB=1; iB<=_Data_var_11->GetSize(); ++iB) {
      float eStat =  bkgdErr_10->GetBinError(iB);
      float Err = TMath::Sqrt(eStat*eStat);
      bkgdErr_10->SetBinError(iB,Err);
      bkgdErr_10->SetBinContent(iB,_QCD->GetBinContent(iB));
    }

    float chi2_10 = 0;
    for (int iB=1; iB<=_Data_var_11->GetSize(); ++iB) {
      float xData = _Data_var_10->GetBinContent(iB);
      float xMC = _QCD->GetBinContent(iB);
      if (xMC>1e-1) {
	float diff2 = (xData-xMC)*(xData-xMC);
	chi2_10 += diff2/xMC;
      }
     
    }
    std::cout << std::endl;
    std::cout << "Chi2 = " << chi2_10 << std::endl;
    std::cout << std::endl;


    TH1D * ratioH_10 = (TH1D*)_QCD->Clone("ratioH_10");
    TH1D * ratioErrH_10 = (TH1D*)bkgdErr_10->Clone("ratioErrH_10");
    ratioH_10->SetMarkerColor(1);
    ratioH_10->SetMarkerStyle(20);
    ratioH_10->SetMarkerSize(1.4);
    ratioH_10->SetLineColor(1);
    ratioH_10->GetXaxis()->SetLabelFont(42);
    ratioH_10->GetXaxis()->SetLabelOffset(0.02);
    ratioH_10->GetXaxis()->SetLabelSize(0.15);
    ratioH_10->GetXaxis()->SetNdivisions(505);
    ratioH_10->GetXaxis()->SetTickLength(0.07);

    ratioH_10->GetYaxis()->SetRangeUser(0.5,2.0);
    ratioH_10->GetYaxis()->SetNdivisions(505);
    ratioH_10->GetYaxis()->SetTitle("Data / MC");
    ratioH_10->GetYaxis()->SetLabelFont(42);
    ratioH_10->GetYaxis()->SetLabelOffset(0.1/5);
    ratioH_10->GetYaxis()->SetLabelSize(0.13);
    ratioH_10->GetYaxis()->SetTitleOffset(0.38);
    ratioH_10->GetYaxis()->SetTickLength(0.02);
    ratioH_10->GetYaxis()->SetTitleSize(0.19);
    ratioH_10->GetYaxis()->SetLabelOffset(0.01);
    ratioH_10->SetMarkerSize(1.2);

    for (int iB=1; iB<=_Data_var_11->GetSize(); ++iB) {
      float x1 = _Data_var_10->GetBinContent(iB);
      float x2 = _QCD->GetBinContent(iB);
      ratioErrH_10->SetBinContent(iB,1.0);
      ratioErrH_10->SetBinError(iB,0.0);
      float xBkg = bkgdErr_10->GetBinContent(iB);
      float errBkg = bkgdErr_10->GetBinError(iB);
      if (xBkg>0) {
	float relErr = errBkg/xBkg;
	ratioErrH_10->SetBinError(iB,relErr);
      }
      if (x1>0 && x2>0) {
	float e1 = sqrt(_QCD->GetBinError(iB)+_Data_var_10->GetBinError(iB));
	float ratio = x1/x2;
	float eratio = e1/x2;
	ratioH_10->SetBinContent(iB,ratio);
	ratioH_10->SetBinError(iB,eratio);
      }
      else {
	ratioH_10->SetBinContent(iB,1000);
      }
    }

    // ------------>Primitives in pad: c1_2
    TPad * c1_2 = new TPad("lower", "pad",0.01,0.11,0.46,0.28);
    c1_2->Draw();
    c1_2->cd();
    c1_2->SetFillColor(0);
    c1_2->SetBorderMode(0);
    c1_2->SetBorderSize(10);
    c1_2->SetGridy();
    c1_2->SetTickx(1);
    c1_2->SetTicky(1);
    c1_2->SetLeftMargin(0.15);
    c1_2->SetRightMargin(0.1);
    c1_2->SetTopMargin(0.1);
    c1_2->SetBottomMargin(0.15);
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

    ratioH_10->Draw("e1");
    ratioH_10->SetTitle("");
    ratioErrH_10->Draw("e2same");

    c1_2->Modified();
    c1_2->RedrawAxis();
    c1->cd();

// ****************************************************************** Corrceted
// ------------>Primitives in pad: c1_3
   TPad *c1_3 = new TPad("c1_3", "c1_3",0.43,0.04,0.86,0.9);
   c1_3->Draw();
   c1_3->cd();
   c1_3->Range(-93.75,-1.992728,531.25,5.643644);
   c1_3->SetFillColor(0);
   c1_3->SetBorderMode(0);
   c1_3->SetBorderSize(10);
   c1_3->SetLogy();
   c1_3->SetTickx(1);
   c1_3->SetTicky(1);
   c1_3->SetLeftMargin(0.15);
   c1_3->SetRightMargin(0.05);
   c1_3->SetTopMargin(0.0);
   c1_3->SetBottomMargin(0.3);
   c1_3->SetFrameFillStyle(0);
   c1_3->SetFrameLineStyle(0);
   c1_3->SetFrameLineWidth(2);
   c1_3->SetFrameBorderMode(0);
   c1_3->SetFrameBorderSize(10);
   c1_3->SetFrameFillStyle(0);
   c1_3->SetFrameLineStyle(0);
   c1_3->SetFrameLineWidth(2);
   c1_3->SetFrameBorderMode(0);
   c1_3->SetFrameBorderSize(10);


  _QCD_SF->Rebin(reb);
  _TTbar_SF->Rebin(reb);
  _VV_SF->Rebin(reb);
  TString YTitle_11 = "Events / xxx GeV";

  _QCD_SF->GetXaxis()->SetLabelFont(42);
  _QCD_SF->GetXaxis()->SetLabelOffset(0.15);
  _QCD_SF->GetXaxis()->SetTitle(XTitle);
  _QCD_SF->GetXaxis()->SetTitleSize(0.05);
  _QCD_SF->GetXaxis()->SetTitleOffset(3.3);
  _QCD_SF->GetXaxis()->SetTitleFont(42);
  _QCD_SF->GetYaxis()->SetNdivisions(505);
  _QCD_SF->GetYaxis()->SetLabelFont(42);
  _QCD_SF->GetYaxis()->SetLabelSize(0.04);
  _QCD_SF->GetYaxis()->SetLabelOffset(0.01);
//  _QCD_SF->GetYaxis()->SetTitle(YTitle_11);
  _QCD_SF->GetYaxis()->SetTitleOffset(1.5);
  _QCD_SF->GetYaxis()->SetTitleSize(0.05);
  _QCD_SF->GetYaxis()->SetTitleFont(42);
  _QCD_SF->SetTitle("");


//  _QCD_SF->Scale(1);

  _QCD_SF->Draw("HIST");
  _QCD_SF->SetFillColor(43);
  _QCD_SF->SetLineWidth(1);
  _QCD_SF->SetLineStyle(1);
  _QCD_SF->SetMaximum(ymax);
  _QCD_SF->SetMinimum(ymin);

  _TTbar_SF->SetFillColor(kOrange);
  _TTbar_SF->SetLineWidth(1);
  _TTbar_SF->SetLineStyle(1);
  _TTbar_SF->SetMaximum(ymax);
  _TTbar_SF->SetMinimum(ymin);
  _TTbar_SF->Draw("HIST same");

  _VV_SF->SetFillColor(kGreen);
  _VV_SF->SetLineWidth(1);
  _VV_SF->SetLineStyle(1);
  _VV_SF->SetMaximum(ymax);
  _VV_SF->SetMinimum(ymin);
  _VV_SF->Draw("HIST same");

  _Data_var_11->Draw("E same");
  _Data_var_11->SetMarkerColor(1);
  _Data_var_11->SetMarkerStyle(20);
  _Data_var_11->SetMarkerSize(1.1);
  _Data_var_11->SetLineColor(1);
  _Data_var_11->SetLineWidth(3);
  _Data_var_11->SetMaximum(ymax);
  _Data_var_11->SetMinimum(ymin);

   c1_3->Modified();
   c1->cd();


// ************************* Errors  ratioHis
    TH1D * bkgdErr_11 = (TH1D*)_QCD_SF->Clone("bkgdErr");
    bkgdErr_11->SetFillStyle(3013);
    bkgdErr_11->SetFillColor(1);
    bkgdErr_11->SetMarkerStyle(21);
    bkgdErr_11->SetMarkerSize(0);

    for (int iB=1; iB<=_Data_var_11->GetSize(); ++iB) {
      float eStat =  bkgdErr_11->GetBinError(iB);
      float Err = TMath::Sqrt(eStat*eStat);
      bkgdErr_11->SetBinError(iB,Err);
      bkgdErr_11->SetBinContent(iB,_QCD_SF->GetBinContent(iB));
    }

    float chi2_11 = 0;
    for (int iB=1; iB<=_Data_var_11->GetSize(); ++iB) {
      float xData = _Data_var_11->GetBinContent(iB);
      float xMC = _QCD_SF->GetBinContent(iB);
      if (xMC>1e-1) {
	float diff2 = (xData-xMC)*(xData-xMC);
	chi2_11 += diff2/xMC;
      }
     
    }
    std::cout << std::endl;
    std::cout << "Chi2 = " << chi2_11 << std::endl;
    std::cout << std::endl;


    TH1D * ratioH_11 = (TH1D*)_QCD_SF->Clone("ratioH_11");
    TH1D * ratioErrH_11 = (TH1D*)bkgdErr_11->Clone("ratioErrH_11");
    ratioH_11->SetMarkerColor(1);
    ratioH_11->SetMarkerStyle(20);
    ratioH_11->SetMarkerSize(1.4);
    ratioH_11->SetLineColor(1);
    ratioH_11->GetXaxis()->SetLabelFont(42);
    ratioH_11->GetXaxis()->SetLabelOffset(0.02);
    ratioH_11->GetXaxis()->SetLabelSize(0.15);
    ratioH_11->GetXaxis()->SetNdivisions(505);
    ratioH_11->GetXaxis()->SetTickLength(0.07);

    ratioH_11->GetYaxis()->SetRangeUser(0.5,2.0);
    ratioH_11->GetYaxis()->SetNdivisions(505);
//    ratioH_11->GetYaxis()->SetTitle("Data / MC");
    ratioH_11->GetYaxis()->SetLabelFont(42);
    ratioH_11->GetYaxis()->SetLabelOffset(0.1/5);
    ratioH_11->GetYaxis()->SetLabelSize(0.13);
    ratioH_11->GetYaxis()->SetTitleOffset(0.38);
    ratioH_11->GetYaxis()->SetTickLength(0.02);
    ratioH_11->GetYaxis()->SetTitleSize(0.19);
    ratioH_11->GetYaxis()->SetLabelOffset(0.01);
    ratioH_11->SetMarkerSize(1.2);

//    for (int iB=1; iB<=_Data_var_11->GetSize(); ++iB) {
    for (int iB=1; iB<=_Data_var_11->GetSize(); ++iB) {
      float x1 = _Data_var_11->GetBinContent(iB);
      float x2 = _QCD_SF->GetBinContent(iB);
      ratioErrH_11->SetBinContent(iB,1.0);
      ratioErrH_11->SetBinError(iB,0.0);
      float xBkg = bkgdErr_11->GetBinContent(iB);
      float errBkg = bkgdErr_11->GetBinError(iB);
      if (xBkg>0) {
	float relErr = errBkg/xBkg;
	ratioErrH_11->SetBinError(iB,relErr);
      }
      if (x1>0 && x2>0) {
	float e1 = sqrt(_QCD->GetBinError(iB)+_Data_var_10->GetBinError(iB));
	float ratio = x1/x2;
	float eratio = e1/x2;
	ratioH_11->SetBinContent(iB,ratio);
	ratioH_11->SetBinError(iB,eratio);
      }
      else {
	ratioH_11->SetBinContent(iB,1000);
      }
    }

    // ------------>Primitives in pad: c1_4
    TPad * c1_4 = new TPad("lower", "pad",0.43,0.11,0.86,0.28);
    c1_4->Draw();
    c1_4->cd();
    c1_4->SetFillColor(0);
    c1_4->SetBorderMode(0);
    c1_4->SetBorderSize(10);
    c1_4->SetGridy();
    c1_4->SetTickx(1);
    c1_4->SetTicky(1);
    c1_4->SetLeftMargin(0.15);
    c1_4->SetRightMargin(0.05);
    c1_4->SetTopMargin(0.1);
    c1_4->SetBottomMargin(0.15);
    c1_4->SetFrameFillStyle(0);
    c1_4->SetFrameLineStyle(0);
    c1_4->SetFrameLineWidth(2);
    c1_4->SetFrameBorderMode(0);
    c1_4->SetFrameBorderSize(10);
    c1_4->SetFrameFillStyle(0);
    c1_4->SetFrameLineStyle(0);
    c1_4->SetFrameLineWidth(2);
    c1_4->SetFrameBorderMode(0);
    c1_4->SetFrameBorderSize(10);

    ratioH_11->Draw("e1");
    ratioH_11->SetTitle("");
    ratioErrH_11->Draw("e2same");

    c1_4->Modified();
    c1_4->RedrawAxis();
    c1->cd();


// ******************************* Legend
    TLegend * leg = new TLegend(0.86,0.4,0.96,0.75);
    leg = new TLegend(0.86,0.7,0.96,0.90);
    SetLegendStyle(leg);
    leg->SetTextSize(0.03);
    leg->AddEntry(_Data_var_11,"Data","p");
    leg->AddEntry(_QCD_SF,"QCD","f");
    leg->AddEntry(_TTbar_SF,"TTbar","f");
    leg->AddEntry(_VV_SF,"VV","f");
    leg->AddEntry(ratioErrH_11,"MC stat. unc.","f");
    leg->Draw("same");

    TPaveText *t1 = new TPaveText(0.16,0.92,0.3,0.97);
    t1->SetFillColor(0);
    t1->SetTextColor(1);
    t1->SetLineWidth(2);
    t1->SetTextFont(22);
    t1->SetTextSize(0.032);
    t1->AddText("#it{                                        CMS 2022 EFG preliminary ~ 260 pb^{-1}}");
    t1->Draw("same");

    TPaveText *t11 = new TPaveText(0.2,0.88,0.3,9);
    t11->SetFillColor(0);
    t11->SetTextColor(1);
    t11->SetLineWidth(2);
    t11->SetTextFont(22);
    t11->SetTextSize(0.035);
    t11->AddText("p_{T} SFs OFF");
//    t11->Draw("same");

    TPaveText *t2 = new TPaveText(0.6,0.92,0.7,0.97);
    t2->SetFillColor(0);
    t2->SetTextColor(1);
    t2->SetLineWidth(2);
    t2->SetTextFont(22);
    t2->SetTextSize(0.032);
    t2->AddText("#it{                                        CMS 2022 EFG preliminary ~ 260 pb^{-1}}");
    t2->Draw("same");

    TPaveText *t22 = new TPaveText(0.6,0.88,0.7,0.9);
    t22->SetFillColor(0);
    t22->SetTextColor(1);
    t22->SetLineWidth(2);
    t22->SetTextFont(22);
    t22->SetTextSize(0.035);
    t22->AddText("p_{T} SFs ON");
//    t22->Draw("same");


    TImage *img = TImage::Create();
    img->FromPad(c1);
    img->WriteImage("Plots/"+variable+".png");
//    Convas_1->Close();

}