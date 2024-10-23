#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"
#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void Make_Histo()
{

  int reb =5;
  TString Era = "";
//  Era = "PreEE/";
//  Era = "PostEE/";

  TString JET = "FatJet1";
//  TString JET = "HBTG";
//  TString VAR = "_pt";
//  TString VAR = "_eta"; reb=2;
//  TString VAR = "_MassSD"; reb=10;
//  TString VAR = "_Tau3OverTau2";
//  TString VAR = "_Xbb";
//  TString VAR = "_Xbb_Legacy";
  TString VAR = "_Xbb_Legacy_AN";

  TString variable = JET+VAR;
//  variable="MET";
//  variable="Dr_JmaxL";
//  variable="Dr_LFJ";
//  variable="Dr_J1FJ";
//  variable="Dr_J2FJ";
//  variable="Lep1_Pt";

  TString XTitle = variable;
// XTitle = "FJet - #tau_{3 }/#tau_{2}";
// XTitle = "FJet-P_{T}";
// XTitle = "FJet-M_{SD}";
 XTitle = "FJet - X_{bb}";
// XTitle = "FJet - T_{Xbb}";
// XTitle = "FJet - #eta";
// XTitle = "Lepton - P_{T}";

  gStyle->SetOptFile(0);
  gStyle->SetOptStat("m");
  SetStyle();

  TCanvas *c1 = new TCanvas("c1", "c1",100,50,850,750);
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

// Files

/*
  // OLD
  TFile* f_Data   = new TFile("Data/"+Era+"Histograms_Data.root");
  TFile* f_QCD    = new TFile("MC/"+Era+"Histograms_QCD.root");
  TFile* f_VV     = new TFile("MC/"+Era+"Histograms_VV.root");
  TFile* f_VJ     = new TFile("MC/"+Era+"Histograms_VJ.root");
  TFile* f_TTbar  = new TFile("MC/"+Era+"Histograms_TTbar.root");
//  TFile* f_TTbar  = new TFile("MC/"+Era+"Histograms_TTbar_T3T2_Corr.root");
//  TFile* f_TTbar  = new TFile("MC/"+Era+"Histograms_TTbar_All_Corr.root");
*/

/*
  TFile* f_Data   = new TFile("NTuples/PostEE/Histograms_Data.root");
  TFile* f_QCD    = new TFile("NTuples/PostEE/Histograms_QCD.root");
  TFile* f_VV     = new TFile("NTuples/PostEE/Histograms_VV.root");
  TFile* f_VJ     = new TFile("NTuples/PostEE/Histograms_VJ.root");
  TFile* f_WtoLNu_4Jets     = new TFile("NTuples/PostEE/Histograms_WtoLNu_4Jets.root");

  TFile* f_TTbar  = new TFile("NTuples/PostEE/Histograms_TTbar.root");
*/


  TFile* f_Data   = new TFile("NTuples/Full_2023/Histograms_Data.root");
  TFile* f_QCD    = new TFile("NTuples/Full_2023/Histograms_QCD.root");
  TFile* f_VV     = new TFile("NTuples/Full_2023/Histograms_VV.root");
  TFile* f_VJ     = new TFile("NTuples/Full_2023/Histograms_VJ.root");
  TFile* f_WtoLNu_4Jets     = new TFile("NTuples/Full_2023/Histograms_WtoLNu_4Jets.root");

//  TFile* f_TTbar  = new TFile("NTuples/Full_2023/Histograms_TTbar.root");
//  TFile* f_TTbar  = new TFile("NTuples/Full_2023/Histograms_TTbar_Corr_Tua3toTau2.root");
  TFile* f_TTbar  = new TFile("NTuples/Full_2023/Histograms_TTbar_Corr_Tua3toTau2and_Xbb.root");



// ****************************************************************** Histograms and K-factors determination

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


  _Data_var->Rebin(reb);
  _VJ_var->Rebin(reb);
  _VV_var->Rebin(reb);
  _QCD_var->Rebin(reb);
  _TTbar_var->Rebin(reb);

//  _VV_var->Add(_VJ_var);
//  _QCD_var->Add(_VV_var);
//  _TTbar_var->Add(_QCD_var);

  _VV_var->Add(_QCD_var);
  _VJ_var->Add(_VV_var);
  _TTbar_var->Add(_VJ_var);

// ****************************************************************** Uncorrceted
// ------------>Primitives in pad: c1_1
   TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.04,0.75,0.9);
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
   c1_1->SetRightMargin(0.05);
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


  TString YTitle_10 = "Events";

  _TTbar_var->GetXaxis()->SetLabelFont(42);
  _TTbar_var->GetXaxis()->SetLabelOffset(0.15);
  _TTbar_var->GetXaxis()->SetTitle(XTitle);
  _TTbar_var->GetXaxis()->SetTitleSize(0.05);
  _TTbar_var->GetXaxis()->SetTitleOffset(3.3);
  _TTbar_var->GetXaxis()->SetTitleFont(42);
  _TTbar_var->GetYaxis()->SetNdivisions(505);
  _TTbar_var->GetYaxis()->SetLabelFont(42);
  _TTbar_var->GetYaxis()->SetLabelSize(0.04);
  _TTbar_var->GetYaxis()->SetLabelOffset(0.01);
  _TTbar_var->GetYaxis()->SetTitle(YTitle_10);
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

  _VJ_var->SetFillColor(kGreen-3);
  _VJ_var->SetLineWidth(1);
  _VJ_var->SetLineStyle(1);
  _VJ_var->SetMaximum(ymax);
  _VJ_var->SetMinimum(ymin);
  _VJ_var->Draw("HIST same");

  _VV_var->SetFillColor(kGreen);
  _VV_var->SetLineWidth(1);
  _VV_var->SetLineStyle(1);
  _VV_var->SetMaximum(ymax);
  _VV_var->SetMinimum(ymin);
  _VV_var->Draw("HIST same");

  _QCD_var->Draw("HIST same");
  _QCD_var->SetFillColor(43);
  _QCD_var->SetLineWidth(1);
  _QCD_var->SetLineStyle(1);
  _QCD_var->SetMaximum(ymax);
  _QCD_var->SetMinimum(ymin);

  _Data_var->Draw("E same");
  _Data_var->SetMarkerColor(1);
  _Data_var->SetMarkerStyle(20);
  _Data_var->SetMarkerSize(1.1);
  _Data_var->SetLineColor(1);
  _Data_var->SetLineWidth(3);
  _Data_var->SetMaximum(ymax);
  _Data_var->SetMinimum(ymin);

   c1_1->Modified();
   c1->cd();


// ************************* Errors  ratioHis
    TH1D * bkgdErr = (TH1D*)_TTbar_var->Clone("bkgdErr");
    bkgdErr->SetFillStyle(3013);
    bkgdErr->SetFillColor(1);
    bkgdErr->SetMarkerStyle(21);
    bkgdErr->SetMarkerSize(0);

    for (int iB=1; iB<=_Data_var->GetSize(); ++iB) {
      float eStat =  bkgdErr->GetBinError(iB);
      float Err = TMath::Sqrt(eStat*eStat);
      bkgdErr->SetBinError(iB,Err);
      bkgdErr->SetBinContent(iB,_TTbar_var->GetBinContent(iB));
    }

    float chi2 = 0;
    for (int iB=1; iB<=_Data_var->GetSize(); ++iB) {
      float xData = _Data_var->GetBinContent(iB);
      float xMC = _TTbar_var->GetBinContent(iB);
      if (xMC>1e-1) {
	float diff2 = (xData-xMC)*(xData-xMC);
	chi2 += diff2/xMC;
      }
     
    }
    std::cout << std::endl;
    std::cout << "Chi2 = " << chi2 << std::endl;
    std::cout << std::endl;


    TH1D * ratioH = (TH1D*)_TTbar_var->Clone("ratioH");
    TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
    ratioH->SetMarkerColor(1);
    ratioH->SetMarkerStyle(20);
    ratioH->SetMarkerSize(1.4);
    ratioH->SetLineColor(1);
    ratioH->GetXaxis()->SetLabelFont(42);
    ratioH->GetXaxis()->SetLabelOffset(0.02);
    ratioH->GetXaxis()->SetLabelSize(0.15);
    ratioH->GetXaxis()->SetNdivisions(505);
    ratioH->GetXaxis()->SetTickLength(0.07);

    ratioH->GetYaxis()->SetRangeUser(0.2,2.0);
    ratioH->GetYaxis()->SetNdivisions(505);
    ratioH->GetYaxis()->SetTitle("Data / MC");
    ratioH->GetYaxis()->SetLabelFont(42);
    ratioH->GetYaxis()->SetLabelOffset(0.1/5);
    ratioH->GetYaxis()->SetLabelSize(0.13);
    ratioH->GetYaxis()->SetTitleOffset(0.38);
    ratioH->GetYaxis()->SetTickLength(0.02);
    ratioH->GetYaxis()->SetTitleSize(0.19);
    ratioH->GetYaxis()->SetLabelOffset(0.01);
    ratioH->SetMarkerSize(1.2);

    for (int iB=1; iB<=_Data_var->GetSize(); ++iB) {
      float x1 = _Data_var->GetBinContent(iB);
      float x2 = _TTbar_var->GetBinContent(iB);
      ratioErrH->SetBinContent(iB,1.0);
      ratioErrH->SetBinError(iB,0.0);
      float xBkg = bkgdErr->GetBinContent(iB);
      float errBkg = bkgdErr->GetBinError(iB);
      if (xBkg>0) {
	float relErr = errBkg/xBkg;
	ratioErrH->SetBinError(iB,relErr);
      }
      if (x1>0 && x2>0) {
	float e1 = _TTbar_var->GetBinError(iB);
	float ratio = x1/x2;
	float eratio = e1/x2;
	ratioH->SetBinContent(iB,ratio);
	ratioH->SetBinError(iB,eratio);
      }
      else {
	ratioH->SetBinContent(iB,1000);
      }
    }

    // ------------>Primitives in pad: c1_2
    TPad * c1_2 = new TPad("lower", "pad",0.01,0.11,0.75,0.28);
    c1_2->Draw();
    c1_2->cd();
    c1_2->SetFillColor(0);
    c1_2->SetBorderMode(0);
    c1_2->SetBorderSize(10);
    c1_2->SetGridy();
    c1_2->SetTickx(1);
    c1_2->SetTicky(1);
    c1_2->SetLeftMargin(0.15);
    c1_2->SetRightMargin(0.05);
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

    ratioH->Draw("e1");
    ratioH->SetTitle("");
    ratioErrH->Draw("e2same");

    c1_2->Modified();
    c1_2->RedrawAxis();
    c1->cd();

// ******************************* Legend
    TLegend * leg = new TLegend(0.86,0.4,0.96,0.75);
    leg = new TLegend(0.75,0.7,0.95,0.90);
    SetLegendStyle(leg);
    leg->SetTextSize(0.03);
    leg->AddEntry(_Data_var,"Data","p");
    leg->AddEntry(_TTbar_var,"TTbar","f");
    leg->AddEntry(_VJ_var,"V+Jets","f");
    leg->AddEntry(_VV_var,"VV","f");
    leg->AddEntry(_QCD_var,"QCD","f");
    leg->AddEntry(ratioErrH,"MC stat. unc.","f");
    leg->Draw("same");

    TPaveText *t1 = new TPaveText(0.16,0.92,0.3,0.97);
    t1->SetFillColor(0);
    t1->SetTextColor(1);
    t1->SetLineWidth(2);
    t1->SetTextFont(22);
    t1->SetTextSize(0.032);
    t1->AddText("                                                           CMS (2023 C-D EGamma + Muon) ~ 21.7 fb^{-1} #it{preliminary}");
    t1->Draw("same");

    TPaveText *t11 = new TPaveText(0.2,0.88,0.3,9);
    t11->SetFillColor(0);
    t11->SetTextColor(1);
    t11->SetLineWidth(2);
    t11->SetTextFont(22);
    t11->SetTextSize(0.035);
    t11->AddText("p_{T} SFs OFF");
//    t11->Draw("same");


    TImage *img = TImage::Create();
    img->FromPad(c1);
    img->WriteImage("Plots/"+variable+".png");
//    Convas_1->Close();

}