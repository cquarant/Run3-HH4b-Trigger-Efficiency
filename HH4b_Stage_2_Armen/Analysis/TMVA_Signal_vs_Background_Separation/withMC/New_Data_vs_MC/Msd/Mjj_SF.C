#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void Mjj_SF()
{

  int reb =1;

// **********************************************************************
  TFile *f = new TFile("Mjj_SF.root","RECREATE");
  Float_t Bound[9] = {0,35,75,130,200,315,450,700,1000};
  TH1D *_SF   = new TH1D("SF","SF",8,Bound);
  _SF->Rebin(reb);
// **********************************************************************

//  Float_t Bound[8] = {0,40,100,170,300,450,700,1000};
//  TH1D *_My_MC       = new TH1D("My_MC","My_MC",7,Bound);
//  TH1D *_My_Data     = new TH1D("My_Data","My_Data",7,Bound);
  TH1D *_My_MC       = new TH1D("My_MC","My_MC",8,Bound);
  TH1D *_My_Data     = new TH1D("My_Data","My_Data",8,Bound);


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

  TString Tau3OverTau2_CUT = "";
//  Tau3OverTau2_CUT = "Cut_0p5/";
//  Tau3OverTau2_CUT = "Cut_0p46/";

  TString JET = "FatJetHBTG";
//  TString JET = "HBTG";
//  TString VAR = "_pt";
//  TString VAR = "_eta";
  TString VAR = "_MassSD";
//  TString VAR = "_Tau3OverTau2";
//  TString VAR = "PNet_QCD";
//  TString VAR = "PNetHbbvsQCD";
//  TString VAR = "PNetMD_Xbb";
//  TString VAR = "PNetMD_Xbb_Modified_as_in_AN_Xbb";

  TString variable = JET+VAR;
//  variable="MET";
//  variable="Dr_LFJ";
//  variable="Lep1_Pt";
  variable = "PTjj";
  TString XTitle = "P_{T}^{JJ}";

  TFile* f_Data   = new TFile("../Data/"+Tau3OverTau2_CUT+"Histograms_Data.root");

  TFile* f_QCD    = new TFile("../MC/"+Tau3OverTau2_CUT+"Histograms_QCD.root");
//  TFile* f_TTbar  = new TFile("../MC/"+Tau3OverTau2_CUT+"Histograms_TTbar_T3T2_Xbb_Corr.root");
  TFile* f_TTbar  = new TFile("../MC/"+Tau3OverTau2_CUT+"Histograms_TTbar.root");
  TFile* f_VV     = new TFile("../MC/"+Tau3OverTau2_CUT+"Histograms_VV.root");
  TFile* f_VJ     = new TFile("../MC/"+Tau3OverTau2_CUT+"Histograms_VJ.root");

  TH1D * _Data_var   = (TH1D*)f_Data->Get(variable);

  TH1D * _QCD_var    = (TH1D*)f_QCD->Get(variable);
  TH1D * _TTbar_var  = (TH1D*)f_TTbar->Get(variable);
  TH1D * _VV_var     = (TH1D*)f_VV->Get(variable);
  TH1D * _VJ_var     = (TH1D*)f_VJ->Get(variable);


    for (int iB=1; iB<=_Data_var->GetSize(); ++iB) {
//       if(_VJ_var->GetBinContent(iB) > 0.8 * _TTbar_var->GetBinContent(iB)) 
//         _VJ_var->SetBinContent(iB,_VJ_var->GetBinContent(iB)/15);

      _Data_var->SetBinContent(iB,_Data_var->GetBinContent(iB) - _VJ_var->GetBinContent(iB) - _VV_var->GetBinContent(iB) - _QCD_var->GetBinContent(iB));
    }

/*
   for(int i=0; i<_TTbar_var->GetSize();i++)
    for(int ib=1; ib<=12; ib++)
      {
       if((i/100.0 >= Bound[ib-1]) && (i/100.0 < Bound[ib]))
        {
          Int_t bin  = _My_MC->GetXaxis()->FindBin(i/100);
         _My_MC->SetBinContent(ib,_My_MC->GetBinContent(ib) + _TTbar_var->GetBinContent(i));
         _My_Data->SetBinContent(ib,_My_Data->GetBinContent(ib)+_Data_var->GetBinContent(i));
        }
      }
*/

   for(int i=1; i<=_My_MC->GetSize();i++)
     {
      _My_MC->SetBinContent(i,0.0);
      _My_Data->SetBinContent(i,0.0);
     }

   for(int i=1; i<=_TTbar_var->GetSize()-2; i++)
     {
        Int_t bin  = _My_MC->GetXaxis()->FindBin(i*10 - 5);
        double current_val_MC = _My_MC->GetBinContent(bin);
        double current_val_Data = _My_Data->GetBinContent(bin);
        current_val_MC = current_val_MC + _TTbar_var->GetBinContent(i);
        current_val_Data = current_val_Data + _Data_var->GetBinContent(i);
       _My_MC->SetBinContent(bin, current_val_MC);
       _My_Data->SetBinContent(bin, current_val_Data);
      }

//   _My_MC=_TTbar_var;
 //  _My_Data=_Data_var;


  _My_MC->Scale(1.0/_My_MC->Integral());
  _My_Data->Scale(1.0/_My_Data->Integral());

  _My_MC->Rebin(reb);
  _My_Data->Rebin(reb);


// ********************************************************************* Data _vs MC
   _TTbar_var= _My_MC;
   _Data_var= _My_Data;

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
//  _TTbar_var->GetYaxis()->SetTitle("Events");
    _TTbar_var->GetYaxis()->SetTitle("Fraction of events");
  _TTbar_var->GetYaxis()->SetTitleOffset(1.5);
  _TTbar_var->GetYaxis()->SetTitleSize(0.05);
  _TTbar_var->GetYaxis()->SetTitleFont(42);
  _TTbar_var->SetTitle("");
//  _TTbar_var->Scale(Offline_BTG_SF);

  double ymin=10e-4;
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

//  _My_MC->Draw("HIST same");
//  _My_MC->SetLineColor(4);
//  _My_MC->SetLineWidth(4);

//  _My_Data->Draw("HIST same");
//  _My_Data->SetLineColor(3);
//  _My_Data->SetLineWidth(4);

    TLegend * leg = new TLegend(0.7,0.5,0.8,0.75);
    leg = new TLegend(0.58,0.8,0.78,0.88);
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
   _SF_pt->SetMaximum(2.0);
   _SF_pt->SetMinimum(0.2);




    TPad * c1_2 = new TPad("lower", "pad",0.1,0.16,0.9,0.35);
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
/*
    TF1 *Mypol1   = new TF1("Mypol1","pol1",0.0,301);
    _SF_pt->Fit("Mypol1","R");
   TH1D * hint1   = new TH1D("hint1","",185,0.0,300);
   TH1D * hint1_68 = new TH1D("hint1_68","",370,0.0,300);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1_68,0.68);
    hint1   ->SetFillColor(kCyan);
    hint1_68->SetFillColor(kYellow);
  hint1_68 ->Draw("e3same");

    TF1 *Mypol2   = new TF1("Mypol2","pol1",300,1000);
    _SF_pt->Fit("Mypol2","R+");
   TH1D * hint2   = new TH1D("hint2","",185,300,1000);
   TH1D * hint2_68 = new TH1D("hint2_68","",370,300,1000);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2_68,0.68);
    hint2   ->SetFillColor(kCyan);
    hint2_68->SetFillColor(kYellow-3);
    hint2_68 ->Draw("e3same");
*/

    TF1 *Mypol   = new TF1("Mypol","pol3",0.0,1000);
    _SF_pt->Fit("Mypol","R");
   TH1D * hint   = new TH1D("hint","",185,0.0,1000);
   TH1D * hint_68 = new TH1D("hint_68","",370,0.0,1000);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68,0.68);
    hint   ->SetFillColor(kCyan);
    hint_68->SetFillColor(kYellow);
    hint_68 ->Draw("e3same");
  _SF_pt  ->Draw("e1same");


    TLegend * leg1 = new TLegend(0.75,0.6,0.95,0.85);
    SetLegendStyle(leg1);
    leg1->SetTextSize(0.14);
    leg1->AddEntry(Mypol,"pol3","l");
    leg1->AddEntry(hint_68,"1 #sigma","f");
//    leg1->AddEntry(hint1_68,"1 #sigma","f");
//    leg1->AddEntry(hint2_68,"1 #sigma","f");
//    leg1->AddEntry(hint,"2#sigma","f");
    leg1->Draw("same");


    c1_2->Modified();
    c1_2->RedrawAxis();
    c1->cd();

//*********************************************************************
// Eff Histrograms to root

   for(int i=1;i<=_SF_pt->GetNbinsX();i++)
     _SF->SetBinContent(i,_SF_pt->GetBinContent(i));

  f->Write();
//*********************************************************************

}