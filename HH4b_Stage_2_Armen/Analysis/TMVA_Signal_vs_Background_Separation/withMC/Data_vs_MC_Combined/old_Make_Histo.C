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

   TString Var_Name[] = {"PT1","PT2","MSD1","MSD2","ETA1","ETA2","Xbb1", "Xbb2"
                        };

   int Nbin[]         = { 100,  100,  100,  100,   100,   100,   100,   100};
   double Xmin[]      = { 200,  200,  40,  40, -3.0, -3.0, 0.0, 0.0};
   double Xmax[]      = { 1000, 1000, 400, 400,  3.0,  3.0,  1.0,  1.0};

  int nVar = sizeof(Var_Name)/sizeof(Var_Name[0]);

for(int variab = 0; variab < nVar; variab++)
 {
  TString variable = Var_Name[variab];
  int nbin = Nbin[variab];
  double xmin = Xmin[variab];
  double xmax = Xmax[variab];
  int rebin =4;

  gStyle->SetOptFile(0);
  gStyle->SetOptStat("m");
  SetStyle();

  TString hist_name[] = {variable+"_VV",variable+"_VJets",variable+"_tt",variable+"_QCD",variable+"_Signal",variable+"_DATA"};
  int Nhisto = sizeof(hist_name)/sizeof(hist_name[0]);

  TCanvas *Convas_1 = new TCanvas();
  Convas_1->SetRightMargin(0.3);
  Convas_1->SetBottomMargin(0.35);

  TFile* f1 = new TFile(variable+"_Output.root");
  TH1D* histoTot[Nhisto];
  TString histoTot_Name[Nhisto];


  for( int j=0; j<Nhisto; j++)
   {
      Float_t xVal = rebin*(xmax-xmin)/nbin;
      TString xtitle = "Events / ";
      if (variable == "PT1" || variable == "PT2" || variable == "PT3" || variable == "PT4"|| variable == "PT5"|| variable == "PTb1" || variable == "PTb2" || variable == "PTq1" || variable == "PTq2" || variable == "MINVQQ" || variable == "MINVBB" || variable =="Mass_in_BDT" || variable =="E_Rest_20" || variable =="PT_Rest_20")
        {
         xtitle+=Form("%.1f",xVal);
         xtitle+=" GeV";
        }
      else
        xtitle+=Form("%.2f",xVal);
      histoTot_Name[j]="histoTot";
      histoTot_Name[j]+=j;
      histoTot[j] = new TH1D(histoTot_Name[j],histoTot_Name[j],nbin,xmin,xmax);
      TH1D * histo = (TH1D*)f1->Get(hist_name[j]);
      histoTot[j]->GetXaxis()->SetTitle(variable);
      if(variable=="BDT")
        histoTot[j]->GetXaxis()->SetTitle("BDT-out");
      histoTot[j]->GetXaxis()->SetLabelFont(42);
      histoTot[j]->GetXaxis()->SetLabelOffset(0.15);
      histoTot[j]->GetXaxis()->SetTitleOffset(2.9);
      histoTot[j]->GetXaxis()->SetTitleFont(42);
      histoTot[j]->GetYaxis()->SetNdivisions(505);
      histoTot[j]->GetYaxis()->SetTitle(xtitle);
      histoTot[j]->GetYaxis()->SetLabelFont(42);
      histoTot[j]->GetYaxis()->SetLabelSize(0.04);
      histoTot[j]->GetYaxis()->SetLabelOffset(0.01);
      histoTot[j]->GetYaxis()->SetTitleOffset(1.3);
      histoTot[j]->GetYaxis()->SetTitleFont(42);
      histoTot[j]->SetTitle("preliminary    2018   54.5 fb^{-1}");

      if(j==0 || j>3)
       histoTot[j]->Add(histo);
      else
       histoTot[j]->Add(histoTot[j-1],histo);
   }


      for(int m=0;m<Nhisto;m++)
        histoTot[m]->Rebin(rebin);


      double ymin=10e-3;
      double ymax=20*(histoTot[5]->GetMaximum());

//QCD
      histoTot[3]->Draw("HIST same");
      histoTot[3]->SetFillColor(43);
      histoTot[3]->SetLineWidth(1);
      histoTot[3]->SetLineStyle(1);
      histoTot[3]->SetMaximum(ymax);
      histoTot[3]->SetMinimum(ymin);
//TTbar
      histoTot[2]->Draw("HIST same");
      histoTot[2]->SetFillColor(kOrange-2);
      histoTot[2]->SetLineWidth(1);
      histoTot[2]->SetLineStyle(1);
      histoTot[2]->SetMaximum(ymax);
      histoTot[2]->SetMinimum(ymin);
//V-Jets
      histoTot[1]->Draw("HIST same");
      histoTot[1]->SetFillColor(kGreen-2);
      histoTot[1]->SetLineWidth(1);
      histoTot[1]->SetLineStyle(1);
      histoTot[1]->SetMaximum(ymax);
      histoTot[1]->SetMinimum(ymin);
//VV
      histoTot[0]->Draw("HIST same");
      histoTot[0]->SetFillColor(kGreen-7);
      histoTot[0]->SetLineWidth(1);
      histoTot[0]->SetLineStyle(1);
      histoTot[0]->SetMaximum(ymax);
      histoTot[0]->SetMinimum(ymin);
//Signal
      histoTot[4]->SetLineColor(kRed);
      histoTot[4]->SetLineWidth(3);
      histoTot[4]->Draw("HIST same");
      histoTot[4]->SetMaximum(ymax);
      histoTot[4]->SetMinimum(ymin);

//Data
      histoTot[5]->Draw("E same");
      histoTot[5]->SetMarkerColor(1);
      histoTot[5]->SetMarkerStyle(20);
      histoTot[5]->SetMarkerSize(1.1);
      histoTot[5]->SetLineColor(1);
      histoTot[5]->SetLineWidth(3);
      histoTot[5]->SetMaximum(ymax);
      histoTot[5]->SetMinimum(ymin);

     Convas_1->SetLogy();

// ************************* Errors  ratioHis
    TH1D * bkgdErr = (TH1D*)histoTot[3]->Clone("bkgdErr");
    bkgdErr->SetFillStyle(3013);
    bkgdErr->SetFillColor(1);
    bkgdErr->SetMarkerStyle(21);
    bkgdErr->SetMarkerSize(0);

    for (int iB=1; iB<=histoTot[5]->GetSize(); ++iB) {
      float eStat =  bkgdErr->GetBinError(iB);
      float Err = TMath::Sqrt(eStat*eStat);
      bkgdErr->SetBinError(iB,Err);
      bkgdErr->SetBinContent(iB,histoTot[3]->GetBinContent(iB));
    }

    float chi2 = 0;
    for (int iB=1; iB<=histoTot[5]->GetSize(); ++iB) {
      float xData = histoTot[5]->GetBinContent(iB);
      float xMC = histoTot[3]->GetBinContent(iB);
      if (xMC>1e-1) {
	float diff2 = (xData-xMC)*(xData-xMC);
	chi2 += diff2/xMC;
      }
     
    }
    std::cout << std::endl;
    std::cout << "Chi2 = " << chi2 << std::endl;
    std::cout << std::endl;


    TH1D * ratioH = (TH1D*)histoTot[5]->Clone("ratioH");
    TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
    ratioH->SetMarkerColor(1);
    ratioH->SetMarkerStyle(20);
    ratioH->SetMarkerSize(1.4);
    ratioH->SetLineColor(1);
    ratioH->GetXaxis()->SetLabelFont(42);
    ratioH->GetXaxis()->SetLabelOffset(0.02);
    ratioH->GetXaxis()->SetLabelSize(0.15);
    ratioH->GetXaxis()->SetNdivisions(505);
//    ratioH->GetXaxis()->SetTitleSize(0.12);
//    ratioH->GetXaxis()->SetTitleOffset(0.35);
    ratioH->GetXaxis()->SetTitleSize(0.12);
    ratioH->GetXaxis()->SetTickLength(0.07);

//    ratioH->GetYaxis()->SetRangeUser(0.2,1.4);               !!!!!!!!!!!!!!! changed to next line
    ratioH->GetYaxis()->SetRangeUser(0.4,2.0);
    ratioH->GetYaxis()->SetNdivisions(505);
    ratioH->GetYaxis()->SetTitle("Data / MC");
    ratioH->GetYaxis()->SetLabelFont(42);
    ratioH->GetYaxis()->SetLabelOffset(0.1/5);
    ratioH->GetYaxis()->SetLabelSize(0.13);
    ratioH->GetYaxis()->SetTitleOffset(0.35);
    ratioH->GetYaxis()->SetTickLength(0.02);
    ratioH->GetYaxis()->SetTitleSize(0.19);
    ratioH->GetYaxis()->SetLabelOffset(0.01);
    ratioH->SetMarkerSize(1.2);

    for (int iB=1; iB<=histoTot[5]->GetSize(); ++iB) {
      float x1 = histoTot[5]->GetBinContent(iB);
      float x2 = histoTot[3]->GetBinContent(iB);
      ratioErrH->SetBinContent(iB,1.0);
      ratioErrH->SetBinError(iB,0.0);
      float xBkg = bkgdErr->GetBinContent(iB);
      float errBkg = bkgdErr->GetBinError(iB);
      if (xBkg>0) {
	float relErr = errBkg/xBkg;
	ratioErrH->SetBinError(iB,relErr);
      }
      if (x1>0 && x2>0) {
	float e1 = histoTot[5]->GetBinError(iB);
	float ratio = x1/x2;
	float eratio = e1/x2;
	ratioH->SetBinContent(iB,ratio);
	ratioH->SetBinError(iB,eratio);
      }
      else {
	ratioH->SetBinContent(iB,1000);
      }
    }

    // ------------>Primitives in pad: lower
    TPad * lower = new TPad("lower", "pad",0,0.12,1,0.33);
    lower->Draw();
    lower->cd();
    lower->SetFillColor(0);
    lower->SetBorderMode(0);
    lower->SetBorderSize(10);
    lower->SetGridy();
    lower->SetTickx(1);
    lower->SetTicky(1);
    lower->SetLeftMargin(0.15);
    lower->SetRightMargin(0.3);
    lower->SetTopMargin(0.1);
    lower->SetBottomMargin(0.15);
    lower->SetFrameFillStyle(0);
    lower->SetFrameLineStyle(0);
    lower->SetFrameLineWidth(2);
    lower->SetFrameBorderMode(0);
    lower->SetFrameBorderSize(10);
    lower->SetFrameFillStyle(0);
    lower->SetFrameLineStyle(0);
    lower->SetFrameLineWidth(2);
    lower->SetFrameBorderMode(0);
    lower->SetFrameBorderSize(10);

    ratioH->Draw("e1");
    ratioH->SetTitle("");
    ratioErrH->Draw("e2same");

    lower->Modified();
    lower->RedrawAxis();
    Convas_1->cd();
    Convas_1->Modified();
    Convas_1->cd();
    Convas_1->SetSelected(Convas_1);

// ******************************* Legend
    TLegend * leg = new TLegend(0.65,0.45,0.85,0.8);
    leg = new TLegend(0.75,0.5,0.9,0.9);
    SetLegendStyle(leg);
    leg->SetTextSize(0.03);
    leg->AddEntry(histoTot[5],"Data","p");
    leg->AddEntry(histoTot_Name[4],"ggF HHto4b, m(H)=125","l");
    leg->AddEntry(histoTot_Name[3],"QCD","f");
    leg->AddEntry(histoTot_Name[2],"t#bar{t}","f");
    leg->AddEntry(histoTot_Name[1],"V + jets","f");
    leg->AddEntry(histoTot_Name[0],"VV","f");
    leg->AddEntry(ratioErrH,"MC stat. unc.","f");
    leg->Draw("same");

    TImage *img = TImage::Create();
    img->FromPad(Convas_1);
    img->WriteImage(variable+".png");
    Convas_1->Close();

 } // End of variab loop
}