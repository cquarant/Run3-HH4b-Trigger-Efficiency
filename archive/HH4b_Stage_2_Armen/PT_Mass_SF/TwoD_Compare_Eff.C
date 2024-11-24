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

  /* Setups */
  TFile *f = new TFile("PT_Mass_2dSF.root","RECREATE");
  Float_t Lower_m[16]   = {0,5,10,20,30,40,50,60,80,100,120,150,200,250,300,350};
  Float_t Lower_pt[46]  = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,320,340,360,380,400,420,440,460,480,500,550,600,700,800,1000};
  
  // Histograms for Data and MC efficiencies, Bined in (pT, Mass)
  // const char *name, const char *title, Int_t nbinsx, const Double_t *xbins, Int_t nbinsy, const Double_t *ybins
  TH2D *_Eff_Data_ETA0 = new TH2D("Eff_Data_ETA0","Eff_Data_ETA0", 45, Lower_pt, 15, Lower_m);
  TH2D *_Eff_MC_ETA0 = new TH2D("Eff_MC_ETA0","Eff_MC_ETA0", 45, Lower_pt, 15, Lower_m);

  // Rebinning
  int rebPT = 3;  // every 3 bins are merged => pT dimension is reduced by factor of 3
  int rebM = 1;   // no merging
  _Eff_Data_ETA0->Rebin2D(rebPT, rebM);
  _Eff_MC_ETA0->Rebin2D(rebPT, rebM);

  /* Styling */
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  SetStyle();
  gStyle->SetPaintTextFormat("1.2f");

  TString PlotTitle="HLT_AK8PFJet230_SoftDropMass40     (|#eta| < 2.5)";

  TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle,100,52,900,686);
  // c1->SetGridy();
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

  /* File handling */
  // TString erA="PreBPix";
  TString erA="PostBPix";
  TFile* fMC = new TFile("MC/QCD_"+erA+".root");
  TFile* fData = new TFile("Data/Data_"+erA+".root");

  /* Get MC and data efficiencies */
  // TH2D * _MC_Probe  = (TH2D*)fMC->Get("FatJet2_Pt_MassN");
  // TH2D * _MC_ProbeM = (TH2D*)fMC->Get("FatJet2_Pt_MassN_M");
  // TH2D * _Data_Probe  = (TH2D*)fData->Get("FatJet2_Pt_MassN");
  // TH2D * _Data_ProbeM = (TH2D*)fData->Get("FatJet2_Pt_MassN_M");
  
  TH2D * _MC_Probe  = (TH2D*)fMC->Get("FatJet2_Pt_Mass");    // all signal events
  TH2D * _MC_ProbeM = (TH2D*)fMC->Get("FatJet2_Pt_Mass_M");  // selected signal events
  _MC_Probe->Rebin2D(rebPT, rebM);
  _MC_ProbeM->Rebin2D(rebPT, rebM);
  _MC_ProbeM->Sumw2();  // sum of squares of weights (to propagate errors correctly)
  _MC_ProbeM->Divide(_MC_Probe);  // MC efficiency
  
  TH2D * _Data_Probe  = (TH2D*)fData->Get("FatJet2_Pt_Mass");   // all signal events
  TH2D * _Data_ProbeM = (TH2D*)fData->Get("FatJet2_Pt_Mass_M"); // selected signal events
  _Data_Probe ->Rebin2D(rebPT, rebM);
  _Data_ProbeM->Rebin2D(rebPT, rebM);
  _Data_ProbeM->Sumw2();  // sum of squares of weights (to propagate errors correctly)
  _Data_ProbeM->Divide(_Data_Probe);  // Data efficiency

  /* Calculate SF */
  // _Data_ProbeM=_MC_ProbeM;     // For efficiency
  // _MC_ProbeM=_Data_ProbeM;     // For efficiency

  _Data_ProbeM->SetMarkerStyle(23);
  // _Data_ProbeM->SetLineColor(2);
  // _Data_ProbeM->SetMarkerColor(2);
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

  _Data_ProbeM->GetZaxis()->SetTitle("Efficiency ");
  // _Data_ProbeM->GetZaxis()->SetTitle("Data / MC ");
  // _Data_ProbeM->SetTitle("2.3 < #eta < 4.7");

  _Data_ProbeM->Divide(_MC_ProbeM); 

  _Data_ProbeM->SetMaximum(1.8);
  // _Data_ProbeM->SetMaximum(1.0);
  _Data_ProbeM->SetMinimum(0.0);
  // _Data_ProbeM->Draw("colz");
  _Data_ProbeM->Draw("colz texte");

   c1->Modified();
   c1->cd();
   c1->SaveAs("Mass_PT/Eff_PT_MassSD_"+erA+".png");

  // Save the histograms
  for (int i = 0; i <=_Data_ProbeM->GetNbinsX(); i++) {
    for (int j = 0; j<=_Data_ProbeM->GetNbinsY(); j++) {
      _Eff_Data_ETA0->SetBinContent(i, j,_Data_ProbeM->GetBinContent(i, j));
      _Eff_Data_ETA0->SetBinError(i, j,_Data_ProbeM->GetBinError(i, j));
      _Eff_MC_ETA0->SetBinContent(i, j,_MC_ProbeM->GetBinContent(i, j));
      _Eff_MC_ETA0->SetBinError(i, j,_MC_ProbeM->GetBinError(i, j));
    }
  }
  f->Write();

}