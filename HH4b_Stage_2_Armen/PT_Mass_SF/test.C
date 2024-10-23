#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void test()
{

  TFile* fEff = new TFile("PT_Mass_2dSF.root");
  TH2D * _Eff  = (TH2D*)fEff->Get("Eff_Data_ETA0");
   
  double FatJet2_pt = 320;
  double FatJet2_MassSD = 73;
  Int_t bin_PT_2     = _Eff->GetXaxis()->FindBin(FatJet2_pt);
  Int_t bin_Mass_2   = _Eff->GetYaxis()->FindBin(FatJet2_MassSD);
  cout   <<_Eff->GetBinContent(bin_PT_2,bin_Mass_2)<<endl;


   for(int i=_Eff->GetNbinsY();i>0;i--)
    {
     for(int j=1;j<=_Eff->GetNbinsX();j++)
      cout<<_Eff->GetBinContent(j,i)<<"   ";
     cout<<endl;
    }

}