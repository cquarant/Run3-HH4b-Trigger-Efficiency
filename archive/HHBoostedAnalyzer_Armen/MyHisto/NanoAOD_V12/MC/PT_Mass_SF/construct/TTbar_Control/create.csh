#!/bin/tcsh
set input = $1
cat includes.txt  >  Making_Histo_$input.C
echo 'void Making_Histo_'$input'()' >> Making_Histo_$input.C
echo '{' >> Making_Histo_$input.C
echo ' TFile *f = new TFile("Histograms_'$input'.root","RECREATE");' >> Making_Histo_$input.C
echo '' >> Making_Histo_$input.C
cat histos.txt    >> Making_Histo_$input.C
#echo 'TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/MC/QCD/postEE/'$input'.root");' >> Making_Histo_$input.C
#echo 'TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/MC/QCD/postEE/TTbar_Control/'$input'.root");' >> Making_Histo_$input.C
echo 'TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/MC/postEE_SFs/'$input'.root");' >> Making_Histo_$input.C
echo '' >> Making_Histo_$input.C
cat body.txt      >> Making_Histo_$input.C
echo '   weight = (weight/SumGenWeights)*XSec_'$input'*Lumi;' >> Making_Histo_$input.C
cat fill.txt      >> Making_Histo_$input.C
