#!/bin/tcsh
set input = $1
cat includes.txt  >  Making_Histo_$input.C
echo 'void Making_Histo_'$input'()' >> Making_Histo_$input.C
echo '{' >> Making_Histo_$input.C
cat JEC_JER_PostBPix.txt >> Making_Histo_$input.C
echo ' TFile *f = new TFile("Histograms_'$input'.root","RECREATE");' >> Making_Histo_$input.C
echo '' >> Making_Histo_$input.C
cat histos.txt    >> Making_Histo_$input.C
echo 'TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/NewMC/PostBPix/'$input'.root");' >> Making_Histo_$input.C
echo '' >> Making_Histo_$input.C
cat body.txt      >> Making_Histo_$input.C
echo '   weight = (weight/SumGenWeights)*XSec_'$input'*Lumi;' >> Making_Histo_$input.C
cat fill.txt      >> Making_Histo_$input.C
