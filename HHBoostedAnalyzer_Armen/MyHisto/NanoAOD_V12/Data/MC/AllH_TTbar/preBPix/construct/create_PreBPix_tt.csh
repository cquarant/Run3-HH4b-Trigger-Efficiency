#!/bin/tcsh
set input = $1
cat includes_tt.txt  >  Making_Tree_$input.C
echo 'void Making_Tree_'$input'()' >> Making_Tree_$input.C
echo '{' >> Making_Tree_$input.C
cat JEC_JER_PreBPix.txt >> Making_Tree_$input.C
echo ' TFile *f = new TFile("Histograms_'$input'.root","RECREATE");' >> Making_Tree_$input.C
echo '' >> Making_Tree_$input.C
cat histos.txt    >> Making_Tree_$input.C
echo 'TFile *f1 = new TFile("/eos/home-t/tumasyan/HHTo4B/Data_2023/Legacy/PreBPix/'$input'.root");' >> Making_Tree_$input.C
echo '' >> Making_Tree_$input.C
cat body.txt      >> Making_Tree_$input.C
echo '   weight = (weight/SumGenWeights)*XSec_'$input'*Lumi;' >> Making_Tree_$input.C
cat fill.txt      >> Making_Tree_$input.C
