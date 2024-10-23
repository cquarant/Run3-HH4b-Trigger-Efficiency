cp ~/public/x509up_u28489 /tmp/

input="/afs/cern.ch/work/t/tumasyan/HHTo4B/2023/CMSSW_13_1_0/src/HHBoostedAnalyzer/MyJOBs_2023/NanoAOD_V12/MC/datasets.txt"
while IFS= read -r line
do

  echo "$line"
  for i in $(dasgoclient --query="file dataset=$line"); do 
  echo "'"$i"'" >>RootFile;
  done

done < "$input"



