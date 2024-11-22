set -e;
# check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}/output
TMP_DIR=${SCRIPT_DIR}/tmp
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TMP_DIR}


era_tag="2023"
# TTBar
ttbar_type="TTtoLNu2Q"
channels=("Muon" "EGamma" "Leptonic")

for channel in ${channels[@]}; do
    sample_dir="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"
    sample_path=${sample_dir}/${era_tag}/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8.root

    output_path=${TMP_DIR}/Histograms_MC_TTtoLNu2Q_${era_tag}_${channel}.root

    pu_path=${PROJ_ROOT}/pileups/pu_${era_tag}.txt
    param_path=${PROJ_ROOT}/parameters/parameters_${era_tag}.txt
    jec_path=${PROJ_ROOT}/JECs
    jec_path_ak4=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi.txt
    jec_path_ak8=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK8PFPuppi.txt

    sf_path=${PROJ_ROOT}/Trigger_SFs/TTbar/PT_Mass_2dSF_2023.root

    root -l -b -q "histo_mc.cpp(\"${ttbar_type}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${pu_path}\", \"${sf_path}\", \"${param_path}\", \"${jec_path_ak4}\", \"${jec_path_ak8}\")"
done


# QCD
channel="QCD"
HT_BINS=("100to200" "200to400" "400to600" "600to800" "800to1000" "1000to1200" "1200to1500" "1500to2000" "2000toInf")
for ht_bin in ${HT_BINS[@]}; do
    sample_dir="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"
    if [ "${ht_bin}" == "2000toInf" ]; then
        sample_path=${sample_dir}/${era_tag}/QCD-4Jets_HT-2000.root
    else
        sample_path=${sample_dir}/${era_tag}/QCD-4Jets_HT-${ht_bin}.root
    fi

    output_path=${TMP_DIR}/Histograms_MC_QCD-4Jets_HT-${ht_bin}_${era_tag}.root

    pu_path=${PROJ_ROOT}/pileups/pu_${era_tag}.txt
    param_path=${PROJ_ROOT}/parameters/parameters_${era_tag}.txt
    jec_path=${PROJ_ROOT}/JECs
    jec_path_ak4=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi.txt
    jec_path_ak8=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK8PFPuppi.txt

    sf_path=${PROJ_ROOT}/Trigger_SFs/QCD/PT_Mass_2dSF_2023.root

    qcd_type="QCD_HT_${ht_bin}"

    root -l -b -q "histo_mc.cpp(\"${qcd_type}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${pu_path}\", \"${sf_path}\", \"${param_path}\", \"${jec_path_ak4}\", \"${jec_path_ak8}\")"
done

# Top
cp ${TMP_DIR}/Histograms_MC_TTtoLNu2Q_${era_tag}_Leptonic.root ${OUTPUT_DIR}/Histograms_${era_tag}_MC_top.root

# QCD
qcd_path=${OUTPUT_DIR}/Histograms_${era_tag}_MC_QCD.root
echo "Hadd into one file: ${target_path}"
hadd -f ${qcd_path} ${TMP_DIR}/Histograms_MC_QCD-4Jets_HT-*_${era_tag}.root

# Combined
output_path=${OUTPUT_DIR}/Histograms_${era_tag}_MC.root
ttbar_path=${OUTPUT_DIR}/Histograms_MC_TTtoLNu2Q_${era_tag}.root
echo "Hadd into one file: ${output_path}"
hadd -f ${output_path} ${ttbar_path} ${qcd_path}