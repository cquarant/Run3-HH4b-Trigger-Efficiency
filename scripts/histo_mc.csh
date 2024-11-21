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


ERA_TAG="2023"
# TTBar
ttbar_type="TTtoLNu2Q"
channels=("Muon" "EGamma" "Leptonic")

for channel in ${channels[@]}; do
    sample_dir="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"
    sample_path=${sample_dir}/${ERA_TAG}/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8.root

    output_path=${TMP_DIR}/Histograms_MC_TTtoLNu2Q_${ERA_TAG}_${channel}.root

    pu_path=${PROJ_ROOT}/pileups/pu_${ERA_TAG}.txt
    param_path=${PROJ_ROOT}/parameters/parameters_${ERA_TAG}.txt
    jec_path=${PROJ_ROOT}/JECs
    jec_path_ak4=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi.txt
    jec_path_ak8=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK8PFPuppi.txt

    sf_path=${PROJ_ROOT}/Trigger_SFs/PT_Mass_SF_TTbar/PT_Mass_2dSF_2023.root

    root -l -b -q "histo_mc.cpp(\"${ttbar_type}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${pu_path}\", \"${sf_path}\", \"${param_path}\", \"${jec_path_ak4}\", \"${jec_path_ak8}\")"
done
cp ${TMP_DIR}/Histograms_MC_TTtoLNu2Q_${ERA_TAG}_Leptonic.root ${OUTPUT_DIR}/Histograms_MC_TTtoLNu2Q_${ERA_TAG}_Leptonic.root

# QCD
channel="QCD"
HT_BINS=("100to200" "200to400" "400to600" "600to800" "600to1000" "1000to1200" "1200to1500" "1500to2000" "2000toInf")
for ht_bin in ${HT_BINS[@]}; do
    sample_dir="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"
    sample_path=${sample_dir}/${ERA_TAG}/QCD-4Jets_HT-${ht_bin}.root

    output_path=${TMP_DIR}/Histograms_MC_QCD-4Jets_HT-${ht_bin}_${ERA_TAG}.root

    pu_path=${PROJ_ROOT}/pileups/pu_${ERA_TAG}.txt
    param_path=${PROJ_ROOT}/parameters/parameters_${ERA_TAG}.txt
    jec_path=${PROJ_ROOT}/JECs
    jec_path_ak4=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi.txt
    jec_path_ak8=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK8PFPuppi.txt

    sf_path=${PROJ_ROOT}/Trigger_SFs/PT_Mass_SF_TTbar/PT_Mass_2dSF_2023.root

    qcd_type="QCD_HT_${ht_bin}"

    root -l -b -q "histo_mc.cpp(\"${qcd_type}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${pu_path}\", \"${sf_path}\", \"${param_path}\", \"${jec_path_ak4}\", \"${jec_path_ak8}\")"
done

# Hadd
qcd_path=${OUTPUT_DIR}/Histograms_MC_QCD_${ERA_TAG}.root
echo "Hadd into one file: ${target_path}"
hadd -f ${qcd_path} ${TMP_DIR}/Histograms_MC_QCD-4Jets_HT-*_${ERA_TAG}.root

# Hadd final
output_path=${OUTPUT_DIR}/Histograms_MC_${ERA_TAG}.root
ttbar_path=${OUTPUT_DIR}/Histograms_MC_TTtoLNu2Q_${ERA_TAG}_Leptonic.root
echo "Hadd into one file: ${output_path}"
hadd -f ${output_path} ${ttbar_path} ${qcd_path}