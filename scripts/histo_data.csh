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

# 2022: 2022C, 2022D
# 2022EE: 2022E, 2022F, 2022G
# 2023: 2023C
# 2023BPix: 2023D

era_tag="2023"
ttbar_type="TTtoLNu2Q"
channels=("JetMET" "Muon" "EGamma")
run_tag="2023C"

sample_dir="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"

for channel in "${channels[@]}"; do
    echo "channel: ${channel}"
    versions=("v1" "v2" "v3")

    for version in "${versions[@]}"; do
        echo "version: ${version}"
        sample_path=${sample_dir}/${era_tag}/${channel}_${run_tag}_${version}.root
        output_path=${TMP_DIR}/Histograms_Data_${channel}_${run_tag}_${version}.root
        jec_path=${PROJ_ROOT}/JECs
        jec_path_L2Relative=${jec_path}/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2Relative_AK8PFPuppi.txt;
        jec_path_L2L3Residual=${jec_path}/Summer23Prompt23_RunCv123_V1_DATA/Summer23Prompt23_RunCv123_V1_DATA_L2L3Residual_AK8PFPuppi.txt;

        root -l -b -q "histo_data.cpp(\"${run_tag}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${jec_path_L2Relative}\", \"${jec_path_L2L3Residual}\")"
    done

    version="v4"
    echo "version: ${version}"
    sample_path=${sample_dir}/${era_tag}/${channel}_${run_tag}_${version}.root
    output_path=${TMP_DIR}/Histograms_Data_${channel}_${run_tag}_${version}.root
    jec_path=${PROJ_ROOT}/JECs
    jec_path_L2Relative=${jec_path}/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2Relative_AK8PFPuppi.txt;
    jec_path_L2L3Residual=${jec_path}/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_L2L3Residual_AK8PFPuppi.txt;

    root -l -b -q "histo_data.cpp(\"${run_tag}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${jec_path_L2Relative}\", \"${jec_path_L2L3Residual}\")"

    # Hadd
    target_path=${TMP_DIR}/Histograms_Data_${era_tag}_${channel}.root
    echo "Hadd into one file: ${target_path}"
    hadd -f ${target_path} ${TMP_DIR}/Histograms_Data_${channel}_${run_tag}_v*.root
done

# Hadd into leptonic channel
target_path=${OUTPUT_DIR}/Histograms_Data_${era_tag}_Leptonic.root
echo "Hadd into one file: ${target_path}"
hadd -f ${target_path} ${TMP_DIR}/Histograms_Data_${era_tag}_Muon.root ${TMP_DIR}/Histograms_Data_${era_tag}_EGamma.root

# Hadd into one file
target_path=${OUTPUT_DIR}/Histograms_Data_${era_tag}.root
echo "Hadd into one file: ${target_path}"
hadd -f ${target_path} ${TMP_DIR}/Histograms_Data_${era_tag}_Leptonic.root ${TMP_DIR}/Histograms_Data_${era_tag}_JetMET.root