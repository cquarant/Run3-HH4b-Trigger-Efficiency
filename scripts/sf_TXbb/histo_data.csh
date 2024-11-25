#!/bin/bash -e

# Check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

# Directory setup
PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}/output
TMP_DIR=${SCRIPT_DIR}/tmp
mkdir -p ${OUTPUT_DIR} ${TMP_DIR}

# Configuration
declare -A era_runs=(
    ["2022"]="2022C 2022D"
    ["2022EE"]="2022E 2022F 2022G"
    ["2023"]="2023C"
    ["2023BPix"]="2023D"
)

# JEC configurations
declare -A jec_configs=(
    ["2022"]="Summer22_22Sep2023_RunCD_V2_DATA"
    ["2022EE_RunE"]="Summer22EE_22Sep2023_RunE_V2_DATA"
    ["2022EE_RunF"]="Summer22EE_22Sep2023_RunF_V2_DATA"
    ["2022EE_RunG"]="Summer22EE_22Sep2023_RunG_V2_DATA"
    ["2023_v123"]="Summer23Prompt23_RunCv123_V1_DATA"
    ["2023_v4"]="Summer23Prompt23_RunCv4_V1_DATA"
    ["2023BPix"]="Summer23BPixPrompt23_RunD_V1_DATA"
)

channels=("JetMET" "Muon" "EGamma")
sample_dir="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"
jec_base="${PROJ_ROOT}/JECs"

process_single_file() {
    local era=$1
    local channel=$2
    local version=$3
    local run_tag=$4
    local jec_config=$5

    if [ "$era" == "2023BPix" ] || [ "$era" == "2022EE" ]; then
        # 2023BPix was in 2023, 2022EE was in 2022
        local year=${era:0:4}
        local sample_path="${sample_dir}/${year}/${channel}_${run_tag}${version}.root"
    else
        local sample_path="${sample_dir}/${era}/${channel}_${run_tag}.root"
    fi
    
    local output_path="${TMP_DIR}/Histograms_${era}_data_${channel}_${run_tag}${version}.root"
    
    # Set appropriate JEC paths
    local jec_dir="${jec_base}/${jec_config}"
    local jec_path_L2Relative="${jec_dir}/${jec_config}_L2Relative_AK8PFPuppi.txt"
    local jec_path_L2L3Residual="${jec_dir}/${jec_config}_L2L3Residual_AK8PFPuppi.txt"

    echo "Processing ${era} ${channel} ${run_tag}${version}"
    root -l -b -q "histo_data.cpp(\"${run_tag}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${jec_path_L2Relative}\", \"${jec_path_L2L3Residual}\")"
}

process_era() {
    local era=$1
    echo "Processing era: ${era}"
    
    case ${era} in
        "2023")
            for run_tag in ${era_runs[${era}]}; do
                for channel in "${channels[@]}"; do
                    # Process v1-v3 with first JEC set
                    for version in v{1..3}; do
                        process_single_file "$era" "$channel" "$version" "$run_tag" "${jec_configs["2023_v123"]}"
                    done
                    # Process v4 with its specific JEC set
                    process_single_file "$era" "$channel" "v4" "$run_tag" "${jec_configs["2023_v4"]}"
                done
            done
            ;;
            
        "2023BPix")
            for run_tag in ${era_runs[${era}]}; do
                for channel in "${channels[@]}"; do
                    for version in v{1..2}; do
                        process_single_file "$era" "$channel" "$version" "$run_tag" "${jec_configs["2023BPix"]}"
                    done
                done
            done
            ;;
            
        "2022")
            for run_tag in ${era_runs[${era}]}; do
                for channel in "${channels[@]}"; do
                    process_single_file "$era" "$channel" "" "$run_tag" "${jec_configs["2022"]}"
                done
            done
            ;;
            
        "2022EE")
            # Treat _RunE/_RunF/_RunG like version tags
            for run_tag in ${era_runs[${era}]}; do
                for channel in "${channels[@]}"; do
                    local run_period=${run_tag:4:1}  # Extract E, F, or G
                    local version="_Run${run_period}"
                    local jec_key="2022EE_Run${run_period}"
                    process_single_file "$era" "$channel" "$version" "$run_tag" "${jec_configs[$jec_key]}"
                done
            done
            ;;
    esac
    
    # Combine all runs for each channel
    for channel in "${channels[@]}"; do
        channel_output="${OUTPUT_DIR}/Histograms_${era}_data_${channel}.root"
        echo "Combining runs for ${channel} into: ${channel_output}"
        hadd -f "${channel_output}" ${TMP_DIR}/Histograms_${era}_data_${channel}_*.root
    done
    
    # Copy JetMET channel to QCD
    cp ${OUTPUT_DIR}/Histograms_${era}_data_JetMET.root ${OUTPUT_DIR}/Histograms_${era}_data_QCD.root
    
    # Combine leptonic channels (Muon + EGamma)
    leptonic_output="${OUTPUT_DIR}/Histograms_${era}_data_TTbar.root"
    echo "Combining leptonic channels into: ${leptonic_output}"
    hadd -f "${leptonic_output}" \
        ${OUTPUT_DIR}/Histograms_${era}_data_Muon.root \
        ${OUTPUT_DIR}/Histograms_${era}_data_EGamma.root
    
    # Final combination (all channels)
    final_output="${OUTPUT_DIR}/Histograms_${era}_data.root"
    echo "Creating final combined output: ${final_output}"
    hadd -f "${final_output}" \
        ${OUTPUT_DIR}/Histograms_${era}_data_JetMET.root \
        ${OUTPUT_DIR}/Histograms_${era}_data_Muon.root \
        ${OUTPUT_DIR}/Histograms_${era}_data_EGamma.root
}

# Process each era
declare -a eras=("2022" "2022EE" "2023" "2023BPix")
for era in "${eras[@]}"; do
    process_era "$era"
done