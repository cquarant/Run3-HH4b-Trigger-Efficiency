#!/bin/bash -e

# Check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    exit 1
fi

# Directory setup
PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR_ERAS=${SCRIPT_DIR}/trees/eras
OUTPUT_DIR_YEARS=${SCRIPT_DIR}/trees/years
TMP_DIR=${SCRIPT_DIR}/tmp/trees/data
mkdir -p ${OUTPUT_DIR_ERAS} ${OUTPUT_DIR_YEARS} ${TMP_DIR}

declare -A YEAR_DICT=(
    ["2022"]="2022 2022EE"
    ["2023"]="2023 2023BPix"
)

# Configuration
declare -A era_runs=(
    ["2022"]="2022C 2022D"
    ["2022EE"]="2022E 2022F 2022G"
    ["2023"]="2023C"
    ["2023BPix"]="2023D"
)

# JEC configurations
declare -A jec_configs=(
    ["2022_RunC"]="Summer22_22Sep2023_RunCD_V2_DATA"
    ["2022_RunD"]="Summer22_22Sep2023_RunCD_V2_DATA"
    ["2022EE_RunE"]="Summer22EE_22Sep2023_RunE_V2_DATA"
    ["2022EE_RunF"]="Summer22EE_22Sep2023_RunF_V2_DATA"
    ["2022EE_RunG"]="Summer22EE_22Sep2023_RunG_V2_DATA"
    ["2023_v123"]="Summer23Prompt23_RunCv123_V1_DATA"
    ["2023_v4"]="Summer23Prompt23_RunCv4_V1_DATA"
    ["2023BPix"]="Summer23BPixPrompt23_RunD_V1_DATA"
)

CHANNELS=("JetMET")
SAMPLE_DIR="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"
jec_base="${PROJ_ROOT}/JECs"

process_single_file() {
    local era=$1
    local channel=$2
    local version=$3
    local run_tag=$4
    
    # Get the year from era
    local year=${era:0:4}
    
    # Extract run period (C/D/E/F/G) from run_tag
    local run_period=${run_tag: -1}

    # Determine JEC config key based on era and version/run
    local jec_key
    if [ "$era" == "2023" ]; then
        if [ "$version" == "v4" ]; then
            jec_key="2023_v4"
        else
            jec_key="2023_v123"
        fi
    elif [ "$era" == "2023BPix" ]; then
        jec_key="2023BPix"
    else
        jec_key="${era}_Run${run_period}"
    fi
    
    local jec_config="${jec_configs[$jec_key]}"

    # Handle input file version suffix for different eras
    local input_version=""
    if [ "$era" == "2023" ] || [ "$era" == "2023BPix" ]; then
        input_version="_${version}"
    fi

    # For output files, add Run prefix for 2022/2022EE
    local output_tag
    if [ "$era" == "2022" ] || [ "$era" == "2022EE" ]; then
        output_tag="Run${run_tag}"
    else
        output_tag="${run_tag}${input_version}"
    fi

    # Set up paths
    local sample_path="${SAMPLE_DIR}/${year}/${channel}_${run_tag}${input_version}.root"
    local output_path="${TMP_DIR}/Histograms_${era}_data_${channel}_${output_tag}.root"
    
    # Set appropriate JEC paths
    local jec_dir="${jec_base}/${jec_config}"
    local jec_path_L2Relative="${jec_dir}/${jec_config}_L2Relative_AK8PFPuppi.txt"
    local jec_path_L2L3Residual="${jec_dir}/${jec_config}_L2L3Residual_AK8PFPuppi.txt"

    echo "Processing ${era} ${channel} ${run_tag}${input_version} with ${jec_config} JECs"
    echo "Input file: ${sample_path}"
    root -l -b -q "tree_data.cpp(\"${year}\", \"${run_tag}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${jec_path_L2Relative}\", \"${jec_path_L2L3Residual}\")"
}

process_era() {
    local era=$1
    echo "Processing era: ${era}"
    
    case ${era} in
        "2023")
            for run_tag in ${era_runs[${era}]}; do
                for channel in "${CHANNELS[@]}"; do
                    # Process v1-v3 with first JEC set
                    for version in v1 v2 v3; do
                        process_single_file "$era" "$channel" "$version" "$run_tag"
                    done
                    # Process v4 with its specific JEC set
                    process_single_file "$era" "$channel" "v4" "$run_tag"
                done
            done
            ;;
            
        "2023BPix")
            for run_tag in ${era_runs[${era}]}; do
                for channel in "${CHANNELS[@]}"; do
                    for version in v1 v2; do
                        process_single_file "$era" "$channel" "$version" "$run_tag"
                    done
                done
            done
            ;;
            
        "2022"|"2022EE")
            for run_tag in ${era_runs[${era}]}; do
                for channel in "${CHANNELS[@]}"; do
                    process_single_file "$era" "$channel" "" "$run_tag"
                done
            done
            ;;
    esac
    
    # Combine all runs for each channel
    for channel in "${CHANNELS[@]}"; do
        channel_output="${TMP_DIR}/Histograms_${era}_data_${channel}.root"
        echo "Combining runs for ${channel} into: ${channel_output}"
        hadd -f "${channel_output}" ${TMP_DIR}/Histograms_${era}_data_${channel}_*.root
    done
    
    # Final combination (all CHANNELS)
    final_output="${OUTPUT_DIR_ERAS}/Histograms_${era}_data.root"
    cp ${TMP_DIR}/Histograms_${era}_data_JetMET.root ${final_output}
    echo "Creating final output: ${final_output}"
}

combine_years() {
    # hadd eras into years using YEAR_DICT
    for year in "${!YEAR_DICT[@]}"; do
        year_files=()
        for era in ${YEAR_DICT[$year]}; do
            year_files+=("${OUTPUT_DIR_ERAS}/Histograms_${era}_data.root")
        done
        
        # Combine all eras for the year
        year_output="${OUTPUT_DIR_YEARS}/Histograms_${year}_data.root"
        echo "Combining eras for ${year} into: ${year_output}"
        hadd -f "${year_output}" "${year_files[@]}"
    done
}

if [ $# -ne 1 ]; then
    # process all eras
    eras=("2022" "2022EE" "2023" "2023BPix")
    for era in "${eras[@]}"; do
        process_era ${era} &
    done
    wait
    echo "All data processing completed!"

    echo "Combining all years..."
    combine_years
    echo "All years combined successfully!"
else
    # process single era
    process_era $1
fi
