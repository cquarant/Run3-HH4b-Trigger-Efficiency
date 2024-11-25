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
SAMPLE_DIR="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"

# Era Configuration
declare -A era_runs=(
    ["2022"]="2022C 2022D"
    ["2022EE"]="2022E 2022F 2022G"
    ["2023"]="2023C"
    ["2023BPix"]="2023D"
)

# JEC configurations
declare -A jec_configs=(
    ["2023"]="Summer23Prompt23_V1_MC"
    ["2023BPix"]="Summer23BPixPrompt23_V1_MC"
    ["2022"]="Summer22_22Sep2023_V2_MC"
    ["2022EE"]="Summer22EE_22Sep2023_V2_MC"
)

get_era_paths() {
    local era=$1
    local process=$2  # "ttbar" or "qcd"
    
    case ${era} in
        "2023")
            pu_path="${PROJ_ROOT}/pileups/pu_2023C.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_${era}.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            ;;
        "2023BPix")
            # Use same structure as 2023 but with BPix specific paths
            pu_path="${PROJ_ROOT}/pileups/pu_2023D.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_2023.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            ;;
        "2022")
            pu_path="${PROJ_ROOT}/pileups/pu_2022CD.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_${era}.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            ;;
        "2022EE")
            pu_path="${PROJ_ROOT}/pileups/pu_2022EFG.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_${era}.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            ;;
        *)
            echo "Error: Invalid era ${era}"
            return 1
            ;;
    esac
    
    # Set JEC paths
    jec_path_ak4="${jec_base}/${jec_configs[${era}]}_L2Relative_AK4PFPuppi.txt"
    jec_path_ak8="${jec_base}/${jec_configs[${era}]}_L2Relative_AK8PFPuppi.txt"

    # Trigger efficiency paths (as a function of mass and pt)
    if [ "$process" == "ttbar" ]; then
        sf_path="${PROJ_ROOT}/scripts/sf_mass_pt/output/efficiency_mass_pt_2023_TTbar.root"
    else
        sf_path="${PROJ_ROOT}/scripts/sf_mass_pt/output/efficiency_mass_pt_2023_QCD.root"
    fi
}

process_ttbar() {
    local era=$1
    echo "Processing TTBar samples for era: ${era}..."
    
    # Get paths specific to TTBar processing
    get_era_paths ${era} "ttbar" || return 1
    
    local channels=("Muon" "EGamma" "Leptonic")
    local ttbar_type="TTtoLNu2Q"
    
    # Find TTBar sample with wildcard
    local ttbar_file=$(ls ${SAMPLE_DIR}/${era}/TTtoLNu2Q*.root 2>/dev/null | head -n1)
    
    if [ -z "${ttbar_file}" ]; then
        echo "Warning: No TTBar file found for era ${era}"
        return 1
    fi
    
    for channel in "${channels[@]}"; do
        local output_path=${TMP_DIR}/Histograms_${era}_MC_TTtoLNu2Q_${channel}.root
        
        echo "Processing TTBar for channel: ${channel}"
        root -l -b -q "histo_mc.cpp(\"${ttbar_type}\", \"${channel}\", \"${ttbar_file}\", \
            \"${output_path}\", \"${pu_path}\", \"${sf_path}\", \"${param_path}\", \
            \"${jec_path_ak4}\", \"${jec_path_ak8}\")"
    done
    
    cp ${TMP_DIR}/Histograms_${era}_MC_TTtoLNu2Q_Leptonic.root \
        ${OUTPUT_DIR}/Histograms_${era}_MC_TTbar.root
}

process_qcd() {
    local era=$1
    echo "Processing QCD samples for era: ${era}..."
    
    # Get paths specific to QCD processing
    get_era_paths ${era} "qcd" || return 1
    
    local channel="QCD"
    local HT_BINS=("100to200" "200to400" "400to600" "600to800" "800to1000" "1000to1200" "1200to1500" "1500to2000" "2000toInf")
    
    for ht_bin in "${HT_BINS[@]}"; do
        # Handle special case for 2000toInf
        if [ "${ht_bin}" == "2000toInf" ]; then
            search_pattern="QCD-4Jets_HT-200*.root"
        else
            search_pattern="QCD-4Jets_HT-${ht_bin}*.root"
        fi
        
        # Find the exact file name
        local sample_file=$(ls ${SAMPLE_DIR}/${era}/${search_pattern} 2>/dev/null | head -n1)
        
        if [ -z "${sample_file}" ]; then
            echo "Warning: No matching file found for HT bin ${ht_bin} in ${era}"
            continue
        fi

        echo "Processing QCD HT bin: ${ht_bin}"
        echo "Sample file: ${sample_file}"

        local output_path=${TMP_DIR}/Histograms_${era}_MC_QCD-4Jets_HT-${ht_bin}.root
        local qcd_type="QCD_HT_${ht_bin}"

        root -l -b -q "histo_mc.cpp(\"${qcd_type}\", \"${channel}\", \"${sample_file}\", \
            \"${output_path}\", \"${pu_path}\", \"${sf_path}\", \"${param_path}\", \
            \"${jec_path_ak4}\", \"${jec_path_ak8}\")"
    done
    
    # Combine all QCD samples
    local qcd_path=${OUTPUT_DIR}/Histograms_${era}_MC_QCD.root
    echo "Combining QCD samples for era ${era}"
    hadd -f ${qcd_path} ${TMP_DIR}/Histograms_${era}_MC_QCD-4Jets_HT-*.root
}

combine_final_output() {
    local era=$1
    local final_output="${OUTPUT_DIR}/Histograms_${era}_MC.root"
    local ttbar_input="${OUTPUT_DIR}/Histograms_${era}_MC_TTbar.root"
    local qcd_input="${OUTPUT_DIR}/Histograms_${era}_MC_QCD.root"
    
    echo "Creating final combined output for era ${era}: ${final_output}"
    hadd -f ${final_output} ${ttbar_input} ${qcd_input}
}

process_era() {
    local era=$1
    echo "Processing era: ${era}"

    # Process QCD samples
    process_qcd ${era}
    
    # Process TTBar samples
    process_ttbar ${era}

    # Combine all MC samples
    combine_final_output ${era}
    
    echo "MC processing for era ${era} completed successfully!"
}

main() {
    # declare -a eras=("2022" "2022EE" "2023" "2023BPix")
    declare -a eras=("2023" "2023BPix")
    for era in "${eras[@]}"; do
        echo "Starting processing for era: ${era}"
        if process_era ${era}; then
            echo "Successfully processed era ${era}"
        else
            echo "Warning: Some processing steps failed for era ${era}"
        fi
    done
    echo "All MC processing completed!"
}

# Execute main function
main
