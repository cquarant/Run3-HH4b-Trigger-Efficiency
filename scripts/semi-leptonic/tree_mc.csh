#!/bin/bash -e

# Check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

# Directory setup
PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}/trees
TMP_DIR=${SCRIPT_DIR}/tmp/trees/mc
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

# JER configurations
declare -A jer_configs=(
    ["2023"]="Summer23Prompt23_RunCv1234_JRV1_MC"
    ["2023BPix"]="Summer23BPixPrompt23_RunD_JRV1_MC"
    ["2022"]="Summer22_22Sep2023_JRV1_MC"
    ["2022EE"]="Summer22EE_22Sep2023_JRV1_MC"
)

get_era_paths() {
    local era=$1
    
    case ${era} in
        "2023")
            pu_path="${PROJ_ROOT}/pileups/pu_2023C.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_${era}.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            jer_base="${PROJ_ROOT}/JERs/${jer_configs[${era}]}"
            ;;
        "2023BPix")
            pu_path="${PROJ_ROOT}/pileups/pu_2023D.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_${era}.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            jer_base="${PROJ_ROOT}/JERs/${jer_configs[${era}]}"
            ;;
        "2022")
            pu_path="${PROJ_ROOT}/pileups/pu_2022CD.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_${era}.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            jer_base="${PROJ_ROOT}/JERs/${jer_configs[${era}]}"
            ;;
        "2022EE")
            pu_path="${PROJ_ROOT}/pileups/pu_2022EFG.txt"
            param_path="${PROJ_ROOT}/parameters/parameters_${era}.txt"
            jec_base="${PROJ_ROOT}/JECs/${jec_configs[${era}]}"
            jer_base="${PROJ_ROOT}/JERs/${jer_configs[${era}]}"
            ;;
        *)
            echo "Error: Invalid era ${era}"
            exit 1
            ;;
    esac
    
    # Set JEC path
    jec_path="${jec_base}/${jec_configs[${era}]}_L2Relative_AK8PFPuppi.txt"

    # Set JER paths
    jer_path="${jer_base}/${jer_configs[${era}]}_PtResolution_AK8PFPuppi.txt"
    jer_path_sf="${jer_base}/${jer_configs[${era}]}_SF_AK8PFPuppi.txt"
}

process_DYto2L() {
    local era=$1
    local year=${era:0:4}
    
    get_era_paths ${era} || exit 1
    
    JETS=(0 1 2)
    for jet in "${JETS[@]}"; do
        # DYto2L-2Jets_MLL-50_0J.root
        local input_file=$(ls ${SAMPLE_DIR}/${era}/DYto2L-2Jets_MLL-50_${jet}J*.root 2>/dev/null | head -n1)
        
        if [ -z "${input_file}" ]; then
            echo "Warning: No DY file found for era ${era}"
            exit 1
        fi
        local sample_type="DYto2L_2Jets_MLL_50_${jet}J"
        local output_path="${TMP_DIR}/Histograms_${era}_MC_DYto2L_2Jets_MLL_50_${jet}J.root"
        
        root -l -b -q "tree_mc.cpp(\"${year}\", \"${sample_type}\", \"${input_file}\", \
            \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
            \"${jer_path}\", \"${jer_path_sf}\")"
    done
}


process_QCD() {
    local era=$1
    local year=${era:0:4}
    
    get_era_paths ${era} || exit 1
    
    local HT_BINS=("200to400" "400to600" "600to800" "800to1000" "1000to1200" "1200to1500" "1500to2000" "2000toInf")
    
    for ht_bin in "${HT_BINS[@]}"; do
        if [ "${ht_bin}" == "2000toInf" ]; then
            search_pattern="QCD-4Jets_HT-2000*.root"
        else
            search_pattern="QCD-4Jets_HT-${ht_bin}*.root"
        fi

        local input_file=$(ls ${SAMPLE_DIR}/${era}/${search_pattern} 2>/dev/null | head -n1)
        
        if [ -z "${input_file}" ]; then
            echo "Warning: No QCD file found for era ${era}"
            exit 1
        fi
        local sample_type="QCD_HT_${ht_bin}"
        local output_path=${TMP_DIR}/Histograms_${era}_MC_QCD-4Jets_HT-${ht_bin}.root

        root -l -b -q "tree_mc.cpp(\"${year}\", \"${sample_type}\", \"${input_file}\", \
            \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
            \"${jer_path}\", \"${jer_path_sf}\")"
    done

    local output_path=${OUTPUT_DIR}/Histograms_${era}_MC_QCD.root
    hadd -f ${output_path} ${TMP_DIR}/Histograms_${era}_MC_QCD-4Jets_HT-*.root
}

process_TTbar() {
    local era=$1
    local year=${era:0:4}
    
    get_era_paths ${era} || exit 1

    local TTbar_types=("TTtoLNu2Q" "TTto2L2Nu" "TTto4Q")

    for ttbar_type in "${TTbar_types[@]}"; do
        local input_file=$(ls ${SAMPLE_DIR}/${era}/${ttbar_type}*.root 2>/dev/null | head -n1)
        
        if [ -z "${input_file}" ]; then
            echo "Warning: No TTBar file found for era ${era}"
            exit 1
        fi
        local output_path="${TMP_DIR}/Histograms_${era}_MC_${ttbar_type}.root"
        
        root -l -b -q "tree_mc.cpp(\"${year}\", \"${ttbar_type}\", \"${input_file}\", \
            \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
            \"${jer_path}\", \"${jer_path_sf}\")"
    done
    
    local output_path="${OUTPUT_DIR}/Histograms_${era}_MC_TTbar.root"
    hadd -f ${output_path} ${TMP_DIR}/Histograms_${era}_MC_TTto*.root
}

process_VV() {
    local era=$1
    local year=${era:0:4}
    
    get_era_paths ${era} || exit 1

    local types=("WW" "WZ" "ZZ")
    if [ "${era}" == "2023BPix" ]; then
        types=("WW" "WZ")
    fi
    local hadd_inputs=""

    for type in "${types[@]}"; do

        local input_file=$(ls ${SAMPLE_DIR}/${era}/${type}.root 2>/dev/null | head -n1)
        
        if [ -z "${input_file}" ]; then
            echo "Warning: No VV file found for era ${era}"
            continue
        fi
        local output_path="${TMP_DIR}/Histograms_${era}_MC_${type}.root"
        
        root -l -b -q "tree_mc.cpp(\"${year}\", \"${type}\", \"${input_file}\", \
            \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
            \"${jer_path}\", \"${jer_path_sf}\")"
        hadd_inputs="${hadd_inputs} ${output_path}"
    done

    local output_path="${OUTPUT_DIR}/Histograms_${era}_MC_VV.root"
    hadd -f ${output_path} ${hadd_inputs}
}

process_Wto2Q() {
    local era=$1
    local year=${era:0:4}
    
    get_era_paths ${era} || exit 1

    PT_BINS=("100to200" "200to400" "400to600" "600")
    JETS=(1 2)
    for pt_bin in "${PT_BINS[@]}"; do
        for jet in "${JETS[@]}"; do
            local input_file=$(ls ${SAMPLE_DIR}/${era}/Wto2Q-2Jets_PTQQ-${pt_bin}_${jet}J*.root 2>/dev/null | head -n1)
            if [ -z "${input_file}" ]; then
                echo "Warning: No Wto2Q file found for era ${era}"
                exit 1
            fi
            local sample_type="Wto2Q_2Jets_PTQQ_${pt_bin}_${jet}J"
            local output_path="${TMP_DIR}/Histograms_${era}_MC_Wto2Q_2Jets_PTQQ_${pt_bin}_${jet}J.root"
            
            root -l -b -q "tree_mc.cpp(\"${year}\", \"${sample_type}\", \"${input_file}\", \
                \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
                \"${jer_path}\", \"${jer_path_sf}\")"
        done
    done
}

process_WtoLNu() {
    local era=$1
    local year=${era:0:4}
    
    get_era_paths ${era} || exit 1
    JETS=(0 1 2)
    if [ ${era} == "2023BPix" ]; then
        JETS=(0 1)  # TODO: Fix 2 for 2023BPix (corrupted file)
    fi
    for jet in "${JETS[@]}"; do
        local input_file=$(ls ${SAMPLE_DIR}/${era}/WtoLNu-2Jets_${jet}J.root 2>/dev/null | head -n1)
        
        if [ -z "${input_file}" ]; then
            echo "Warning: No WtoLNu file found for era ${era}"
            exit 1
        fi
        local sample_type="WtoLNu_2Jets_${jet}J"
        local output_path="${TMP_DIR}/Histograms_${era}_MC_WtoLNu_2Jets_${jet}J.root"
        
        root -l -b -q "tree_mc.cpp(\"${year}\", \"${sample_type}\", \"${input_file}\", \
            \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
            \"${jer_path}\", \"${jer_path_sf}\")"
    done
    output="${OUTPUT_DIR}/Histograms_${era}_WtoLNu.root"
    hadd -f ${output} ${TMP_DIR}/Histograms_${era}_MC_WtoLNu_*.root
}

process_Zto2Q() {
    local era=$1
    local year=${era:0:4}
    
    get_era_paths ${era} || exit 1
    PT_BINS=("100to200" "200to400" "400to600" "600")
    JETS=(1 2)
    for pt_bin in "${PT_BINS[@]}"; do
        for jet in "${JETS[@]}"; do
            local input_file=$(ls ${SAMPLE_DIR}/${era}/Zto2Q-2Jets_PTQQ-${pt_bin}_${jet}J*.root 2>/dev/null | head -n1)
            
            if [ -z "${input_file}" ]; then
                echo "Warning: No Zto2Q file found for era ${era}"
                exit 1
            fi
            local sample_type="Zto2Q_2Jets_PTQQ_${pt_bin}_${jet}J"
            local output_path="${TMP_DIR}/Histograms_${era}_MC_Zto2Q_2Jets_PTQQ_${pt_bin}_${jet}J.root"
            
            root -l -b -q "tree_mc.cpp(\"${year}\", \"${sample_type}\", \"${input_file}\", \
                \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
                \"${jer_path}\", \"${jer_path_sf}\")"
        done
    done
}

process_VJ() {
    local era=$1
    local output_path="${OUTPUT_DIR}/Histograms_${era}_MC_VJ.root"
    hadd -f ${output_path} ${TMP_DIR}/Histograms_${era}_MC_Wto2Q_2Jets_PTQQ_*.root ${TMP_DIR}/Histograms_${era}_MC_WtoLNu*.root ${TMP_DIR}/Histograms_${era}_MC_Zto2Q_2Jets_PTQQ_*.root ${TMP_DIR}/Histograms_${era}_MC_DYto2L_2Jets_MLL_50_*.root
}

process_ttHto2B() {
    local era=$1
    local year=${era:0:4}
    local output_path="${OUTPUT_DIR}/Histograms_${era}_MC_ttHto2B.root"
    local input_file=$(ls ${SAMPLE_DIR}/${era}/ttHto2B*.root 2>/dev/null | head -n1)

    if [ -z "${input_file}" ]; then
        echo "Warning: No ttHto2B file found for era ${era}"
        exit 1
    fi
    local sample_type="ttHto2B_M_125"
    local output_path="${TMP_DIR}/Histograms_${era}_MC_ttHto2B.root"
    
    root -l -b -q "tree_mc.cpp(\"${year}\", \"${sample_type}\", \"${input_file}\", \
        \"${output_path}\", \"${pu_path}\", \"${param_path}\", \"${jec_path}\", \
        \"${jer_path}\", \"${jer_path_sf}\")"

}

process_era() {
    local era=$1
    echo "Processing era: ${era}"

    process_QCD ${era}
    process_TTbar ${era}
    process_VV ${era}
    process_Wto2Q ${era}
    process_Zto2Q ${era}
    process_WtoLNu ${era}
    process_DYto2L ${era}
    process_VJ ${era}

    echo "MC processing for era ${era} completed successfully!"
}

# check args
if [ $# -ne 1 ]; then
    # process all eras
    eras=("2022" "2022EE" "2023" "2023BPix")
    for era in "${eras[@]}"; do
        process_era ${era} &
    done
    wait
    echo "All MC processing completed!"
else
    # process single era
    process_era $1
fi