#!/bin/bash -e

# Check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

# Directory setup
PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
INPUT_DIR=${SCRIPT_DIR}/trees
OUTPUT_DIR=${SCRIPT_DIR}/trees/merged
mkdir -p ${OUTPUT_DIR}

MC_TYPES=("QCD" "TTbar" "VJ" "VV")

process_year() {
    year=$1

    case ${year} in
        2022)
            eras=("2022" "2022EE")
            ;;
        2023)
            eras=("2023" "2023BPix")
            ;;
        *)
            echo "Year ${year} not supported."
            return
            ;;
    esac

    # MC
    for mc_type in ${MC_TYPES[@]}; do
        output_path=${OUTPUT_DIR}/Histograms_${year}_MC_${mc_type}.root
        input_paths=""
        for era in ${eras[@]}; do
            input_paths="${input_paths} ${INPUT_DIR}/Histograms_${era}_MC_${mc_type}.root"
        done
        hadd -f ${output_path} ${input_paths}
    done

    # Data
    output_path=${OUTPUT_DIR}/Histograms_${year}_data.root
    input_paths=""
    for era in ${eras[@]}; do
        input_paths="${input_paths} ${INPUT_DIR}/Histograms_${era}_data.root"
    done
    hadd -f ${output_path} ${input_paths}
}

# check args
if [ $# -ne 1 ]; then
    # process all eras
    eras=("2022" "2023")
    for era in "${eras[@]}"; do
        process_year ${era} &
    done
    wait
    echo "All years processing completed!"
else
    # process single era
    process_year $1
fi