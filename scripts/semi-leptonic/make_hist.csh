#!/bin/bash -e

# Check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    exit 1
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
TREE_DIR=${SCRIPT_DIR}/trees
HIST_DIR=${SCRIPT_DIR}/hists
mkdir -p ${HIST_DIR}

MC_CHANNELS=("QCD" "TTbar" "VJ" "VV")

process_data() {
    era=$1;

    input_path="${TREE_DIR}/Histograms_${era}_data.root"
    output_path="${HIST_DIR}/Histograms_${era}_data.root"

    # Make histograms
    root -l -b -q "${SCRIPT_DIR}/make_hist.cpp(\"${input_path}\", \"${output_path}\")"
}

process_mc() {
    era=$1;

    for channel in "${MC_CHANNELS[@]}"; do
        input_path="${TREE_DIR}/Histograms_${era}_MC_${channel}.root"
        output_path="${HIST_DIR}/Histograms_${era}_MC_${channel}.root"

        # Make histograms
        root -l -b -q "${SCRIPT_DIR}/make_hist.cpp(\"${input_path}\", \"${output_path}\")"
    done
}


if [ $# -ne 1 ]; then
    # process all eras
    eras=("2022" "2022EE" "2023" "2023BPix")
    for era in "${eras[@]}"; do
        process_data ${era}
        process_mc ${era}
    done
    wait
    echo "All data processing completed!"
else
    # process single era
    process_data $1
    process_mc $1
fi