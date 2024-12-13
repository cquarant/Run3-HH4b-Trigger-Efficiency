#!/bin/bash -e

# Check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

# Directory setup
PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
HIST_DIR=${SCRIPT_DIR}/hists
SF_DIR=${SCRIPT_DIR}/SFs
TREE_DIR=${SCRIPT_DIR}/trees
mkdir -p ${SF_DIR}
cd ${SCRIPT_DIR};

get_paths() {
    era=$1;

    path_data="${HIST_DIR}/Histograms_${era}_data.root"
    path_QCD="${HIST_DIR}/Histograms_${era}_MC_QCD.root"
    path_VV="${HIST_DIR}/Histograms_${era}_MC_VV.root"
    path_VJ="${HIST_DIR}/Histograms_${era}_MC_VJ.root"
    path_TTbar="${HIST_DIR}/Histograms_${era}_MC_TTbar.root"
    path_sf_tau32="${SF_DIR}/SF_tau3overtau2_${era}.root"
    path_sf_TXbb="${SF_DIR}/SF_TXbb_${era}.root"
}

process_era() {
    era=$1;

    get_paths $era
    root -l -b -q "SF_tau32.cpp(\"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${path_sf_tau32}\")"
    
    # histogram file paths
    input_path="${TREE_DIR}/Histograms_${era}_MC_TTbar.root"
    output_path="${HIST_DIR}/Histograms_${era}_MC_TTbar_tau32.root"

    root -l -b -q "make_hist.cpp(\"${input_path}\", \"${output_path}\", \"${path_sf_tau32}\")"

    # TXbb SF after tau32
    path_TTbar="${HIST_DIR}/Histograms_${era}_MC_TTbar_tau32.root"
    root -l -b -q "SF_TXbb.cpp(\"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${path_sf_TXbb}\")"

    input_path="${TREE_DIR}/Histograms_${era}_MC_TTbar.root"
    output_path="${HIST_DIR}/Histograms_${era}_MC_TTbar_tau32_TXbb.root"
    root -l -b -q "make_hist.cpp(\"${input_path}\", \"${output_path}\", \"${path_sf_tau32}\", \"${path_sf_TXbb}\")"
    
}

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