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

MC_CHANNELS=("QCD" "TTbar" "VJ" "VV")

process_data() {
    era=$1;

    local input_path="${TREE_DIR}/Histograms_${era}_data.root"
    local output_path="${HIST_DIR}/Histograms_${era}_data.root"

    # Make histograms
    root -l -b -q "${SCRIPT_DIR}/make_hist.cpp(\"${input_path}\", \"${output_path}\")"
}

process_mc() {
    era=$1;

    for channel in "${MC_CHANNELS[@]}"; do
        local input_path="${TREE_DIR}/Histograms_${era}_MC_${channel}.root"
        local output_path="${HIST_DIR}/Histograms_${era}_MC_${channel}.root"

        # Make histograms
        root -l -b -q "${SCRIPT_DIR}/make_hist.cpp(\"${input_path}\", \"${output_path}\")"
    done
}

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

    # prepare histograms from trees
    process_data ${era}
    process_mc ${era}

    # derive SFs for tau3/tau2
    root -l -b -q "SF_tau32.cpp(\"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${path_sf_tau32}\")"
    # make histograms with tau32 SF to TTbar
    input_path="${TREE_DIR}/Histograms_${era}_MC_TTbar.root"
    output_path="${HIST_DIR}/Histograms_${era}_MC_TTbar_tau32.root"
    root -l -b -q "make_hist.cpp(\"${input_path}\", \"${output_path}\", \"${path_sf_tau32}\")"

    # derive SFs for TXbb after applying tau32 SF
    path_TTbar="${HIST_DIR}/Histograms_${era}_MC_TTbar_tau32.root"  # after applying tau32 SF
    root -l -b -q "SF_TXbb.cpp(\"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${path_sf_TXbb}\")"
    # make histograms with tau32 and TXbb SF to TTbar
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