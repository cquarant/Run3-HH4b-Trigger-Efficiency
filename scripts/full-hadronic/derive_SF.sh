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
SF_DIR=${PROJ_ROOT}/TTbar_SFs
TREE_DIR=${SCRIPT_DIR}/trees/years
mkdir -p ${SF_DIR} ${HIST_DIR}
cd ${SCRIPT_DIR};

MC_CHANNELS=("QCD" "TTbar" "VJ" "VV" "ttHto2B")

process_data() {
    year=$1;

    local input_path="${TREE_DIR}/Histograms_${year}_data.root"
    local output_path="${HIST_DIR}/Histograms_${year}_data.root"

    # Make histograms
    root -l -b -q "${SCRIPT_DIR}/make_hist.cpp(\"${input_path}\", \"${output_path}\")"
}

process_mc() {
    year=$1;

    for channel in "${MC_CHANNELS[@]}"; do
        local input_path="${TREE_DIR}/Histograms_${year}_MC_${channel}.root"
        local output_path="${HIST_DIR}/Histograms_${year}_MC_${channel}.root"

        # Make histograms
        root -l -b -q "${SCRIPT_DIR}/make_hist.cpp(\"${input_path}\", \"${output_path}\")"
    done
}

get_paths() {
    year=$1;

    path_data="${HIST_DIR}/Histograms_${year}_data.root"
    path_QCD="${HIST_DIR}/Histograms_${year}_MC_QCD.root"
    path_VV="${HIST_DIR}/Histograms_${year}_MC_VV.root"
    path_VJ="${HIST_DIR}/Histograms_${year}_MC_VJ.root"
    path_TTbar="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb.root"
    path_ttHto2B="${HIST_DIR}/Histograms_${year}_MC_ttHto2B.root"
    path_sf_tau32="${SF_DIR}/SF_tau3overtau2_${year}.root"
    path_sf_TXbb="${SF_DIR}/SF_TXbb_${year}.root"
    path_sf_pTjj="${SF_DIR}/SF_pTjj_${year}.root"
}

process_year() {
    year=$1;
    get_paths $year

    # prepare histograms from trees
    process_data ${year}
    process_mc ${year}

    # derive prescale factors
    python3 ${SCRIPT_DIR}/kfact.py --year ${year} --ttbar-corr "tau32_TXbb"

    # derive SFs for pTjj
    root -l -b -q "SF_pTjj.cpp(\"${year}\", \"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${path_ttHto2B}\", \"${path_sf_pTjj}\")"
    # make histograms with pTjj SF to TTbar
    input_path="${TREE_DIR}/Histograms_${year}_MC_TTbar.root"
    output_path="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb_pTjj.root"
    root -l -b -q "make_hist.cpp(\"${input_path}\", \"${output_path}\", \"${path_sf_tau32}\", \"${path_sf_TXbb}\", \"${path_sf_pTjj}\")"
    # plot
    echo "Plotting..."
    ${SCRIPT_DIR}/plot.csh ${year}
    
}

if [ $# -ne 1 ]; then
    # process all years
    # years=("2022" "2022EE" "2023" "2023BPix")
    years=("2022" "2023")
    for year in "${years[@]}"; do
        process_year ${year}
    done
    wait
    echo "All MC processing completed!"
else
    # process single year
    process_year $1
fi