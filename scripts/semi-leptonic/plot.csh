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
PLOT_DIR=${SCRIPT_DIR}/plots
mkdir -p ${PLOT_DIR}
cd ${SCRIPT_DIR};

# Variables and their labels
declare -A LABEL_DICT
VARS=("lep1_pt" "FatJet1_pt" "FatJet1_MassSD" "FatJet1_GloParT_XbbVsQCD" "FatJet1_Tau3OverTau2")
LABEL_DICT["lep1_pt"]="p_{T}^{lep} [GeV]"
LABEL_DICT["FatJet1_pt"]="p_{T} [GeV]"
LABEL_DICT["FatJet1_MassSD"]="m_{SD} [GeV]"
LABEL_DICT["FatJet1_GloParT_XbbVsQCD"]="T_{Xbb}"
LABEL_DICT["FatJet1_Tau3OverTau2"]="#tau_{3}/#tau_{2}"

get_paths() {
    era=$1;

    path_data="${HIST_DIR}/Histograms_${era}_data.root"
    path_QCD="${HIST_DIR}/Histograms_${era}_MC_QCD.root"
    path_VV="${HIST_DIR}/Histograms_${era}_MC_VV.root"
    path_VJ="${HIST_DIR}/Histograms_${era}_MC_VJ.root"
    path_TTbar="${HIST_DIR}/Histograms_${era}_MC_TTbar.root"
}

process_era() {
    era=$1;

    get_paths ${era}
    
    # Loop over all variables
    for var in "${VARS[@]}"; do
        var_label="${LABEL_DICT[$var]}"
        output_path="${PLOT_DIR}/var_${era}_${var}.pdf"
        echo "Processing ${var} for era ${era}"
        root -l -b -q "plot.cpp(\"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${output_path}\", \"${var}\", \"${var_label}\")"

        # tau32 correction
        output_path="${PLOT_DIR}/var_${era}_${var}_tau32corr.pdf"
        path_TTbar="${HIST_DIR}/Histograms_${era}_MC_TTbar_tau32.root"
        echo "Processing ${var} with tau32 correction for era ${era}"
        root -l -b -q "plot.cpp(\"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${output_path}\", \"${var}\", \"${var_label}\")"

        # tau32+TXbb correction
        output_path="${PLOT_DIR}/var_${era}_${var}_tau32TXbbcorr.pdf"
        path_TTbar="${HIST_DIR}/Histograms_${era}_MC_TTbar_tau32_TXbb.root"
        echo "Processing ${var} with tau32+TXbb correction for era ${era}"
        root -l -b -q "plot.cpp(\"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${output_path}\", \"${var}\", \"${var_label}\")"
    done
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