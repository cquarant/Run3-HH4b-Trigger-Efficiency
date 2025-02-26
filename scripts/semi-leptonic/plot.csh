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
cd ${SCRIPT_DIR};
# BIN="0.8,0.9,0.94,0.97,0.99,1.0"
# PLOT_DIR=${SCRIPT_DIR}/plots_sfBDTBin
BIN="0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1.0"
PLOT_DIR=${SCRIPT_DIR}/plots_sfBDTFinerBin
# BIN="0.990,0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998,0.999,1.0"
# PLOT_DIR=${SCRIPT_DIR}/plots_sfBDTEvenFinerBin
mkdir -p ${PLOT_DIR}

# Variables and their labels
declare -A LABEL_DICT
VARS=("lep1_pt" "FatJet1_pt" "FatJet1_MassSD" "FatJet1_GloParT_XbbVsQCD" "FatJet1_ParticleNetLegacy_XbbVsQCD" "FatJet1_Tau3OverTau2")
LABEL_DICT["lep1_pt"]="Lepton p_{T} [GeV]"
LABEL_DICT["FatJet1_pt"]="FatJet p_{T} [GeV]"
LABEL_DICT["FatJet1_MassSD"]="FatJet m_{SD} [GeV]"
LABEL_DICT["FatJet1_GloParT_XbbVsQCD"]="FatJet T_{Xbb}"
LABEL_DICT["FatJet1_ParticleNetLegacy_XbbVsQCD"]="FatJet T_{Xbb}"
LABEL_DICT["FatJet1_Tau3OverTau2"]="FatJet #tau_{3}/#tau_{2}"

get_paths() {
    year=$1;

    path_data="${HIST_DIR}/Histograms_${year}_data.root"
    path_QCD="${HIST_DIR}/Histograms_${year}_MC_QCD.root"
    path_VV="${HIST_DIR}/Histograms_${year}_MC_VV.root"
    path_VJ="${HIST_DIR}/Histograms_${year}_MC_VJ.root"
    path_TTbar="${HIST_DIR}/Histograms_${year}_MC_TTbar.root"
    path_ttHto2B="${HIST_DIR}/Histograms_${year}_MC_ttHto2B.root"
}

process_year() {
    year=$1;

    get_paths ${year}
    
    # Loop over all variables
    for var in "${VARS[@]}"; do
        var_label="${LABEL_DICT[$var]}"
        output_path="${PLOT_DIR}/var_${year}_${var}.pdf"
        echo "Processing ${var} for year ${year}"
        root -l -b -q "plot.cpp(\"${year}\", \"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar}\", \"${path_ttHto2B}\", \"${output_path}\", \"${var}\", \"${var_label}\", \"${BIN}\")"

        # tau32 correction
        output_path="${PLOT_DIR}/var_${year}_${var}_tau32corr.pdf"
        path_TTbar_tau32="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32.root"
        echo "Processing ${var} with tau32 correction for year ${year}"
        root -l -b -q "plot.cpp(\"${year}\", \"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar_tau32}\", \"${path_ttHto2B}\", \"${output_path}\", \"${var}\", \"${var_label}\", \"${BIN}\")"

        # tau32+TXbb correction
        output_path="${PLOT_DIR}/var_${year}_${var}_tau32TXbbcorr.pdf"
        path_TTbar_tau32_TXbb="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb.root"
        echo "Processing ${var} with tau32+TXbb correction for year ${year}"
        root -l -b -q "plot.cpp(\"${year}\", \"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar_tau32_TXbb}\", \"${path_ttHto2B}\", \"${output_path}\", \"${var}\", \"${var_label}\", \"${BIN}\")"
    done
}

if [ $# -ne 1 ]; then
    # process all years
    # years=("2022" "2022EE" "2023" "2023BPix")
    years=("2022" "2023")
    for year in "${years[@]}"; do
        process_year ${year} &
    done
    wait
    echo "All MC processing completed!"
else
    # process single year
    process_year $1
fi