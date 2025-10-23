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
cd ${SCRIPT_DIR};
BIN="0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0"
PLOT_DIR=${SCRIPT_DIR}/plots
mkdir -p ${PLOT_DIR}

# Variables and their labels
declare -A LABEL_DICT
# VARS=("FatJet1_Tau3OverTau2" "FatJet1_GloParT_XbbVsQCD" "FatJet1_GloParT_MassVis" "FatJet2_Tau3OverTau2" "FatJet2_GloParT_XbbVsQCD" "FatJet2_GloParT_MassVis" "FatJet1_pt" "FatJet1_MassSD" "FatJet2_pt" "FatJet2_MassSD")
VARS=("pTjj")
LABEL_DICT["FatJet1_pt"]="FatJet p_{T} [GeV]"
LABEL_DICT["FatJet1_MassSD"]="FatJet m_{SD} [GeV]"
LABEL_DICT["FatJet1_GloParT_XbbVsQCD"]="FatJet T_{Xbb}"
LABEL_DICT["FatJet1_GloParT_MassVis"]="FatJet M_{GloParT} [GeV]"
LABEL_DICT["FatJet1_Tau3OverTau2"]="FatJet #tau_{3}/#tau_{2}"

LABEL_DICT["FatJet2_pt"]="FatJet p_{T} [GeV]"
LABEL_DICT["FatJet2_MassSD"]="FatJet m_{SD} [GeV]"
LABEL_DICT["FatJet2_GloParT_XbbVsQCD"]="FatJet T_{Xbb}"
LABEL_DICT["FatJet2_GloParT_MassVis"]="FatJet M_{GloParT} [GeV]"
LABEL_DICT["FatJet2_Tau3OverTau2"]="FatJet #tau_{3}/#tau_{2}"

LABEL_DICT["pTjj"]="p_{T}(jj) [GeV]"

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

        # tau32+TXbb+ptjj correction
        output_path="${PLOT_DIR}/var_${year}_${var}_tau32TXbbpTjjcorr.pdf"
        path_TTbar_tau32_TXbb_ptjj="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb_pTjj.root"
        echo "Processing ${var} with tau32+TXbb+ptjj correction for year ${year}"
        root -l -b -q "plot.cpp(\"${year}\", \"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${path_TTbar_tau32_TXbb_ptjj}\", \"${path_ttHto2B}\", \"${output_path}\", \"${var}\", \"${var_label}\", \"${BIN}\")"
    done
}

if [ $# -ne 1 ]; then
    # process all years
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