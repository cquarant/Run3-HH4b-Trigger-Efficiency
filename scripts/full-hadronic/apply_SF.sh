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

# Plot variables
BIN="0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0"
PLOT_DIR=${SCRIPT_DIR}/plots
mkdir -p ${PLOT_DIR}

# Variables and their labels
declare -A LABEL_DICT
VARS=("pTjj" "FatJet1_Tau3OverTau2" "FatJet1_GloParT_XbbVsQCD" "FatJet1_GloParT_MassVis" "FatJet2_Tau3OverTau2" "FatJet2_GloParT_XbbVsQCD" "FatJet2_GloParT_MassVis" "FatJet1_pt" "FatJet1_MassSD" "FatJet2_pt" "FatJet2_MassSD")

LABEL_DICT["pTjj"]="p_{T}^{JJ} [GeV]"

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

get_paths() {
    year=$1;

    path_data="${HIST_DIR}/Histograms_${year}_data.root"
    path_QCD="${HIST_DIR}/Histograms_${year}_MC_QCD.root"
    path_VV="${HIST_DIR}/Histograms_${year}_MC_VV.root"
    path_VJ="${HIST_DIR}/Histograms_${year}_MC_VJ.root"
    path_TTbar="${HIST_DIR}/Histograms_${year}_MC_TTbar.root"
    path_ttHto2B="${HIST_DIR}/Histograms_${year}_MC_ttHto2B.root"

    path_sf_tau32="${SF_DIR}/SF_tau3overtau2_${year}.root"
    path_sf_TXbb="${SF_DIR}/SF_TXbb_${year}.root"
    path_sf_pTjj="${SF_DIR}/SF_pTjj_${year}.root"
}

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

plot_year() {
    year=$1
    corr_flag=$2  # New parameter: tau32, tau32_TXbb, tau32_TXbb_pTjj

    path_data="${HIST_DIR}/Histograms_${year}_data.root"
    path_QCD="${HIST_DIR}/Histograms_${year}_MC_QCD.root"
    path_VV="${HIST_DIR}/Histograms_${year}_MC_VV.root"
    path_VJ="${HIST_DIR}/Histograms_${year}_MC_VJ.root"
    path_TTbar="${HIST_DIR}/Histograms_${year}_MC_TTbar.root"
    path_ttHto2B="${HIST_DIR}/Histograms_${year}_MC_ttHto2B.root"
    
    # Determine which TTbar path to use based on corr_flag
    ttbar_path="${path_TTbar}"
    suffix=""
    
    if [ "${corr_flag}" == "tau32" ]; then
        ttbar_path="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32.root"
        suffix="_tau32corr"
    elif [ "${corr_flag}" == "tau32_TXbb" ]; then
        ttbar_path="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb.root"
        suffix="_tau32TXbbcorr"
    elif [ "${corr_flag}" == "tau32_TXbb_pTjj" ]; then
        ttbar_path="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb_pTjj.root"
        suffix="_tau32TXbbpTjjcorr"
    else 
        ttbar_path="${HIST_DIR}/Histograms_${year}_MC_TTbar.root"
        suffix=""
    fi
    
    # Loop over all variables
    for var in "${VARS[@]}"; do
        local var_label="${LABEL_DICT[$var]}"
        local output_path="${PLOT_DIR}/var_${year}_${var}${suffix}.pdf"
        echo "Processing ${var} for year ${year} with ${corr_flag:-no} correction"
        root -l -b -q "plot.cpp(\"${year}\", \"${path_data}\", \"${path_QCD}\", \"${path_VV}\", \"${path_VJ}\", \"${ttbar_path}\", \"${path_ttHto2B}\", \"${output_path}\", \"${var}\", \"${var_label}\", \"${BIN}\")"
    done
}

process_year() {
    year=$1;
    get_paths $year

    # prepare histograms from trees
    echo "Processing data and MC for year ${year}..."
    process_data ${year}
    process_mc ${year}
    wait
    echo "Data and MC processing completed for year ${year}!"

    # derive prescale factors
    python3 ${SCRIPT_DIR}/kfact.py --year ${year}
    plot_year ${year}

    # make histograms with tau32 SF to TTbar
    input_path="${TREE_DIR}/Histograms_${year}_MC_TTbar.root"
    output_path="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32.root"
    root -l -b -q "make_hist.cpp(\"${input_path}\", \"${output_path}\", \"${path_sf_tau32}\")"
    
    # derive prescale factors for tau32+TXbb
    python3 ${SCRIPT_DIR}/kfact.py --year ${year} --ttbar-corr "tau32"
    plot_year ${year} "tau32"
    
    # make histograms with tau32 and TXbb SF to TTbar
    input_path="${TREE_DIR}/Histograms_${year}_MC_TTbar.root"
    output_path="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb.root"
    root -l -b -q "make_hist.cpp(\"${input_path}\", \"${output_path}\", \"${path_sf_tau32}\", \"${path_sf_TXbb}\")"

    # derive prescale factors for tau32+TXbb
    python3 ${SCRIPT_DIR}/kfact.py --year ${year} --ttbar-corr "tau32_TXbb"
    plot_year ${year} "tau32_TXbb"

    # make histograms with tau32+TXbb+pTjj SF to TTbar
    input_path="${TREE_DIR}/Histograms_${year}_MC_TTbar.root"
    output_path="${HIST_DIR}/Histograms_${year}_MC_TTbar_tau32_TXbb_pTjj.root"
    root -l -b -q "make_hist.cpp(\"${input_path}\", \"${output_path}\", \"${path_sf_tau32}\", \"${path_sf_TXbb}\", \"${path_sf_pTjj}\")"

    # derive prescale factors for tau32+TXbb+pTjj
    python3 ${SCRIPT_DIR}/kfact.py --year ${year} --ttbar-corr "tau32_TXbb_pTjj"
    plot_year ${year} "tau32_TXbb_pTjj"
    
}

if [ $# -ne 1 ]; then
    # process all years
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