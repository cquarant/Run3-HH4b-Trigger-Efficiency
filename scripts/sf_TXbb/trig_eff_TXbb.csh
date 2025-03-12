# check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}/output
TMP_DIR=${SCRIPT_DIR}/tmp
# FIG_DIR=${SCRIPT_DIR}/figures
# BINS="0.8,0.9,0.94,0.97,0.99,1.0"
BINS="0.8,0.825,0.85,0.875,0.9,0.92,0.94,0.97,0.99,1.0"
FIG_DIR=${SCRIPT_DIR}/figures_zoomed
mkdir -p ${OUTPUT_DIR}
mkdir -p ${FIG_DIR}

DATA_TYPES=("QCD" "TTbar")
# DATA_TYPES=("TTbar" "QCD")
TAGGER_NAMES=("GloParT" "PNetLegacy")

process_era() {
    local era_tag=$1
    echo "Era: ${era_tag}"
    
    for tagger_name in ${TAGGER_NAMES[@]}; do
        for data_type in ${DATA_TYPES[@]}; do
            echo "Type: ${data_type}"
            hist_mc_path="${OUTPUT_DIR}/Histograms_${era_tag}_MC_${data_type}.root"
            hist_data_path="${OUTPUT_DIR}/Histograms_${era_tag}_data_${data_type}.root"
            output_path="${OUTPUT_DIR}/efficiency_${tagger_name}_${era_tag}_${data_type}.root"
            figure_path="${FIG_DIR}/efficiency_${tagger_name}_${era_tag}_${data_type}.pdf"
            root -l -b -q "trig_eff_TXbb.cpp(\"${hist_mc_path}\", \"${hist_data_path}\", \"${tagger_name}\", \"${output_path}\", \"${figure_path}\", \"${BINS}\")"
        done
    done
}

if [ $# -ne 1 ]; then
    # process all eras
    era_tags=("2022" "2022EE" "2023" "2023BPix")
    for era_tag in ${era_tags[@]}; do
        process_era ${era_tag}
    done
    wait
    echo "All done"
else
    # process single era
    process_era $1
fi