# check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}/output
TMP_DIR=${SCRIPT_DIR}/tmp
FIG_DIR=${SCRIPT_DIR}/figures
FIG_TMP_DIR=${SCRIPT_DIR}/figures
mkdir -p ${OUTPUT_DIR}
mkdir -p ${FIG_DIR}
mkdir -p ${FIG_TMP_DIR}

DATA_TYPES=("QCD" "TTbar")
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
            root -l -b -q "trig_eff_TXbb.cpp(\"${hist_mc_path}\", \"${hist_data_path}\", \"${tagger_name}\", \"${output_path}\", \"${figure_path}\")"
        done

        # together
        hist_mc_path="${OUTPUT_DIR}/Histograms_${era_tag}_MC.root"
        hist_data_path="${OUTPUT_DIR}/Histograms_${era_tag}_data.root"
        output_path="${OUTPUT_DIR}/efficiency_${tagger_name}_${era_tag}.root"
        figure_path="${FIG_DIR}/efficiency_${tagger_name}_${era_tag}.pdf"
        root -l -b -q "trig_eff_TXbb.cpp(\"${hist_mc_path}\", \"${hist_data_path}\", \"${tagger_name}\", \"${output_path}\", \"${figure_path}\")"
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