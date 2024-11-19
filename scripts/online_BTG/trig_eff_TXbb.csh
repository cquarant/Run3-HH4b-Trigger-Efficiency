# check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}/output
FIG_DIR=${SCRIPT_DIR}/figures
mkdir -p ${OUTPUT_DIR}
mkdir -p ${FIG_DIR}

era_tag="2023"
channels=("Muon" "EGamma" "Leptonic")

for channel in ${channels[@]}; do
    echo "channel: ${channel}"
    hist_mc_path="${OUTPUT_DIR}/Histograms_MC_TTtoLNu2Q_${era_tag}_${channel}.root"
    hist_data_path="${OUTPUT_DIR}/Histograms_Data_${era_tag}_${channel}.root"
    output_path="${OUTPUT_DIR}/efficiency_TXbb_${era_tag}_${channel}.root"
    figure_mc_path="${FIG_DIR}/efficiency_TXbb_mc_${era_tag}_${channel}.pdf"
    figure_data_path="${FIG_DIR}/efficiency_TXbb_data_${era_tag}_${channel}.pdf"
    root -l -b -q "trig_eff_TXbb.cpp(\"${hist_mc_path}\", \"${hist_data_path}\", \"${output_path}\", \"${figure_mc_path}\", \"${figure_data_path}\")"
done