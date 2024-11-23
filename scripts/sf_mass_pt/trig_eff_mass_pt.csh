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
mkdir -p ${OUTPUT_DIR}
mkdir -p ${FIG_DIR}

era_tags=("2023" "2023BPix")
for era_tag in ${era_tags[@]}; do
    echo "Era: ${era_tag}"
    data_types=("QCD" "TTbar")
    for data_type in ${data_types[@]}; do
        echo "Type: ${data_type}"
        hist_mc_path="${OUTPUT_DIR}/Histograms_${era_tag}_MC_${data_type}.root"
        hist_data_path="${OUTPUT_DIR}/Histograms_${era_tag}_data_${data_type}.root"
        output_path="${OUTPUT_DIR}/efficiency_mass_pt_${era_tag}_${data_type}.root"
        figure_mc_path="${FIG_DIR}/efficiency_mass_pt_${era_tag}_MC_${data_type}.pdf"
        figure_data_path="${FIG_DIR}/efficiency_mass_pt_${era_tag}_data_${data_type}.pdf"
        figure_sf_path="${FIG_DIR}/efficiency_mass_pt_${era_tag}_SF_${data_type}.pdf"
        root -l -b -q "trig_eff_mass_pt.cpp(\"${hist_mc_path}\", \"${hist_data_path}\", \"${output_path}\", \"${figure_mc_path}\", \"${figure_data_path}\", \"${figure_sf_path}\")"
    done

    # together
    hist_mc_path="${OUTPUT_DIR}/Histograms_${era_tag}_MC.root"
    hist_data_path="${OUTPUT_DIR}/Histograms_${era_tag}_data.root"
    output_path="${OUTPUT_DIR}/efficiency_mass_pt_${era_tag}.root"
    figure_mc_path="${FIG_DIR}/efficiency_mass_pt_${era_tag}_MC.pdf"
    figure_data_path="${FIG_DIR}/efficiency_mass_pt_${era_tag}_data.pdf"
    figure_sf_path="${FIG_DIR}/efficiency_mass_pt_${era_tag}_SF.pdf"
    root -l -b -q "trig_eff_mass_pt.cpp(\"${hist_mc_path}\", \"${hist_data_path}\", \"${output_path}\", \"${figure_mc_path}\", \"${figure_data_path}\", \"${figure_sf_path}\")"
done