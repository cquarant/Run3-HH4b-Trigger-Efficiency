# check if CMSSW_BASE is set
if [ -z "${CMSSW_BASE}" ]; then
    echo "CMSSW_BASE is not set. Please set it first."
    return
fi

PROJ_ROOT="${CMSSW_BASE}/src"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}/output
mkdir -p ${OUTPUT_DIR}


ERA_TAG="2023"
ttbar_type="TTtoLNu2Q"
channels=("Muon" "EGamma" "Leptonic")

for channel in ${channels[@]}; do
    sample_dir="/eos/uscms/store/group/lpcdihiggsboost/sixie/analyzer/HHTo4BNtupler/ArmenVersion/nano/run3/combined"
    sample_path=${sample_dir}/${ERA_TAG}/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8.root

    output_path=${OUTPUT_DIR}/Histograms_MC_TTtoLNu2Q_${ERA_TAG}_${channel}.root

    pu_path=${PROJ_ROOT}/pileups/pu_${ERA_TAG}.txt
    param_path=${PROJ_ROOT}/parameters/parameters_${ERA_TAG}.txt
    jec_path=${PROJ_ROOT}/JECs
    jec_path_ak4=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi.txt
    jec_path_ak8=${jec_path}/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_L2Relative_AK8PFPuppi.txt

    sf_path=${PROJ_ROOT}/Trigger_SFs/PT_Mass_SF_TTbar/PT_Mass_2dSF_2023.root


    root -l -b -q "histo_mc.cpp(\"${ttbar_type}\", \"${channel}\", \"${sample_path}\", \"${output_path}\", \"${pu_path}\", \"${sf_path}\", \"${param_path}\", \"${jec_path_ak4}\", \"${jec_path_ak8}\")"
done
