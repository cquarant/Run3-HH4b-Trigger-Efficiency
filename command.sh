# error if base not set up
if [ -z "$CMSSW_BASE" ]; then
    echo "CMSSW_BASE not set up. Run cmsenv in a CMSSW release."
    return 1
fi

set -xe

############### Step 1: Trigger Efficiency Histograms ###############
proj_dir=$CMSSW_BASE/src/HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12
# Online_BTG_Legacy (MC)
cd $proj_dir/MC/Online_BTG_Legacy/preBPix; ./run_jobs.csh
cd $proj_dir/MC/Online_BTG_Legacy/postBPix; ./run_jobs.csh

# Online_BTG_TTbar_Legacy (MC + Data)
cd $proj_dir/MC/Online_BTG_TTbar_Legacy; ./run_jobs.csh
cd $proj_dir/Data/Online_BTG_TTbar_Legacy; ./run_jobs.csh

# PT_Mass_SF (MC + Data)
cd $proj_dir/MC/PT_Mass_SF/preBPix; ./run_jobs.csh
cd $proj_dir/MC/PT_Mass_SF/postBPix; ./run_jobs.csh
cd $proj_dir/Data/PT_Mass_SF; ./run_jobs.csh

# PT_Mass_SF_TTbar (MC)
cd $proj_dir/MC/PT_Mass_SF_TTbar; ./run_jobs.csh