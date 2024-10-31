# Trigger Efficiency Measurement for Run-3 Boosted $HH \to b \bar{b} b \bar{b}$ Analysis at the CMS Experiment
- Make histogram: `HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/<Data or MC>/SemiL_TTbar`
    - `Online_BTG_TTbar_Legacy`
    - `PT_Mass_SF`
- Extract the corrections and make the plots: `HH4b_Stage_2_Armen`

# Relevant Code Generation Files
## Data


## MC
- `Online_BTG_Legacy`
    - [preBPix](HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/MC/Online_BTG_Legacy/preBPix/code_gen.py)
    - [postBPix](HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/MC/Online_BTG_Legacy/postBPix/code_gen.py)
- [`Online_BTG_TTbar_Legacy`](HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/MC/Online_BTG_TTbar_Legacy/code_gen.py)
- `PT_Mass_SF`
    - [preBPix](HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/MC/PT_Mass_SF/preBPix/code_gen.py)
    - [postBPix](HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/MC/PT_Mass_SF/postBPix/code_gen.py)
- [`PT_Mass_SF_TTbar`](HHBoostedAnalyzer_Armen/MyHisto/NanoAOD_V12/MC/PT_Mass_SF_TTbar/code_gen.py)