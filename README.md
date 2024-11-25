# Trigger Efficiency Measurement for Run-3 Boosted $HH \to b \bar{b} b \bar{b}$ Analysis at the CMS Experiment
## Steps
1. Derive the trigger efficiencies and scale factors as a function of $(m_\mathrm{SD}, p_\mathrm{T})$: [scripts/sf_mass_pt](scripts/sf_mass_pt)
    - Run `scripts/sf_mass_pt/histo_data.csh` and `scripts/sf_mass_pt/histo_mc.csh`
    - Run `trig_eff_mass_pt.csh`
2. Derive the trigger efficiencies and scale factors as a function of $(m_\mathrm{SD}, p_\mathrm{T}, \eta)$: [scripts/sf_TXbb](scripts/sf_TXbb)
    - Run `scripts/sf_mass_pt/scripts/sf_TXbb/histo_data.csh` and `scripts/sf_TXbb/histo_mc.csh`
    - Run `scripts/sf_mass_pt/trig_eff_TXbb.csh`


## References
- [Run-3 Analysis Twiki](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis)
    - Get the luminosity and era tags there; needed for [`parameters/`](parameters/)
- [How to work with files for Good Luminosity Sections in JSON format](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile)
    - Run 3 data EOS: `/eos/user/c/cmsdqm/www/CAF/certification/`
    - [Run 3 data web](https://cms-service-dqmdc.web.cern.ch/CAF/certification/) 
    - Good Lumi List in [`GoodLumiList/`](GoodLumiList/) (`*Golden.json`)
        - [2022](https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/)
        - [2023](https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/)
- [Jet Energy Resolution and Corrections Analysis Recommendation](https://cms-jerc.web.cern.ch/Recommendations/)
    - [txt files](https://github.com/cms-jet/JECDatabase/tree/master/textFiles)