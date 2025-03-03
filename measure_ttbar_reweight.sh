#!/bin/bash


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"/scripts;
cd $SCRIPT_DIR; cmsenv;

# SF(mass, pt)
cd $SCRIPT_DIR/sf_mass_pt; ./histo_mc_ttbar_reweight.csh; ./trig_eff_mass_pt_ttbar_reweight.csh;

# SF(TXbb)
cd $SCRIPT_DIR/sf_TXbb; ./histo_mc_ttbar_reweight.csh; ./trig_eff_TXbb_ttbar_reweight.csh;

