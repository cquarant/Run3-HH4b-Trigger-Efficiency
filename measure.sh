#!/bin/bash


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"/scripts;
cd $SCRIPT_DIR; cmsenv;

# SF(mass, pt)
cd $SCRIPT_DIR/sf_mass_pt; ./histo_data.csh; ./histo_mc.csh; ./trig_eff_mass_pt.csh;

# SF(TXbb)
cd $SCRIPT_DIR/sf_TXbb; ./histo_data.csh; ./histo_mc.csh; ./trig_eff_TXbb.csh;

# ttbar background correction
cd $SCRIPT_DIR/semi-leptonic; ./tree_data.csh; ./tree_mc.csh; ./merge_trees.sh; ./derive_SF.sh;
