#!/bin/bash -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REGIONS=("QCD" "TTbar")

for region in ${REGIONS[@]}; do
    python $SCRIPT_DIR/trigger_weight.py --region ${region}
done
