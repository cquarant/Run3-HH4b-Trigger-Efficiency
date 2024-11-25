#!/bin/bash

eras=("2022C" "2022D" "2022E" "2022F" "2022G" "2023C" "2023D")
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

for era in ${eras[@]}; do
    echo "Processing ${era}"
    input=${SCRIPT_DIR}/jsons/GoodLumiList_Map_${era}.json
    output=${SCRIPT_DIR}/GoodLumiList_Map_${era}.txt
    python3 ${SCRIPT_DIR}/process_json.py --input ${input} --output ${output}
done