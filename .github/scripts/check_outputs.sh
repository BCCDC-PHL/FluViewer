#!/usr/bin/env bash

set -e -o pipefail

source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate check-outputs

.github/scripts/check_outputs.py \
    --kkuchinski-outdir .github/data/test_output/fluviewer-kkuchinski \
    --bccdc-phl-outdir .github/data/test_output/fluviewer-bccdc-phl \
    -o artifacts/check_outputs_results.csv
