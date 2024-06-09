#!/usr/bin/env bash

set -e -o pipefail

.github/scripts/check_outputs.py \
    --analysis-outdir-kkuchinski .github/data/test_output/fluviewer-kkuchinski \
    --analysis-outdir-bccdc-phl .github/data/test_output/fluviewer-bccdc-phl \
    -o artifacts/check_outputs_results.csv
