#!/bin/bash

.github/scripts/check_outputs.py \
    --analysis-outdir-upstream .github/data/test_output/KevinKuchinski-FluViewer-output \
    --analysis-outdir-origin .github/data/test_output/BCCDC-PHL-FluViewer-output \
    --outdir artifacts
