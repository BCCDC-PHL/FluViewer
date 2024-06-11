#!/usr/bin/env bash

.github/scripts/check_outputs.py \
    --analysis-outdir-upstream .github/data/test_output/KevinKuchinski-FluViewer-output \
    --analysis-outdir-origin .github/data/test_output/BCCDC-PHL-FluViewer-output \
    --outdir artifacts

column -ts ',' artifacts/check_outputs_summary.csv

while read -r test_name test_result; do
    if [ "$test_result" == "FAIL" ]; then
	echo "Test $test_name failed"
	exit 1
    fi
done < tail -n+2 artifacts/check_outputs_summary.csv
