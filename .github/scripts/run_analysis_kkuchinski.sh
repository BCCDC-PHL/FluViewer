#!/bin/bash

set -eo pipefail

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda activate fluviewer-kkuchinski

# Check for a sign that we're in the GitHub Actions environment.
# Prevents these settings from being applied in other environments.
if [ -z "${GITHUB_ACTIONS}" ]; then 
    echo "Not running in GitHub Actions environment."
    num_threads=16
else
    echo "Running in GitHub Actions environment."
    num_threads=2
fi
echo "Number of threads used for analysis: ${num_threads}"

database_version="v0.1.8"

mkdir -p .github/data/test_output/fluviewer-kkuchinski

while IFS=, read -r sample_id assembly; do
    echo "Analyzing sample ${sample_id}"
    rm -rf ./${sample_id}
    FluViewer \
	-T ${num_threads} \
	-g \
	-f .github/data/fastq/${sample_id}_R1.fastq.gz \
	-r .github/data/fastq/${sample_id}_R2.fastq.gz \
	-d .github/data/fluviewer_db-${database_version}/FluViewer_db.fa \
	-n ${sample_id}

    mv ${sample_id} .github/data/test_output/fluviewer-kkuchinski/${sample_id}
    echo "Finished analyzing sample ${sample_id}"

done < .github/data/reads_to_simulate.csv
