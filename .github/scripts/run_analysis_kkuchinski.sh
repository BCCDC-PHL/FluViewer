#!/bin/bash

set -eo pipefail

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda activate fluviewer-kkuchinski

# Check for a sign that we're in the GitHub Actions environment.
# Prevents these settings from being applied in other environments.
if [ -z "${GITHUB_ACTIONS}" ]; then 
    echo "Not running in GitHub Actions environment."
else
    echo "Running in GitHub Actions environment."
fi

database_version="v0.1.8"

mkdir -p .github/data/test_output_kkuchinski

while IFS=, read -r sample_id assembly; do
    echo ${sample_id}
    FluViewer \
	-f .github/data/fastq/${sample_id}_R1.fastq.gz \
	-r .github/data/fastq/${sample_id}_R2.fastq.gz \
	-d .github/data/fluviewer_db/fluviewer_db \
	-n ${sample_id}

done < .github/data/reads_to_simulate.csv    
    