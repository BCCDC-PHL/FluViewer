#!/bin/bash

set -eo pipefail

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda activate fluviewer-bccdc-phl

# Check for a sign that we're in the GitHub Actions environment.
# Prevents these settings from being applied in other environments.
if [ -z "${GITHUB_ACTIONS}" ]; then 
    echo "Not running in GitHub Actions environment."
else
    echo "Running in GitHub Actions environment."
fi

database_version="v0.1.8"

mkdir -p .github/data/test_output/fluviewer-bccdc-phl

while IFS=, read -r sample_id assembly; do
    echo ${sample_id}
    fluviewer \
	--forward-reads .github/data/fastq/${sample_id}_R1.fastq.gz \
	--reverse-reads .github/data/fastq/${sample_id}_R2.fastq.gz \
	--database .github/data/fluviewer_db/fluviewer_db \
	--outdir .github/data/test_output/fluviewer-bccdc-phl \
	--output-name ${sample_id}

done < .github/data/reads_to_simulate.csv    
    
