#!/bin/bash


source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate sra-tools

mkdir -p .github/data/fastq

pushd .github/data/fastq

while IFS=',' read -r sample_id; do
    fasterq-dump ${sample_id}
    gzip *.fastq
done < .github/data/sra_samples_to_download.csv

popd

