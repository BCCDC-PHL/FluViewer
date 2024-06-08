#!/bin/bash


source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate art

mkdir -p .github/data/fastq

while IFS=',' read -r sample_id assembly; do
    art_illumina \
	--paired \
	--in ${assembly} \
	--fcov 500 \
	--len 150 \
	--mflen 400 \
	--sdev 100 \
	--rndSeed 42 \
	--qShift 0 \
	--qShift2 0 \
	--out .github/data/fastq/${sample_id}_R

    rm -f .github/data/fastq/${sample_id}_R1.aln
    rm -f .github/data/fastq/${sample_id}_R2.aln

    mv .github/data/fastq/${sample_id}_R1.fq .github/data/fastq/${sample_id}_R1.fastq
    mv .github/data/fastq/${sample_id}_R2.fq .github/data/fastq/${sample_id}_R2.fastq

    gzip -f .github/data/fastq/${sample_id}_R1.fastq
    gzip -f .github/data/fastq/${sample_id}_R2.fastq

done < .github/data/reads_to_simulate.csv

