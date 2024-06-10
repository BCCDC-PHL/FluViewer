#!/bin/bash

# Run the same steps as the CI pipeline, but with the local environment

set -e

rm -rf .github/data/test_output/*

# Install dependencies
art_env_dir="${HOME}/.conda/envs/art"
if [ ! -d ${art_env_dir} ]; then
    .github/scripts/create_art_environment.sh
fi

.github/scripts/download_assemblies.sh

.github/scripts/simulate_reads.sh

.github/scripts/download_fluviewer_db.sh

fluviewer_kkuchinski_env_dir="${HOME}/.conda/envs/fluviewer-kkuchinski"
if [ ! -d ${fluviewer_kkuchinski_env_dir} ]; then
    .github/scripts/install_fluviewer_kkuchinski.sh
fi

.github/scripts/run_analysis_upstream.sh

fluviewer_bccdc_phl_env_dir="${HOME}/.conda/envs/fluviewer-bccdc-phl"
if [ ! -d ${fluviewer_bccdc_phl_env_dir} ]; then
    .github/scripts/install_fluviewer_bccdc-phl.sh
fi

.github/scripts/run_analysis_origin.sh

.github/scripts/check_outputs.sh


# Cleanup

rm -rf .github/data/assemblies/*
rm -rf .github/data/fastq/*
rm -rf .github/data/fluviewer_db-v*
