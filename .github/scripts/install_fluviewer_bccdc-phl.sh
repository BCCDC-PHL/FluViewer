#!/bin/bash

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda env create -f .github/environments/fluviewer-bccdc-phl.yml -n fluviewer-bccdc-phl

conda activate fluviewer-bccdc-phl

pip install .
