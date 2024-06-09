#!/bin/bash

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda env create -f environment.yaml -n fluviewer-bccdc-phl

conda activate fluviewer-bccdc-phl

pip install .
