#!/bin/bash

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda env create -f environment.yml -n fluviewer-bccdc-phl

conda activate fluviewer-bccdc-phl

pip install .
