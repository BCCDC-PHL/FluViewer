#!/bin/bash

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

fluviewer_version="0.1.11"

conda env create -f environment.yaml -n fluviewer-kkuchinski

conda activate fluviewer-kkuchinski

pip install fluviewer==${fluviewer_version}
