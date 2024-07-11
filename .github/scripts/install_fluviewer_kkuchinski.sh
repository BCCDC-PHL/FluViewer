#!/bin/bash

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda env create -f .github/environments/fluviewer-kkuchinski.yaml
