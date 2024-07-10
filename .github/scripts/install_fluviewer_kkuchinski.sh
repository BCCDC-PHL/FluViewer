#!/bin/bash

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

fluviewer_version="0.1.11"

conda env create -f .github/environments/fluviewer-kkuchinski.yml -n fluviewer-kkuchinski
