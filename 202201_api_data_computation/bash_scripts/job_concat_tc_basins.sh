#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
. ~/venv/climada_dev/bin/activate
bsub -R "rusage[mem=20000]" python3 ../python_scripts/tc_concat_basins.py
