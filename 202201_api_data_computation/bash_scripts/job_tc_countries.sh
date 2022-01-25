#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
conda activate climada_env

bsub -W 20:00 -R "rusage[mem=20000]" python3 python_scripts/make_tc_countries.py