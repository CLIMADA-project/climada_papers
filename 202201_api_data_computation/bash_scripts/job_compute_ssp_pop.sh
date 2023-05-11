#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
. ~/venv/climada_dev/bin/activate
bsub -W 4:00 -R "rusage[mem=20000]" python3 ../python_scripts/make_ssp_pop.py