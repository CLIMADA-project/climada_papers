#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=20000

. ~/venv/climada_dev/bin/activate

python3 ../python_scripts/tc_concat_basins.py
