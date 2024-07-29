#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=80000


. ~/venv/climada_dev/bin/activate

python3 ../python_scripts/compute_yearsets_warming.py
