#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=30000

. ~/venv/climada_dev/bin/activate
python3 ../python_scripts/save_100y_event.py