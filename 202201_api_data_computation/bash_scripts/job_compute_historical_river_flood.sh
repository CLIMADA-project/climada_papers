#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=60000

. ~/venv/climada_dev/bin/activate

year=1980
year2=2010
scenario=hist

python3 ../python_scripts/compute_river_flood.py $year $year2 $scenario
