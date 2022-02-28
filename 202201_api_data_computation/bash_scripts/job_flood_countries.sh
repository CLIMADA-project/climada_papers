#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
conda activate climada_env

for year in 2010 2030 2050 2070
do
for rcp in rcp26 rcp60 rcp85
do
year2=$((year + 20))
bsub -W 20:00 -R "rusage[mem=20000]" python3 ../python_scripts/compute_flood_countries.py $year $year2 $rcp
done
done