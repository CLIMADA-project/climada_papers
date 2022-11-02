#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
. ~/venv/climada_dev/bin/activate
sleep_time=$((($RANDOM % 5) + 5))

echo $sleep_time

sleep $sleep_time

year=1980
year2=2000
scenario=hist

bsub -W 4:00 -R "rusage[mem=40000]" python3 ../python_scripts/compute_river_flood.py $year $year2 $scenario ''
