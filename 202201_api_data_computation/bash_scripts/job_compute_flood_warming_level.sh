#!/bin/bash
. ~/.bashrc
. ~/venv/climada_dev/bin/activate

for warming_level in 1 2
do
sleep_time=$((($RANDOM % 5) + 5))

echo $sleep_time

sleep $sleep_time

bsub -W 4:00 -R "rusage[mem=40000]" python3 ../python_scripts/compute_river_flood_degree_warming.py $warming_level ''
done