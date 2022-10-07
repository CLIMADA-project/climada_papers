#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
conda activate climada_env

for year in 2005 2030 2050 2070
do
year2=$((year + 20))
for rcp in rcp26 rcp60 rcp85
do
        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        bsub -W 4:00 -R "rusage[mem=40000]" python3 ../python_scripts/compute_river_flood.py $year $year2 $rcp ''
done
done