#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
conda activate climada_env

for basin in "NI" "SI" "NA" "SP" "WP" "SA" "EP"

do

        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        echo $basin

        bsub -n 1 -W 20:00 -R "rusage[mem=100000]" python3 ../python_scripts/compute_tc_genesis_basin.py $basin 10  ''

done

for basin in "SI" "NI" "NA" "SP" "WP" "SA" "EP"

do

        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        echo $basin

#        bsub -W 60:00 -R "rusage[mem=100000]" python3 ../python_scripts/compute_tc_genesis_basin.py $basin 50  ''

done