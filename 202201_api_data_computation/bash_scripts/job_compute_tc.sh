#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
conda activate climada_env

for basin in "SI" "NI" "NA" "SP" "WP" "SA" "EP"

do

        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        echo $basin

        bsub -R -W 10:00 "rusage[mem=20000]" python3 ../python_scripts/compute_tc_genesis_basin.py $basin 10  ''

done

for basin in "SI" "NI" "NA" "SP" "WP" "SA" "EP"

do

        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        echo $basin

        bsub -R "rusage[mem=20000]" python3 ../python_scripts/compute_tc_genesis_basin.py $basin 50  ''

done