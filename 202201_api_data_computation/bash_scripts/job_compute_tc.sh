#!/bin/bash
PYTHON_SCRIPT=$1
. ~/.bashrc
. ~/venv/climada_dev/bin/activate
for basin in "WP" #"SA" "SI" #"EP" #"NI" "SI" "NA" "SP"

do

        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        echo $basin

        sbatch -n 1 --cpus-per-task=1 --time=20:00:00 --mem-per-cpu=150000 --wrap="python3 ../python_scripts/compute_tc_genesis_basin.py $basin 10"

done
