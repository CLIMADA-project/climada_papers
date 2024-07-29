#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=20000

. ~/venv/climada_dev/bin/activate

for basin in "NI" "SI" "NA" "SP" "WP" "SA" "EP"

do

        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        echo $basin

       python3 ../python_scripts/compute_tc_climate_change.py $basin 10  ''

done
