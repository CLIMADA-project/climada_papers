#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10000

for basin in "NI" "SI" "NA" "SP" "WP" "SA" "EP"
do

        sleep_time=$((($RANDOM % 5) + 5))

        echo $sleep_time

        sleep $sleep_time

        echo $basin

        python3 ../python_scripts/compute_tc_tracks.py $basin 10  ''

done
