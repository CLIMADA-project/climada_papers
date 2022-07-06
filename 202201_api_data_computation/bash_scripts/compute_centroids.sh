PYTHON_SCRIPT=$1
. ~/.bashrc
conda activate climada_env

bsub -W 20:00 -n 20 python3 ../python_scripts/compute_centroids.py