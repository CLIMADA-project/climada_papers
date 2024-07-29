PYTHON_SCRIPT=$1
. ~/.bashrc
. ~/venv/climada_dev/bin/activate
bsub -W 60:00 -n 1 python3 ../python_scripts/compute_centroids.py