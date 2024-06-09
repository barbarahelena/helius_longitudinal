#!/bin/bash
#SBATCH -c 96
#SBATCH --mem=35G
#SBATCH --time=10:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast_adj.py \
    -name hypertension \
    -path /projects/0/prjs0784/helius_longitudinal/hypertension \
    -x class \
    -cv 3 \
    -test 0.4 \
    -n 30 \
    -tune \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/grid_tuning.json