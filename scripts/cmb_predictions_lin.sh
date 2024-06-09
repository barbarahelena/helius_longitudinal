#!/bin/bash
#SBATCH -c 96
#SBATCH --mem=35G
#SBATCH --time=10:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast_adj.py \
    -name sbp \
    -path /projects/0/prjs0784/helius_longitudinal/sbp \
    -x reg \
    -n 10 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb_2.json
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast_adj.py \
    -name dbp \
    -path /projects/0/prjs0784/helius_longitudinal/dbp \
    -x reg \
    -n 10 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb_2.json
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast_adj.py \
    -name hba1c \
    -path /projects/0/prjs0784/helius_longitudinal/hba1c \
    -x reg \
    -n 10 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb_2.json
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast_adj.py \
    -name ldl \
    -path /projects/0/prjs0784/helius_longitudinal/ldl \
    -x reg \
    -n 10 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb_2.json