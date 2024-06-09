#!/bin/bash
#SBATCH -c 96
#SBATCH --mem=68G
#SBATCH --time=5:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
    -name diabetes \
    -path /projects/0/prjs0784/helius_longitudinal/diabetes \
    -x class \
    -test 0.4 \
    -n 20 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
    -name hypertension \
    -path /projects/0/prjs0784/helius_longitudinal/hypertension \
    -x class \
    -test 0.4 \
    -n 20 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json 
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
    -name metsyn \
    -path /projects/0/prjs0784/helius_longitudinal/metsyn \
    -x class \
    -test 0.4 \
    -n 20 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json 
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
    -name lld \
    -path /projects/0/prjs0784/helius_longitudinal/lld \
    -x class \
    -test 0.4 \
    -n 20 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json
    
# python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
#     -name diabetes \
#     -path /projects/0/prjs0784/helius_longitudinal/diabetes \
#     -x class \
#     -n 200 \
#     -rand_seed 4321 \
#     -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json \
#     -permute
# python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
#     -name hypertension \
#     -path /projects/0/prjs0784/helius_longitudinal/hypertension \
#     -x class \
#     -n 200 \
#     -rand_seed 4321 \
#     -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json \
#     -permute
# python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
#     -name metsyn \
#     -path /projects/0/prjs0784/helius_longitudinal/metsyn \
#     -x class \
#     -n 200 \
#     -rand_seed 4321 \
#     -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json \
#     -permute
# python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
#     -name lld \
#     -path /projects/0/prjs0784/helius_longitudinal/lld \
#     -x class \
#     -n 200 \
#     -rand_seed 4321 \
#     -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json \
#     -permute