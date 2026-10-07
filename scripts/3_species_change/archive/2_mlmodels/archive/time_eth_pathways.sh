#!/bin/bash
#SBATCH -c 96
#SBATCH --mem=68G
#SBATCH --time=5:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
# python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
#     -name timepoints \
#     -path /projects/0/prjs0784/helius_longitudinal/timepoint_pathways \
#     -x class \
#     -n 200 \
#     -rand_seed 4321 \
#     -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json
# python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
#     -name dutch_timepoints \
#     -path /projects/0/prjs0784/helius_longitudinal/timepoint_dutch_pathways \
#     -x class \
#     -n 200 \
#     -rand_seed 4321 \
#     -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json 
# python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
#     -name sas_timepoints \
#     -path /projects/0/prjs0784/helius_longitudinal/timepoint_sas_pathways \
#     -x class \
#     -n 200 \
#     -rand_seed 4321 \
#     -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json 
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
    -name eth_baseline \
    -path /projects/0/prjs0784/helius_longitudinal/eth_base_pathways \
    -x class \
    -n 200 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json
python /projects/0/prjs0784/helius_longitudinal/scripts/XGBeast.py \
    -name eth_followup \
    -path /projects/0/prjs0784/helius_longitudinal/eth_fu_pathways \
    -x class \
    -n 200 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/helius_longitudinal/scripts/param_grid_mb.json