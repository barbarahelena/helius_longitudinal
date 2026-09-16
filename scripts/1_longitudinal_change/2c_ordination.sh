#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate rpackages
Rscript /projects/0/prjs0784/helius_longitudinal/scripts/1_longitudinal_change/2c_ordination.R
