#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate rpackages
Rscript /projects/0/prjs0784/helius_longitudinal/scripts/ordination_calc.R