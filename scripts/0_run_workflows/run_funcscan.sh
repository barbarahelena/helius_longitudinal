#!/bin/bash
#SBATCH --cpus-per-task=2       # Number of cores
#SBATCH --mem=8G                # Memory per node
#SBATCH --time=2:00:00         # Runtime limit (HH:MM:SS)
#SBATCH --job-name=funcscan       # Job name
#SBATCH --output=funcscan_%j.out  # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run pipelines/funcscan \
  --input funcscan/bin_samplesheet.csv \
  --outdir funcscan_results \
  -profile snellius \
  --run_arg_screening \
  --arg_skip_abricate \
  --arg_skip_rgi \
  --arg_skip_deeparg \
  --arg_skip_fargene \
  --arg_skip_argnorm