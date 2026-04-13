#!/bin/bash
#SBATCH --cpus-per-task=2          # Number of cores
#SBATCH --mem=8G                   # Memory per node
#SBATCH --time=20:00:00             # Runtime limit (HH:MM:SS)
#SBATCH --job-name=amrfinder       # Job name
#SBATCH --output=amrfinder_%j.out  # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run pipelines/amrfinderflow \
  --input amrfinderflow/bin_samplesheet.csv \
  --fastqs amrfinderflow/fastqs_sheet.csv \
  --outdir amrfinderflow_results \
  -profile snellius \
  -resume